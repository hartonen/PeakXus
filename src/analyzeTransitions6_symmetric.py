# Copyright (c) Tuomo Hartonen, 2015-2016
#
# THIS PROGRAM COMES WITH ABSOLUTELY NO WARRANTY!
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 2 as 
# published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

from ReadContainer6 import ReadContainer
from operator import itemgetter
from scipy import stats
from time import clock
from math import factorial
from math import log
from random import randint

import multiprocessing as mp
import numpy as np
import csv

class tPoints:

#CONSTRUCTOR
    def __init__(self,w,l_limit,u_limit,mindist,s):
        self.w = w #2w+1 is the width of region around the 5'-end of a read at transition point
        self.l_limit = l_limit #i-l_limit is the lowest possible index for 5'-end of + reads when i is the 5'-end of the - reads
        self.u_limit = u_limit #i-u_limit is the highest possible index for 5'-end of the + reads
        self.mindist = mindist #minimum distance between peaks
        self.s = s #sliding window size used to calculate the distribution of 5'-ends

        self.tPoints = {}
        self.num_tPoints = 0
        self.numwindows = {}
        self.read_counts_plus = {}
        self.read_counts_minus = {}
        self.numtests = 0
        self.umi_count = 0
        #key = chromosome name
        #value = numpy array, [c_(0),...,c_(P-1),index,read count,p-value,cluster id]

        #this is just for testing
        self.transitions = []
    #end __init__

    def getUMIcount(self): return self.umi_count
    def getNumTests(self): return self.numtests
    def calcTransitions_counts(self,chrom,reads,allpairs,nproc):

        #calculating the densities of 5'-ends on both strands 
        self.read_counts_plus[chrom] = reads.getAllReadHisto(chrom,'++',self.s,nproc)
        self.read_counts_minus[chrom] = reads.getAllReadHisto(chrom,'--',self.s,nproc)

        #initializing the array that holds the top 1000 peaks per chromosome
        self.tPoints[chrom] = np.zeros((1001,3))
        #Format of the self.tPoints[chrom]-array:
        #[j=left edge of binding site,i=t-point coordinate(right edge of binding site),total 5'-end count]
        latest_i = -200
        for i in range(0,len(self.read_counts_plus[chrom])-100,1):
            current = sum(self.read_counts_plus[chrom][i:i+101])+sum(self.read_counts_minus[chrom][i:i+101])
            
            #checking if the read count is within the top 1000
            if current<self.tPoints[chrom][-1,-1]: continue

            #checking if the current point overlaps with the latest added
            if i-self.tPoints[chrom][latest_i,0]<=100:
                if current>self.tPoints[chrom][latest_i,2]:
                    self.tPoints[chrom][latest_i,0] = i
                    self.tPoints[chrom][latest_i,1] = i+101
                    self.tPoints[chrom][latest_i,2] = current
                    self.tPoints[chrom] = self.tPoints[chrom][self.tPoints[chrom][:,2].argsort()][::-1]

                    latest_i = np.where(self.tPoints[chrom][:,0]==i)[0][0]
            else:
                #adding the current point to the list of top 1000
                self.tPoints[chrom][-1,0] = i
                self.tPoints[chrom][-1,1] = i+101
                self.tPoints[chrom][-1,2] = current
                self.tPoints[chrom] = self.tPoints[chrom][self.tPoints[chrom][:,2].argsort()][::-1]
            
                latest_i = np.where(self.tPoints[chrom][:,0]==i)[0][0]
    #end
        
    def calcTransitions(self,chrom,reads,allpairs,nproc,allreads):
        #chrom = chromosome name
        #reads = ReadContainer-object

        tpoint_count = 0
        N,n = reads.getChromSize(chrom)

        read_counts = reads.getReadHisto(chrom,'-',1,1)
        counts_set = set([(i,read_counts[i]) for i in range(0,len(read_counts)) if read_counts[i]!=0])

        #getAllReadHisto returns histograms that are not filtered for duplicates with UMIs

        if allreads==1: self.read_counts_plus[chrom] = reads.getAllReadHisto(chrom,'++',self.s,nproc)
        else: self.read_counts_plus[chrom] = reads.getReadHisto(chrom,'++',self.s,nproc)
        if allreads==1: self.read_counts_minus[chrom] = reads.getAllReadHisto(chrom,'--',self.s,nproc)
        else: self.read_counts_minus[chrom] = reads.getReadHisto(chrom,'--',self.s,nproc)

        self.umi_count = np.sum(self.read_counts_plus[chrom])+np.sum(self.read_counts_minus[chrom])
        print "UMI count="+str(self.umi_count)+" | ",

        self.numwindows[chrom] = len(read_counts)
        #creating a dictionary where
        #key=index, value=count
        nonzeros = {}
        for s in counts_set: nonzeros[s[0]] = s[1]

        #now all successive index pairs +,- are considered transition points. The exact transition point coo
        #rdinate is the coordinate exactly at halfway between + and -

        #NOTE! For now on we do not consider the beginning and end of the chromosome

        keys = sorted(nonzeros.keys())
        self.tPoints[chrom] = np.zeros((len(keys)/2,3+2*(2+self.l_limit+2*self.w)+2))
        #Format of the self.tPoints[chrom]-array:
        #[j=left edge of binding site,i=t-point coordinate(right edge of binding site),size of the binding site region,c(j-w),...,c(i+w),possible 0's if region is smaller than max size,p-value,peak score]

        for k in range(1,len(keys)):
            #i = current genomic index of the t-point candidate
            i = keys[k]

            if i>N-self.w-self.l_limit/2: break
            if i<self.l_limit+self.w: continue

            #if current read count >0, this is not a t-point
            if nonzeros[i]>0: continue
            #if last read count <0, this is not a t-point
            if nonzeros[keys[k-1]]<0: continue

            #this is just for testing
            self.transitions.append(i)

            #checking if the t-point array is full
            if tpoint_count+(self.l_limit+2*self.w)**2>=len(self.tPoints[chrom]): self.tPoints[chrom] = np.concatenate((self.tPoints[chrom],np.zeros((len(keys)/4,3+2*(2+self.l_limit+2*self.w)+2))))

            #finding position j between i-l_limit and i-u_limit, where c_+(j) is highest
            #if there are no + reads, the t-point candidate is skipped
            if self.read_counts_plus[chrom][i-self.l_limit/2:i].max()<1: continue
            
            #setting i to the middle point between keys[k] and keys[k-1]
            i = (keys[k-1]+i)/2

            if allpairs==1:
                all_j = [j for j in range(i-self.l_limit/2,i-1) if self.read_counts_plus[chrom][j]>0]
                all_i = [i for i in range(i+1,i+self.l_limit/2) if self.read_counts_minus[chrom][i]>0]
            else:
                #the following is used only if we take only the j with highest positive read count and i with highest negative read count into account
                max_indices_plus = np.where(self.read_counts_plus[chrom][i-self.l_limit/2:i][::-1]==max(self.read_counts_plus[chrom][i-self.l_limit/2:i]))[0]
                max_indices_minus = np.where(self.read_counts_minus[chrom][i:i+self.l_limit/2]==max(self.read_counts_minus[chrom][i:i+self.l_limit/2]))[0]

                #if len(max_indices_plus)>1 or len(max_indices_minus)>1: continue
                aux_j = max_indices_plus[0]#np.random.choice(max_indices_plus)
                aux_i = max_indices_minus[0]
                
                #if aux_i+aux_j>=32 and aux_i+aux_j<=34:
                #print str(aux_j)+"\t",
                #print str(aux_i)+"\t",

                aux_j = i-aux_j-1#aux_j+i-(self.l_limit/2)
                #aux_j = max(max_indices_plus)+i-(self.l_limit/2)


                #aux_i = max_indices_minus[0]#np.random.choice(max_indices_minus)
            
                aux_i = aux_i+i
                #print aux_i-aux_j
                #aux_i = min(max_indices_minus)+i+1

                #aux_j = len(self.read_counts_plus[chrom][i-self.l_limit/2:i-1])-self.read_counts_plus[chrom][i-self.l_limit/2:i-1][::-1].argmax()-1+i-(self.l_limit/2+1)
                #aux_i = self.read_counts_minus[chrom][i+1:i+self.l_limit/2].argmax()+i+1


                #all_i = [ind+i+1 for ind in max_indices_minus]
                #all_j = [ind+i-(self.l_limit/2) for ind in max_indices_plus]
                #if np.amax(self.read_counts_plus[chrom][i-self.l_limit/2:i-1])>2 and np.amax(self.read_counts_minus[chrom][i+1:i+self.l_limit/2])>0:
                all_i = [aux_i]
                all_j = [aux_j]
            #adding the region
            for i in all_i:
                for j in all_j:
                    self.tPoints[chrom][tpoint_count,0] = j
                    self.tPoints[chrom][tpoint_count,1] = i
                    self.tPoints[chrom][tpoint_count,2] = 2*(i-j+2*self.w) #2*N = twice the length of the binding site region. i-j+2*w is the size of the binding site region
                    self.tPoints[chrom][tpoint_count,3:3+(i-j+2*self.w)] = self.read_counts_plus[chrom][j-self.w:i+self.w] #positive strand reads
                    self.tPoints[chrom][tpoint_count,3+(i-j+2*self.w):3+2*((i-j+2*self.w))] = self.read_counts_minus[chrom][j-self.w:i+self.w] #negative strand reads

                    tpoint_count += 1
    
        #whole genome analyzed for transition points
        #deleting the unused part of the t-point array
        self.tPoints[chrom] = self.tPoints[chrom][:tpoint_count,:]
        self.num_tPoints += tpoint_count
        
    #end calcTransitions

    def testBackground(self,chrom,p_threshold,yates,nproc,pseudo,allowoverlap):
        #p_threshold = p-value threshold for rejecting the null hypothesis
        #that observed read counts come from Poisson distribution
        #pseudo = pseudocount added to each bin when calculating read distributions
        #for significance testing

        #testing all t-point-regions against the background

        #Testing is done by looking at the difference between positions of reads
        #that point away from the peak summit (background) against the ones that
        #point towards the peak summit (signal)

        pool = mp.Pool(processes=nproc)

        res = [pool.apply_async(backgroundModelPositions,args=(self.tPoints[chrom][k,3:3+int(self.tPoints[chrom][k,2]/2)],self.tPoints[chrom][k,3+int(self.tPoints[chrom][k,2]/2):3+2*int(self.tPoints[chrom][k,2]/2)],p_threshold,(self.tPoints[chrom][k,2]/2),yates,pseudo)) for k in range(0,len(self.tPoints[chrom]))]

        #gathering the results from parallel tests
        res = [p.get() for p in res]
        for k in range(0,len(res)):
            if res[k][0] == True: self.tPoints[chrom][k,-1] = -1
            else:
                current_start = int(self.tPoints[chrom][k,0])
                current_end = int(self.tPoints[chrom][k,1])
                current_middle = int(current_start+current_end)/2
                if yates==0:
                    #G-test
                    self.tPoints[chrom][k,-2] = res[k][1]
                    self.tPoints[chrom][k,-1] = res[k][2]*np.sum(self.read_counts_plus[chrom][current_start:current_middle+1])+np.sum(self.read_counts_minus[chrom][current_middle:current_end+1])-np.sum(self.read_counts_plus[chrom][current_middle:current_end+1])-np.sum(self.read_counts_minus[chrom][current_start:current_middle+1])
                else:
                    #KS-test
                    self.tPoints[chrom][k,-2] = res[k][1]
                    self.tPoints[chrom][k,-1] = res[k][2]
        pool.close()
        pool.terminate()
        pool.join()

        print "num t-points before testing: ",
        print len(self.tPoints[chrom]),
        self.numtests = len(self.tPoints[chrom])
        #removing the random regions
        self.tPoints[chrom] = self.tPoints[chrom][np.where(self.tPoints[chrom][:,-1]>-1)[0],:]
        print " num significant t-points: ",
        print len(self.tPoints[chrom]),

        #removing overlapping regions
        if len(self.tPoints[chrom])>0 and allowoverlap==0: self.delOverlap(allowoverlap)
    #end testBackground

    def delOverlap(self,allowoverlap):
        #this function checks the tPoints-dictionary for overlapping peak regions
        #when regions overlap, only the one with the highest score is accepted

        for c in self.tPoints:
            for i in range(1,len(self.tPoints[c])):
                if self.tPoints[c][i,-2]==-1: continue
                #now looking back for all regions that overlap with the current one
                #bases the binding site candidate i covers
                current_site = set(range(int(self.tPoints[c][i,0]),int(self.tPoints[c][i,1])+1))
                current_start = self.tPoints[c][i,0]
                current_end = self.tPoints[c][i,1]
                current_middle = (current_end+current_start)/2
                current_S = self.tPoints[c][i,-1]

                smallest_p_index = i #index of the largest score
                j = i
                deleted = []
                while j>0:
                    j -= 1
                    if self.tPoints[c][j,-2]==-1: continue
                    #keeping track of all the indices that are to be deleted
                    candidate_start = int(self.tPoints[c][j,0])
                    candidate_end = int(self.tPoints[c][j,1])
                    candidate_middle = (candidate_end+candidate_start)/2
                    candidate_S = self.tPoints[c][j,-1]

                    if len(current_site.intersection(set(range(candidate_start,candidate_end+1))))>0:
                        #this means the sites do overlap
                        if (candidate_S>current_S):
                            deleted.append(smallest_p_index)
                            smallest_p_index = j
                        else: deleted.append(j)
                    else:
                        #this means that there were no overlapping peak candidates
                        break
                #now marking the overlapping peaks with a high p-value as deleted

                if allowoverlap==0:
                    for d in deleted: self.tPoints[c][d,-2] = -1
            
            #deleting the high p-value overlapping peak candidates from chromosome c
            if allowoverlap==0: self.tPoints[c] = self.tPoints[c][np.where(self.tPoints[c][:,-2]>-1)[0],:]
            print ", num t-points after removing overlap: ",
            print len(self.tPoints[c])
    #end delOverlap

    def getTRegions(self,chrom): return self.tPoints[chrom]

    def getTransitions(self,chrom): return self.tPoints[chrom][:,self.P]

    def getNumwindows(self,chrom): return self.numwindows[chrom]
    
    def getNumtPoints(self,chrom): return self.num_tPoints

    def getNumTests(self): return self.numtests

    def saveTransitions(self,filename):

        with open(filename,'wb') as csvfile:
            w = csv.writer(csvfile,delimiter='\t')
            for t in self.transitions: w.writerow([t])

#end class

def backgroundModelPositions(observed_plus,observed_minus,p_threshold,N,yates,pseudo):
    #observed_plus = observed reads on positive strand (histogram of 5' end counts per position)
    #observed_minus = observed reads on negative srand
    #p_threshold = p-value threshold for rejecting the null hypothesis
    #that observed read counts come from Poisson distribution
    #N = length of peak region = i-j+2*w

    #we can use either G-test or chi-squared test
    #+ strand ->
    #           |
    #   (1)     |      (2)
    #----------T|F--------
    #   (4)     |      (3)
    #           |
    #            <- - strand
    #Here we test if positions of reads at regions (1) and (3) come from different
    #distribution than positions of reads at (2) and (4)
    signal_pos = []
    bg_pos = []
    #calculating the signal and background read positions
    for i in range(0,int(N)):
        if i<N/2:
            signal_pos += [i for j in range(0,int(observed_plus[i]))]
            bg_pos += [i for j in range(0,int(observed_minus[i]))]
        else:
            signal_pos += [i-int(N/2) for j in range(0,int(observed_minus[i]))]
            bg_pos += [i-int(N/2) for j in range(0,int(observed_plus[i]))]

    #computing histograms
    if yates<2:
        bins = [i for i in range(0,int(N/2)+1)] #maximum distance from summit is N/2! 
        contig = np.array([np.histogram(signal_pos,bins=bins)[0],np.histogram(bg_pos,bins=bins)[0]])

    #now we can calculate the test statistics
    if yates==1: p,G = chi2Yates(contig[0],contig[1])
    elif yates==0: p,G = gtest(contig[0],contig[1],pseudo)
    else: p,G = kstest_laplace(signal_pos,bg_pos,pseudo,p_threshold,N)

    if p>p_threshold:
        #Reads CAN be explained by the random background
        return (True,p,G)
    else: return (False,p,G)
#end backgroundModelPositions        


def backgroundModelCombined(observed_plus,observed_minus,p_threshold,N,N_avg):
    #THIS FUNCTION IS NOT USED!
    #observed_plus = observed reads on positive strand
    #observed_minus = observed reads on negative srand
    #p_threshold = p-value threshold for rejecting the null hypothesis
    #that observed read counts come from Poisson distribution
    #N = length of peak region
    #N_avg = number of test averaged per peak region

    #we can use either G-test or chi-squared test

    #separating the reads pointing towards the peak from reads pointing
    #away from the peak
    if sum(observed_plus)<2 and sum(observed_minus)<2: return (True, 1.0)
    N = int(N)
    reads_towards_plus = np.zeros(N)
    reads_away_plus = np.zeros(N)
    reads_towards_minus = np.zeros(N)
    reads_away_minus = np.zeros(N)

    for i in range(0,N):
        if i<N/2:
            reads_towards_plus[i] += observed_plus[i]
            reads_away_minus[i] += observed_minus[i]
        else:
            reads_away_plus[i] += observed_plus[i]
            reads_towards_minus[i] +=  observed_minus[i]

    x_toward = sum(reads_towards_plus)+sum(reads_towards_minus)
    x_away = sum(reads_away_plus)+sum(reads_away_minus)
    mean_moved = x_toward-x_away

    #Now drawing the moved reads from "towards-reads" N_avg times with the mean number of
    #moved reads being x_toward-x_away (and uniformly distributed)

    distances_signal = []
    distances_bg = []

    for i in range(0,N_avg):
        if mean_moved>0: r1 = mean_moved#randint(0,mean_moved)
        reads_signal_plus = reads_towards_plus
        reads_signal_minus = reads_towards_minus
        if mean_moved>=(sum(reads_signal_plus)+sum(reads_signal_minus)-1): return (True, 1.0)
        reads_bg_plus = reads_away_plus
        reads_bg_minus = reads_away_minus
        if mean_moved>0:
            moved = 0
            while moved<r1:
                #trying to move a read from signal to background
                r2 = randint(0,N-1)
                if r2<N/2:
                    #this means we are on the left side of the peak
                    if reads_signal_plus[r2]>0:
                        #checking if there is a read to move at this position
                        reads_signal_plus[r2] -= 1
                        reads_bg_plus[r2] += 1
                        moved += 1
                else:
                    if reads_signal_minus[r2]>0:
                        reads_signal_minus[r2] -= 1
                        reads_bg_minus[r2] += 1
                        moved += 1

        #now all the r1 reads are moved and we can calculate the histograms

        for j in range(0,N):
            for k in range(0,N):
                if j==k: continue
                distances_signal += [abs(j-k) for n in range(0,int(reads_signal_plus[j]*reads_signal_minus[k]))]
                distances_bg += [abs(j-k) for n in range(0,int(reads_bg_plus[j]*reads_bg_minus[k]))]

    #now we can calculate the test statistics for aggregated results

    #first we create the same bin edges for both histograms
    if len(distances_signal)<1 or len(distances_bg)<1: return (True, 1.0)
    bin_range = (0,max([max(distances_signal),max(distances_bg)]))
    nbins = 10

    p = gtest(np.histogram(distances_signal,bins=nbins,range=bin_range)[0]/float(N_avg),np.histogram(distances_bg,bins=nbins,range=bin_range)[0]/float(N_avg))
    #p = self.gtest([float(h)/N_avg for h in hist_signal],[float(h)/N_avg for h in hist_bg])

    if p>=p_threshold:
        #Reads CAN be explained by the random background
        return (True,p)
    else: return (False,p)
#end backgroundModelCombined

def gtest(O,E,alpha):
    #performs the goodness-of-fit G-test
    #O = non-normalized frequency distribution of observed events
    #E = non-normalized frequency distribution of expected events
    #O and E have to be of same length
    #alpha = pseudocount
    
    E = E[:]+alpha
    O = O[:]+alpha
    
    dof = len(O)-1
    G = 2*sum([float(O[i])*log(float(O[i])/float(E[i])) for i in range(0,len(O)) if (O[i]>0 and E[i]>0)])

    p = 1-stats.chi2.cdf(G,df=dof)

    return p, G
#end gtest

def kstest(O,E):
    #performs the two sample KS-test
    #O = list of positions of signal reads
    #E = list of positions of background reads
    #O and E need not to have same length

    KS,p = stats.ks_2samp(O,E)

    return p, KS

def kstest_laplace(O,E,alpha,p_thres,N):
    #performs the two sample KS-test with Laplace smoothing
    #O = list of positions of signal reads
    #E = list of positions of background reads
    #alpha = smoothing parameter
    #smoothing:
    # y' = y+alpha/(N+alpha*y)

    #calculating histograms for O and E
    bins = range(0,int(N/2)+1)
    O_hist = np.histogram(O,bins=bins)[0]
    E_hist = np.histogram(E,bins=bins)[0]
    #smoothing
    O_hist = (O_hist[:]+alpha)/(len(O)+alpha*len(O_hist))
    E_hist = (E_hist[:]+alpha)/(len(E)+alpha*len(E_hist))
    #calculating the adjusted sample sizes after smoothing
    nO = np.rint(len(O)+alpha*len(O_hist))
    nE = np.rint(len(E)+alpha*len(E_hist))
    #cumulative distribution functions
    O_cum = np.cumsum(O_hist)
    E_cum = np.cumsum(E_hist)
    max_diff = max(np.abs(O_hist-E_hist))
    #critical value of KS-distribution:
    D = get_KS_2samp_critical(nO,nE,p_thres)
    if max_diff>D:
        #this means that the null hypothesis that the distributions are identical is rejected
        #at level p_thres
        return p_thres, max_diff
    else:
        return 1.0, max_diff
    
def get_KS_2samp_critical(n1,n2,alpha):
    #returns the critical value of Kolmogorov-Smirnov distribution
    #in two sample case for sample sizes n1 and n2(smaller)
    #alpha = confidence level (0.001,0.05,0.1)

    #This is the table of critical values (http://www.soest.hawaii.edu/wessel/courses/gg313/Critical_KS.pdf)
    #key = (n2,n1)
    #value = (critical value for alpha=0.01,critical value for alpha=0.05) (if value=-1, cannot reject null hypothesis)

    C = {(2,8):(10.0,1.0),(2,9):(10.0,1.0),(2,10):(10.0,1.0),(2,11):(10.0,1.0),(2,12):(10.0,1.0),(3,5):(10.0,1.0),(3,6):(10.0,1.0),(3,7):(10.0,1.0),(3,8):(1.0,21.0/24.0),(3,9):(1.0,24.0/27.0),(3,10):(1.0,27.0/30.0),(3,11):(1.0,30.0/33.0),(3,12):(1.0,30.0/36.0),(4,4):(10.0,1.0),(4,5):(10.0,1.0),(4,6):(10.0,20.0/24.0),(4,7):(1.0,24.0/28.0),(4,8):(1.0,28.0/32.0),(4,9):(32.0/36.0,28.0/36.0),(4,10):(36.0/40.0,30.0/40.0),(4,11):(40.0/44.0,33.0/44.0),(4,12):(44.0/48.0,36.0/48.0),(5,6):(1.0,24.0/30.0),(5,7):(1.0,30.0/35.0),(5,8):(35.0/40.0,30.0/40.0),(5,9):(40.0/45.0,35.0/45.0),(5,10):(45.0/50.0,40.0/50.0),(5,11):(45.0/55.0,39.0/55.0),(5,12):(50.0/60.0,43.0/60.0),(6,6):(1.0,30.0/36.0),(6,7):(36.0/42.0,30.0/42.0),(6,8):(40.0/48.0,34.0/48.0),(6,9):(45.0/54.0,39.0/54.0),(6,10):(48.0/60.0,40.0/60.0),(6,11):(54.0/66.0,43.0/66.0),(6,12):(60.0/72.0,48.0/72.0),(7,7):(42.0/49.0,42.0/49.0),(7,8):(48.0/56.0,46.0/56.0),(7,9):(49.0/63.0,42.0/63.0),(7,10):(53.0/70.0,46.0/70.0),(7,11):(59.0/77.0,48.0/77.0),(7,12):(60.0/84.0,53.0/84.0),(8,8):(56.0/64.0,48.0/64.0),(8,9):(55.0/72.0,46.0/72.0),(8,10):(60.0/80.0,48.0/80.0),(8,11):(64.0/88.0,53.0/88.0),(8,12):(68.0/96.0,60.0/96.0),(9,9):(63.0/81.0,54.0/81.0),(9,10):(70.0/90.0,53.0/90.0),(9,11):(70.0/99.0,59.0/99.0),(9,12):(75.0/108.0,63.0/108.0),(10,10):(80.0/100.0,70.0/100.0),(10,11):(77.0/110.0,60.0/110.0),(10,12):(80.0/120.0,66.0/120.0),(11,11):(88.0/121.0,77.0/121.0),(11,12):(86.0/132.0,72.0/132.0),(12,12):(84.0/144.0,96.0/144.0)}

    #coefficient for calculating critical value when n1>12
    c_alpha = {0.01:1.63,0.05:1.73}

    if n1<n2:
        aux = n2
        n2 = n1
        n1 = aux

    if n2<1: n2 = 1
    
    if n1>12: D = c_alpha[alpha]*np.sqrt((n1+n2)/(n1*n2))
    elif (n2,n1) not in C: D = 10.0
    elif alpha==0.01: D = C[(n2,n1)][0]
    else: D = C[(n2,n1)][1]
    return D
        
def chi2Yates(O,E):
    #performs the chi2 goodness-of-fit test with Yates' correction
    #O = non-normalized frequency distribution of observed events
    #E = non-normalized frequency distribution of expected events
    #O and E have to be of same length

    dof = len(O)-1
    X2 = sum((abs(O[i]-E[i])-0.5)**2/E[i] for i in range(0,len(O)))
    return 1-stats.chi2.cdf(X2,df=dof), X2



