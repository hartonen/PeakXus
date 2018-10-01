#!/usr/bin/env python
# Copyright (c) Tuomo Hartonen, 2015-2016
#
# THIS PROGRAM COMES WITH ABSOLUTELY NO WARRANTY!
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 2 as 
# published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

import argparse
import csv
import pysam

import matplotlib
#Using the Agg backend for plotting so that no X-server is required
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from operator import itemgetter
from os import system

from ReadContainer6 import ReadContainer
import clusterRegions6 as cr
#Functions for plotting fancy results from peak calling

from plotting import plotHeatMap
from plotting import plot1d
from plotting import plotHistogram

def fancyResults():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    #MANDATORY PARAMETERS
    parser.add_argument("peakfile",help="Full path to an igv-file containing called peaks (peaks must be presorted in the desired order!).",type=str)
    parser.add_argument("readfile",help="Full path to the bam-file used as input for peak calling.",type=str)
    parser.add_argument("outdir",help="Full path to output directory.",type=str)
    
    #OPTIONAL PARAMETERS
    #if these are not given, defaults are used
    parser.add_argument("--motiffile",help="Full path to a gff-file containing the high-affinity recognition sequence (HARS) locations, if not provided, part of the results are not output.",type=str,default=None)
    parser.add_argument("-N","--numpeaks",help="Number of top peaks analyzed.",type=int,default=1000)
    parser.add_argument("-c","--comparison",help="If 0 (=default), calculating distance between middle of the peak and middle of the motif, if 1, calculating distance between right end of the motif and right end of the peak.",type=int,choices=[0,1],default=0)
    parser.add_argument("-m","--maxdist",help="Maximum distance between peak and motif (default=10).",type=int,default=10)
    parser.add_argument("-w","--peakwidth",help="Width of the region plotted around the peaks in heatmaps (default=200bps).",type=int,default=200)
    parser.add_argument("-e","--expname",help="Experiment name for plot titles.",type=str,default="PeakXus")
    parser.add_argument("-n","--nproc",help="Number of parallel processes used in clustering (default=1).",type=int,default=1)
    parser.add_argument("-k",help="Number of clusters (default=4).",type=int,default=4)
    parser.add_argument("-p","--npass",help="Number of times clustering is run (default=1).",type=int,default=1)
    parser.add_argument("-t","--test",help="Statistical test used, 0=G-test, 2=KS-test (default=0).",type=int,choices=[0,2],default=0)
    parser.add_argument("-H","--skipheatmap",help="1, if plotting heatmaps is omitted,0 (=default) otherwise.",type=int,choices=[0,1],default=0)
    parser.add_argument("-S","--smoothing",help="Window size used in calculating read distributions (default=1).",type=int,default=1)
    parser.add_argument("-ms","--motif_start",help="Start position of the motif, if middle of the motif is at 0",type=int,default=-8)
    parser.add_argument("-ml","--motif_len",help="Motif length.",type=int,default=17)
    parser.add_argument("--truereads",help="1, if only reads pointing towards the peak summit are plotted, 0 if all reads are plotted (default=0).",type=int,choices=[0,1],default=0)
    parser.add_argument("--pseudo",help="Pseudocount value for g-test (default=1).",type=float,default=1.0)
    parser.add_argument("--ctest",help="g if g-test is used as distance for clustering, euc if euclidian distance.",type=str,choices=['g','euc'],default='g')
    parser.add_argument("-u","--UMIs",help="Full path to the file containing true UMIs (one barcode per row).",type=str,default=None)
    parser.add_argument("-l","--umilen",help="Length of the used UMI-labels (default=5).",type=int,default=5)
    parser.add_argument("--skipclustering",help="If 1, clustering is skipped, otherwise 0 (default).",type=int,choices=[0,1],default=0)
    parser.add_argument("-s",help="Peak file column used as the peak score, numbering starting from 0 (default=6).",type=int,default=6)


    args = parser.parse_args()
    N = args.numpeaks
    c = args.comparison
    m = args.maxdist
    w = args.peakwidth
    h = args.skipheatmap
    S = args.smoothing

    matplotlib.rcParams.update({'font.size': 26})
    matplotlib.rcParams.update({'font.weight': 'bold'})

    #################
    #READING IN DATA#
    #################

    #reading in peaks
    peaks = []
    chroms = set()
    count = 0
    with open(args.peakfile,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter='\t')
        for row in r:
            if row[0].count('chromosome')>0: continue
            peaks.append(row)
            if row[0] not in chroms: chroms.add(row[0])
            count += 1
            if count>=3*N: break
    topPeaks = peaks[:N]

    #reading in motifs
    motif_size = 0
    motifs = {}
    #dictionary structure:
    #key = chromosome name
    #value = [[list of + strand motifs],[list of - strand motifs]]
    if args.motiffile!=None:
        with open(args.motiffile,'rb') as csvfile:
            r = csv.reader(csvfile,delimiter='\t')
            for row in r:
                chrom = row[0]
                if chrom not in chroms: continue
                start = int(float(row[3]))
                end = int(float(row[4]))
                strand = row[5]
                if c==0: loc = (end-start)/2+start #middle position of HARS
                elif c==1: loc = end #right edge of HARS
                elif c==2: loc = start #left edge of HARS

                motif_size = end-start
                if chrom not in motifs:
                    if strand=='+': motifs[chrom] = [set([loc]),set()]
                    else: motifs[chrom] = [set(),set([loc])]
                else:
                    if strand=='+': motifs[chrom][0].add(loc)
                    else: motifs[chrom][1].add(loc)

    #reading in sequences
    samfile = pysam.Samfile(args.readfile,'rb')
    if args.UMIs!=None:
        #reading in UMIs
        umis = set()
        with open(args.UMIs,'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            for row in r: umis.add(row[0])

    #########################
    #HISTOGRAM OF PEAK SIZES#
    #########################

    #peak size (total UMI-count pointing towards the peak summit) is column 4, peak score is column number 6 (starting from 0)
    sizes = []
    scores = []
    for p in peaks:
        sizes.append(int(float(p[4])))
        scores.append(float(p[args.s]))

    #plotting sizes
    histo = np.histogram(np.array(sizes),bins=range(0,max(sizes)+1))
    plot1d(histo[1][:-1],[histo[0]],args.outdir+"peak_size_histo.png",xlabel="UMIs pointing towards peak summit",title=args.expname)#,yscale='log',Nyticks=None,xscale='log')

    #plotting scores
    histo = np.histogram(np.array(scores),bins=1000)

    plot1d(histo[1][:-1],[histo[0]],args.outdir+"peak_score_histo.png",xlabel="Peak score",title=args.expname,xscale='log')#,yscale='log',Nyticks=None,xscale='log')

    ########################
    #HEATMAP OF TOP N PEAKS#
    ########################
    
    #creating matrix of the read densities of top peaks
    E = np.zeros((2*N,w+1))
    ind = 0
    for p in topPeaks:
        chrom = p[0]
        start = int(float(p[1]))
        end = int(float(p[2]))
        summit = (end-start)/2+start
        if args.UMIs!=None:
            used_umis = []
            #List format:
            #(fiveprime end location,strand,UMI-label)
            #This is used to assure that each UMI is counted only once
        if summit-w/2-100<0: continue
        for read in samfile.fetch(chrom,summit-w/2-100,summit+w/2+1+100):
            if read.is_unmapped: continue
            strand = '+'
            if read.is_reverse: strand = '-'
            
            #determine five prime end
            if strand == '+': fiveprime = read.reference_start+1#read.aend-read.rlen+1
            else: fiveprime = read.reference_end#read.aend

            if args.UMIs!=None:
                index = [i for i,v in enumerate(read.tags) if v[0].count('BC')>0]
                umi = read.tags[index[0]][1][:args.umilen]
                umi = umi.upper()
                if umi not in umis: continue
                if (fiveprime,strand,umi) not in used_umis: used_umis.append((fiveprime,strand,umi))
                else: continue

            if fiveprime in range(summit-w/2,summit+w/2+1):
                if strand=='+':
                    if args.truereads==1:
                        #only reads starting from the left side of the summit are considered
                        if fiveprime<summit: E[ind,fiveprime-(summit-w/2)] += 1.0
                    else: E[ind,fiveprime-(summit-w/2)] += 1.0
                else:
                    if args.truereads==1:
                        #only reads starting from the right side of the summit are considered
                        if fiveprime>summit: E[ind+1,fiveprime-(summit-w/2)] -= 1.0
                    else: E[ind+1,fiveprime-(summit-w/2)] -= 1.0
        
        ind += 2
    #plotting

    #Smoothing the read 5'-end count distribution using np.convolve if S>1
    if S>1:
        for i in range(0,N): E[i,:] = np.convolve(E[i,:],np.ones(S),'same')/float(S)

    ##########################
    #CLUSTERING THE TOP PEAKS#
    ##########################

    if args.UMIs!=None: ylabel = "Average UMI 5'-end count"
    else: ylabel = "Average read 5'-end count"
    x = np.array([i for i in range(-w/2,w/2+1)])
    y = np.array([i for i in range(0,int(np.shape(E)[0]))])

    if args.skipclustering<1:

        #The read 5'-end count distributions areound the top peaks are clustered using the k-medoids algorithm from Biopython


        E_clustered,centroids,clusters = cr.clusterRegions(E[range(0,len(E),2),:],np.abs(E[range(1,len(E),2),:]),args.k,args.npass,args.pseudo,args.nproc,args.ctest)

        #Plotting the cluster centroids
        ind = 0
        for c in centroids:
            ind += 1
            plot1d(x,[c[1],c[2]],args.outdir+"Centroid"+str(ind)+"_top_"+str(N)+".png",xlabel="Distance from peak summit",ylabel=ylabel,colors=['r','b'],title="Centroid="+str(ind)+", "+str(len(np.where(clusters==c[0])[0]))+" peaks")

        #Plotting a heatmap of each cluster
        ind = 0
        for c in set(clusters):
            ind += 1
            E_aux = E_clustered[np.where(E_clustered[:,-1]==c)]
            plotHeatMap(E_aux[:,:-1],x,np.array([i for i in range(0,int(np.shape(E_aux[:,:-1])[0]))]),args.outdir+"Heatmap_top_"+str(N)+"cluster"+str(ind)+".png",xlabel="Distance from peak summit [bps]",ylabel="Cluster"+str(ind),xticks=[-w/2,0,w/2])

    plotHeatMap(E,x,y,args.outdir+"Heatmap_top_"+str(N)+".png",xlabel="Distance from peak summit [bps]",xticks=[-w/2,0,w/2],title=args.expname)

    #############################################
    #AVERAGE READ DENSITY PROFILE OF TOP N PEAKS#
    #############################################

    plot1d(x,[np.mean(np.abs(E[range(0,len(E),2)]),axis=0),np.mean(np.abs(E[range(1,len(E),2)]),axis=0)],args.outdir+"Avg_top_"+str(N)+".png",colors=['r','b'],xlabel="Distance from peak summit",ylabel=ylabel,xticks=[-w/2,0,w/2],title=args.expname,errorbar=[np.std(E[range(0,len(E),2)],axis=0),np.std(E[range(0,len(E),1)],axis=0)])  

    if args.motiffile==None: return

    ###########################################
    #HEATMAP OF TOP N PEAKS WITH A MOTIF MATCH#
    ###########################################

    #selecting top N peaks overlapping with an HARS
    topPeaks = []
    dists_from_summit = []
    strands = []
    hitcount = 0
    allowed_dists = [i for i in range(0,m+1)]
    matching_motifs = []

    #A new set of top peaks is created as not all the top 1000 peaks overlap with an HARS
    for peak in peaks:
        if hitcount>=N: break
        chrom = peak[0]
        if chrom not in motifs: continue
        peak_start = int(float(peak[1]))
        peak_end = int(float(peak[2]))
        if args.comparison==0: peak_loc = (peak_end-peak_start)/2+peak_start
        elif args.comparison==1: peak_loc = peak_end
        elif args.comparison==2: peak_loc = peak_start
        

        for d in allowed_dists:
            if (peak_loc+d in motifs[chrom][0]) or (peak_loc-d in motifs[chrom][0]) or (peak_loc+d in motifs[chrom][1]) or (peak_loc-d in motifs[chrom][1]):
                #This means we have a hit, a peak and an HARS overlap
                hitcount += 1
                topPeaks.append(peak)
                #removing the used HARS
                if peak_loc+d in motifs[chrom][0]:
                    dists_from_summit.append(d)
                    motifs[chrom][0].remove(peak_loc+d)
                    strands.append('+')
                    matching_motifs.append(peak_loc+d)
                elif peak_loc+d in motifs[chrom][1]:
                    dists_from_summit.append(d)
                    motifs[chrom][1].remove(peak_loc+d)
                    strands.append('-')
                    matching_motifs.append(peak_loc+d)
                elif peak_loc-d in motifs[chrom][0]:
                    dists_from_summit.append(-d)
                    motifs[chrom][0].remove(peak_loc-d)
                    strands.append('+')
                    matching_motifs.append(peak_loc-d)
                else:
                    dists_from_summit.append(-d)
                    motifs[chrom][1].remove(peak_loc-d)
                    strands.append('-')
                    matching_motifs.append(peak_loc-d)
                break

    #creating matrix of the read densities of top peaks
    E = np.zeros((2*N,w+1))
    ind = 0

    for p in range(0,len(topPeaks)):
        chrom = topPeaks[p][0]
        start = int(float(topPeaks[p][1]))
        end = int(float(topPeaks[p][2]))
        summit = (end-start)/2+start
        loc = summit+dists_from_summit[p]

        if args.UMIs!=None: used_umis = []

        for read in samfile.fetch(chrom,loc-w/2-100,loc+2/2+1+100):
            if read.is_unmapped: continue
            strand = '+'
            if read.is_reverse: strand = '-'
            
            #determine five prime end
            if strand == '+': fiveprime = read.reference_start+1
            else: fiveprime = read.reference_end

            if args.UMIs!=None:
                index = [i for i,v in enumerate(read.tags) if v[0].count('BC')>0]
                umi = read.tags[index[0]][1][:args.umilen]
                umi = umi.upper()
                if umi not in umis: continue
                if (fiveprime,strand,umi) not in used_umis: used_umis.append((fiveprime,strand,umi))
                else: continue

            if fiveprime in range(loc-w/2,loc+w/2+1):
                if strand=='+':
                    if args.truereads==1:
                        #true reads are on the left side of the motif center
                        if fiveprime<loc: E[ind,fiveprime-(loc-w/2)] += 1.0
                    else: E[ind,fiveprime-(loc-w/2)] += 1.0
                else:
                    if args.truereads==1:
                        #true reads are on the right side of the motif center
                        if fiveprime>loc: E[ind+1,fiveprime-(loc-w/2)] -= 1.0
                    else: E[ind+1,fiveprime-(loc-w/2)] -= 1.0
        ind += 2
    #plotting

    #Smoothing using np.convolve
    if S>1:
        for i in range(0,N): E[i,:] = np.convolve(E[i,:],np.ones(S),'same')/float(S)

    #calculating the average peak width
    widths = []
    for i in range(0,2*N,2):
        d = 0.0
        left = E[i,0:w/2][::-1]
        right = E[i+1,w/2+1:]
        if len(np.where(left>0)[0])>0: d += np.where(left>0)[0][0]
        else: continue
        if len(np.where(right<0)[0])>0: d += np.where(right<0)[0][0]
        else: continue
        widths.append(d)

    x = np.array([i for i in range(-w/2,w/2+1)])
    y = np.array([i for i in range(0,np.shape(E)[0])])
    plotHeatMap(E,x,y,args.outdir+"Heatmap_motif_match_top_"+str(N)+".png",xlabel="Distance from motif center [bps]",xticks=[-w/2,0,w/2],yticks=[])

    #############################################
    #CLUSTERING THE TOP PEAKS WITH A MOTIF MATCH#
    #############################################

    if args.skipclustering<0:
        E_clustered,centroids,clusters = cr.clusterRegions(E[range(0,len(E),2),:],np.abs(E[range(1,len(E),2),:]),args.k,args.npass,args.pseudo,args.nproc,args.ctest)

        #Plotting the cluster centroids
        x = np.array([i for i in range(-w/2,w/2+1)])
        y = np.array([i for i in range(0,int(np.shape(E)[0]))])
        ind = 0
        for c in centroids:
            ind += 1
            plot1d(x,[c[1],c[2]],args.outdir+"Centroid"+str(ind)+"_top_"+str(N)+"motif_match.png",xlabel="Distance from peak summit",ylabel=ylabel,colors=['r','b'],title="Centroid="+str(ind)+", "+str(len(np.where(clusters==c[0])[0]))+" peaks")

        #Plotting a heatmap of each cluster
        ind = 0
        for c in set(clusters):
            ind += 1
            E_aux = E_clustered[np.where(E_clustered[:,-1]==c)]
            plotHeatMap(E_aux[:,:-1],x,np.array([i for i in range(0,int(np.shape(E_aux[:,:-1])[0]))]),args.outdir+"Heatmap_top_"+str(N)+"cluster"+str(ind)+"motif_match.png",xlabel="Distance from peak summit [bps]",ylabel="Cluster"+str(ind),xticks=[-w/2,0,w/2],title=args.expname)
    

    ################################################################
    #AVERAGE READ DENSITY PROFILE OF TOP N PEAKS WITH A MOTIF MATCH#
    ################################################################

    plot1d(x,[np.mean(np.abs(E[range(0,len(E),2)]),axis=0),np.mean(np.abs(E[range(1,len(E),2)]),axis=0)],args.outdir+"Avg_motif_match_top_"+str(N)+".png",title=args.expname,colors=['r','b'],xlabel="Distance from motif center",ylabel=ylabel,ylim=0,errorbar=[np.std(E[range(0,len(E),2)],axis=0),np.std(E[range(0,len(E),1)],axis=0)])   

    ##########################################################################
    #AVERAGE READ DENSITY PROFILE OF TOP N PEAKS WITH AN ORIENTED MOTIF MATCH#
    ##########################################################################

    #for motifs on + strand: left=start, right=end
    #for motifs on - strand: left=end, right=start

    #creating matrix of the read densities of top peaks
    #list matching_motifs contains the middle positions of motifs and
    #variable motif_size the width of motif

    E = np.zeros((2*N,w+1))
    ind = 0
    dists_from_start = []
    dists_from_end = []
    for p in range(0,len(topPeaks)):
        chrom = topPeaks[p][0]
        start = int(float(topPeaks[p][1]))
        end = int(float(topPeaks[p][2]))
        summit = (end-start)/2+start
        motif_loc = matching_motifs[p]

        if args.UMIs!=None: used_umis = []

        if strands[p]=='+':
            loc = summit+dists_from_summit[p]
            #start of motif is paired with start of peak
            dists_from_start.append(abs(start-(motif_loc-motif_size/2)))
            dists_from_end.append(abs(end-(motif_loc+motif_size/2)))
        else:
            loc = summit-dists_from_summit[p]
            #start of motif is paired with end of peak
            dists_from_start.append(abs(end-(motif_loc+motif_size/2)))
            dists_from_end.append(abs(start-(motif_loc-motif_size/2)))

        for read in samfile.fetch(chrom,motif_loc-w/2-100,motif_loc+2/2+1+100):
            if read.is_unmapped: continue
            strand = '+'
            if read.is_reverse: strand = '-'
            
            #determine five prime end
            if strand == '+': fiveprime = read.reference_start+1
            else: fiveprime = read.reference_end

            if args.UMIs!=None:
                index = [i for i,v in enumerate(read.tags) if v[0].count('BC')>0]
                umi = read.tags[index[0]][1][:args.umilen]
                umi = umi.upper()
                if umi not in umis: continue
                if (fiveprime,strand,umi) not in used_umis: used_umis.append((fiveprime,strand,umi))
                else: continue

            if fiveprime in range(motif_loc-w/2,motif_loc+w/2+1):
                if strand=='+':
                    if args.truereads==1:
                        #true reads are on the left side of the motif center
                        if fiveprime<motif_loc: E[ind,fiveprime-(motif_loc-w/2)] += 1.0
                    else: E[ind,fiveprime-(motif_loc-w/2)] += 1.0
                else:
                    if args.truereads==1:
                        #true reads are on the right side of the motif center
                        if fiveprime>motif_loc: E[ind+1,fiveprime-(motif_loc-w/2)] -= 1.0
                    else: E[ind+1,fiveprime-(motif_loc-w/2)] -= 1.0

        #now orienting the read density with respect to motif direction
        if strands[p]=='-':
            aux = -1*E[ind+1,:][::-1]
            E[ind+1,:] = -1*E[ind,:][::-1]
            E[ind,:] = aux

        ind += 2
            
    #plotting

    #Smoothing using np.convolve
    if S>1:
        for i in range(0,N): E[i,:] = np.convolve(E[i,:],np.ones(S),'same')/float(S)

    #calculating average peak size and its standard deviation
    std_peak_size = np.std(np.array(widths))
    avg_peak_size = np.mean(np.array(widths))
    med_peak_size = np.median(np.array(widths))

    x = np.array([i for i in range(-w/2,w/2+1)])
    y = np.array([i for i in range(0,np.shape(E)[0])])
    
    plotHeatMap(E,x,y,args.outdir+"Heatmap_motif_match_oriented_top"+str(N)+".png",xlabel="Distance from motif center [bps]",xticks=[-w/2,0,w/2+1],yticks=[],title=args.expname)

    #######################################################
    #CLUSTERING THE TOP PEAKS WITH AN ORIENTED MOTIF MATCH#
    #######################################################

    if args.skipclustering<0:
        E_clustered,centroids,clusters = cr.clusterRegions(E[range(0,len(E),2),:],np.abs(E[range(1,len(E),2),:]),args.k,args.npass,args.pseudo,args.nproc,args.ctest)

        #Plotting the cluster centroids
        x = np.array([i for i in range(-w/2,w/2+1)])
        y = np.array([i for i in range(0,int(np.shape(E)[0]))])
        ind = 0
        for c in centroids:
            ind += 1
            plot1d(x,[c[1],c[2]],args.outdir+"Centroid"+str(ind)+"_top_"+str(N)+"_oriented_motif_match.png",xlabel="Distance from peak summit",ylabel=ylabel,colors=['r','b'],title="Centroid="+str(ind)+", "+str(len(np.where(clusters==c[0])[0]))+" peaks")
        #Plotting a heatmap of each cluster
        ind = 0
        for c in set(clusters):
            ind += 1
            E_aux = E_clustered[np.where(E_clustered[:,-1]==c)]
            plotHeatMap(E_aux[:,:-1],x,np.array([i for i in range(0,int(np.shape(E_aux[:,:-1])[0]))]),args.outdir+"Heatmap_top_"+str(N)+"cluster"+str(ind)+"_oriented_motif_match.png",xlabel="Distance from peak summit (bp's)",ylabel="Cluster"+str(ind),xticks=[-w/2,0,w/2],title=args.expname)  

    ##########################################################################
    #AVERAGE READ DENSITY PROFILE OF TOP N PEAKS WITH AN ORIENTED MOTIF MATCH#
    ##########################################################################

    auxtitle = "Mean peak size="+str(avg_peak_size)+", std="+str(std_peak_size)+", median="+str(med_peak_size)
    plot1d(x,[np.mean(np.abs(E[range(0,len(E),2)]),axis=0),np.mean(np.abs(E[range(1,len(E),2)]),axis=0)],args.outdir+"Avg_motif_match_top_oriented_"+str(N)+".png",colors=['r','b'],xlabel="Distance from motif center",ylabel=ylabel,title=args.expname,m_left=args.motif_start,m_len=args.motif_len,opacity=0.3,errorbar=[np.std(E[range(0,len(E),2)],axis=0),np.std(E[range(0,len(E),1)],axis=0)])   

#end

fancyResults()
