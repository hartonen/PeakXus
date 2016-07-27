# Copyright (c) Tuomo Hartonen, 2015-2016
#
# THIS PROGRAM COMES WITH ABSOLUTELY NO WARRANTY!
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 2 as 
# published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

from Bio.Cluster import *
import numpy as np
import multiprocessing as mp
from scipy.spatial.distance import euclidean

import analyzeTransitions6_symmetric as at

def clusterRegions(tRegions_plus,tRegions_minus,k,npass,pseudo,nproc,test):
    #clusters the rows in tRegions using k-medoids clustering from biopython
    #tRegions = numpy array, rows are the read densities (+ reads minus - reads)
    #k = number of clusters
    #npass = number of times clustering is run starting from different seeds
    #pseudo = pseudocount for g-test
    #pairwise p-values are used in the distance matrix
    #test = euc if euclidean distance, g if g-test

    N = np.shape(tRegions_plus)[0]
    L = len(tRegions_plus[0,:])
    #calculating the distance matrix
    dm = np.zeros(N*(N-1)/2)
    inds = []
    for i in range(0,N-1): inds += [(i,j) for j in range(i+1,N)]

    if nproc<2:
        dm = [distance([tRegions_plus[i[0]]+tRegions_minus[i[0]]],[tRegions_plus[i[1]]+tRegions_minus[i[1]]],pseudo,test) for i in inds]
    else:
        pool = mp.Pool(processes=nproc)
        res = [pool.apply_async(distance,args=([tRegions_plus[i[0]]+tRegions_minus[i[0]]],[tRegions_plus[i[1]]+tRegions_minus[i[1]]],pseudo,test)) for i in inds]
        #gathering the results from parallel tests
        dm = np.array([p.get() for p in res])

        pool.close()
        pool.terminate()
        pool.join()

    #k-medoids clustering
    clusters,error,nfound = kmedoids(dm,nclusters=k,npass=npass,initialid=None)
    
    tRegions = np.zeros((N,L+1))
    ind = 0
    for i in range(0,len(clusters)):
        #adding the cluster id to the last element of the corresponding histogram
        tRegions[ind,:-1] = tRegions_plus[i,:]-tRegions_minus[i,:]
        tRegions[ind,-1] = clusters[i]
        ind += 1
    #calculating the cluster centroids
    IDs = set(clusters)
    centroids = []

    for ID in IDs:
        centroids.append([int(ID)])
        centroids[-1].append(np.mean(tRegions_plus[np.where(clusters==ID)[0]],axis=0))
        centroids[-1].append(np.mean(tRegions_minus[np.where(clusters==ID)[0]],axis=0))


    return tRegions,centroids,clusters

def distance(p1,p2,pseudo,test):
    #p1 = [+,+,+,+,-,-,-,-]_1
    #p2 = [+,+,+,+,-,-,-,-]_2
    
    p1 = p1[0]
    p2 = p2[0]
    #calculating the distance histograms for peaks pw1 and p2
    N = len(p1)
    p1_pos = []
    p2_pos = []
    for i in range(0,N):
        if i<N/2:
            p1_pos += [i for j in range(0,int(p1[i]))]
            p2_pos += [i for j in range(0,int(p2[i]))]
        else:
            #print p1[i]
            p1_pos += [i-int(N/2) for j in range(0,int(p1[i]))]
            p2_pos += [i-int(N/2) for j in range(0,int(p2[i]))]

    bins = [i for i in range(0,int(N/2)+1)] #maximum distance from summit is N/2!
    contig = np.array([np.histogram(p1_pos,bins=bins)[0],np.histogram(p2_pos,bins=bins)[0]])
    #pi = [[histogram of + strand reads],[histogram of - strand reads]]
    #histogram of - strand reads is reversed here to put the origo at the
    #5' end of the region for both strands
    #p1 = p1[0]+p1[1][::-1]
    #p2 = p2[0]+p2[1][::-1]
    if test == 'g': return at.gtest(contig[0],contig[1],pseudo)[1]
    elif test == 'euc': return euclidean(contig[0],contig[1])

#end

def e_distance(p1,p2):
    #p1 = [+,+,+,+,-,-,-,-]_1
    #p2 = [+,+,+,+,-,-,-,-]_2

    return euclidean(p1,p2)
    
