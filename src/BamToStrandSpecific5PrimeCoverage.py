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
import string
import csv
import pysam
import multiprocessing as mp
import numpy as np

import matplotlib
#using the Agg-backend for plotting so that X-server is not required
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from plotting import plot1d

def BamToStrandSpecific5PrimeCoverage():

    matplotlib.rcParams.update({'font.size': 26})
    matplotlib.rcParams.update({'font.weight': 'bold'})

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    
    #MANDATORY PARAMETERS
    parser.add_argument("infile",help="Full path to the input bam-file.",type=str)
    parser.add_argument("outdir",help="Full path to the output directory.",type=str)
    parser.add_argument("chromsizes",help="Full path to a file with chromosome names and sizes.",type=str)

    parser.add_argument("-s","--strand",help="Which strand reads are output.",type=str,choices=['+','-'],default='+')
    parser.add_argument("-u","--umifile",help="Full path to file containing used UMIs.",type=str,default=None)
    parser.add_argument("-n","--nproc",help="Number of parallel processes.",type=int,default=1)
    parser.add_argument("-l","--umilen",help="Length of the UMI (starting from the 5'-end).",type=int,default=5)

    args = parser.parse_args()

    #Reading in the reference UMI-labels from file
    if args.umifile!=None:
        ref_umis = set()
        with open(args.umifile,'rb') as csvfile:
            r = csv.reader(csvfile,delimiter='\t')
            for line in r: ref_umis.add(line[0].upper())

    #Reading in chromosome sizes:
    chroms = []
    with open(args.chromsizes,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter="\t")
        for row in r: chroms.append([row[0],int(row[1])])
    histo = [np.zeros(len(ref_umis)),np.zeros(len(ref_umis)+1)]

    #Calculating the 5'-end coverage for each chromosome
    if args.nproc<2:
        for c in chroms:
            aux_histo = calcCoverageFast(c[0],args.strand,args.infile,args.outdir+c[0]+"_"+args.strand+".igv",ref_umis,args.umilen)
            histo[0] += aux_histo[0]
            histo[1] = aux_histo[1]
    else:
        pool = mp.Pool(processes=args.nproc)
        res = [pool.apply_async(calcCoverageFast,args=(c[0],args.strand,args.infile,args.outdir+c[0]+"_"+args.strand+".igv",ref_umis,args.umilen)) for c in chroms]
        res = [p.get() for p in res]
        for aux_histo in res:
            histo[0] += aux_histo[0]
            histo[1] = aux_histo[1]

        pool.close()
        pool.terminate()
        pool.join()

    #Plotting the histogram of UMI-counts per position
    
    index = np.where(histo[0][1:]>0)[0][-1]

    if index<len(ref_umis)-10: index += 10
    if args.strand=='+': colors = ['r']
    else: colors = ['b']
    legend = [args.strand+" strand"]
    plot1d(histo[1][0:index],[histo[0][0:index]],args.outdir+"/UMI_histo_"+args.strand+".png",colors=colors,xlabel="Number of UMIs per position",title="UMI-count histogram",yscale='log',Nyticks=None,Nxticks=6,legend=legend)

    print "total number of UMIs="+str(np.sum(histo[0][0:index]))

#end

def calcCoverageFast(chrom,strand,infile,outfile,umis,umilen):

    samfile = pysam.Samfile(infile,'rb')
    iter = samfile.fetch()

    counts = {}
    #key = position
    #value = set of unique umis

    for read in iter:
        if read.is_unmapped: continue

        readchrom = samfile.getrname(read.rname)
        if readchrom!=chrom: continue
        readstrand = '+'
        if read.is_reverse: readstrand = '-'
        if readstrand!=strand: continue

        if readstrand == '+': fiveprime = read.reference_start+1
        else: fiveprime = read.reference_end
        
        
        #reading in the UMI
        if umis!=None:
            index = [i for i,v in enumerate(read.tags) if v[0].count('BC')>0]
            umi = read.tags[index[0]][1]
            umi = umi.upper()
            umi = umi[0:umilen]
        
            if umi not in umis: continue

        if fiveprime not in counts:
            if umis!=None: counts[fiveprime] = set([umi])
            else: counts[fiveprime] = [1]
        else:
            if umis!=None: counts[fiveprime].add(umi)
            else: counts[fiveprime].append(1)


    #strand specific coverage calculated, saving to file
    with open(outfile,'wb') as csvfile:
        w = csv.writer(csvfile,delimiter="\t")
        w.writerow(["chromosome","start","end","id","signal"])
        for pos in sorted(counts.keys()): w.writerow([chrom,pos-1,pos-1,pos,len(counts[pos])])

    #calculating the UMI-histogram
    count_list = []
    for key in counts: count_list.append(len(counts[key]))
    histo = np.histogram(np.array(count_list),bins=range(0,len(umis)+1))
    
    return histo

BamToStrandSpecific5PrimeCoverage()
