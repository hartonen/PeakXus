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

import csv
import argparse
import pysam

import numpy as np
from scipy import special
from scipy.stats import binom_test
import multiprocessing as mp

from operator import itemgetter
from time import time

def testASB():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    #MANDATORY PARAMETERS
    parser.add_argument("peakfile",help="Full path to an igv-file containing the called peaks.",type=str,nargs=1)
    parser.add_argument("readfile",help="Full path to the bam-file used as input for peak calling.",type=str,nargs=1)
    parser.add_argument("umifile",help="Full path to the file containing the UMI-labels for the input bam-file each label on its own row.",type=str,nargs=1)
    parser.add_argument("snpfile",help="Full path to a vcf-file containing the SNPs. Either one file containing all SNPs or one file for each chromosome.",type=str,nargs='+')

    #OPTIONAL PARAMETERS
    parser.add_argument("-m",help="Maximum distance for a SNP from a peak edge for it to be considered to overlap with a peak (default=5 bps).",type=int,default=5)
    parser.add_argument("-n",help="Number of parallel processes used (default=1).",type=int,default=1)
    parser.add_argument("--nofilter",help="If 1, locations with 0 reads on either of the alleles are included to the output, if 0 these locations are filtered out (default=0).",type=int,choices=[0,1],default=0)
    parser.add_argument("-l",help="Length of the UMI-labels used (default=5 bps).",type=int,default=5)

    args = parser.parse_args()


    #################
    #READING IN DATA#
    #################

    peaks = []
    chroms = set()
    with open(args.peakfile[0],'rb') as csvfile:
        r = csv.reader(csvfile,delimiter='\t')
        for row in r:
            if row[0].count('chromosome')>0: continue
            peaks.append(row)
            if row[0] not in chroms: chroms.add(row[0])

    #creating a dictionary of peaks
    peakdict = {}
    #Structure of the dictionary:
    #key = chromosome
    #value = [start,end,id]
    for p in peaks:
        chrom = p[0]
        if chrom not in peakdict: peakdict[chrom] = [[int(float(p[1])),int(float(p[2])),p[3],p[4]]]
        else: peakdict[chrom].append([int(float(p[1])),int(float(p[2])),p[3],p[4]])

    #making sure that the peaks are sorted by the start coordinate
    for chrom in peakdict: peakdict[chrom].sort(key=lambda x: x[1])
    umis = set()
    with open(args.umifile[0],'rb') as csvfile:
        r = csv.reader(csvfile,delimiter='\t')
        for row in r: umis.add(row[0])

    tic = time()

    overlapsites = {}
    #Structure of the dictionary of sites where a SNP overlaps with a peak
    #key = chrom
    #value = [[pos1,snp_id1,peak_id1,REF_snp,ALT_snp],...]
    #analyzing SNP overlap with peaks

    if args.n<2:
        for snpfile in args.snpfile:
            with open(snpfile,'rb') as csvfile:
                r = csv.reader(csvfile,delimiter='\t')
                for row in r:
                    chrom = row[0]
                    if chrom not in chroms: continue
                    elif chrom not in overlapsites: overlapsites[chrom] = []
                    pos = int(row[1])
                    snp_id = row[2]
                    REF_snp = row[3]
                    ALT_snp = row[4]

                    #testing proximity of peaks
                    for p in peakdict[chrom]:
                        start = int(float(p[0]))
                        end = int(float(p[1]))
                        if (start-args.m)<=pos and (end+args.m)>=pos:
                            #peak and SNP do overlap
                            overlapsites[chrom].append([pos,snp_id,p[2],REF_snp,ALT_snp])
                            break
                        if (end+args.m)>pos: break

    else:

        #getting the list of chromosomes for which we have the SNPs
        chromlist = []
        for snpfile in args.snpfile:
            with open(snpfile,'rb') as csvfile:
                r = csv.reader(csvfile,delimiter='\t')
                for row in r:
                    chrom = row[0]
                    chromlist.append(chrom)
                    break

        pool = mp.Pool(processes=args.n)
        res = [pool.apply_async(testOverlap,args=(args.snpfile[i],peakdict[chromlist[i]],chromlist[i],args.m)) for i in range(0,len(chromlist))]
        res = [p.get() for p in res]
        for r in res: overlapsites[r[0]] = r[1]

        pool.close()
        pool.terminate()
        pool.join()

    #now all overlapping sites have been selected
    #next acquiring the sequences under the peaks

    toc = time()
    #print "Analyzing overlap took "+str(toc-tic)+" s"
    results = [] # results[i] = [chrom,pos,snp_id,peak_id,ref allele,alt allele,p-value] 
    #reading in sequences
    samfile = pysam.AlignmentFile(args.readfile[0],'rb')
    #assuming read length is constant
    
    skip = True
    for chrom in overlapsites.keys():
        for site in overlapsites[chrom]:
            pos = site[0]
            REF_snp = site[3]
            ALT_snp = site[4]
            used_umis = [] #used_umis[i] = (strand,fiveprime,umi)
            letters = {REF_snp:0.0,ALT_snp:0.0} #key = letter, value = count
            letters_allreads = {REF_snp:0.0,ALT_snp:0.0} #same statistics for all reads i.e. no umis
            for pileupcolumn in samfile.pileup(chrom,pos,pos+1):
                #now looping through all reads that overlap with position pos
                if pos-1!=pileupcolumn.pos: continue
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del: continue
                    index = [i for i,v in enumerate(pileupread.alignment.tags) if v[0].count('BC')>0]
                    umi = pileupread.alignment.tags[index[0]][1][:args.l]
                    umi = umi.upper()

                    #each umi is counted only once per each strand
                    strand = '+'
                    if pileupread.alignment.is_reverse: strand = '-'

                    #determine five prime end
                    if strand == '+': fiveprime = pileupread.alignment.aend-pileupread.alignment.rlen+1
                    else: fiveprime = pileupread.alignment.aend

                    l = pileupread.alignment.query_sequence[pileupread.query_position]
                    if l!=REF_snp and l!=ALT_snp: continue
                    #this stores read counts without umis for reference
                    letters_allreads[l] += 1.0
                    
                    #filtering artefacts with umis
                    if ((umi in umis) and ((strand,fiveprime,umi) not in used_umis)):
                        letters[l] += 1.0
                        used_umis.append((strand,fiveprime,umi))

            #now all reads overlapping with pos are analyzed
            #deleting SNPs with 0 hits to either allele
            if args.nofilter!=1:
                if letters[REF_snp]<1 or letters[ALT_snp]<1: continue

            results.append([chrom]+site+[letters[REF_snp],letters[ALT_snp],letters_allreads[REF_snp],letters_allreads[ALT_snp]])

    #printing results
    print "#chrom\tlocation\tSNP_id\tpeak_location\tREF_allele\tALT_allele\tN_REF\tN_ALT\tallreads_REF\tallreads_ALT"
    for res in results:
        for r in res: print str(r)+"\t",
        print ""


def testOverlap(snpfile,peaks,chrom,m):
    #snpfile = path to the file containing the SNPs
    #peaks = dictionary containing the called peaks
    #chrom = the chromosome to be analyzed
    #m = maximum distance from peak edge to snp
    overlapsites = [chrom,[]]

    with open(snpfile,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter='\t')
        last_peak_index = 0 #index where to start searching for an overlapping peak
        for row in r:
            chrom = row[0]
            pos = int(row[1]) #SNP position
            if pos<(peaks[0][0]-m): continue
            snp_id = row[2]
            REF_snp = row[3]
            ALT_snp = row[4]
                
            #testing proximity of peaks
            for i in range (last_peak_index,len(peaks)):
                start = peaks[i][0] #peak start position
                end = peaks[i][1] #peak end position

                if (start-m)<=pos and (end+m)>=pos:
                    #peak and SNP do overlap
                    overlapsites[1].append([pos,snp_id,peaks[i][2],REF_snp,ALT_snp])
                    last_peak_index = i
                if (end+m)>pos: break

    return overlapsites

testASB()



