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
import sys
import csv
import numpy as np

from time import time
from os import listdir
from glob import glob

from ReadContainer6 import ReadContainer
from analyzeTransitions6_symmetric import tPoints

def peakC():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    #MANDATORY PARAMETERS
    parser.add_argument("infile",help="Full path to the input bam-file.",type=str)
    parser.add_argument("outdir",help="Full path to the output directory.",type=str)
    parser.add_argument("chromnames",help="Full path to a file containing chromosome names and sizes each chromosome on its own line, name and size separated by tab.",type=str)

    #OPTIONAL PARAMETERS
    #if these are not given, defaults are used

    parser.add_argument("-m","--mindist",help="minimum distance between two peaks (default=30)",type=int,default=30)
    parser.add_argument("-w", "--winsize", help="width of the window around 5'-end (default=5)", type=int, default=5)
    parser.add_argument("-l_l", "--l_limit", help="l_limit+w is the largest allowed width for a peak (default=60)", type=int, default=60)
    parser.add_argument("-u_l", "--u_limit", help="u_limit is the smallest allowed peak width (default=5)", type=int,default=5)
    parser.add_argument("-b", "--b_threshold", help="Threshold p-value for background model. Default is 0.05", type=float, default=0.05)
    parser.add_argument("-r", "--saveregions", help="If 1, read densities for all peak regions are output into a file, if 0, not (default=0).",type=int, choices=[0,1], default=0)
    parser.add_argument("-u", "--readumis", type=str, help="Full path to the file containing the UMI-sequences (1 per row, default=None).",default=None)
    parser.add_argument("-v", "--verbosity", help="verbosity level, 1=verbose on, 0=verbose off", type=int, choices=[0,1], default=1)
    parser.add_argument("-t", "--timing", help="Timing on(1)/off(0)", type=int, choices=[0,1], default=1)
    parser.add_argument("-y","--yates",help="1, if chi2-test with Yates' correction is used, 0 meaning G-test is used. 2=KS-test (default=0)",type=int,choices=[0,1,2],default=0)
    parser.add_argument("-d","--allpairs",help="1, if all positive peaks are considered for each candidate site, 0 if only the highest (default=0)",type=int,choices=[0,1],default=0)
    parser.add_argument("-s","--sliding_window",help="Width of a sliding window used in averaging density of 5' ends of reads. This should be an odd number (default=1).",type=int,default=1)
    parser.add_argument("-n","--nproc",help="Number of parallel processes used in calculating significance of peaks (default=1).",type=int,default=1)
    parser.add_argument("-o","--onlysignal",help="Only report the top 1000 regions with highest read count for each chromosome if o=1 (default=0)",type=int,choices=[0,1],default=0)
    parser.add_argument("-p","--pseudocount",help="Pseudocount added to each position when testing significance of candidate binding sites (default=5.0).",type=float,default=5.0)
    parser.add_argument("-a","--allowoverlap",help="1 if overlapping peaks are allowed (0 otherwise, which is the default)",type=int,choices=[0,1],default=0)
    parser.add_argument("-l","--umilen",help="Length of UMI-label (default=5).",type=int,default=5)
    parser.add_argument("-A","--allreads",help="If 1, all reads are always used in scoring of peak candidates, if 0 (default), UMIs are used if possible.",type=int,choices=[0,1],default=0)

    args = parser.parse_args()

    #################################
    #Testing input, setting defaults#
    #################################


    w = args.winsize
    l_limit = args.l_limit
    u_limit = args.u_limit
    b = args.b_threshold
    mindist = args.mindist
    s = args.sliding_window

    #counter for total number of test conducted
    test_counter = 0
    umi_counter = 0
    read_counter = 0
    peak_counter = 0
    save_umi = args.readumis

    #testing that the input files open
    #log is written to args.outdir[0]/log.txt
    if args.verbosity==1:
        f = open(args.outdir+"log.txt",'a')
        f.write("Testing that input files exist...")
        f.close()

    try:
        #testing that all the individual files open
        f = open(args.infile,'r')
        f.close()
        chroms = []
        with open(args.chromnames,'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            for line in r: chroms.append(line[0])

        #testing that the output directory is writable
        f = open(args.outdir+"all_transition_points.igv",'w')
        f.close()
        if args.verbosity==1:
            f = open(args.outdir+"log.txt",'a')
            f.write("done!")
            f.close()

    except Exception as detail:
        if args.verbosity==1:
            sys.stderr.write("Input file doesn't exist or output directory is not writable!")
            print detail
        sys.exit(1)

    
    for chrom in chroms:

        if args.verbosity==1:
            f = open(args.outdir+"log.txt",'a')
            f.write("\n\nProcessing "+args.infile+"\n")
            f.close()

        #################
        #Reading in data#
        #################

        if args.timing==1:
            times = []
            start = time()

        if args.verbosity==1:
            f = open(args.outdir+"log.txt",'a')
            f.write("Reading in data... ")
            f.close()

        reads = ReadContainer()

        #if readSam returns False there are not reads mapping to chrom and we can proceed to next chromosome
        if not reads.readSam(args.infile,save_umi,args.outdir,args.umilen,args.allreads,chrom): continue
        reads.sortReads()

        read_counter += reads.getReadCount()

        if args.timing==1:
            end = time()
            times.append(end-start)
            if args.verbosity==1:
                f = open(args.outdir+"log.txt",'a')
                f.write("Reading data took "+str(end-start)+" s\n")
                f.close()

        #Analysing one chromosome at a time saves memory

        if args.verbosity==1: print chrom+": ",

        #######################
        #Creating peak regions#
        #######################

        if args.verbosity==1:
            f = open(args.outdir+"log.txt",'a')
            f.write("Creating peak regions... ")
            f.close()
        if args.timing==1: start = time()

        peak_regions = tPoints(w,l_limit,u_limit,mindist,s)
        if args.onlysignal==1: peak_regions.calcTransitions_counts(chrom,reads,args.allpairs,args.nproc)
        else: peak_regions.calcTransitions(chrom,reads,args.allpairs,args.nproc,args.allreads)

        umi_counter += peak_regions.getUMIcount()

        if args.timing==1:
            end = time()
            times.append(end-start)
            if args.verbosity==1:
                f = open(args.outdir+"log.txt",'a')
                f.write("Creating peak regions took "+str(end-start)+" s\n")
                f.close()

        #######################
        #Estimating background#
        #######################

        if b>0:

            if args.verbosity==1:
                f = open(args.outdir+"log.txt",'a')
                f.write("Estimating background... ")
                f.close()
            if args.timing==1: start = time()

            if args.onlysignal==0:
                peak_regions.testBackground(chrom,b,args.yates,args.nproc,args.pseudocount,args.allowoverlap)
                test_counter += peak_regions.getNumTests()

            if args.timing==1:
                end = time()
                times.append(end-start)
                if args.verbosity==1:
                    f = open(args.outdir+"log.txt",'a')
                    f.write("Estimating background took "+str(end-start)+" s\n")
                    f.close()
        tRegions = peak_regions.getTRegions(chrom)

        peak_counter += tRegions.shape[0]
        #####################
        #Saving called peaks#
        #####################

        if args.saveregions==1: regionfile = args.outdir+"region_w="+str(w)+"_l_limit="+str(l_limit)+"_k="+str(k)+"_u_limit="+str(u_limit)+"_"+chrom+".csv"
        else: regionfile = None

        f = open(args.outdir+chrom+"_transition_points.igv",'w')
        if args.onlysignal==1:
            f.write("chromosome\tstart\tend\tid\tsignal\n")
            for t in range(0,len(tRegions)): f.write(chrom+"\t"+str(tRegions[t,0])+"\t"+str(tRegions[t,1])+"\t"+str(tRegions[t,0])+"-"+str(tRegions[t,1])+"\t"+str(tRegions[t,2])+"\n")
        else:
            f.write("chromosome\tstart\tend\tid\tsignal\tp-value\tscore\n")
            for t in range(0,len(tRegions)): f.write(chrom+"\t"+str(tRegions[t,0])+"\t"+str(tRegions[t,1])+"\t"+str(tRegions[t,0])+"-"+str(tRegions[t,1])+"\t"+str(sum(tRegions[t,3:-1]))+"\t"+str(tRegions[t,-2])+"\t"+str(tRegions[t,-1])+"\n")
        f.close()

        if regionfile != None:

            Regions = np.zeros((tRegions.shape[0],2*(l_limit+2*w)+1))
            for r in range(0,len(tRegions)): Regions[r,:] = tRegions[r,2:-1]
            np.savetxt(regionfile,Regions,delimiter=',')

        #analysis done for one chromosome

    #saving a short summary-file of the run
    f = open(args.outdir+"PeakXus_summary.txt",'wb')
    f.write("Input file:")
    f.write("    "+args.infile+"\n")
    f.write("\n")
    f.write("Total number of reads                  = "+str(read_counter)+"\n")
    f.write("Total number of UMIs                   = "+str(umi_counter)+"\n")
    f.write("reads/UMIs                             = "+str(float(read_counter)/float(umi_counter))+"\n")
    f.write("Total number of peaks (no FDR cut-off) = "+str(peak_counter))
    f.close()

    #saving the total number of tests conducted
    f = open(args.outdir+"numtests.txt",'w')
    f.write(str(test_counter))
    f.close()

#end

peakC()
