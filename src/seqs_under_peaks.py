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
import string

from operator import itemgetter
from glob import glob
from Bio import SeqIO
from random import randint

def seqs_under_peaks():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    
    #MANDATORY PARAMETERS
    parser.add_argument("genome",help="Full path to the reference genome file in fasta-format.",type=str)
    parser.add_argument("peaks",help="Full path to the file containing peaks in igv-format.",type=str)
    parser.add_argument("outname",help="Full path for the output file.",type=str)
    
    #OPTIONAL PARAMETERS
    parser.add_argument("-N","--numpeaks",help="How many of the top peaks are analyzed (default=100).",type=int,default=100)
    parser.add_argument("-w","--peakwidth",help="Width of the sequence around the peak summit that is saved (w+1+w), default=25.",type=int,default=25)
    parser.add_argument("-b","--background",help="1, if a same sized random set of sequences is output, 0 (=default) otherwise.",type=int,default=0,choices=[0,1])
    parser.add_argument("-v","--varlen",help="If yes, region sizes are taken from the input bed-file and -w is overridden, if no (=default), argument -w is used to determine region size.", type=str,choices=['no','yes'],default='no')

    args = parser.parse_args()

    N = args.numpeaks
    w = args.peakwidth

    peaks = []
    #Structure of the list of lists:
    #[[chr,peak location,peak width],...]

    #reading in the peaks
    with open(args.peaks,'r') as peakfile:
        r = csv.reader(peakfile,delimiter='\t')
        for row in r:
            if row[0].count('chromosome')>0: continue
            chrom = row[0]
            start = int(float(row[1]))
            end = int(float(row[2]))
            if args.varlen=='no':
                if end-start>1: loc = start+(end-start)/2
                else: loc = start
                peaks.append([chrom,loc,end-start])
            else: peaks.append([chrom,start,end])

    #selecting the required number of top peaks
    peaks = peaks[:N]

    sequences = []
    names = []
    if args.background==1: background = []

    #sorting found peaks according to chromosome
    top_peaks = {}
    #Structure of the dictionary:
    #key = chromosome
    #value = [[loc1,len1],[loc2,len2],...]
    for i in range(0,len(peaks)):
        if peaks[i][0] not in top_peaks: top_peaks[peaks[i][0]] = [[peaks[i][1],peaks[i][2]]]
        else: top_peaks[peaks[i][0]].append([peaks[i][1],peaks[i][2]])

    #retrieving the underlying sequences from reference genome
    handle = open(args.genome,'rU')
    for record in SeqIO.parse(handle,"fasta"):
        chrom = record.id
        if chrom not in top_peaks: continue
        seq = record.seq
        #saving the sequence for each peak in this chromosome
        for peak in top_peaks[chrom]:
            if args.varlen=='no': sequences.append(seq[peak[0]-w:peak[0]+w])
            else: sequences.append(seq[peak[0]:peak[1]])
            names.append(chrom+"_"+str(peak[0]))
            if args.background==1:
                if args.varlen=='no':
                    r = randint(0,len(seq)-2*w)
                    background.append(seq[r:r+2*w])
                else:
                    r = randint(0,len(seq)-(peak[1]-peak[0]))
                    background.append(seq[r:r+(peak[1]-peak[0])])
    handle.close()

    f = open(args.outname,'w')
    for s in range(0,len(sequences)):
        f.write(">"+str(s)+"_"+str(names[s])+"\n"+str(sequences[s]).upper()+"\n")
    f.close()

    if args.background==1:
        f = open(args.outname+".bg",'w')
        for s in range(0,len(background)):
            f.write(">"+str(s)+"_"+str(names[s])+"\n"+str(background[s]).upper()+"\n")
        f.close()
#end


seqs_under_peaks()
