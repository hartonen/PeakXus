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


def removeOverlappingPeaks():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()

    # MANDATORY PARAMETERS
    parser.add_argument("peakfile",help="Full path to the input igv-file.", type=str)
    parser.add_argument("blacklist",help="Full path to the bed-file containing the blacklisted genomic regions.",type=str)
    parser.add_argument("outfile",help="Full path to the output directory.", type=str)

    args = parser.parse_args()

    #Reading in the blacklisted regions
    blacklist = {}
    #key = chromosome
    #value = [[start_coordinate1,end_coordinate1],[start_coordinate2,end_coordinate2],...]
    with open(args.blacklist,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter='\t')
        for row in r:
            chrom = row[0]
            #print(row)
            if chrom not in blacklist: blacklist[chrom] = [[int(float(row[1])),int(float(row[2]))]]
            else: blacklist[chrom].append([int(float(row[1])),int(float(row[2]))])

    #removing the peaks that overlap with a blacklisted region
    with open(args.outfile,'wb') as outfile:
        w = csv.writer(outfile,delimiter='\t')
        with open(args.peakfile,'rb') as peakfile:
            r = csv.reader(peakfile,delimiter='\t')
            for row in r:
                if row[0].count('chrom')>1: w.writerow(row)
                else:
                    chrom = row[0]
                    if chrom in blacklist:
                        start = int(float(row[1]))
                        end = int(float(row[2]))
                        saving = True
                        for b in blacklist[chrom]:
                            if (start>=b[0] and start<=b[1]) or (end>=b[0] and end<=b[1]):
                                #there is overlap and the peak is removed
                                saving = False
                                break
                        if saving: w.writerow(row)
                    else: w.writerow(row)
#end

removeOverlappingPeaks()
