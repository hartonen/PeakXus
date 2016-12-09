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
import numpy as np
from operator import itemgetter


def FDR():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    
    #MANDATORY PARAMETERS
    parser.add_argument("peakfile",help="Full path to the file containing called peaks in igv-format.",type=str)
    parser.add_argument("outfile",help="Full path to output file (.igv).",type=str)

    #OPTIONAL PARAMETERS
    #if these are not given, defaults are used
    parser.add_argument("-t","--test_type",help="BH=Benjamini-Hochberg (default), BY=Benjamini-Hochberg-Yekutieli.",type=str,choices=['BH','BY'],default='BH')
    parser.add_argument("-m","--num_tests",help="Total number of tests conducted. If not given, the input file is assumed to contain all the test results.",type=int,default=None)
    args = parser.parse_args()

    data = [] # = [[row1],[row2],...]

    with open(args.peakfile,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter="\t")
        for row in r:
            if row[0]=="chromosome": continue
            for i in range(1,len(row)):
                if i==3: continue
                row[i] = float(row[i])
            data.append(row)

    #sorting the list according to p-value (secondary by score)
    data = sorted(data,key=lambda x: (x[5],-x[6]))

    if args.num_tests==None: m = len(data)
    else: m = args.num_tests

    #calculating the Benjamini-Hochberg corrected p-values
    #https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini.E2.80.93Hochberg.E2.80.93Yekutieli_procedure
    if args.test_type=='BY': c = sum([1.0/i for i in range(1,m)])
    else: c = 1.0
    inds = np.ones(len(data))
    for i in range(0,len(data)): inds[i] = (i+1)*inds[i]
    ps = c*m*np.array([d[-2] for d in data])/inds
    ps = np.minimum.accumulate(ps[::-1])
    ps = np.minimum(ps,np.ones(len(data)))
    ps_rev = ps[::-1]
    for i in range(len(ps_rev)): data[i].append(ps_rev[i])
    
    #saving the results as an igv-file sorted by chromosome and location
    data = sorted(data,key=lambda x: (x[0],x[1]))

    with open(args.outfile,'wb') as csvfile:
        w = csv.writer(csvfile,delimiter="\t")
        w.writerow(["chromosome","start","end","id","signal","p-value","score","FDR"])
        for row in data: w.writerow(row)
#FDR ends

FDR()
