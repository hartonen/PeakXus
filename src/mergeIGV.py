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

from glob import glob

def mergeIGV():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    
    #MANDATORY PARAMETERS
    parser.add_argument("inpath",help="Full path to the input igv-files. Can include wildcards/regular expressions.",type=str,nargs='+')
    parser.add_argument("outname",help="Full path to the output igv-file.",type=str,nargs=1)

    args = parser.parse_args()

    first = True
    data = []

    for f in args.inpath:
        with open(f,'rb') as csvfile:
            r = csv.reader(csvfile,delimiter='\t')
            for line in r:
                if line[0]=='chromosome':
                    if not first: continue
                
                data.append(line)
                first = False

    with open(args.outname[0],'wb') as csvfile:
        w = csv.writer(csvfile,delimiter="\t")
        for line in data: w.writerow(line)

#end

mergeIGV()
