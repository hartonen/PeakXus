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

import matplotlib
#Using the Agg-backend so that X-server is not needed
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def motifsPerPeak():


    matplotlib.rcParams.update({'font.size': 26})
    matplotlib.rcParams.update({'font.weight': 'bold'})

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    
    #MANDATORY PARAMETERS
    parser.add_argument("peakfiles",help="Full path(s) to the peak-file(s) in igv-format, sorted in the desired order.",type=str,nargs='+')
    parser.add_argument("outdir",help="Full path for the output directory.",type=str,nargs=1)
    parser.add_argument("motiffile",help="Full path to the file containing the HARS locations in gff-format.",type=str,nargs=1)
    
    #OPTIONAL PARAMETERS
    #if these are not given, defaults are used
    parser.add_argument("-m",help="Maximum distance from peak summit to motif center (default=20).",type=int,default=20)
    parser.add_argument("-o",help="If 1, each motif is only counted matching to one peak, if 0, there is no limit (default=0).",type=int,default=0)
    parser.add_argument('-t','--title',help="Title to the plot.",type=str,default=None)
    parser.add_argument('-l','--legend',help="Legend text to the plot.",type=str,default=None,nargs='+')

    args = parser.parse_args()

    matplotlib.rcParams.update({'font.size': 18})
    styles = ['-r','-b','-g','-m','-k','-c','-y','--r','--b','--g','--m','--k','--c','--y','-.r','-.b','-.g','-.m','-.k','-.c','-.y',':r',':b',':g',':m',':k',':c',':y','-*r','-*b','-*g','-*m','-*k','-*c','-*y','-Dr','-Db','-Dg','-Dm','-Dk','-Dc','-Dy','-or','-ob','-og','-om','-ok','-oc','-oy']
    largest_x = 0
    largest_y = 0

    for i in range(0,len(args.peakfiles)):

        peaks = []
        motifs = {}
        #Dictionary structure:
        #key = chrom
        #value = set of bases covered by motifs
        peakcount = 0

        #reading in peaks
        with open(args.peakfiles[i],'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            for row in r:
                if row[0].count("chrom")>0: continue
                peakcount += 1
                peaks.append(row)

        #largest_x is the largest value needed on the x-axis
        if peakcount>largest_x: largest_x = peakcount

        #reading in motifs
        with open(args.motiffile[0],'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            for row in r:
                chrom = row[0]
                start = int(float(row[3]))
                end = int(float(row[4]))
                if chrom not in motifs: motifs[chrom] = set([(start+end)/2])
                else: motifs[chrom].add((start+end)/2)

        #analyzing overlap
        x_ax = [x for x in range(0,peakcount,10)] #x-axis values
        y_ax = [0.0] #y-axis values

        with open(args.outdir[0]+"peaks_with_pwm_hit_m="+str(args.m)+"_o="+str(args.o)+".igv",'wb') as csvfile:
            w = csv.writer(csvfile,delimiter="\t")

            for j in range(1,len(x_ax)):
                hits = 0.0
                for peak in peaks[x_ax[j-1]:x_ax[j]]:
                    chrom = peak[0]
                    if chrom not in motifs: continue
                    start = int(float(peak[1]))
                    end = int(float(peak[2]))
                    summit = (start+end)/2
                    for k in range(summit-args.m,summit+args.m):
                        if k in motifs[chrom]:
                            hits += 1.0
                            #saving the peak to the cluster with the pwm-hit
                            w.writerow(peak)
                            #checking if the motif-location needs to be deleted
                            if args.o==1: motifs[chrom].discard(k)
                            break
                #hits calculated for x_ax[j], adding the corresponding y-axis value
                y_ax.append(y_ax[-1]+hits)
        #overlap calculated for args.peakfiles[i]
        #saving results to .tsv-file
        f = open(args.outdir[0]+"motifsPerPeak_m="+str(args.m)+"_o="+str(args.o)+"_"+args.legend[i]+".tsv",'wb')
        for j in range(0,len(x_ax)): f.write(str(x_ax[j])+"\t"+str(y_ax[j])+"\n")
        f.close()
        #plotting
        if y_ax[-1]>largest_y: largest_y = y_ax[-1]
        if args.legend!=None: plt.plot(x_ax,y_ax,styles[i],linewidth=4,label=args.legend[i])
        else: plt.plot(x_ax,y_ax,styles[i],linewidth=4)

    #all input files analyzed
    plt.xlim((0,largest_x))
    plt.ylim((0,largest_y))
    plt.xlabel('Best scoring peaks',fontsize=26,weight='bold')
    plt.ylabel('Number of motifs found',fontsize=26,weight='bold')
    if args.legend!=None: plt.legend(loc=4)
    if args.title!=None: plt.suptitle(args.title,y=0.93)
    plt.tight_layout()
    plt.savefig(args.outdir[0]+"motifsPerPeak_m="+str(args.m)+"o="+str(args.o)+".png")
    plt.clf()
#end


motifsPerPeak()
