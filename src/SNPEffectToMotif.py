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
import math
import matplotlib
#Using the Agg backend for plotting so that no X-server is required
matplotlib.use('Agg')

import numpy as np

from matplotlib import pyplot as plt
from Bio import SeqIO
from math import log
from scipy.stats import pearsonr
from scipy.stats import binom_test
from scipy import special
from math import factorial

def SNPEffectToMotif():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    #MANDATORY PARAMETERS
    parser.add_argument("asbfile",help="Full path to asb-analysis results file (i.e. the output from testASB.py).",type=str)
    parser.add_argument("freqmatrix",help="Full path to the frequency matrix file, i.e. the file containing the PWM.",type=str)
    parser.add_argument("motiffile",help="Full path to file containing the high affinity recognition sequence locations in gff-format.",type=str)
    parser.add_argument("refgenome",help="Full path to the reference genome in fasta-format.",type=str)
    parser.add_argument("outdir",help="Full path to output directory.",type=str)

    #OPTIONAL PARAMETERS
    parser.add_argument("-p",help="Pseudocount for generating the PFM.",type=float,default=0.01)
    parser.add_argument("-t",help="P-value threshold for ASB calls (default=0.01).",type=float,default=0.01)
    parser.add_argument("-AR",help="Full path to a file containing the local genomic allelic ratios (gARs) for all input SNPs. If not provided, the gAR is assumed to be 0.5 for all locations.",type=str,default=None)
    parser.add_argument("-m",help="Maximum distance for a SNP from a peak edge for it to be considered to overlap with a peak (default=5 bps).",type=int,default=5)
    parser.add_argument("--reads",help="1 if calculating the p-values for pure read counts, 0 if usign UMI-counts (0=default).",type=int,choices=[0,1],default=0)
    parser.add_argument("--legend",help="Legend to the figures.",type=str,default=None)
    parser.add_argument("--nocolor",help="1 if otput is not color coded, when 0, the color coding is such that SNPs where the alterante allele introduces an extra CG are colored red and SNPs where the reference allele has an extra CG are colored yellow (0=default).",type=int,choices=[0,1],default=1)
    parser.add_argument("--imprintfile",help="This can be used to provide an additional bed-file containing a list of imprinting control regions whose overlap is tested against the ASB significant SNPs.(default=None)",type=str,default=None)
    parser.add_argument("--test",help="1=binomial test, 2=Audic-Claverie test (2=default).",type=int,choices=[1,2],default=2)
    args = parser.parse_args()
    m = args.m

    matplotlib.rcParams.update({'font.size':20})
    matplotlib.rcParams.update({'font.weight':'bold'})

    #creating the pwm
    pwm = readMatrix(args.freqmatrix,args.p,[0.25,0.25,0.25,0.25])
    pwm_rc = pwm[::-1,::-1]
    s = pwm.shape[1] #length of the pwm
    #reading in asb-analysis results
    asb,asb_count = readASB(args.asbfile,args.t,args.reads,args.AR,args.test,args.outdir+"ASB_p-values.tsv")
    #reading in motifs
    motifs = readMotifs(args.motiffile,asb.keys())
    #This is a dictionary with dictionaries as values, structure:
    #key = chromosome
    #    key = +/-
    #    value = 5'-ends of motifs 

    #now deleting all asb sites that do not overlap with a high affinity motif
    newsites = {}
    nearbycount = 0
    significant_nearbycount = 0
    nearbytest = []
    for i in range(1,s/2+1):
        nearbytest.append(-1*i)
        nearbytest.append(s+i)

    for c in asb.keys():
        for site in asb[c]:
            found = False
            snp = int(site[0])
            peak_start = int(float(site[2].split('-')[0]))
            peak_end = int(float(site[2].split('-')[1]))
            #determining the position index in a high affinity motif for snp
            for i in range(0,s):
                if (snp not in range(peak_start-m,peak_end+m+1)): break
                if (snp-i in motifs[c]['+']) and (snp in range(peak_start-m,peak_end+m+1)):
                    #motif is on + strand
                    #the matching motif starts at position snp-i
                    #this means that snp is at index i
                    auxsite = list(site)
                    auxsite.append(i)
                    auxsite.append('+')
                    if c not in newsites: newsites[c] = [auxsite]
                    else: newsites[c].append(auxsite)
                    found = True
                elif snp-i in motifs[c]['-'] and (snp in range(peak_start-m,peak_end+m+1)):
                    #motif is on - strand
                    auxsite = list(site)
                    auxsite.append(i)
                    auxsite.append('-')
                    if c not in newsites: newsites[c] = [auxsite]
                    else: newsites[c].append(auxsite)
                    found = True
            if not found:                
                for i in nearbytest:
                    if snp-i in motifs[c]['+']:
                        nearbycount += 1
                        break
                    elif snp-i in motifs[c]['-']:
                        nearbycount += 1
                        break

    #retrieving the underlying sequences from reference genome
    handle = open(args.refgenome,'rU')
    sitecount = 0
    for record in SeqIO.parse(handle,"fasta"):
        chrom = record.id
        if chrom not in newsites.keys(): continue
        seq = record.seq
        #saving the sequence for each motif in this chromosome
        for i in range(0,len(newsites[chrom])):
            sitecount += 1
            sequence = str(seq[int(newsites[chrom][i][0])-newsites[chrom][i][-2]-1:int(newsites[chrom][i][0])-newsites[chrom][i][-2]+s-1]).upper()
            newsites[chrom][i].append(sequence)
    handle.close()
    
    #calculating scores for motifs with ref and alt allele for each remaining site
    scores = np.zeros((sitecount,2))
    xaxis = np.zeros(sitecount) #affinity change
    yaxis_rar = np.zeros(sitecount) #reference allele ratio N_alt/N_ref
    yaxis = np.zeros(sitecount) #asb-score = -log(pval)
    higher_all_count = np.zeros(sitecount) #higher allele count for each SNP location
    positions = np.zeros(sitecount)

    colors = []
    #blue if reference allele produces binding site with CG
    #red if alternate allele produces binding site with CG
    #black otherwise
    if args.imprintfile!=None:
        imprinted = []
        #True if SNP overlaps with imprinted region, False otherwise
        imprintRegions = {}
        #key = chromosome
        #value = [(start1,end2),(start2,end2),...]
        with open(args.imprintfile,'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            for row in r:
                chrom = row[0]
                if chrom not in imprintRegions: imprintRegions[chrom] = [(int(row[1]),int(row[2]))]
                else: imprintRegions[chrom].append((int(row[1]),int(row[2])))

    i = 0
    for c in newsites:
        for j in range(0,len(newsites[c])):
            REF_seq = newsites[c][j][-1]
            ind = newsites[c][j][-3] #SNP location when sequence is on + strand
            strand = newsites[c][j][-2]
	    if strand=='-': ind = s-1-ind
            snp_pos = newsites[c][j][0]
            REF = newsites[c][j][3]
            ALT = newsites[c][j][4]
            REF_seq = REF_seq[:ind]+REF+REF_seq[ind+1:]
            ALT_seq = newsites[c][j][-1]
            ALT_seq = ALT_seq[:ind]+ALT+ALT_seq[ind+1:]
            if args.reads==0:
                N_ref = float(newsites[c][j][6])
                N_alt = float(newsites[c][j][7])
            else:
                N_ref = float(newsites[c][j][8])
                N_alt = float(newsites[c][j][9])
                
            hit = False
            if strand=='+':
                scores[i,:] = np.array([calcScore(pwm,REF_seq),calcScore(pwm,ALT_seq)])
                if ind<(len(REF_seq)-1):
                    if REF=='C' and REF_seq[ind+1]=='G':
                        if args.nocolor==1: colors.append('k')
                        else: colors.append('y')
                        hit = True
                    if (not hit) and ALT=='C' and ALT_seq[ind+1]=='G':
                        if args.nocolor==1: colors.append('k')
                        else: colors.append('r')
                        hit = True
                if ind>0:
                    if (not hit) and REF_seq[ind-1]=='C' and REF=='G':
                        if args.nocolor==1: colors.append('k')
                        else: colors.append('y')
                        hit = True
                    if (not hit) and ALT_seq[ind-1]=='C' and ALT=='G':
                        if args.nocolor==1: colors.append('k')
                        else: colors.append('r')
                        hit = True
            else:               
                scores[i,:] = np.array([calcScore(pwm_rc,REF_seq),calcScore(pwm_rc,ALT_seq)])
                if ind<(len(REF_seq)-1):
                    if REF=='C' and REF_seq[ind+1]=='G':
                        if args.nocolor==1: colors.append('k')
                        else: colors.append('y')
                        hit = True
                    if (not hit) and ALT=='C' and ALT_seq[ind+1]=='G':
                        if args.nocolor==1: colors.append('k')
                        else: colors.append('r')
                        hit = True
                if ind>0:
                    if (not hit) and REF_seq[ind-1]=='C' and REF=='G':
                        if args.nocolor==1: colors.append('k')
                        else: colors.append('y')
                        hit = True
                    if (not hit) and ALT_seq[ind-1]=='C' and ALT=='G':
                        if args.nocolor==1: colors.append('k')
                        else: colors.append('r')
                        hit = True

            if not hit: colors.append('k')
            yaxis[i] = -log(float(newsites[c][j][5]))
            yaxis_rar[i] = N_ref/(N_alt+N_ref)
            higher_all_count[i] = max([N_alt,N_ref])

            #testing if SNP overlaps an imprint control region if imprintfile is provided
            if args.imprintfile!=None:
                found = False
                if c in imprintRegions:
                    for site in imprintRegions[c]:
                        if int(snp_pos)>=site[0] and int(snp_pos)<=site[1]:
                            found = True
                            break
                imprinted.append(found)

            positions[i] = ind
            xaxis[i] = (scores[i,0])-(scores[i,1])
            i += 1

    #calculating the pearson correlation coefficient between -log(p) and affinity change
    r_pearson,p_pearson = pearsonr(xaxis,yaxis)

    #plotting
    sc = plt.scatter(xaxis,yaxis,s=np.pi*10*np.ones(sitecount),c=colors,label=args.legend)
    t = [0,0.2,0.4,0.6,0.8,1.0]

    #if imprintfile is provided, plotting green circles around SNPs at ICRs
    if args.imprintfile!=None:
        x_imprinted = []
        y_imprinted = []
        for i in range(0,len(imprinted)):
            if imprinted[i]:
                x_imprinted.append(xaxis[i])
                y_imprinted.append(yaxis[i])
        sc2 = plt.scatter(x_imprinted,y_imprinted,s=np.pi*30*np.ones(sitecount),facecolors='none',edgecolors='g',linewidth=2)
    plt.ylabel("-log(p-value)",fontsize=20,weight='bold')
    plt.xlabel("Affinity change (S_ref-S_alt)",fontsize=20,weight='bold')
    if args.legend!=None: plt.legend(loc=2)
    plt.suptitle("r = "+str(round(r_pearson,2))+", p-value = "+str('{:.2e}'.format(p_pearson)))
    plt.tight_layout()
    plt.savefig(args.outdir+"pval_vs_affinity_change.png")
    plt.clf()

    #plotting the reference allele ratio vs affinity change

    #calculating the pearson correlation coefficient between reference allele ratio and affinity change
    r_pearson,p_pearson = pearsonr(xaxis,yaxis_rar)
    
    sc = plt.scatter(xaxis,yaxis_rar,s=np.pi*10*np.ones(sitecount),c=colors,label=args.legend)
    t = [0,0.2,0.4,0.6,0.8,1.0]

    #if imprintfile is provided, plotting green circles around SNPs at ICRs
    if args.imprintfile!=None:
        x_imprinted = []
        y_imprinted = []
        for i in range(0,len(imprinted)):
            if imprinted[i]:
                x_imprinted.append(xaxis[i])
                y_imprinted.append(yaxis_rar[i])
        sc2 = plt.scatter(x_imprinted,y_imprinted,s=np.pi*30*np.ones(sitecount),facecolors='none',edgecolors='g',linewidth=2)

    plt.ylabel("Reference allele ratio",fontsize=20,weight='bold')
    plt.xlabel("Affinity change",fontsize=20,weight='bold')
    if args.legend!=None:
        leg = plt.legend(loc=2,fancybox=True)
        leg.get_frame().set_alpha(0.5)

    plt.suptitle("r = "+str(round(r_pearson,2))+", p-value = "+str('{:.2e}'.format(p_pearson)))
    plt.tight_layout()
    plt.savefig(args.outdir+"refallele_ratio_vs_affinity_change.png")
    plt.clf()

    createHTML(asb_count,sitecount,nearbycount,args)

    print "Total number of SNPs overlapping with an HARS = "+str(sitecount)
    print "Number of SNPs within the flanking "+str(2*s/2)+" bps of an HARS-site = "+str(nearbycount)
#end

def createHTML(asb_count,SNP_overlap_HARS,SNP_flanking_HARS,args):
    #asb_count = total number of peaks overlapping with a SNP
    #SNP_overlap_HARS = total number of peaks overlapping with a SNP and an HARS
    #SNP_flanking_HARS = total number of SNPs flanking an HARS

    f = open(args.outdir+"ASB_results.html",'wb')

    ###########################
    #The required declarations#
    ###########################
    f.write("<!DOCTYPE html>\n<html>\n<head>")

    ############
    #Page style#
    ############

    f.write('<style>')
    f.write('h1 {')
    f.write('color: #ccccff;')
    f.write('font-family: verdana;')
    f.write('font-size: 400%;')
    f.write('}')

    f.write('h2 {')
    f.write('color: white;')
    f.write('font-family: verdana;')
    f.write('font-size: 300%;')
    f.write('}')

    f.write('h3 {')
    f.write('color: white;')
    f.write('font-family: verdana;')
    f.write('font-size: 200%;')
    f.write('}')

    f.write('table.td {')
    f.write('height: 30px;')
    f.write('horizontal-align: center;')
    f.write('}')

    f.write('table.td.div {')
    f.write('height: 100%;')
    f.write('}')
    
    f.write('div.caption {')
    f.write('color: black;')
    f.write('align: center;')
    f.write('font-family: Times;')
    f.write('font-size: 100%;')
    f.write('font-style: italic;')
    f.write('}')

    f.write('p {')
    f.write('color: black;')
    f.write('font-family: Times;')
    f.write('font-size: 160%;')
    f.write('}')

    f.write('p.white {')
    f.write('color: white;')
    f.write('font-family: verdana;')
    f.write('font-style: bold')
    f.write('font-size: 160%;')
    f.write('}')

    f.write('p.small {')
    f.write('color: black;')
    f.write('font-family: Times;')
    f.write('font-size: 100%;')
    f.write('}')
    f.write('</style>\n<body>')


    
    #########
    #Heading#
    #########
    f.write('<div style="background-color:white; color:#ccccff; padding:20px;">')
    f.write("<h1>ASB-ANALYSIS RESULTS</h1>")
    f.write('</div>')

    #########
    #Results#
    #########
    f.write('<table style="width:100%">')
    f.write('<tr>')
    f.write("<td><p>Total number of peaks overlapping a SNP="+str(asb_count)+"</p></td>")
    f.write("</tr>")
    f.write("<tr>")
    f.write("<td><p>Total number of SNPs overlapping a peak and an HARS="+str(SNP_overlap_HARS)+"</p></td>")
    f.write("</tr>")
    f.write("</tr></table>")

    f.write('<table style="width:100%">')
    f.write("<tr>")
    f.write('<td><div class="image"><img src="'+args.outdir+'refallele_ratio_vs_affinity_change.png" alt="Reference allele ratio vs affinity change" style="width:800px;height:600px;"></div></td>') 
    f.write('<td><div class="image"><img src="'+args.outdir+'pval_vs_affinity_change.png" alt="Peak size histogram" style="width:800px;height:600px;"></div></td>')
    f.write("</tr>")
    f.write('<tr>')
    f.write('<td><div class="caption">Reference allele ratio as a function of binding sequence affinity change. Y-axis shows reference allele ratio while x-axis is the affinity change (reference minus alternate sequence affinity). SNPs with p-value > '+str(args.t)+' are filtered out. Red dots represent SNPs where alternate allele creates a CG-site to the sequence but reference does not. Yellow dots represent SNPs where sequence with reference allele has an extra CG. Other SNPs are black. Value of Pearson correlation coefficient (r) along with the corresponding p-value are shown above the figure.</div></td>')
    f.write('<td><div class="caption">-log(p-value) as a function of binding sequence affinity change. Y-axis shows negative logarithm of the p-value while x-axis is the affinity change (reference minus alternate sequence affinity). SNPs with p-value > '+str(args.t)+' are filtered out. Red dots represent SNPs where alternate allele creates a CG-site to the sequence but reference does not. Yellow dots represent SNPs where sequence with reference allele has an extra CG. Other SNPs are black. Value of Pearson correlation coefficient (r) along with the corresponding p-value are shown above the figure.</div></td>')
    f.write('</tr>')
    f.write("</table>")

    ##################
    #Input parameters#
    ##################

    f.write('<div style="background-color:#ccccff; color:white; padding:20px;">')
    f.write('<h3>INPUT PARAMETERS FOR SNPEffectToMotif.py</h3>')
    f.write('</div>')
    arguments = []

    args_dict = vars(args)
    for k in args_dict: arguments.append(str(k)+"="+str(args_dict[k]))
    f.write('<p class="small">\n')
    for a in arguments: f.write(a+'<br>\n')
    f.write('</p>')

    f.write('<div style="background-color:#ccccff; color:white; padding:20px;">')
    f.write('<p class="white">If you use PeakXus in your work, please cite:<br>Hartonen T, Sahu B, Dave K, Kivioja T and Taipale J, PeakXus: Comprehensive Transcription Factor Binding Site Discovery From ChIP-Nexus and ChIP-exo Experiments, to be published.</p>')
    f.write('</div>')

    #end
    f.write("</body>\n</html>")
    f.close()


def readMotifs(path,chroms):
    #path = path to motif file
    #chroms = chromosomes used in analysis

    motifs = {}
    for c in chroms: motifs[c] = {'+':set(),'-':set()}

    #saving the 5'-end index 
    with open(path,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter='\t')
        for row in r:
            chrom = row[0]
            if chrom in chroms:
                strand = row[6]
                start = int(row[3])
                
                motifs[chrom][strand].add(start)

    return motifs

def readASB(path,threshold,reads,AR,test,outfile):
    #AR = path to file containing the allelic ratios
    #outfile = path to a tsv-file where the ASB-results are saved
    counter = 0
    #reading in the allelic ratios
    if AR!=None:
        ARs = {}
        #key = chromosome
        #    value/key = {(id1,pos1),...}
        #        value = AR1
        with open(AR,'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            for row in r:
                chrom = row[0]
                id = row[2]
                pos = int(row[1])
                refcount = float(row[3])
                altcount = float(row[4])
                if chrom not in ARs: ARs[chrom] = {(id,pos):(refcount,altcount)}
                else: ARs[chrom][(id,pos)] = (refcount,altcount)
    else: ARs = None

    asb = {}
    #key = chrom
    #value = list of sites on chromosome chrom

    with open(path,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter='\t')
        with open(outfile,'wb') as out:
            w = csv.writer(out,delimiter='\t')
            w.writerow(["#chrom","location","SNP_id","peak_location","REF_allele","ALT_allele","N_REF","N_ALT","allreads_REF","allreads_ALT","p-value"])
            for row in r:
                if row[0].count('#')>0: continue
                #filtering out variations that aren't snps
                counter += 1
                if len(row[4])>1 or len(row[5])>1: continue
                #filtering out SNPs where there are zero UMI hits
                if float(row[-3])+float(row[-4])<1: continue
                chrom = row[0]
                pos = int(row[1])
                id = row[2]
                if row[-1]=='': row = row[1:-1]
                else: row = row[1:]
                if ARs!=None:
                    if (id,pos) not in ARs[chrom]: continue
                    wgs_refcount = ARs[chrom][(id,pos)][0]
                    wgs_altcount = ARs[chrom][(id,pos)][1]
                    g = float(wgs_refcount)/(float(wgs_refcount)+float(wgs_altcount))
                else: g = 0.5

                if reads==1:
                    #this means we calculate p-values for pure read counts
                    refreads = float(row[-2])
                    altreads = float(row[-1])
                    #if refreads<3 or altreads<3: continue #this is a test
                    #test if there is too few hits to rarer allele
                    if test==1: pval = binom_test(refreads,refreads+altreads,g)
                    else:
                        #winflat test
                        #x = reference allele hits in wgs
                        #N1= total num of reads in wgs
                        #y = reference allele hits in chip
                        #N2= total num of reads in chip
                        pval = winflat(wgs_refcount,wgs_refcount+wgs_altcount,refreads,refreads+altreads)
                    if pval==0: pval = 0.0000000000000001
                    row = row[0:5]+[pval]+row[6:]
                else:
                    #this means p-value is calculated from umi counts
                    refumis = float(row[-4])
                    altumis = float(row[-3])
                    #if refumis<3 or altumis<3: continue #this is a test
                    #test if there is too few hits to rarer allele
                    if test==1: pval = binom_test(refumis,refumis+altumis,g)
                    else:
                        #winflat test
                        #x = reference allele hits in wgs
                        #N1= total num of reads in wgs
                        #y = reference allele hits in chip
                        #N2= total num of reads in chip
                        pval = winflat(wgs_refcount,wgs_refcount+wgs_altcount,refumis,refumis+altumis)
                    if pval==0: pval = 0.0000000000000001
                    row = row[0:5]+[pval]+row[6:]

                if threshold<0:
                    aux = abs(threshold)
                    if pval<=aux: continue
                else:
                    if pval>threshold: continue
                if chrom not in asb: asb[chrom] = [row]
                else: asb[chrom].append(row)
                w.writerow(row)

    return asb,counter
#end

def readMatrix(path,p,bgfreqs):
    #path = path to frequency matrix file
    #p = pseudocount
    #bgfreqs = background frequencies = [p(A),p(C),p(G),p(T)]

    pwm = []
    with open(path,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter='\t')
        for row in r: pwm.append([float(i) for i in row])

    pwm = np.array(pwm)
    #applying the pseudocount, normalizing and calculating scores
    for i in range(0,pwm.shape[1]):
        pwm[:,i] += np.sum(pwm[:,i])*p
        pwm[:,i] /= np.sum(pwm[:,i])
        pwm[:,i] = np.log(pwm[:,i]/bgfreqs)

    return pwm
#end

def calcScore(pwm,motif):
    score = 100.0
    for i in range(0,len(motif)):
        l = motif[i]
        if l=='A': score += pwm[0,i]
        elif l=='C': score += pwm[1,i]
        elif l=='G': score += pwm[2,i]
        else: score += pwm[3,i]
    return score
#end

def revComp(seq):
    #returns the reverse complement sequence of seq
    #seq = seq[::-1]
    newseq = ""
    for s in seq[::-1]:
        if s=='A': newseq += 'T'
        elif s=='T': newseq += 'A'
        elif s=='G': newseq += 'C'
        elif s=='C': newseq += 'G'
    return newseq

def binomm(k,n): return special.binom(float(n),float(k))*0.5**float(n)

def winflat(x,N1,y,N2):
    #Audric & Claverie Genome Research 1997 equation 2
    #p(y|x)=(N2/N1)^y*((x+y)!/(x!y!(1+N2/N1)^(x+y+1)))
    #-> ln(p(y|x)) = y*ln(N2/N1)+ln((x+y)!)-ln(x!)-ln(y!)-(x+y+1)*ln(1+N2/N1)
    #                            =xy_fac    =x_fac =y_fac
    #x = float(x)
    #y = float(y)
    N1 = float(N1)
    N2 = float(N2)
    
    if N1<1: return 1.0
    #if numbers are too large for calculating factorials, we use the Stirling approximation
    #ln(n) = nln(n)-n
    try: x_fac = math.log(factorial(x))
    except OverflowError: x_fac = x*math.log(x)-x
    try: y_fac = math.log(factorial(y))
    except OverflowError: y_fac = y*math.log(y)-y
    try: xy_fac = math.log(factorial(x+y))
    except OverflowError: xy_fac = (x+y)*math.log(x+y)-(x+y)

    try:
        logp = y*math.log(N2/N1)+xy_fac-x_fac-y_fac-(x+y+1)*math.log(1+(N2/N1))
        return math.exp(logp)
    except OverflowError as e:
        print e
        return 1.0

SNPEffectToMotif()

