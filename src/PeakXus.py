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
from os import system
from sys import exit
import csv
import textwrap

import multiprocessing as mp

def PeakXus():

    ########################
    #COMMAND LINE ARGUMENTS#
    ########################

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description=textwrap.dedent('''\

    ***  *****   *   *  * *   * *   *  ****
    *  * *      * *  * *   * *  *   * * 
    ***  *****  ***  **     *   *   *  ****
    *    *     *   * * *   * *  *   *      *
    *    ***** *   * *  * *   *  ***  ****

    Copyright (c) Tuomo Hartonen, 2015-2016

    THIS PROGRAM COMES WITH ABSOLUTELY NO WARRANTY!
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License version 2 as 
    published by the Free Software Foundation.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see http://www.gnu.org/licenses/.
    
    ----------
    |OVERVIEW|
    ----------

    This script performs the following functionality if a fastq-file is given as an input:
    
      1) Trimming adapters from ends of reads using Cutadapt-tool.
      2) Aligning the reads in the input fastq-file to a specified reference genome using bwa aln-algorithm.
      3) Peak calling using PeakXus.
      4-6) Plotting different kind of figures for graphical inspection of the peak calling results.
      7) Mining the best peaks for common sequence motifs with MEME.
    
    Steps 1-2 are skipped if the input is already in bam-format. Steps 4-6 are optional as part of the output figures are produced only if a list of high-affinity recognition sequence (HARS) locations is given as an input (see --matrixhits option below). Summary of the results is shown in PeakXus_results.html. All reported peaks can be found from all_transition_points.igv in igv-format.

    Below is a detailed description of all the possible input parameters for different stages of the pipeline.

    ----------------
    |USAGE EXAMPLES|
    ----------------

    Majority of the optional input parameter values seldom need tweaking. Below are listed some example
    calls of the pipeline for conducting different tasks. Input file locations are assumed to be:

    ref_genome_path/genome.idx = reference genome index
    ref_genome_path/wg.fasta   = reference genome in fasta-format
    exp_path/exp.fastq         = reads from ChIP-Nexus experiment in fastq-format
    exp_path/exp.bam           = aligned reads from ChIP-Nexus experiment in bam-format
    exp_path/UMIs.txt          = UMI-labels used in the experiment
    exp_path/chroms.txt        = chromosome names in the ChIP-Nexus experiment
    motif_path/matrix_hits.gff = hits of the corresponding TF-pwm to the used reference genome in gff-format

    It is important to give the correct length of UMI-labels for PeakXus! In the following examples it is assumed that the UMI-label length is 5 bps.

      i) PEAK CALLING ONLY FROM A BAM-FILE
      --------------------------------------

        without UMIs:
        PeakXus.py exp_path/exp.bam out/ exp_path/chroms.txt

        with UMIs:
        PeakXus.py exp_path/exp.bam out/ exp_path/chroms.txt --UMIs exp_path/UMIs.txt --l3 5

      ii) PEAK CALLING ONLY FROM A FASTQ-FILE
      ---------------------------------------

        * aligning with 8 parallel processes
        PeakXus.py exp_path/exp.fastq out/ exp_path/chroms.txt --genome2 ref_genome_path/genome.idx --t2 8

        * UMIs of length 5 used, no parallel processes in aligning
          - five first bases of each read are saved to the BC-field of the bam-file
        PeakXus.py exp_path/exp.fastq out/ exp_path/chroms.txt --genome2 ref_genome_path/genome.idx --UMIs exp_path/UMIs.txt --B2 5 --l3 5

      iii) PEAK CALLING FROM A BAM-FILE, PLOTTING GRAPHICAL RESULTS
      -------------------------------------------------------------

        * assuming the high-affinity recognition sequence width is 8 bps
        PeakXus.py exp_path/exp.bam out/ exp_path/chroms.txt --matrixhits motif_path/matrix_hits.gff --ml4 8 --ms4 -4 --l3 5

      MORE EXAMPLES AND DETAILS: https://github.com/hartonen/PeakXus/wiki

    ----------------------------------
    |DETAILED LIST OF INPUT ARGUMENTS|
    ----------------------------------
    '''))
    
    #INPUT FILES, MANDATORY
    parser.add_argument("fastq",help="Full path to the input fastq/bam-file. If input is a bam-file, alignment-step is skipped and the pipeline starts straight from peak calling. File type is determined based on the file name extension.",type=str)
    parser.add_argument("outdir",help="Full path to output directory. This has to exist and be writable!.",type=str)
    parser.add_argument("chromnames",help="Full path to a file containing chromosome names and sizes each chromosome on its own line, name and size separated by tab.",type=str)
    
    #INPUT FILES, OPTIONAL
    parser.add_argument("--UMIs",help="Full path to a file containing all UMI-labels. The file should have two tab-separated columns, first is the UMI name (e.g. BC1) and second is the actual label sequence (e.g. ATTAG). If this argument is not provided, UMIs are not used in the analysis.",type=str,default=None)
    parser.add_argument("--matrixhits",help="Full path to a gff-file containing the locations of the HARS sites (High-Affinity Recognition Sequence) in the appropropriate reference genome. If this file is not provided, part of the results are not calculated.",type=str,default=None)
    parser.add_argument("--blacklist",help="Full path to a bed-file containing the start and end coordinates of blacklisted genomic regions. Peaks overlapping these regions are filtered out from the results if blacklist is given.",type=str,default=None)

    #0 FASTQC, input parameters
    group0 = parser.add_argument_group('0. FASTQC','These arguments are passed to the FastQC-tool.')
    group0.add_argument('--skip0',help='If 1, FastQC-analysis will be skipped (default=0).',type=int,choices=[0,1],default=0)

    #1 CUTADAPT, input parameters
    group1 = parser.add_argument_group('1. CUTADAPT','These arguments are passed to the Cutadapt-tool for trimming adapters from the\nends of reads from the input fastq-file. This option is not available for bam-files.')
    group1.add_argument("--adapters1",help="Adapter-sequences trimmed from the ends of reads using Cutadatp. If not provided, this step is skipped.",type=str,nargs='+',default=None)

    #2 BWA ALN, input parameters
    group2 = parser.add_argument_group('2. BWA','These arguments are passed to the bwa aln-algorithm for aligning the reads to\nthe chosen reference genome.')
    group2.add_argument("--genome2",help="Full path to the reference genome index. This is a mandatory parameter for fastq-input!",type=str,default=None)
    group2.add_argument("--B2",help="Length of UMI-labels in bps.",type=int,default=None,nargs=1)
    group2.add_argument("--q2",help="Quality threshold for base calls (default=20).",type=int,default=20)
    group2.add_argument("--t2",help="Number of parallel processes used by the aligner (default=1).",type=int,default=1)

    #3 PEAK CALLING, input parameters
    group3 = parser.add_argument_group('3. PEAKXUS','These arguments are passed to the PeakXus-peak caller.')
    group3.add_argument("--w3", help="Width of the window around the 5'-end of a candidate peak. (default=5).", type=int, default=5)
    group3.add_argument("--l_l3", help="l_limit is the largest allowed width for a peak (default=60).", type=int, default=60)
    group3.add_argument("--l3", help="Length of the UMI-labels used (default=5).",type=int,default=5)
    group3.add_argument("--b3",help="Threshold p-value for background model (default=0.05).", type=float, default=0.05)
    group3.add_argument("--s3",help="Width of a sliding window used in averaging density of 5' ends of reads. This should be an odd number (default=1).",type=int,default=1)
    group3.add_argument("--n3",help="Number of parallel processes used in calculating significance of peaks (default=1).",type=int,default=1)
    group3.add_argument("--p3",help="Pseudocount added to each position when testing significance of candidate peaks (default=1.0).",type=float,default=1.0)

    #4 READ DENSITIES AROUND PEAKS, input parameters
    group4 = parser.add_argument_group('4. READ COUNT PROFILE FIGURES','These arguments are passed to a script that plots figures detailing the read\ncount profiles around called peaks. Some of the figures are not produced if\n--matrixhits argument is not given.')
    group4.add_argument("--N4",help="Number of top peaks considered (default=1000).",type=int,default=1000)
    group4.add_argument("--m4",help="Maximum distance between peak and motif for them to be considered to overlap (default=10).",type=int,default=10)
    group4.add_argument("--w4",help="Width of the plotted region around the peak in heatmaps (default=100 bps).",type=int,default=100)
    group4.add_argument("--e4",help="Experiment name for plot titles.",type=str,default="PeakXus")
    group4.add_argument("--n4",help="Number of parallel processes used in clustering the peaks (default=1).",type=int,default=1)
    group4.add_argument("--k4",help="Number of clusters for k-means clustering (default=4).",type=int,default=4)
    group4.add_argument("--S4",help="Window size used in calculating read distributions (default=1).",type=int,default=1)
    group4.add_argument("--ms4",help="Start position of the HARS when coordinate of the center of the HARS is 0 (default=-8).",type=int,default=-8)
    group4.add_argument("--ml4",help="Motif length (default=17).",type=int,default=17)
    
    #5 BINDING MOTIF HITS PER PEAK, input parameters
    group5 = parser.add_argument_group('5. ANALYSIS OF OVERLAP BETWEEN PEAKS AND HARS-SITES','These arguments are passed to a script that calculates the number of peaks that\noverlap with a TF-specific high-affinity binding sequence. These results are not\ncalculated if --matrixhits argument is not given.')
    group5.add_argument("--m5",help="Maximum distance from peak summit to motif center. (default=20).",type=int,default=20)
    group5.add_argument("--o5",help="If 1, each motif is only counted matching to one peak, if 0, there is no limit (default=1).",type=int,default=1)
    group5.add_argument('--title5',help="Figure title.",type=str,default="PeakXus")

    #6 PEAK DISTANCE FROM MOTIF, input parameters
    group6 = parser.add_argument_group('6. ANALYSIS OF PEAK LOCALIZATION ACCURACY','These arguments are passed to a script that calculates the distances of called\npeaks from TF-specific high-affinity binding sequences. These results are not\ncalculated if --matrixhits argument is not given.')
    group6.add_argument("--l6","--legend",help="Figure legend.",type=str,default="PeakXus")
    group6.add_argument("--m6",help="Maximum distance between peak and motif, i.e. the end point of the x-axis (default=100).",type=int,default=100)
    group6.add_argument("--N6",help="Number of top peaks considered (default=1000).",type=int,default=1000)
    group6.add_argument("--center6",help="If 1, only center to center distances calculated, otherwise also left edge to left edge, and right edge to right edge plots are produced (default=0).",type=int,default=0)

    parser.add_argument("-v","--verbosity",help="1 (default) if progress is printed to screen, 0 otherwise.",type=int,default=1)

    #7 MEME, input parameters
    group7 = parser.add_argument_group('7. MEME','These arguments are passed to MEME for motif enrichment analysis.')
    group7.add_argument("--wg7",help="Full path to reference genome file in fasta-format, this is mandatory for MEME analysis to be run!",type=str,default=None)
    group7.add_argument("--minw7",help="Minimun motif width (default=4).",type=int,default=4)
    group7.add_argument("--maxw7",help="Maximum motif width (default=50).",type=int,default=50)
    group7.add_argument("--p7",help="Use parallel version with <p> (default=1) processors.",type=int,default=1)

    a = parser.parse_args()

    #log-file is written into a.outdir[0]/log.txt
    f = open(a.outdir+"log.txt",'a')
    f.write("INPUT PARAMETERS:\n")
    f.write(str(a))
    f.close()

    if a.verbosity==1: print "Testing input parameters...",
    if not testinput(a): exit(0)
    if a.verbosity==1: print "succesful!"

    chroms = []
    with open(a.chromnames,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter="\t")
        for row in r: chroms.append(row[0])

    #0) RUNNING FASTQC FOR QUALITY CHECK OF THE INPUT FILES

    if a.skip0==0:
        if a.verbosity==1: print "Running FastQC...",
        system("mkdir "+str(a.outdir)+"/FastQC/")
        command = "fastqc "+a.fastq
        system(command+" -q -o "+str(a.outdir)+"/FastQC")
        if a.verbosity==1: print "done!"

    #This defines the number of parallel samtools-processes used
    nproc = max([a.t2,a.n3,a.p7])

    if a.fastq[-6:]==".fastq":
        #1) TRIMMING ADAPTERS FROM THE ENDS OF READS

        aln_inname = a.fastq
        if a.adapters1!=None:
            if a.verbosity==1: print "Trimming adapters...",
            ca_call = "~/.local/bin/cutadapt -m 22 -O 4 -e 0.2 "
            for adapter in a.adapters1: ca_call += "-a "+adapter+" "
            ca_call += a.fastq
            ca_call += " > "
            ca_call += a.fastq[:-6]+"_trimmed.fastq"
            system(ca_call)
            aln_inname = a.fastq[:-6]+"_trimmed.fastq"
            if a.verbosity==1: print "succesful!"

        #2) ALIGNING
        if a.verbosity==1: print "Aligning...",
        
        full_bam = a.outdir+aln_inname.split('/')[-1][:-6]+"_sorted_filtered.bam"

        if a.B2==None:
            #This means no UMIs are used
            aln_call = "bwa_align.py "+aln_inname+" "+a.genome2+" "+a.outdir+" -t "+str(a.t2)+" -q "+str(a.q2)
            system(aln_call)

        else:
            #This means UMIs are used
            aln_call = "bwa_align.py "+aln_inname+" "+a.genome2+" "+a.outdir+" -t "+str(a.t2)+" -q "+str(a.q2)+" -B "+str(a.B2[0])
            system(aln_call)
            
        #Using UMI-labels of different lengts is deprecated
        #else:
        #    for B in a.B2:
        #        aln_call = "bwa_align.py "+aln_inname+" "+a.genome2+" "+a.outdir[0]+"BC"+str(B)+"_ -t #"+str(a.t2)+" -q "+str(a.q2)+" -B "+str(B)
        #        system(aln_call)
        #    system("samtools merge "+full_bam+" "+a.outdir[0]+"BC*_"+aln_inname.split('/')[-1][:-6]+"_sorted_filtered.bam")
        #    system("samtools index "+full_bam)
        if a.verbosity==1: print "succesful!"

    else:
        #This means that the input file is already aligned and in bam-format
        full_bam = a.fastq
            
    #3 PEAK CALLING

    #creating a UMI-input file for peak calling and later applications (only one column,the actual UMI-label).
    if a.UMIs!=None:
        with open(a.UMIs,'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            with open(a.outdir+"UMI.bc",'wb') as outfile:
                w = csv.writer(outfile,delimiter="\t")
                for row in r: w.writerow([row[1]])

    #Peak calling is performed by calling the PeakC6.py-script
    if a.verbosity==1: print "Calling peaks..."
    numtests = 0
        
    pc_call ="peakC6.py "+full_bam+" "+a.outdir+" "+a.chromnames

    pc_call += " -w "+str(a.w3)+" -l_l "+str(a.l_l3)+" -l "+str(a.l3)+" -b "+str(a.b3)+" -s "+str(a.s3)+" -n "+str(a.n3)+" -p "+str(a.p3)
    if a.UMIs!=None: pc_call += " -u "+a.outdir+"UMI.bc"

    if a.verbosity==1: print pc_call
    system(pc_call)
    with open(a.outdir+"numtests.txt",'r') as csvfile:
        r = csv.reader(csvfile,delimiter="\t")
        for row in r:
            numtests += int(float(row[0]))
            break

    #Merging the results from different chromosomes
    system("mergeIGV.py "+a.outdir+"chr*.igv "+a.outdir+"all_transition_points_nofdr.igv")
    #Removing the auxiliary files that have the peak calling peaks in separate files for different chromosomes
    system("rm "+a.outdir+"chr*.igv")

    #if a bed-file containing blacklisted genomic regions is provided, peaks overlapping the blacklisted regions are filtered out
    if a.blacklist!=None:
        system("removeOverlappingPeaks.py "+a.outdir+"all_transition_points_nofdr.igv "+a.blacklist+" "+a.outdir+"all_transitions_points_filtered.igv")
        system("mv "+a.outdir+"all_transition_points_nofdr.igv "+a.outdir+"all_transition_points_nofdr_nofilter.igv")
        system("mv "+a.outdir+"all_transitions_points_filtered.igv "+a.outdir+"all_transition_points_nofdr.igv")

    #Calculating the false discovery rates with the Benjamini-Hochberg procedure
    system("FDR.py "+a.outdir+"all_transition_points_nofdr.igv "+a.outdir+"all_transition_points.igv -m "+str(numtests))

    #Sorting the transition points according to the peak score
    system("sort -k7,7gr "+a.outdir+"all_transition_points.igv > "+a.outdir+"all_transition_points_sorted.igv")

    if a.verbosity==1: print "succesful!"

    #Creating igv-files that show the 5'-end count distributions of UMIs/reads
    if a.UMIs!=None:
        nproc = max([a.t2,a.n3,a.p7])
        for s in ['+','-']:
            system("BamToStrandSpecific5PrimeCoverage.py "+full_bam+" "+a.outdir+" "+a.chromnames+" -s "+s+" -u "+a.outdir+"UMI.bc -n "+str(nproc)+" -l "+str(a.l3))
            #merging the files and deleting auxiliary files
            system("mergeIGV.py "+a.outdir+"chr*_"+s+".igv "+a.outdir+"5primeCoverage_"+s+".igv")
            system("rm "+a.outdir+"chr*_"+s+".igv")

    #4 READ COUNT DENSITIES AROUND PEAKS

    if a.verbosity==1: print "Plotting results...",

        
    if a.matrixhits!=None: fr_call = "fancyResults6_lowMem.py "+a.outdir+"all_transition_points_sorted.igv "+full_bam+" "+a.outdir+" --motiffile "+a.matrixhits+" -N "+str(a.N4)+" -m "+str(a.m4)+" -w "+str(a.w4)+" -e "+str(a.e4)+" -n "+str(a.n4)+" -k "+str(a.k4)+" -S "+str(a.S4)+" -ms "+str(a.ms4)+" -ml "+str(a.ml4)+" -l "+str(a.l3)
    else: fr_call = "fancyResults6_lowMem.py "+a.outdir+"all_transition_points_sorted.igv "+full_bam+" "+a.outdir+" -N "+str(a.N4)+" -m "+str(a.m4)+" -w "+str(a.w4)+" -e "+str(a.e4)+" -n "+str(a.n4)+" -k "+str(a.k4)+" -S "+str(a.S4)+" -ms "+str(a.ms4)+" -ml "+str(a.ml4)+" -l "+str(a.l3)
    if a.UMIs!=None: fr_call = fr_call+" -u "+a.outdir+"/UMI.bc"
    system(fr_call)

    #5 BINDING MOTIF HITS PER PEAK
    if a.matrixhits!=None:
        #This is only calculated if the list of high-affinity recognition sequence locations is given as input.
        mh_call = "motifsPerPeak.py "+a.outdir+"all_transition_points_sorted.igv "+a.outdir+" "+a.matrixhits+" -m "+str(a.m5)+" -o "+str(a.o5)+" -t "+a.title5+" -l PeakXus"
        system(mh_call)

    #6 DISTANCE FROM MOTIF
    if a.matrixhits!=None:
        #This is only calculated if the list of high-affinity recognition sequence locations is given as input.
        if a.N6==None: cd_call = "cumul_peak_dist_from_motif.py "+a.matrixhits+" "+a.outdir+"dist_from_motif_m="+str(a.m6)+"_N=all.png "+a.outdir+"all_transition_points_sorted.igv -p g -l "+a.l6+" -m "+str(a.m6)+" --center "+str(a.center6)
        else: cd_call = "cumul_peak_dist_from_motif.py "+a.matrixhits+" "+a.outdir+"dist_from_motif_m="+str(a.m6)+"_N="+str(a.N6)+".png "+a.outdir+"all_transition_points_sorted.igv -p g -l "+a.l6+" -m "+str(a.m6)+" -N "+str(a.N6)+" --center "+str(a.center6)
        system(cd_call)
    if a.verbosity==1: print "succesful!"

    #7 MEME
    if a.wg7!=None:
        #Creating fasta-input for MEME
        system("seqs_under_peaks.py "+a.wg7+" "+a.outdir+"all_transition_points_sorted.igv "+a.outdir+"top"+str(a.N4)+"_peaks.fasta -N "+str(a.N4)+" -w 50")
        #Calling MEME
        system("meme "+a.outdir+"top"+str(a.N4)+"_peaks.fasta -oc "+a.outdir+"/meme_results/ -minw "+str(a.minw7)+" -maxw "+str(a.maxw7)+" -p "+str(a.p7)+" -dna -mod anr -nmotifs 5 -maxsites "+str(a.N4)+" -revcomp -nostatus")

    #CREATING THE HTML-OUTPUT

    html_call = "cd "+a.outdir+" && createHTMLoutput.py PeakXus_results.html ./ "+str(a.N4)+" --umifigs 1"
    if a.matrixhits!=None: html_call += " --motifFigs 1"
    system(html_call)

    if a.verbosity!=None: print "PeakXus succesfully terminated!"

#end

def samtools(call_string,verbosity):
    
    if verbosity==1: print call_string
    system(call_string)
    return 1

def testinput(a):
    #this function test that the most crucial input parameters are sensible
    
    #1 Reference genome
    if a.genome2==None and a.fastq[-5:]=='fastq':
        print "Reference genome index not given! Terminating..."
        return False
         
    if a.wg7==None:
        print "** Reference genome file in fasta-format not given, some of the graphical output is not produced."
    else:
        try:
            f = open(a.wg7,'rb')
            f.close()
        except Exception:
            print e
            print "Reference genome file cannot be opened! Terminating..."
            return False

    #2 fastq/bam input
    try:
        f = open(a.fastq,'rb')
        f.close()

    except Exception as e:
        print e
        print "Input file cannot be opened! Terminating..."
        return False

    #3 UMI-file
    if a.UMIs!=None:
        try:
            aux = []
            with open(a.UMIs,'rb') as csvfile:
                r = csv.reader(csvfile,delimiter="\t")
                for row in r: aux = [row[0],row[1]]
        except Exception as e:
            print e
            print "Error with the UMI-file! Terminating..."
            return False
        
    #4 chromosome names-file
    try:
        aux = ""
        with open(a.chromnames,'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            for row in r: aux = row[0]

    except Exception as e:
        print e
        print "Problem with chromnames-file! Terminating..."
        return False

    #5 matrix hits-file
    if a.matrixhits!=None:
        try:
            f = open(a.matrixhits,'rb')
            f.close()
        except Exception as e:
            print e
            print "matrixhits-file could not be opened! Terminating..."
            return False
    else:
        print "** List of binding matrix hits to genome not provided, some of the graphical results are not printed."

    return True
        
            

PeakXus()
