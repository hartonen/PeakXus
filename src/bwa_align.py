#!/usr/bin/env python
#
# THIS PROGRAM COMES WITH ABSOLUTELY NO WARRANTY!
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 2 as 
# published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

import os
import argparse

def bwa_align():

    
    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    
    #MANDATORY PARAMETERS
    parser.add_argument("infile",help="Full path to the input fastq-file.",type=str)
    parser.add_argument("hgfile",help="Full path to the reference human genome.",type=str)
    parser.add_argument("outdir",help="Full path to the output directory.",type=str)

    #OPTIONAL PARAMETERS
    parser.add_argument("-t","--numcores",help="Number of parallel processes used by bwa (default=1).",type=int,default=1)
    parser.add_argument("-b","--bamfile",help="1 (default) if BAM-file is created, 0 if not.",type=int,default=1)
    parser.add_argument("-B","--barcode",help="Length of the UMI-labels (default=0, no UMIs used).",type=int,default=0)
    parser.add_argument("-q","--quality",help="Quality limit for the alignments (default=20).",type=int,default=20)
    
    args = parser.parse_args()

    inname = args.infile.split('/')[-1][:-6]
    #Aligning the reads with bwa aln
    print "bwa aln -t "+str(args.numcores)+" -B "+str(args.barcode)+" -q "+str(args.quality)+" "+args.hgfile+" "+args.infile+" > "+args.outdir+inname+".bwa"
    os.system("bwa aln -t "+str(args.numcores)+" -B "+str(args.barcode)+" -q "+str(args.quality)+" "+args.hgfile+" "+args.infile+" > "+args.outdir+inname+".bwa")

    #Converting the bwa-file to a sam-file
    print "bwa samse "+args.hgfile+" "+args.outdir+inname+".bwa "+args.infile+" > "+args.outdir+inname+".sam"
    os.system("bwa samse "+args.hgfile+" "+args.outdir+inname+".bwa "+args.infile+" > "+args.outdir+inname+".sam")

    if args.bamfile==1:
        #Converting the sam-file to a sorted bam-file
        print "samtools view -bS -q "+str(args.quality)+" "+args.outdir+inname+".sam | samtools sort - "+args.outdir+inname+"_sorted"
        os.system("samtools view -bS -q "+str(args.quality)+" "+args.outdir+inname+".sam | samtools sort - "+args.outdir+inname+"_sorted")

        #indexing
        print "samtools index "+args.outdir+inname+"_sorted.bam"
        os.system("samtools index "+args.outdir+inname+"_sorted.bam")

        #removing unmapped reads
        print "samtools view -h -F 4 -b "+args.outdir+inname+"_sorted.bam > "+args.outdir+inname+"_sorted_filtered.bam"
        os.system("samtools view -h -F 4 -b "+args.outdir+inname+"_sorted.bam > "+args.outdir+inname+"_sorted_filtered.bam")
        
        print "samtools index "+args.outdir+inname+"_sorted_filtered.bam"
        os.system("samtools index "+args.outdir+inname+"_sorted_filtered.bam")

        #removing the .sam-file
        os.system("rm "+args.outdir+inname+".sam")
#end

bwa_align()
