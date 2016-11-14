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

def createHTMLoutput():

    ########################
    #Command line arguments#
    ########################

    parser = argparse.ArgumentParser()

    #MANDATORY PARAMETERS
    parser.add_argument("outfile",help="Full path to the output HTML-file.",type=str)
    parser.add_argument("indir",help="Full path to the directory containing the PeakXus-results.",type=str)
    parser.add_argument("N",help="Total number of peaks analyzed.",type=int)

    #OPTIONAL PARAMETERS
    parser.add_argument("--motifFigs",help="If 1, also the motif-match images were produced by PeakXus and are embedded to the html-page, otherwise 0 (default).",type=int,choices=[0,1],default=0)
    parser.add_argument("--umifigs",help="If 1, the UMI-histograms were produced by PeakXus and are embedded to the html-page, otherwise 0 (default).",type=int,choices=[0,1],default=0)

    args = parser.parse_args()

    f = open(args.outfile,'wb')

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
    f.write("<h1>PEAKXUS RESULTS</h1>")
    f.write('</div>')

    ########################
    #Summary of the results#
    ########################
    #reading in the summary results from PeakXus output
    totreads = 0.0
    totumis = 0.0
    totpeaks = 0.0
    with open(args.indir+"PeakXus_summary.txt",'rb') as csvfile:
        r = csv.reader(csvfile,delimiter=" ")
        for row in r:
            if len(row)<4: continue
            if row[0].count("Total")>0:
                if row[3].count("reads")>0: totreads = float(row[-1])
                elif row[3].count("UMIs")>0: totumis = float(row[-1])
                elif row[3].count("peaks")>0: totpeaks = float(row[-1])

    f.write('<div style="background-color:#ccccff; color:white; padding:20px;">')
    f.write("<h2>SUMMARY</h2>")
    f.write('</div>')
    
    f.write('<table style="width:100%">')
    f.write('<tr>')
    f.write("<td><p>Total number of peaks (no FDR-cutoff)</p></td>")
    f.write("<td><p>=</p></td>")
    f.write("<td><p>"+str(totpeaks)+"</p></td>")
    f.write("</tr>")
    f.write("<tr>")
    f.write("<td><p>Total number of reads</p></td>")
    f.write("<td><p>=</p></td>")
    f.write("<td><p>"+str(totreads)+"</p></td>")
    f.write("</tr>")
    f.write("<tr>")
    f.write("<td><p>Total number of UMIs</p></td>")
    f.write("<td><p>=</p></td>")
    f.write("<td><p>"+str(totumis)+"</p></td>")
    f.write("</tr>")
    f.write("<tr>")
    f.write("<td><p>Reads/UMIs</p></td>")
    f.write("<td><p>=</p></td>")
    f.write("<td><p>"+str(totreads/totumis)+"</p></td>")
    f.write("</tr></table>")

    fastqcname = glob("FastQC/*.html")
    if len(fastqcname)==0: fastqcname = "none.html"
    else: fastqcname = fastqcname[0]
    f.write('<p><a href="'+fastqcname+'">Quality control results from FastQC for the input library.</a></p>')

    f.write('<table style="width:100%">')
    f.write("<tr>")
    f.write("<td><p>The best MEME-motif:</p></td>")
    f.write('<td><img src="meme_results/logo1.png" alt="Motif 1"></td>') 
    f.write('<td><img src="meme_results/logo_rc1.png" alt="Motif 1 rc"></td>')
    f.write("</tr>")
    f.write("</table>")

    f.write('<p><a href="meme_results/meme.html">Full MEME-results</a></p>')
    
    ######################
    #Peak size statistics#
    ######################
    #This is always produced
    f.write('<div style="background-color:#ccccff; color:white; padding:20px;">')
    f.write('<h2>PEAK SIZE STATISTICS</h2>')
    f.write('</div>')
    f.write('<table style="width:100%">')
    f.write("<tr>")
    f.write('<td><div class="image"><img src="peak_score_histo.png" alt="Peak score histogram" style="width:800px;height:600px;"></div></td>') 
    f.write('<td><div class="image"><img src="peak_size_histo.png" alt="Peak size histogram" style="width:800px;height:600px;"></div></td>')
    f.write("</tr>")
    f.write('<tr>')
    f.write('<td><div class="caption">Histogram of peak scores for all reported peaks.</div></td>')
    f.write('<td><div class="caption">Histogram of unique read 5\'-end counts within peak borders (Reads pointing towards peak summit minus reads pointing away from peak summit). If UMIs were not used in the analysis, raw read counts are reported.</div></td>')
    f.write('</tr>')
    f.write("</table>")

    ######################
    #UMI-count histograms#
    ######################
    #This is optional and works only if UMI-histograms are created by PeakXus

    if args.umifigs==1:
        f.write('<div style="background-color:#ccccff; color:white; padding:20px;">')
        f.write('<h2>HISTOGRAM OF UMI-COUNTS</h2>')
        f.write('</div>')
        f.write('<table style="width:100%">')
        f.write("<tr>")
        f.write('<td><img src="UMI_histo_+.png" alt="+ strand UMI-count histogram" style="width:800px;height:600px;"></td>') 
        f.write('<td><img src="UMI_histo_-.png" alt="- strand UMI-count histogram" style="width:800px;height:600px;"></td>')
        f.write("</tr>")
        f.write('<tr>')
        f.write('<td><div class="caption">Histogram of UMI-counts per position on the sense strand.</div></td>')
        f.write('<td><div class="caption">Histogram of UMI-counts per position on the antisense strand.</div></td>')
        f.write('</tr>')
        f.write("</table>")

    #########################
    #Images of the top peaks#
    #########################
    #Images of top N peaks are always produced
    
    f.write('<div style="background-color:#ccccff; color:white; padding:20px;">')
    f.write('<h2>BREAKDOWN OF TOP '+str(args.N)+' PEAKS</h2>')
    f.write('</div>')
    f.write('<table style="width:100%">')
    f.write("<tr>")
    f.write('<td><img src="Avg_top_'+str(args.N)+'.png" alt="Top 1000 peaks" style="width:800px;height:600px;"></td>') 
    f.write('<td><img src="Heatmap_top_'+str(args.N)+'.png" alt="Top 1000 peaks" style="width:800px;height:600px;"></td>')
    f.write("</tr>")
    f.write('<tr>')
    f.write('<td><div class="caption">Average read 5\'-end count around the top '+str(args.N)+' peaks. Red curve corresponds to read 5\'-end counts on the sense strand, and blue on the antisense strand, respectively.</div></td>')
    f.write('<td><div class="caption">Heatmap presentation of the same peaks shown on the left. Red color marks the sense strand read 5\'-end count, and blue the antisense strand, respectively.</div></td>')
    f.write('</tr>')
    if args.motifFigs==1:
        f.write("<tr>")
        f.write('<td><img src="Avg_motif_match_top_oriented_'+str(args.N)+'.png" alt="Top 1000 peaks with a motif match" style="width:800px;height:600px;"></td>') 
        f.write('<td><img src="Heatmap_motif_match_oriented_top'+str(args.N)+'.png" alt="Top 1000 peaks with a motif match" style="width:800px;height:600px;"></td>')
        f.write("</tr>")
        f.write('<tr>')
        f.write('<td><div class="caption">Average read 5\'-end count around the top '+str(args.N)+' peaks with an underlying High-Affinity Recognition Sequence (HARS) hit. The read counts are oriented according to the HARS-match so that the red peak corresponds to the 5\'-side of the motif. The green shading indicates the width of the HARS. The red and blue shadings correspond to half of the standard deviation over the top '+str(args.N)+' peaks.</div></td>')
        f.write('<td><div class="caption">Heatmap presentation of the same peaks shown on the left. The read counts are oriented according to the HARS-match so that the red peak corresponds to the 5\'-side of the motif.</div></td>')
        f.write('</tr>')
    f.write("</table>")

    #####################
    #Images of all peaks#
    #####################
    #These are optional and require a list of PWM-hits as input for PeakXus
    if args.motifFigs==1:
        motifsperpeakname = glob("motifsPerPeak*.png")[0]
        distfrommotifname = glob("dist_from_motif*.png")[0]
        f.write('<div style="background-color:#ccccff; color:white; padding:20px;">')
        f.write('<h2>PEAK POSITIONING RELATIVE TO PWM-HITS (HARS-SITES)</h2>')
        f.write('</div>')
        f.write('<table style="width:100%">')
        f.write("<tr>")
        f.write('<td><img src="'+motifsperpeakname+'" alt="PWM-hits per peak" style="width:800px;height:600px;"></td>') 
        f.write('<td><img src="'+distfrommotifname+'" alt="Distance between peaks and PWM-hits" style="width:800px;height:600px;"></td>')
        f.write("</tr>")
        f.write('<tr>')
        f.write('<td><div class="caption">Number of HARS-sites overlapping with a peak as a function of the number of peaks analyzed. Peaks are ranked according to peak score.</div></td>')
        f.write('<td><div class="caption">Cumulative fraction of the top '+str(args.N)+' peaks that are within a distance indicated on the x-axis from the nearest HARS. The topmost subfigure shows the distances between left border of the peak and left edge of the HARS, the subfigure in the middle the distances between the summit of the peak and the center of the HARS and the bottom subfigure the distance between the right border of the peak and the right edge of the HARS.</div></td>')
        f.write('</tr>')
        f.write("</table>")

    ##################
    #Input parameters#
    ##################

    f.write('<div style="background-color:#ccccff; color:white; padding:20px;">')
    f.write('<h3>INPUT PARAMETERS FOR PEAKXUS</h3>')
    f.write('</div>')
    arguments = ""
    g = open(args.indir+"/log.txt",'rb')
    for row in g:
        if row.count('Namespace(')>0:
            arguments = row[10:-1]
    arguments = arguments.split(', ')
    f.write('<p class="small">\n')
    for a in arguments: f.write(a+'<br>\n')
    f.write('</p>')


    f.write('<div style="background-color:#ccccff; color:white; padding:20px;">')
    f.write('<p class="white">If you use PeakXus in your work, please cite:<br>Hartonen, T., Sahu, B., Dave, K., Kivioja, T., & Taipale, J. (2016). PeakXus: comprehensive transcription factor binding site discovery from ChIP-Nexus and ChIP-Exo experiments. Bioinformatics, 32(17), i629-i638.</p>')
    f.write('</div>')

    #end
    f.write("</body>\n</html>")
    f.close()

#end

createHTMLoutput()
