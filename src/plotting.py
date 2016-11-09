# Copyright (c) Tuomo Hartonen, 2015-2016
#
# THIS PROGRAM COMES WITH ABSOLUTELY NO WARRANTY!
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 2 as 
# published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

import matplotlib
#Using the Agg-backend so that X-server is not needed
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm

def plotHeatMap(E,x,y,outname,color='seismic',title=None,xlabel=None,ylabel=None,xticks=None,xticklabels=None,yticks=None,yticklabels=None,figsize=None,resolution=900):
    #plots a heatmap using pcolormesh from matplotlib

    fig = plt.figure()
    gs = gridspec.GridSpec(100,100,top=0.94,bottom=0.11,left=0.04,right=0.9)
    ax = fig.add_subplot(gs[:,0:91])
    axC = fig.add_subplot(gs[:,93:99])

    #transforming the heatmap values for better visualization
    E_new = np.zeros(shape=(E.shape[0],E.shape[1]))
    ind = 0

    while ind<E.shape[0]:
        E_new[ind,:] = np.log(E[ind,:]+1)
        E_new[ind+1,:] = -1*np.log(abs(E[ind+1,:])+1)
        ind += 2
    quadmesh = ax.pcolormesh(x,y,E_new,cmap=cm.seismic,vmin=-1*np.abs(E_new).max(),vmax=np.abs(E_new).max())
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(quadmesh,ax=ax,cax=axC)
    if title!=None: ax.set_title(title,fontsize=26,weight='bold')
    if ylabel!=None: ax.set_ylabel(ylabel,fontsize=26,weight='bold')
    if xlabel!=None: ax.set_xlabel(xlabel,fontsize=26,weight='bold')
    if xticks!=None: ax.set_xticks(xticks)
    if yticks!=None: ax.set_yticks(yticks)
    else: ax.set_yticks([])
    if figsize==None: fig.set_size_inches((10,7.5))

    plt.savefig(outname,dpi=resolution)
    plt.clf()
    plt.close()

def plot1d(x,y,outname,legend=None,colors=['r','b','g','c','m','k'],title=None,xlabel=None,xticks=None,xticklabels=None,ylabel=None,yticks=None,yticklabels=None,figsize=None,resolution=150,opacity=1,m_left=None,m_len=None,errorbar=None,yscale=None,xscale=None,Nyticks=3,Nxticks=None,xlim=None,ylim=None,xrotation=None):
    #plots a normal 1d-plot using plot from matplotlib

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    for i in range(0,len(y)):
        if legend==None: ax.plot(x,y[i],color=colors[i],linewidth=4)
        else: ax.plot(x,y[i],color=colors[i],linewidth=4,label=legend[i])
    maxy = max([max(y[i]) for i in range(0,len(y))])
    #error bars
    if errorbar!=None:
        for i in range(0,len(errorbar)): ax.fill_between(x,y[i]-errorbar[i]/4.0,y[i]+errorbar[i]/4.0,alpha=0.2,edgecolor=colors[i],facecolor=colors[i])

    if opacity<1:
        x = []
        i = 0
        while i<(m_len-1):
            x.append(m_left+i)
            i += 1
        y = [maxy for i in x]
        ax.bar(x,y,width=1,color='g',alpha=opacity,edgecolor=None,linewidth=0)

    if title!=None: ax.set_title(title)
    if ylabel!=None: ax.set_ylabel(ylabel,fontsize=26,weight='bold')
    if xlabel!=None: ax.set_xlabel(xlabel,fontsize=26,weight='bold')
    if xticks!=None: ax.set_xticks(xticks)
    if yticks!=None: ax.set_yticks(yticks)
    if figsize!=None: fig.set_size_inches((8,6))
    if legend!=None: ax.legend(loc='upper right')
    if yscale!=None: ax.set_yscale(yscale)
    if xscale!=None: ax.set_xscale(xscale)
    if ylim!=None: ax.axis([min(x), max(x), 0, max([max(y[i]) for i in range(0,len(y))])])

    if Nyticks!=None: plt.locator_params(axis='y',nbins=Nyticks,tight=True,trim=False)
    if Nxticks!=None: plt.locator_params(axis='x',nbins=Nxticks,tight=True,trim=False)
    plt.tight_layout()

    plt.savefig(outname,dpi=resolution)
    plt.clf()
    plt.close()

def plotHistogram(data,bins,outname,color='r',legend=None,xlabel=None,ylabel=None,title=None,figsize=None,resolution=None):
    #Plots a histogram using bar from matplotlib

    hist,edges = np.histogram(data,bins=bins,density=False)
    plt.bar(edges[:-1],hist,color=color,width=edges[1]-edges[0],label=legend+"\nmean="+str(np.mean(np.array(data)))+"\nmedian="+str(np.median(np.array(data))))
    if xlabel!=None: plt.xlabel(xlabel)
    if ylabel!=None: plt.ylabel(ylabel)
    plt.margins(0.05)
    plt.subplots_adjust(bottom=0.08,top=0.87)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.17),ncol=1, fancybox=True, shadow=True)
    
    plt.savefig(outname)
    plt.clf()
    plt.close()
