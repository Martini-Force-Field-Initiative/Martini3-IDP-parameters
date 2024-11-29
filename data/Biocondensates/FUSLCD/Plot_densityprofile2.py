#! /usr/bin/python

import sys
import os
import numpy as np
import argparse
import pandas as pd
import seaborn as sns
#sns.set(color_codes=True)
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager



def plot_rg(filename1,figout):
    font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
    font_prop = font_manager.FontProperties(fname=font_path, size=20)
    tick_prop = font_manager.FontProperties(fname=font_path, size=16)
    legend_prop = font_manager.FontProperties(fname=font_path, size=15)     
    os.system("sed -i 's/^@/#/g' %s " %filename1)
    fig, ax = plt.subplots(figsize=(8,3.8))
    data1 = np.loadtxt(filename1)
    position_z = data1[:,0]
    factor=17168/20340
    densitys_protein = data1[:,1]*factor

        
    cm = plt.cm.get_cmap('tab20')
    ax.plot(position_z,densitys_protein,color='green',alpha=0.8,linewidth=2.5,label='Protein')
    ax.plot(position_z,data1[:,2],color='grey',alpha=0.8,linewidth=2.5,label='Water')
    
    
    ax.set_xlabel('Z (nm)',fontproperties=font_prop)
    ax.set_ylabel('Density (mg/ml)',fontproperties=font_prop)
    #ax.set_title('FUSLCD Density profle along Z Axis',fontproperties=font_prop)
    plt.ylim(0,1050)
    plt.legend(prop=legend_prop)
    #leg=plt.legend(loc=1, labelspacing=0.1, prop=leg_prop, scatterpoints=1, markerscale=1, numpoints=1,handlelength=1.5)
    #leg.get_frame().set_linewidth(0.0)
    #leg.get_frame().set_alpha(0.1)
    plt.xticks(fontproperties=tick_prop)
    plt.yticks(fontproperties=tick_prop)
    plt.savefig(figout,dpi=600,bbox_inches='tight')
    plt.show()


plot_rg('density-8000_32000ns.xvg','density2.png')

