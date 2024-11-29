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
    fig, ax = plt.subplots(figsize=(8,4.7))
    data1 = np.loadtxt(filename1)
    position_z = data1[:,0]
    factor=17168/20340
    densitys = data1[:,1:]*factor
    
    data2 = np.loadtxt('density_Water.xvg')
    density_Water = data2[:,1:]     
    
    cm = plt.cm.get_cmap('tab20')
    ax.plot(position_z,densitys[:,0],color=cm.colors[0],label='0-4 us',alpha=0.8)
    ax.plot(position_z,densitys[:,1],color=cm.colors[2],label='4-8 us',alpha=0.8)
    ax.plot(position_z,densitys[:,2],color=cm.colors[4],label='8-12 us',alpha=0.8)
    ax.plot(position_z,densitys[:,3],color=cm.colors[6],label='12-16 us',alpha=0.8)
    ax.plot(position_z,densitys[:,4],color=cm.colors[8],label='16-20 us',alpha=0.8)
    ax.plot(position_z,densitys[:,5],color=cm.colors[10],label='20-24 us',alpha=0.8)
    ax.plot(position_z,densitys[:,6],color=cm.colors[12],label='24-28 us',alpha=0.8)
    ax.plot(position_z,densitys[:,7],color=cm.colors[14],label='28-32 us',alpha=0.8)
    
    ax.plot(position_z,density_Water[:,0],color=cm.colors[0],alpha=0.8,linestyle='--')
    ax.plot(position_z,density_Water[:,1],color=cm.colors[2],alpha=0.8,linestyle='--')
    ax.plot(position_z,density_Water[:,2],color=cm.colors[4],alpha=0.8,linestyle='--')
    ax.plot(position_z,density_Water[:,3],color=cm.colors[6],alpha=0.8,linestyle='--')
    ax.plot(position_z,density_Water[:,4],color=cm.colors[8],alpha=0.8,linestyle='--')
    ax.plot(position_z,density_Water[:,5],color=cm.colors[10],alpha=0.8,linestyle='--')
    ax.plot(position_z,density_Water[:,6],color=cm.colors[12],alpha=0.8,linestyle='--')
    ax.plot(position_z,density_Water[:,7],color=cm.colors[14],alpha=0.8,linestyle='--')      
    
    ax.set_xlabel('Z Axis position (nm)',fontproperties=font_prop)
    ax.set_ylabel('Mass Concentration (mg/ml)',fontproperties=font_prop)
    ax.set_title('FUSLCD Density profle along Z Axis',fontproperties=font_prop)
    plt.ylim(0,1000)
    plt.legend(prop=legend_prop)
    #leg=plt.legend(loc=1, labelspacing=0.1, prop=leg_prop, scatterpoints=1, markerscale=1, numpoints=1,handlelength=1.5)
    #leg.get_frame().set_linewidth(0.0)
    #leg.get_frame().set_alpha(0.1)
    plt.xticks(fontproperties=tick_prop)
    plt.yticks(fontproperties=tick_prop) 
    plt.savefig(figout,dpi=600,bbox_inches='tight')
    plt.show()


plot_rg('density-collect-32000ns-interva4000ns.xvg','density.png')

