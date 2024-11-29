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
    factor=(20863*15+12203*18)/(26154*15+13986*18)
    densitys = data1[:,1:]*factor

    data2 = np.loadtxt('density_Water.xvg')
    density_Water = data2[:,1:] 
    
    cm = plt.cm.get_cmap('tab20')
    #ax.plot(position_z,densitys[:,0],color=cm.colors[0],label='0-3 us',alpha=0.8)
    #ax.plot(position_z,densitys[:,1],color=cm.colors[2],label='3-6 us',alpha=0.8)
    ax.plot(position_z,(densitys[:,2]+densitys[:,1]+densitys[:,0])/3,color=cm.colors[0],label='0-9 us',alpha=0.8)
    ax.plot(position_z,densitys[:,3],color=cm.colors[2],label='9-12 us',alpha=0.8)
    ax.plot(position_z,densitys[:,4],color=cm.colors[4],label='12-15 us',alpha=0.8)
    ax.plot(position_z,densitys[:,5],color=cm.colors[6],label='15-18 us',alpha=0.8)
    ax.plot(position_z,densitys[:,6],color=cm.colors[12],label='18-21 us',alpha=0.8)
    ax.plot(position_z,densitys[:,7],color=cm.colors[14],label='21-24 us',alpha=0.8)
    
    #ax.plot(position_z,density_Water[:,0],color=cm.colors[0],alpha=0.8,linestyle='--')
    #ax.plot(position_z,density_Water[:,1],color=cm.colors[2],alpha=0.8,linestyle='--')
    ax.plot(position_z,(density_Water[:,2]+density_Water[:,1]+density_Water[:,0])/3,color=cm.colors[4],alpha=0.8,linestyle='--')
    ax.plot(position_z,density_Water[:,3],color=cm.colors[6],alpha=0.8,linestyle='--')
    ax.plot(position_z,density_Water[:,4],color=cm.colors[8],alpha=0.8,linestyle='--')
    ax.plot(position_z,density_Water[:,5],color=cm.colors[10],alpha=0.8,linestyle='--')
    ax.plot(position_z,density_Water[:,6],color=cm.colors[12],alpha=0.8,linestyle='--')
    ax.plot(position_z,density_Water[:,7],color=cm.colors[14],alpha=0.8,linestyle='--')      
    
    ax.set_xlabel('Z Axis position (nm)',fontproperties=font_prop)
    ax.set_ylabel('Mass Concentration (mg/ml)',fontproperties=font_prop)
    ax.set_title('ProA+H1 Density profle along Z Axis',fontproperties=font_prop)
    plt.ylim(0,1000)
    plt.legend(prop=legend_prop)
    plt.xticks(fontproperties=tick_prop)
    plt.yticks(fontproperties=tick_prop)
    plt.savefig(figout,dpi=600,bbox_inches='tight')
    plt.show()


plot_rg('density-collect-24000ns-interval3000ns.xvg','density.png')

