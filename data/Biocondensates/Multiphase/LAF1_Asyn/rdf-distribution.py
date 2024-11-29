#! /usr/bin/python
# -*- coding: utf-8 -*-
# calculate gyrate distribution in different trajectory interval to analyze convergence and conformation ensemble
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
from statistics import mean


filename1='rdf.xvg'
filename2='rdf2.xvg'
os.system("sed -i 's/^@/#/g' %s " %filename1)
os.system("sed -i 's/^@/#/g' %s " %filename2)
fig, ax = plt.subplots(figsize=(8,4.7))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
tick_prop = font_manager.FontProperties(fname=font_path, size=16)
legend_prop = font_manager.FontProperties(fname=font_path, size=15)    
data1 = np.loadtxt(filename1)
rdf1 = data1[:,1]
print(mean(rdf1))
data2 = np.loadtxt(filename2)
rdf2 = data2[:,1]
print(mean(rdf2))
ax.plot(rdf1,label='LAF1-LAF1',alpha=1,color='#3d5a9c',linewidth=0.8)
ax.plot(rdf2,label='LAF1-Asyn',alpha=1,color='#b03d56',linewidth=0.8)

#ax.spines[['right', 'top']].set_visible(False)      
ax.set_ylabel('g(r)',fontproperties=font_prop)
ax.set_xlabel('Distance (nm)',fontproperties=font_prop)
plt.xlim(0,4000)
ax.set_xticks([0,1000,2000,3000,4000])
ax.set_xticklabels([0,2,4,6,8])
plt.xticks(fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)    
plt.legend(prop=legend_prop)
plt.savefig('RDF.png',dpi=600,bbox_inches='tight')
plt.show()

