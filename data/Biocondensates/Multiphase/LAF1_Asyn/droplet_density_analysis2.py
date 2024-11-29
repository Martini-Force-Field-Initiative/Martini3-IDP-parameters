#! /usr/bin/python

import sys
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis import *
import math
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

traj=Universe('production_clean.gro','production_5500to10000ns.xtc')
#traj=Universe('production_clean.gro')
LAF1=traj.select_atoms('index 0-5069 ')
Asyn=traj.select_atoms('index 5070-8609 ')
condensate=traj.select_atoms('index 0-5069')

hist_Asyn_collect=np.zeros(49)
hist_LAF1_collect=np.zeros(49)
density_Asyn=np.zeros(49)
density_LAF1=np.zeros(49)
for index, ts in enumerate(traj.trajectory): 
    dist_Asyn = contacts.distance_array(condensate.center_of_geometry(), Asyn.positions)
    dist_LAF1 = contacts.distance_array(condensate.center_of_geometry(), LAF1.positions)
    #bindwidth=0.2nm
    hist_Asyn,edge=np.histogram(dist_Asyn,bins=np.arange(0,100,2))
    hist_LAF1,edge=np.histogram(dist_LAF1,bins=np.arange(0,100,2))
    hist_Asyn_collect=hist_Asyn_collect+hist_Asyn
    hist_LAF1_collect=hist_LAF1_collect+hist_LAF1
frames=len(traj.trajectory)
hist_Asyn_ave=hist_Asyn_collect/frames
hist_LAF1_ave=hist_LAF1_collect/frames

bins=[int((edge[i]+edge[i+1])/2) for i in range(49)]
print(bins)
bins=np.array(bins)
density_Asyn=[hist_Asyn_ave[i]/(4*np.pi*bins[i]**2) for i in range(49)]
density_LAF1=[hist_LAF1_ave[i]/(4*np.pi*bins[i]**2) for i in range(49)]

np.savetxt('droplet_density2_Asyn.xvg',hist_Asyn_ave,delimiter='	')
np.savetxt('droplet_density2_LAF1.xvg',hist_LAF1_ave,delimiter='	')

fig, ax = plt.subplots(figsize=(4.5,4.5))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
tick_prop = font_manager.FontProperties(fname=font_path, size=16)
legend_prop = font_manager.FontProperties(fname=font_path, size=15)  
ax.plot(bins[5:]/10,density_Asyn[5:],label='Asyn',alpha=1,color='green')
ax.plot(bins[5:]/10,density_LAF1[5:],label='LAF1',alpha=1,color='orange')
ax.set_ylabel('Radial Density',fontproperties=font_prop)
ax.set_xlabel('Radial Distance (nm)',fontproperties=font_prop)
plt.xticks(fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)
#plt.legend(prop=legend_prop)
plt.savefig('Droplet_density2.png',dpi=600,bbox_inches='tight')
plt.show()

    
    



