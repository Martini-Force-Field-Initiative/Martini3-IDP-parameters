#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Liguo
"""

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib.pyplot as plt
import os


u = mda.Universe('../production-pbc.gro','../production-pbc.xtc')
ag = u.select_atoms('all')
print(len(ag.atoms))
os.system('mkdir cg_frame')

start = 0
end = -1
#every 2.5ns, to have 4000 frames
step = 5

for index, ts in enumerate(u.trajectory[start:end:step]): 
    ag.write('cg_frame/frame{}.gro'.format(index))
   

    


