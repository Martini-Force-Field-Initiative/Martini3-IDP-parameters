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

stepsize=2.5 #ns
for i in range(4000):
    time=i*stepsize
    u1 = mda.Universe('AA_backmapping/AA_frame{}.gro'.format(i))
    u1_CA=u1.select_atoms('name CA')
    u1_CA.write('AA_backmapping/AA_CA_frame{}.gro'.format(i))
    os.system('gmx trjconv -f AA_backmapping/AA_CA_frame{}.gro -o AA_backmapping/AA_CA_frame{}.xtc -t0 {}'.format(i,i,time))
os.system('gmx trjcat -f AA_backmapping/AA_CA_frame*.xtc -o backward_CA_traj.xtc')
   

    


