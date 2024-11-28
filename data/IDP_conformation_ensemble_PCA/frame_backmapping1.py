#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Liguo
"""

import os
#os.system('mkdir ../AA_backmapping')
bashfile='/home/liguo/Downloads/backward-v5/backward-v5/initram-v5.sh'
for i in range(1000):
    os.system('bash ' + bashfile + ' -f ../cg_frame/frame{}.gro -p topol_ACTR.top -to charmm36 -o ../AA_backmapping/AA_frame{}.gro'.format(i,i))


    


