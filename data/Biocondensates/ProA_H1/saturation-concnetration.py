#! /usr/bin/python

import sys
import os
import numpy as np

import pandas as pd


intervals=[12000,15000,18000,21000,24000]
factor=(20863*15+12203*18)/(26154*15+13986*18)
data=[]
with open('Dilute-phase.txt','w') as g:
    print('Begin\tEnd\tDiluteconcentration',file=g)
    for i in range(4):
        former=intervals[i]
        latter=intervals[i+1]
        print('density interval: {} to {}'.format(former,latter))

        filename='density-{}_{}ns.xvg'.format(former,latter)
        print(filename)
        os.system("sed -i 's/^@/#/g' %s " %filename)
        data_interval = np.loadtxt(filename)
        print(data_interval.shape)
        dilute1=data_interval[(data_interval[:,0]>24)]
        dilute2=data_interval[(data_interval[:,0]<-24)]
        print(dilute1)
        dilute_conc=(np.mean(dilute1[:,1])+np.mean(dilute2[:,1]))*factor/2
        print('{}\t{}\t{}'.format(former,latter,dilute_conc),file=g)



