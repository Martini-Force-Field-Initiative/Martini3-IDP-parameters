#! /usr/bin/python
# -*- coding: utf-8 -*-
# @Date    : 2018-01-30 14:05:51
# @Author  : WDD 
# @Link    : https://github.com/dongdawn
# @Version : v1
import sys
import os
import numpy as np

import pandas as pd


intervals=[4000,8000,12000,16000,20000,24000,28000]
factor=2253/2700
data=[]
with open('Dilute-phase.txt','w') as g:
    print('Begin\tEnd\tDiluteconcentration',file=g)
    for i in range(6):
        former=intervals[i]
        latter=intervals[i+1]
        print('density interval: {} to {}'.format(former,latter))

        filename='density-{}_{}ns.xvg'.format(former,latter)
        print(filename)
        os.system("sed -i 's/^@/#/g' %s " %filename)
        data_interval = np.loadtxt(filename)
        print(data_interval.shape)
        dilute1=data_interval[(data_interval[:,0]>12)]
        dilute2=data_interval[(data_interval[:,0]<-12)]
        print(dilute1)
        dilute_conc=(np.mean(dilute1[:,1])+np.mean(dilute2[:,1]))*factor/2
        print('{}\t{}\t{}'.format(former,latter,dilute_conc),file=g)



