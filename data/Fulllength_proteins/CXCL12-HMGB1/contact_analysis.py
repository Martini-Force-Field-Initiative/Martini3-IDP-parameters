#! /usr/bin/python

import sys
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis import *
import math
import matplotlib.pyplot as plt
from MDAnalysis.analysis import contacts

def contacts_within_cutoff(u, group_a, group_b, radius=6):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries)

#clean not-protein atoms
#convert gro to pdb, add chain information to allow us select each residue in 2 monomers 
results_contact=[]

for i in range(214):
    res_A=i+1
    print('***HMGB1 residue{}, CXCL12***'.format(res_A))
    traj=Universe('production-pbc.pdb','production-pbc.xtc')
    resA=traj.select_atoms('resid {} and segid A'.format(res_A))
    resB=traj.select_atoms(' segid B')
    ca=contacts_within_cutoff(traj, resA, resB, radius=6)
    frames=ca.shape[0]
    print(frames)  
    ca_contact=ca[np.where(ca[:,1]>0)]
    ca_noncontact=ca[np.where(ca[:,1]==0)]
    print(ca_contact.shape[0])
    #print(ca_noncontact.shape[0])
    num_contacts=ca_contact.shape[0]
    occupancy_contacts=num_contacts/frames
    print(occupancy_contacts)
    results_contact.append([res_A,occupancy_contacts])
print(results_contact)
results_contact=np.array(results_contact)    
np.savetxt('contacts_occupancys_residue.xvg',results_contact,delimiter='	')


    
    



