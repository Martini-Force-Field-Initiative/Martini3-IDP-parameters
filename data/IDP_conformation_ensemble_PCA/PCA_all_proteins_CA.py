import pyemma
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import numpy as np
np.bool = np.bool_
import pickle as pkl
import sys

def save_pickle(filename, pickle_obj):
    with open(filename, 'wb') as f:
        pkl.dump(pickle_obj, f)

#Featurize pairwise CA distances based on gro-file
pdb_file = 'AA_backmapping/AA_CA_frame0.gro'
distances_feat = pyemma.coordinates.featurizer(pdb_file)
CA_atoms = distances_feat.select('name CA')
distances_feat.add_distances(CA_atoms , periodic=False)

#Calculate pairwise CA distances for each type of force field
traj_M3IDP = 'backward_CA_traj.xtc'
distances_data_M3IDP = pyemma.coordinates.load(traj_M3IDP, features=distances_feat)
frames_M3IDP=distances_data_M3IDP.shape[0]
print(frames_M3IDP)
distances_data_M3IDP = distances_data_M3IDP[0:]

traj_M3 = '../../../StandardM3/ACTR/backmapping/backward_CA_traj.xtc'
distances_data_M3 = pyemma.coordinates.load(traj_M3, features=distances_feat)
frames_M3=distances_data_M3.shape[0]
print(frames_M3)
distances_data_M3 = distances_data_M3[0:]

traj_ref1 = 'pnas2018b-ACTR-a99SBdisp-protein-merge-skip10.xtc'
distances_data_ref1 = pyemma.coordinates.load(traj_ref1, features=distances_feat)
frames_ref1=distances_data_ref1.shape[0]
print(frames_ref1)
distances_data_ref1 = distances_data_ref1[0:]

traj_ref2 = 'pnas2018b-ACTR-a03ws-protein-merge-skip10.xtc'
distances_data_ref2 = pyemma.coordinates.load(traj_ref2, features=distances_feat)
frames_ref2=distances_data_ref2.shape[0]
print(frames_ref2)
distances_data_ref2 = distances_data_ref2[0:]


#distances_all = np.concatenate((distances_data_M3IDP, distances_data_M3, distances_data_ref1, distances_data_ref2), axis=0)
#print(distances_all.shape)

#pca = pyemma.coordinates.pca(distances_all, dim=2)

#pca_output = pca.get_output()
#print(pca_output)
#save_pickle('PCA_output.pkl', pca_output)

pca_output = pkl.load( open( "./PCA_output.pkl", "rb" ))
pca_concatenated = np.concatenate(pca_output)

index1=frames_M3IDP
index2=frames_M3IDP+frames_M3
index3=frames_M3IDP+frames_M3+frames_ref1

fig, axes = plt.subplots(1, 4, figsize=(22, 4), sharex=True)
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=16)
tick_prop = font_manager.FontProperties(fname=font_path, size=14)
legend_prop = font_manager.FontProperties(fname=font_path, size=12)  
pyemma.plots.plot_density(*pca_concatenated[:index1, :2].T, ax=axes[0], cbar=False,nbins=40, logscale=True)
pyemma.plots.plot_density(*pca_concatenated[index1:index2, :2].T, ax=axes[1], cbar=False,nbins=40, logscale=True)
pyemma.plots.plot_density(*pca_concatenated[index2:index3, :2].T, ax=axes[2], cbar=False,nbins=40, logscale=True)
pyemma.plots.plot_density(*pca_concatenated[index3:, :2].T, ax=axes[3], cbar=True,nbins=40, logscale=True)
axes[0].set_title('Martini3-IDP Rg=3.05nm',fontproperties=font_prop)  
axes[1].set_title('Martini3 Rg=1.67nm',fontproperties=font_prop) 
axes[2].set_title('a99SBdisp Rg=2.13nm',fontproperties=font_prop) 
axes[3].set_title('a03ws Rg=2.22nm',fontproperties=font_prop) 

for ax in axes.flat[0:4]:
    ax.set_xlabel('PC 1',fontproperties=tick_prop)
    ax.set_ylabel('PC 2',fontproperties=tick_prop)

plt.savefig('PCA.png',dpi=600,bbox_inches='tight')
plt.show()
