#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 10:09:52 2025

@author: edm9fd
"""


import numpy as np
from scipy.stats import spearmanr


kernel = np.zeros(shape=(3,3,3)) # 27 voxel kernel
kernel_shape = (3,3,3)

center = [1,1,1]


def spearman(s):           
    center = (s.shape[0]//2,s.shape[1]//2,s.shape[2]//2)
    print(center)
    center_array = selection[center[0],center[1], center[2],:].reshape(selection.shape[3])
    cor_array = []
    for i in range(selection.shape[0]):
        for j in range(selection.shape[1]):
            for k in range(selection.shape[2]):
                selection_array = selection[i,j,k,:].reshape(selection.shape[3])
                res = spearmanr(center_array, selection_array)
                cor_array.append(float(res[0]))
    del cor_array[len(cor_array)//2]
    mean_correlation = np.mean(cor_array)
    return float(mean_correlation)

def localCorrelation(file, kernel_size):
    final_image = np.zeros((file.shape[0],file.shape[1],file.shape[2]))
    for i in range(file.shape[0]+1):
        for j in range(file.shape[1]+1):
            for k in range(file.shape[2]+1):
                center = [i+1,j+1,k+1]
                selection = file[center[0]-1:center[0]+2,
                                     center[1]-1:center[1]+2,
                                     center[2]-1:center[2]+2,
                                     :]
                if np.sum(selection) == 0:
                    next
                else:
                    print(selection)
                    final_image[i,j,k] = spearman(selection)
    return final_image

#%%

random_array = np.random.rand(100,100,10,150)
kernel_shape = (3,3,3)
center = [1,1,1]

final_image = np.zeros((random_array.shape[0],random_array.shape[1],random_array.shape[2]))
for i in range(random_array.shape[0]+1):
    for j in range(random_array.shape[1]+1):
        for k in range(random_array.shape[2]+1):
            center = [i,j,k]
            selection = random_array[center[0]-1:center[0]+2,
                                 center[1]-1:center[1]+2,
                                 center[2]-1:center[2]+2,
                                 :]
            if selection.size == 0:
                next
            else:
                final_image[i,j,k] = spearman(selection)
#%%
import matplotlib.pyplot as plt
import nibabel as nib

mri_data = nib.load('/mnt/arborea/Preclinical/BPA-rat/preproc/rest/cohort8/Rat122/p90/Rat122_rest_final.nii.gz')
data = mri_data.get_fdata()
reho = localCorrelation(data, (3,3,3))




#%%
