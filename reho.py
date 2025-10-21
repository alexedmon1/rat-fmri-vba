#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 14:50:06 2025

@author: edm9fd
"""

import numpy as np
from scipy.stats import rankdata
random_array = np.random.rand(100,100,27,150)

kernel = np.zeros(shape=(3,3,3)) # 27 voxel kernel
kernel_shape = (3,3,3)

center = [1,1,1]

def kendall_w(mat, n, kernel_size):
    center_rank_array = []
    num_obs = kernel_size[0]*kernel_size[1]*kernel_size[2]
    for i in range(n):
        new = mat[:,:,:,i]
        ranked_new = rankdata(new) #R(i) 
        middle = len(ranked_new)//2
        center_value = int(ranked_new[middle]) - int(np.mean(ranked_new))
        center_rank_array.append(center_value) # R
    W = (12 * np.sum(center_rank_array)) / ((n*2)*(num_obs**3 - num_obs)) 
    return W


def ReHo(file, kernel_size):
    final_image = np.zeros((file.shape[0],file.shape[1],file.shape[2]))
    for i in range(file.shape[0]+1):
        for j in range(file.shape[1]+1):
            for k in range(file.shape[2]+1):
                center = [i+1,j+1,k+1]
                selection = file[center[0]-1:center[0]+2,
                                     center[1]-1:center[1]+2,
                                     center[2]-1:center[2]+2,
                                     :]
                if selection.size == 0:
                    next
                else:
                    final_image[i,j,k] = kendall_w(selection, file.shape[3], kernel_size)
    return final_image

#%%
import matplotlib.pyplot as plt
import nibabel as nib

mri_data = nib.load('/mnt/arborea/Preclinical/BPA-rat/preproc/rest/cohort8/Rat122/p90/Rat122_rest_final.nii.gz')
data = mri_data.get_fdata()
reho = ReHo(data, (3,3,3))
