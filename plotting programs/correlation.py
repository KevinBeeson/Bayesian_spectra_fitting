#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 14:31:10 2022

@author: kevin
"""
import numpy as np 
import matplotlib.pyplot as plt 
tmp = np.load("NN_normalized_spectra_fe_alpha_shortest_Blue.npz")
w_array_0 = tmp["w_array_0"]
w_array_1 = tmp["w_array_1"]
w_array_2 = tmp["w_array_2"]
b_array_0 = tmp["b_array_0"]
b_array_1 = tmp["b_array_1"]
b_array_2 = tmp["b_array_2"]
x_min = tmp["x_min"]
x_max = tmp["x_max"]
tmp.close()

multiply=np.matmul(np.matmul(w_array_0.T,w_array_1.T),w_array_2.T)
correlation=np.zeros((8,8))
for i in range(8):
    for j in range(8):
        if i==j:
            correlation[i][j]=1.0
        else:
            correlation[i][j]=np.sum((multiply[i]-multiply[j]))/np.sqrt(np.sum((multiply[i]-multiply[j])**2))/len(multiply[i])
            # correlation[i][j]=np.sum((multiply[i]-multiply[j])/np.sqrt((multiply[i]-multiply[j])**2)/len(multiply[i]))