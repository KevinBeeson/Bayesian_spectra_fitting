#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 10:33:08 2023

@author: kevin
"""
import numpy as np
x='Red'
limits={'Blue':[4705,4908],'Green':[5643,5879],'Red':[6470,6743],'IR':[7577.0,7894.0]}

for x in limits:
    tmp = np.load("NN_normalized_spectra_all_elements_2_"+x+".npz")
    w_array_0 = tmp["w_array_0"]
    w_array_1 = tmp["w_array_1"]
    w_array_2 = tmp["w_array_2"]
    b_array_0 = tmp["b_array_0"]
    b_array_1 = tmp["b_array_1"]
    b_array_2 = tmp["b_array_2"]
    x_min = tmp["x_min"]
    x_max = tmp["x_max"]
    tmp.close()
    NN_coeffs= (w_array_0, w_array_1, w_array_2, b_array_0, b_array_1, b_array_2, x_min, x_max)
    limits={'Blue':[4705,4908],'Green':[5643,5879],'Red':[6470,6743],'IR':[7577.0,7894.0]}
    spacings={'Blue':0.004,'Green':0.005,'Red':0.006,'IR':0.008}
    x_wave_old=np.linspace(limits[x][0],limits[x][1],num=int((limits[x][1]-limits[x][0])/0.004))
    x_wave_new=np.linspace(limits[x][0],limits[x][1],num=int((limits[x][1]-limits[x][0])/spacings[x]))
    
    b_array_2_new=np.interp(x_wave_new, x_wave_old, b_array_2)
    w_array_2_new=np.vstack([np.interp(x_wave_new,x_wave_old,x) for x in w_array_2.T ]).T
    np.savez("NN_normalized_spectra_all_elements_3_"+x+".npz",\
            w_array_0 = w_array_0,\
            w_array_1 = w_array_1,\
            w_array_2 = w_array_2_new,\
            b_array_0 = b_array_0,\
            b_array_1 = b_array_1,\
            b_array_2 = b_array_2_new,\
            x_max=x_max,\
            x_min=x_min,)