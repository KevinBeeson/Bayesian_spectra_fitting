#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 16:30:52 2022

@author: kevin
"""
import numpy as np
import matplotlib.pyplot as plt
data_short=np.load('short_line_spectras.npy',allow_pickle=True)
data_long=np.load('long_line_spectras.npy',allow_pickle=True)
spectras_short=np.array([x[0][0] for x in data_short])
spectras_long=np.array([x[0][0] for x in data_long])

labels=np.array([x[1] for x in data_short])
diff=spectras_long-spectras_short
wave_length=np.linspace(5649, 5873,num=42560)
plt.figure()
plt.plot(wave_length,np.max(abs(diff),axis=0))
plt.plot(wave_length,np.mean(abs(spectras_long)-1,axis=0))

headers=['teff','logg','monh','non_sense','alpha','vmac','vmic','vsini']
for x,x_label in zip(range(7),headers):
    if x!=3:
        plt.figure()
        plt.scatter(labels[:,x],np.max(diff,axis=1),s=0.5)
        plt.xlabel(x_label)
        plt.ylabel('flux_maximum_difference')