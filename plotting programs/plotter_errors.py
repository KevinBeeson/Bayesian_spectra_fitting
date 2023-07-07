#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 08:12:30 2022

@author: kevin
"""

import matplotlib.pyplot as plt
import numpy as np

data=np.load('NN_normalized_spectra_fe_alpha_shortest.npy',allow_pickle=True)
data2=np.load('NN_normalized_spectra_fe_alpha_shorter.npy',allow_pickle=True)
data3=np.load('NN_normalized_spectra_fe_alpha.npy',allow_pickle=True)
ends=[[4718, 4903], [5649, 5873], [6481, 6739],[7590, 7890]]
colours=['Blue','Green','Red','IR']
for value,x in enumerate(colours):
    average=np.mean(abs(data[value]),axis=0)
    average2=np.mean(abs(data2[value]),axis=0)
    average3=np.mean(abs(data3[value]),axis=0)

    step=(ends[value][1]-ends[value][0])/len(average)*10
    ends_temp=[ends[value][0]-step*80,ends[value][1]-step*80]        
    x_line=np.linspace(ends_temp[0],ends_temp[1],num=len(average))
    plt.figure()
    plt.plot(x_line,average,label='shortest')
    plt.plot(x_line,average2,label='medium')
    plt.plot(x_line,average3,label='largest')
    plt.legend(loc='best')
    plt.title(x)
    plt.yscale('log')
