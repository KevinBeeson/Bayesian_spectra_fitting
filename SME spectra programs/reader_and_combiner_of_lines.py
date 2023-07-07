#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 11:15:43 2022

@author: kevin
"""

import numpy as np 
import os
def wave_length_reader(name):
    data=[]
    header=[]
    ending=[]
    with open('line_list/'+name) as f:
        data_temp=[]
        end=False
        for value,line in enumerate(f):
            if line[:9]=="'castelli":
                end=True
            if value<=2:
                header.append(line)
            if value>2 and not end:
                data_temp.append(line)
            # print(line, end='')
            if end:
                ending.append(line)
            if len(data_temp)==4:
                data.append(data_temp)
                data_temp=[]
    data_wavelength=[]
    for x in data:
        data_wavelength.append(float(x[0].split(',')[1]))

    return data,data_wavelength,header,ending
def data_combiner(data_1,wavelength_1,data_2,wavelength_2):
    indentifier=[x[0][:-7] for x in data_1 ]
    data_unique_all=[(x,y) for x,y in zip(data_2,wavelength_2) if not x[0][:-7] in indentifier]
    data_unique=[x[0] for x in data_unique_all]
    wavelength_unique=[x[1] for x in data_unique_all]
    if len(data_unique):
        data_all=np.vstack((data_1,data_unique))
        wave_length_all=np.hstack((wavelength_1,wavelength_unique))
    else:
        data_all=data_1
        wave_length_all=wavelength_1
    arr1inds = wave_length_all.argsort()
    data_all = data_all[arr1inds]
    wave_length_all=wave_length_all[arr1inds]
    return data_all, wave_length_all
band='IR'
names_all=os.listdir('line_list/'+band+'/')

data_all,wavelength_all,header,ending=wave_length_reader(band+'/'+names_all[0])
print(len(data_all))
old_len=len(data_all)
for name in names_all[1:]:
    
    data2,wave_length_range2,heade,ending2r=wave_length_reader(band+'/'+name)

    data_all,wavelength_all=data_combiner(data_all,wavelength_all,data2,wave_length_range2)
    print('added '+str( len(data_all)-old_len))
    old_len=len(data_all)
    location_of_commas=[i for i, c in enumerate(header[0]) if c == ',']
    header[0]=header[0][:location_of_commas[1]+2]+str(len(data_all))+header[0][location_of_commas[2]:]

with open('line_list/'+band+'/combined_line_list_'+band+'.txt','w+') as f:
    for line in header:
        f.write(line)
    for line in data_all:
        for y in line:
            f.write(y)
    for line in ending:
        f.write(line)