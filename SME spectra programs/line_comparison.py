#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 15:11:40 2022

@author: kevin
"""

import numpy as np 
with open('Green_no_molecule.lin','r') as input_file:
    lines_long=input_file.readlines()
with open('Galah_Green_6.lin','r') as input_file:
    lines_short=input_file.readlines()
line_list_long=lines_long[2:]
final_line_list=[]
line_list_short=lines_short[3:]

short_line_list_values=[]
names_short=[]
values=[]
header=lines_long[:2]

for value,x in enumerate(line_list_short[::4]):
    if x.split(',')[0]=="'castelli_ap00k2_T05250G40.krz'":
        break
    name=x.split(',')[0]
    rest=x.split(',')[1:-1]
    rest=[round(float(y),4) for y in rest]
    short_line_list_values.append([name,rest])
short_line_list_values=np.array(short_line_list_values,dtype=object)    
line_list_without_short=[]
line_list_with=[]
for value,x in enumerate(line_list_long[::4]):
    name=x.split(',')[0]
    if name in short_line_list_values[:,0]:
        points=np.where(short_line_list_values[:,0]==name)
        rest=[y for y in short_line_list_values[points]]
        rest=[y[1] for y in rest]
        rest=np.array(rest)
        wavelength=rest[:,0:2]
        wave_to_find=[float(x.split(',')[1]),float(x.split(',')[2])]
        wave_to_find=[round(y,4) for y in wave_to_find]
        if not  wave_to_find in   wavelength:
            line_list_without_short.append(line_list_long[value*4:value*4+4])
        else:
            line_list_with.append(line_list_long[value*4:value*4+4])
    else:
        line_list_without_short.append(line_list_long[value*4:value*4+4])
        
big_string=''
refrences=''
ending=False
for y in line_list_long:
    if y[0]=='*':
        ending=True
    if ending:
        refrences+=y
for x in header:
    big_string+=x
for x in line_list_without_short:

    for y in x:
        big_string+=y
big_string+=refrences
with open('Green_no_over_lap.lin','w') as output:
    # for x in header:
    #     output.write(str(x))
    output.write(big_string)