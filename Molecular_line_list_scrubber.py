#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 15:36:38 2022

@author: kevin
"""

import numpy as np 
solar_abundances={'H': 12.0, 'He': 10.93, 'Li': 1.05, 'Be': 1.38, 'B': 2.7, 'C': 8.43, 'N': 7.83, 'O': 8.69, 
                  'F': 4.56, 'Ne': 7.93, 'Na': 6.24, 'Mg': 7.6, 'Al': 6.45, 'Si': 7.51, 'P': 5.41, 'S': 7.12, 'Cl': 5.5, 'Ar': 6.4, 
                  'K': 5.03, 'Ca': 6.34, 'Sc': 3.15, 'Ti': 4.95, 'V': 3.93, 'Cr': 5.64, 'Mn': 5.43, 'Fe': 7.5, 'Co': 4.99, 'Ni': 6.22, 
                  'Cu': 4.19, 'Zn': 4.56, 'Ga': 3.04, 'Ge': 3.65, 'As': 2.3, 'Se': 3.34, 'Br': 2.54, 'Kr': 3.25, 'Rb': 2.52, 'Sr': 2.87, 
                  'Y': 2.21, 'Zr': 2.58, 'Nb': 1.46, 'Mo': 1.88, 'Tc': np.nan, 'Ru': 1.75, 'Rh': 0.91, 'Pd': 1.57, 'Ag': 0.94, 'Cd': 1.71, 
                  'In': 0.8, 'Sn': 2.04, 'Sb': 1.01, 'Te': 2.18, 'I': 1.55, 'Xe': 2.24, 'Cs': 1.08, 'Ba': 2.18, 'La': 1.1, 'Ce': 1.58, 
                  'Pr': 0.72, 'Nd': 1.42, 'Pm': np.nan, 'Sm': 0.96, 'Eu': 0.52, 'Gd': 1.07, 'Tb': 0.3, 'Dy': 1.1, 'Ho': 0.48, 'Er': 0.92, 
                  'Tm': 0.1, 'Yb': 0.84, 'Lu': 0.1, 'Hf': 0.85, 'Ta': -0.12, 'W': 0.85, 'Re': 0.26, 'Os': 1.4, 'Ir': 1.38, 'Pt': 1.62, 
                  'Au': 0.92, 'Hg': 1.17, 'Tl': 0.9, 'Pb': 1.75, 'Bi': 0.65, 'Po': np.nan, 'At': np.nan, 'Rn': np.nan, 'Fr': np.nan,
                  'Ra': np.nan, 'Ac': np.nan, 'Th': 0.02, 'Pa': np.nan, 'U': -0.54, 'Np': np.nan, 'Pu': np.nan, 'Am': np.nan, 
                  'Cm': np.nan, 'Bk': np.nan, 'Cf': np.nan, 'Es': np.nan}

with open('Galah_IR_full.lin','r') as input_file:
    lines=input_file.readlines()
line_list=lines[2:]
final_line_list=[]
header=lines[:2]

for value,x in enumerate(line_list[::4]):
    name=x.split(',')[0]
    if name[0]!='*':
        charge=int(name.split(' ')[-1][0])
        name=name.split(' ')[0]
        # name=name.replace("'","")
        # if  name in solar_abundances:
        #     if solar_abundances[name]!=np.nan:
        if charge<5:
            final_line_list.append(line_list[value*4:value*4+4])
    else:
        break
big_string=''
refrences=''
ending=False
for y in lines:
    if y[0]=='*':
        ending=True
    if ending:
        refrences+=y
for x in header:
    big_string+=x
for x in final_line_list:

    for y in x:
        big_string+=y
big_string+=refrences
with open('IR_no_molecule.lin','w') as output:
    # for x in header:
    #     output.write(str(x))
    output.write(big_string)