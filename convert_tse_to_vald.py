#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 16:09:42 2022

@author: kevin
"""
import copy
import numpy as np
line_list=np.loadtxt('atomic_lines_very_long.lin',dtype=str,delimiter='\t')
header=['element', 'wave_A', 'wave_nm', 'loggf', 'lower_state_eV',
       'lower_state_cm1', 'lower_j', 'upper_state_eV', 'upper_state_cm1',
       'upper_j', 'upper_g', 'lande_lower', 'lande_upper',
       'spectrum_transition_type', 'turbospectrum_rad', 'rad', 'stark',
       'waals', 'waals_single_gamma_format', 'turbospectrum_fdamp',
       'spectrum_fudge_factor', 'theoretical_depth', 'theoretical_ew',
       'lower_orbital_type', 'upper_orbital_type', 'molecule',
       'spectrum_synthe_isotope', 'ion', 'spectrum_moog_species',
       'turbospectrum_species', 'width_species', 'reference_code',
       'spectrum_support', 'turbospectrum_support', 'moog_support',
       'width_support', 'synthe_support', 'sme_support']
#needed=Spec Ion       WL_air(A)  log gf* E_low(eV) J lo E_up(eV)  J up  lower   upper   
# mean   Rad.   Stark  Waals  depth
limits=[[4700,4910],[5640,5880],[6474,6746],[7730,7900]]
vald_line_list=[x for x in line_list if x[-1]=='T']
blue_line_list=[x for x in vald_line_list if float(x[1])>7730 and float(x[1])<7900]
# blue_line_list=[x for x in blue_line_list if x[0]=='Fe 1']
# header_wanted=[' 4710.00000, 4908.00000,'+ len(blue_line_list)+', 1172314, 1.4 Wavelength region, lines selected, lines processed, Vmicro','Lande factors       Damping parameters  Central'
#         ,'Spec Ion       WL_air(A)  log gf* E_low(eV) J lo E_up(eV)  J up  lower   upper    mean   Rad.   Stark  Waals  depth']
index=[0,1,3,4,6,7,9,11,12,-1,15,16,17,21]
mean=[(float(x[11])+float(x[12]))/2 for  x in blue_line_list]
blue_line_list_mean=[np.append(x,str(round(y,3))) for x,y in zip(blue_line_list,mean)]
header_try=[header[x] for x in index]
vald_blue=[]
for y in blue_line_list_mean:
    line_temp=[y[x] for x in index]
    vald_blue.append(line_temp)      
    
ending="""'castelli_ap00k2_T05250G40.krz',
'H :  0.92','He: -1.11',
'Li:-11.02','Be:-10.72','B : -9.57','C : -3.60','N : -4.20','O : -3.29',
'F : -7.56','Ne: -4.04','Na: -5.79','Mg: -4.54','Al: -5.65','Si: -4.57',
'P : -6.67','S : -4.79','Cl: -6.62','Ar: -5.72','K : -7.00','Ca: -5.76',
'Sc: -8.95','Ti: -7.10','V : -8.12','Cr: -6.45','Mn: -6.73','Fe: -4.62',
'Co: -7.20','Ni: -5.87','Cu: -7.91','Zn: -7.52','Ga: -9.24','Ge: -8.71',
'As: -9.75','Se: -8.71','Br: -9.49','Kr: -8.81','Rb: -9.52','Sr: -9.15',
'Y : -9.88','Zr: -9.52','Nb:-10.70','Mo:-10.20','Tc:-20.08','Ru:-10.28',
'Rh:-11.00','Pd:-10.43','Ag:-11.18','Cd:-10.35','In:-10.46','Sn:-10.12',
'Sb:-11.12','Te: -9.88','I :-10.61','Xe: -9.95','Cs:-10.99','Ba: -9.99',
'La:-10.95','Ce:-10.54','Pr:-11.41','Nd:-10.62','Pm:-20.08','Sm:-11.11',
'Eu:-11.61','Gd:-11.00','Tb:-11.77','Dy:-10.98','Ho:-11.86','Er:-11.19',
'Tm:-12.12','Yb:-11.04','Lu:-12.06','Hf:-11.24','Ta:-12.25','W :-11.01',
'Re:-11.84','Os:-10.67','Ir:-10.77','Pt:-10.32','Au:-11.11','Hg:-10.99',
'Tl:-11.22','Pb:-10.17','Bi:-11.41','Po:-20.08','At:-20.08','Rn:-20.08',
'Fr:-20.08','Ra:-20.08','Ac:-20.08','Th:-12.03','Pa:-20.08','U :-12.62',
'Np:-20.08','Pu:-20.08','Am:-20.08','Cm:-20.08','Bk:-20.08','Cf:-20.08',
'Es:-20.08','END'"""
header_blue=[' 7730.00000, 7900.00000, '+ str(len(vald_blue))+', 1172314, 1.4 Wavelength region, lines selected, lines processed, Vmicro','                                                                     Lande factors       Damping parameters  Central'
                      ,'Spec Ion       WL_air(A)  log gf* E_low(eV) J lo E_up(eV)  J up  lower   upper    mean   Rad.   Stark  Waals  depth']
# for y in vald_blue:
    # y[0]='\''+y[0]+'\''
for y in vald_blue:
    if  not y[0][-1].isdigit():
        y[0]+='1'
        # print(y)
with open('IR_no_NLTE.lin', 'w') as output:
    for x in header_blue:
        output.write(str(x))
        output.write(str('\n'))
    for y in vald_blue:
        line='\''+y[0]+'\''
        line+=',       '+', '.join(y[1:])        
        output.write(line)
        output.write(str('\n'))
        output.write(str('\n'))
        output.write(str('\n'))
        output.write(str('\n'))
    output.write(ending)