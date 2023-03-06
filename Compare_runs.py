#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 09:20:07 2023

@author: kevin
"""

from numba import jit

import numpy as np
import matplotlib.pyplot as plt
import csv
from astropy.io.votable import parse,from_table,writeto
from scipy.stats import kde

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kde
from pathlib import Path
import requests
import re
from os.path import exists
import csv
from scipy.interpolate import interp1d
name='Melotte_22'
votable = parse(name+"_fixed_iron_final_2.xml")
data_slow=votable.get_first_table().to_table(use_names_over_ids=True)

votable = parse(name+"_fast.xml")
data_fast=votable.get_first_table().to_table(use_names_over_ids=True)

parameters_index=['teff','logg','fe_h','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']

parameters_slow=[]
parameters_fast=[]
for param in parameters_index:
    if param in elements:
        to_ask=param+'_Fe'
    else:
        to_ask=param
    parameters_slow.append(np.var(data_slow[to_ask+'_prior']))
    parameters_fast.append(np.var(data_fast[to_ask+'_prior']))
    
fig, ax = plt.subplots(1,1)
x=np.linspace(0,len(parameters_slow),num=len(parameters_slow))
plt.plot(np.array(parameters_slow)/np.array(parameters_fast),label='slow')
# plt.plot(parameters_fast,label='fast')
ax.set_xticks(x)
# Set ticks labels for x-axis
ax.set_xticklabels(parameters_index)
plt.legend(loc='best')