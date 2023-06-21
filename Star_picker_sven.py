#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 15:21:53 2023

@author: kevin
"""


from functools import  partial
from astropy.io.votable import parse
from astropy.io.votable import parse,from_table,writeto

# from The_Payne import spectral_model

from pathlib import Path

from astropy.table import Table,vstack,join
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
# import Galah_tool_py3 as gtools
# import random
all_reduced_data_dr61=fits.getdata('dr6.1.fits',1)
all_reduced_data_dr61=Table(all_reduced_data_dr61)
masks=all_reduced_data_dr61['reduction_flags']==0
all_reduced_data_dr61=all_reduced_data_dr61[masks]
all_reduced_data=fits.getdata('galah_dr4_allspec_220713.fits',1)

all_reduced_data=Table(all_reduced_data)

all_reduced_data=join(all_reduced_data,all_reduced_data_dr61,'sobject_id')

mask=all_reduced_data['flag_sp']==0 
all_reduced_data=all_reduced_data[mask]

mask=all_reduced_data['teff']<4700 
all_reduced_data=all_reduced_data[mask]
mask=all_reduced_data['teff']>4300
all_reduced_data=all_reduced_data[mask]
mask=all_reduced_data['logg']>4.40
all_reduced_data=all_reduced_data[mask]
mask=all_reduced_data['logg']<4.60
all_reduced_data=all_reduced_data[mask]

mask=all_reduced_data['fe_h']<-0.225
all_reduced_data=all_reduced_data[mask]
mask=all_reduced_data['fe_h']>-0.275
all_reduced_data=all_reduced_data[mask]

all_reduced_data.sort('e_fe_h')
#done as the first data peice of data is lacks the IR band
all_reduced_data=all_reduced_data[1:20]

votable=from_table(all_reduced_data)
writeto(votable,"svens_testing_stars.xml")
