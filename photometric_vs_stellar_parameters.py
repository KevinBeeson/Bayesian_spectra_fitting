#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 11:28:41 2022

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

all_reduced_data=fits.getdata('galah_dr4_allspec_220713.fits',1)
all_reduced_data=Table(all_reduced_data)

cluster_name='upper_Scorpius'
votable = parse(cluster_name+"_bp_rp_run.xml")
file_directory = cluster_name+'_reduction/plots/'
Path(file_directory).mkdir(parents=True, exist_ok=True)
data=votable.get_first_table().to_table(use_names_over_ids=True)
# teff_temp=np.mean(data['teff_raw'],axis=1)
# data['teff']=teff_temp
# data['source_id']=np.array(data['source_id'],dtype='str')
cross=join(data,all_reduced_data,keys_left='source_id',keys_right='gaiadr3_source_id',table_names=['photometric', 'spectroscopic'],)
mask=cross['flag_sp']==0
cross=cross[mask]
# mask=cross['teff_photometric']-cross['teff_spectroscopic']<0
# cross=cross[mask]
# cross_smaller=vstack([x for x in cross if (x['teff_1']-x['teff_2']<1000 and x['teff_1']<5000) or x['teff_2']>5000 ])
mask=np.isnan(cross['teff_spectroscopic'])
cross=cross[~mask]

x=cross['teff_spectroscopic']
y=cross['teff_photometric']-cross['teff_spectroscopic']

polynomial_coeff=np.polyfit(x, y, 2)
x_line=np.linspace(min(x)-1000, max(x)+1000)

ynew=np.poly1d(polynomial_coeff)

fig=plt.figure()
plt.scatter(cross['teff_spectroscopic'],cross['teff_photometric']-cross['teff_spectroscopic'],c=cross['logg_spectroscopic'],s=5.0)
plt.plot(x_line,ynew(x_line))
plt.xlabel(r'Spectroscopic $T_{\rm{eff}}$/K')
plt.ylabel(r'Photometric -Spectroscopic $T_{\rm{eff}}$/K')
# plt.xlim((3300,11000))
# plt.ylim((-2600,3700))
# plt.xlim((3700,6500))
# plt.ylim((-1500,1000))

fig.set_size_inches(3.32088003321,3.32088003321/1.61)
plt.tight_layout()

plt.savefig('/home/kevin/Documents/Paper/'+cluster_name+' comparing spectroscopic vs photometric.pdf')

cross['coeff']=np.vstack([polynomial_coeff for x in range(len(cross))])
votable=from_table(cross)
writeto(votable,cluster_name+"_photometric_cross.xml")

plt.figure()
plt.scatter(cross['phot_g_mean_mag'],cross['bp_rp'],c=cross['teff_photometric'])
plt.gca().invert_yaxis()
plt.xlabel(r'BP-RP')
plt.ylabel(r'G')


plt.figure()
plt.scatter(cross['teff_photometric'],cross['logg_photometric'],c='blue')
plt.scatter(cross['teff_spectroscopic'],cross['logg_spectroscopic'],c='red')

plt.gca().invert_yaxis()
plt.gca().invert_xaxis()


# plt.figure()
# plt.scatter(cross['teff_1'],cross['teff_1']-cross['teff_2']-ynew(cross['teff_1']))
# plt.figure()
# plt.scatter(cross['logg'],cross['logg']-cross['logg_r'],c='black')
# plt.plot(x_line,ynew(x_line))
# plt.xlabel(r'Photometric $T_{\rm{eff}}$/K')
# plt.ylabel(r'Photometric -Spectroscopic $T_{\rm{eff}}$/K')

# plt.xlim((3300,11000))
# plt.ylim((-2600,3700))
# plt.savefig('Melotte_22 comparing spectroscopic vs photometric.pdf', type='pdf')