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
%matplotlib ipympl
import numpy as np
import csv
from astropy.io import fits
# import Galah_tool_py3 as gtools
# import random
svens_reduced_data=fits.getdata('galah_dr4_allspec_220713.fits',1)
svens_reduced_data=Table(svens_reduced_data)

all_reduced_data=fits.getdata('dr6.1.fits',1)
all_reduced_data=Table(all_reduced_data)
mask=all_reduced_data['gaia_id']!='None'
all_reduced_data=all_reduced_data[mask]
all_reduced_data['gaia_id']=all_reduced_data['gaia_id'].astype(np.int64)

all_reduced_data_old=fits.getdata('dr6.0.fits',1)
all_reduced_data_old=Table(all_reduced_data_old)
mask=all_reduced_data_old['gaia_id']!='None'
all_reduced_data_old=all_reduced_data_old[mask]
all_reduced_data_old['gaia_id']=all_reduced_data_old['gaia_id'].astype(np.int64)



#stacks all of bp_rp 

cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)
print(cluster_details_all[:,0])
names=cluster_details_all[:,0]
# names[0]=cluster_name
votable = parse(names[0]+"_bp_rp_run.xml")
data_all=votable.get_first_table().to_table(use_names_over_ids=True)
data_all['cluster_name']=[names[0] for x in range(len(data_all))]
plt.figure()
plt.scatter(data_all['phot_g_mean_mag'],data_all['bp_rp'],c=data_all['teff'])
plt.title(names[0])
for name in names[1:]:
    print(name)
    votable = parse(name+"_bp_rp_run.xml")
    data=votable.get_first_table().to_table(use_names_over_ids=True)
    data['cluster_name']=[name for x in range(len(data))]
    data_all=vstack((data_all,data))
    # plt.figure()
    # plt.scatter(data['phot_g_mean_mag'],data['bp_rp'],c=data['teff'])
    # plt.title(name)

cross_sven=join(data_all,svens_reduced_data,keys_left='source_id',keys_right='gaiadr3_source_id',table_names=['photometric', 'spectroscopic'])

cross_reduced=join(data_all,all_reduced_data,keys_left='source_id',keys_right='gaia_id',table_names=['photometric', 'spectroscopic'])
cross_reduced_old=join(data_all,all_reduced_data_old,keys_left='source_id',keys_right='gaia_id',table_names=['photometric', 'spectroscopic'])
nights_reduced=cross_reduced['sobject_id']
nights_to_reduce=[x for x in cross_sven['sobject_id'] if x not in nights_reduced]
nights_to_reduce=np.unique(nights_to_reduce)

unique=join(cross_sven,cross_reduced_old,keys='sobject_id',join_type='outer')
mask=unique['pipeline_version']!='6.0'
unique=unique[mask]

mask=cross_sven['flag_sp']==0
cross_masked=cross_sven[mask]
mask=np.isnan(cross_masked['teff_spectroscopic'])
cross_masked=cross_masked[~mask]

x=cross_masked['teff_spectroscopic']
y=cross_masked['teff_photometric']-cross_masked['teff_spectroscopic']

polynomial_coeff=np.polyfit(x, y, 2)
x_line=np.linspace(min(x)-1000, max(x)+1000)

ynew=np.poly1d(polynomial_coeff)

fig=plt.figure()
plt.scatter(cross_masked['teff_spectroscopic'],cross_masked['teff_photometric']-cross_masked['teff_spectroscopic'],s=1.0,c='black')
plt.plot(x_line,ynew(x_line))
plt.xlabel(r'Spectroscopic $T_{\rm{eff}}$/K')
plt.ylabel(r'$\Delta T_{\rm{eff}}$/K')
plt.xlim((cross_masked['teff_spectroscopic'].min()-100,cross_masked['teff_spectroscopic'].max()+100))
plt.ylim((cross_masked['teff_photometric']-cross_masked['teff_spectroscopic']).min()-100,(cross_masked['teff_photometric']-cross_masked['teff_spectroscopic']).max()+100)
# plt.xlim((3300,11000))
# plt.ylim((-2600,3700))
# plt.xlim((3700,6500))
# plt.ylim((-1500,1000))

fig.set_size_inches(3.32088003321,3.32088003321/1.61)
plt.tight_layout()

plt.savefig('/home/kevin/Documents/Paper/all comparing spectroscopic vs photometric.pdf')

cross_sven['coeff']=np.vstack([polynomial_coeff for x in range(len(cross_sven))])
votable=from_table(cross_sven)
writeto(votable,"open_cluster_photometric_cross.xml")

plt.figure()
plt.scatter(cross_sven['phot_g_mean_mag'],cross_sven['bp_rp'],c=cross_sven['teff_photometric'])
plt.gca().invert_yaxis()
plt.xlabel(r'BP-RP')
plt.ylabel(r'G')


plt.figure()
plt.scatter(cross_sven['teff_photometric'],cross_sven['logg_photometric'],c='blue')
plt.scatter(cross_sven['teff_spectroscopic'],cross_sven['logg_spectroscopic'],c='red')

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