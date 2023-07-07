#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 11:06:17 2022

@author: kevin
"""
from astropy.table import Table,vstack
from astropy.io import fits
import numpy as np 
data=np.load('testing_fit_3.npy',allow_pickle=True)

data_galah=fits.getdata('dr6.0.fits',1)
all_reduced_data=Table(data_galah)

data_galah_small=data_galah[::100]
data_galah_small=Table(data_galah_small[:9000])
error=[x[0] for x in data]
data_galah_small['error_fit']=error
data_clean=[x for x in data_galah_small if np.isnan(x['teff_r'])==False and not x['error_fit'] is None]
data_clean=vstack(data_clean)
data_clean['snr_com']=np.mean(data_clean['snr'],axis=1)
import matplotlib.pyplot as plt
for x in data_clean:
    if abs(x['fe_h_r'])<4 and abs(x['error_fit'])<1e-17:
        print('here')
data_iron=[x for x in data_clean if abs(x['fe_h_r'])<4 and abs(x['error_fit'])<1e5]
data_iron=vstack(data_iron)

plt.figure()
plt.scatter(data_iron['fe_h_r'],data_iron['error_fit'],s=0.5)

data_teff=[x for x in data_clean if abs(x['error_fit'])<1e7]
data_teff=vstack(data_teff)
data_teff['telluric_h2o_com']=np.mean(data_teff['telluric_h2o'],axis=1)
data_teff['telluric_o2_com']=np.nanmean(data_teff['telluric_o2'],axis=1)

data_teff['res_com']=np.min(data_teff['res'],axis=1)


plt.figure()
plt.scatter(data_teff['teff_r'],data_teff['error_fit'],s=0.5)
plt.xlabel('teff_reduction')
plt.ylabel('fit_error')
plt.tight_layout()

plt.figure()
plt.scatter(data_teff['epoch'],data_teff['error_fit'],s=0.5)
plt.ylabel('fit_error')
plt.xlabel('epoch')
plt.tight_layout()
plt.figure()
plt.scatter(data_teff['mean_airmass'],data_teff['error_fit'],s=0.5)
plt.ylabel('fit_error')
plt.xlabel('mean_airmass')
plt.tight_layout()
plt.figure()
plt.scatter(data_teff['snr_com'],data_teff['error_fit'],s=0.5)
plt.ylabel('fit_error')
plt.xlabel('snr_com')

plt.figure()
plt.scatter(data_teff['mean_ra'],data_teff['error_fit'],s=0.5)
plt.ylabel('fit_error')
plt.xlabel('mean_ra')
plt.tight_layout()

plt.figure()
plt.scatter(data_teff['mean_dec'],data_teff['error_fit'],s=0.5)
plt.ylabel('fit_error')
plt.xlabel('mean_dec')
plt.tight_layout()
plt.figure()
plt.scatter(data_teff['mean_zd'],data_teff['error_fit'],s=0.5)
plt.ylabel('fit_error')
plt.xlabel('mean_zd')
plt.tight_layout()
plt.figure()
plt.scatter(data_teff['telluric_h2o_com'],data_teff['error_fit'],s=0.5)
plt.ylabel('fit_error')
plt.xlabel('telluric_h2o_com')
plt.tight_layout()

plt.figure()
plt.scatter(data_teff['telluric_o2_com'],data_teff['error_fit'],s=0.5)
plt.ylabel('fit_error')
plt.xlabel('telluric_o2_com')
plt.tight_layout()

plt.figure()
plt.scatter(data_teff['res_com'],data_teff['error_fit'],s=0.5)
plt.ylabel('fit_error')
plt.xlabel('res_com')
plt.tight_layout()

plt.figure()
plt.scatter(data_teff['mag'],data_teff['error_fit'],s=0.5)
plt.ylabel('fit_error')
plt.xlabel('mag')
plt.tight_layout()


plt.figure()
plt.scatter(data_teff['logg_r'],data_teff['error_fit'],s=0.5)
plt.ylabel('fit_error')
plt.xlabel('logg')


plt.figure()
plt.scatter(data_teff['e_rv_com'],data_teff['error_fit'],s=0.5)
plt.ylabel('fit_error')
plt.xlabel('e_rv')





plt.figure()
plt.hexbin(data_iron['teff_r'],data_iron['logg_r'],data_iron['error_fit'],gridsize=int(1e3))

plt.xlabel('teff')
plt.ylabel('logg_r')
plt.colorbar()

plt.figure()
plt.hexbin(data_iron['teff_r'],data_iron['fe_h_r'],data_iron['error_fit'],gridsize=int(1e3))
plt.xlabel('teff')
plt.ylabel('fe_h')
plt.colorbar()



plt.figure()
plt.hexbin(data_iron['mean_ra'],data_iron['mean_dec'],data_iron['error_fit'],gridsize=int(5e2))
plt.xlabel('mean_ra')
plt.ylabel('mean_dec')
plt.colorbar()

plt.figure()
plt.hexbin(data_iron['alpha_fe_r'],data_iron['teff_r'],data_iron['error_fit'],gridsize=int(5e2))
plt.xlabel('alpha_fe_r')
plt.ylabel('fe_r')
plt.colorbar()