#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 16:09:34 2023

@author: kevin
"""

from scipy.stats import kde
from scipy import integrate
from numba import jit
import numpy as np
from astropy.io.votable import parse
import matplotlib.pyplot as plt

@jit(nopython=True,cache=True)
def normal(x,mu,sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-1/2*((x-mu)/sigma)**2)

cluster_name='Ruprecht_147'
global large_data
votable = parse(cluster_name+"_photometric_cross.xml")
photometric_data=votable.get_first_table().to_table(use_names_over_ids=True)
photometric_data=photometric_data[0]
photometric_probability_density=kde.gaussian_kde([photometric_data['teff_raw'],photometric_data['logg_raw']])
e_teff=np.mean(photometric_data['e_teff_raw'])
e_logg=np.mean(photometric_data['e_logg_raw'])
f= lambda teff_var,logg_var,teff_0,logg_0:photometric_probability_density([teff_var,logg_var])*normal(teff_0,teff_var,e_teff)*normal(logg_0,logg_var,e_logg)
total_int=0
teff_prior=3594 
logg_prior=4.72
teff_line=np.linspace(teff_prior-e_teff*5, teff_prior+e_teff*5,8)
logg_line=np.linspace(logg_prior-e_logg*5, logg_prior+e_logg*5,8)
xi, yi = np.mgrid[teff_prior-e_teff*5:teff_prior+e_teff*5:8*1j, logg_prior-e_logg*5:logg_prior+e_logg*5:8*1j]
zi = photometric_probability_density(np.vstack([xi.flatten(), yi.flatten()]))
plt.contour(xi, yi, zi.reshape(xi.shape),levels=5)

for x in range(len(teff_line)-1):
    for y in range(len(logg_line)-1):
        total_int+=integrate.dblquad(f, logg_line[y], logg_line[y+1], lambda teff_var: teff_line[x], lambda teff_var:teff_line[x+1],args=(teff_prior,logg_prior))[0]
        
