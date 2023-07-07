#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 13:09:41 2022

@author: kevin
"""
from numba import jit
import time
from functools import  partial
from astropy.io.votable import parse
import emcee
# from The_Payne import spectral_model
import scipy
from scipy import signal
from os.path import exists
import subprocess
from pathlib import Path
from datetime import date
import os.path
import logging
from astropy.table import Table,vstack
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
from astropy.io import fits
import copy
from astropy.table import QTable,join
# import Galah_tool_py3 as gtools
import os
# import random
import functools
from multiprocessing.dummy import Pool as ThreadPool 
from scipy.stats import chi2
import warnings
from astropy.table import QTable
from astropy.io.votable import parse,from_table,writeto

cluster_name='Ruprecht_147'
votable = parse(cluster_name+"_photometric_cross.xml")
cross_data=votable.get_first_table().to_table(use_names_over_ids=True)
parameters=['teff','logg','fe_h','vmic','vsini','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
parameters_tags=[]
for x in parameters:
	if x in elements:
		parameters_tags.append(x+'_Fe')
	else:
		parameters_tags.append(x)
sobject_id_existing=[]
new_variables=[]
new_sigmas=[]
convergence=[]

# elem=np.mean(sampler_elem.get_chain(flat=True,discard=sampler_elem.iteration//2),axis=0)
# elem_sigma=np.sqrt(np.var(sampler_elem.get_chain(flat=True,discard=sampler_elem.iteration//2),axis=0))

# new_variables_temp=np.hstack((elem[:5],vrad,elem[5:]))
# new_variables.append(new_variables_temp)
# new_sigma_temp=np.hstack((elem_sigma[:5],vrad_sigma,elem_sigma[5:]))
# new_sigmas.append(new_sigma_temp)
change_ratio=[]
iteration=[]
radial_velocties=[]
for star in cross_data:
    name=str(star['sobject_id'])
    print(name)
    filename = cluster_name+'_reduction/fixed_iron_run_3_no_prior_'+str(name)
    if not (os.path.exists(filename+'_mask_finding_loop.h5') and os.path.exists(filename+'_main_loop.h5') and os.path.exists(filename+'_radial_velocities')):
          print('not synthesized '+ name)
          continue

    # star=cross_data[0]
    radial_velocities_temp=np.load(filename+'_radial_velocities')
    radial_velocties.append(radial_velocities_temp)
    
    sobject_id_existing.append(star['sobject_id'])
    sampler_elem=emcee.backends.HDFBackend(filename+'_main_loop.h5')
    iteration.append(sampler_elem.iteration)
    tags=np.load(filename+'_tags.npy')
    autocorellation=np.round(sampler_elem.get_autocorr_time(tol=0,discard=100)).astype(int)
    data=sampler_elem.get_chain(flat=True,discard=100)
    value=0
    convergence_temp=[]
    new_variables_temp=[]
    new_sigmas_temp=[]
    non_zeros_elements=len(data[0])
    change_ratio_temp=[]
    for param in parameters:
        if param in tags:
            corelation=autocorellation[value]
            correlation_max=len(data)//corelation
            mean_older=np.mean(data[-int(len(data[:,value])*0.2):,value])
            mean_younger=np.mean(data[-int(len(data[:,value])*0.5):-int(len(data)*0.3),value])
            if np.isnan(mean_younger):
                print('oh no')
            change=abs(mean_older-mean_younger)
            standard_deviation=np.sqrt(np.var(data[corelation*4*non_zeros_elements:,value]))
            change_ratio_temp.append(standard_deviation/change)
            if standard_deviation/change>3:
                convergence_temp.append(0)
            else:
                convergence_temp.append(1)
            new_variables_temp.append(np.mean(data[corelation*4*non_zeros_elements:,value]))
            new_sigmas_temp.append(change)
            value+=1
        else:
            change_ratio_temp.append(np.nan)
            convergence_temp.append(2)
            new_variables_temp.append(np.nan)
            new_sigmas_temp.append(np.nan)
    change_ratio.append(change_ratio_temp)
    new_variables.append(new_variables_temp)
    new_sigmas.append(new_sigmas_temp)
    convergence.append(convergence_temp)
change_ratio=np.array(change_ratio)  
change_ratio_no_prior=[np.nanmean(x) for x in change_ratio.T]
sobject_id_existing=np.vstack(sobject_id_existing)
print(sobject_id_existing)
new_sigmas=np.vstack(new_sigmas)
new_variables=np.vstack(new_variables)
convergence=np.vstack(convergence)

parameters_no_prior=[x+'_no_prior' for x in parameters_tags]
tab=QTable(new_variables,
            names=(parameters_no_prior),
            meta={'name': 'first table'})
tab['sobject_id']=sobject_id_existing
for conv,sig,nam in zip(convergence.T,new_sigmas.T,parameters_no_prior):
    tab['e_'+nam]=sig
    tab['flag_'+nam]=conv
tab['radial_velocities_no_prior']=radial_velocties
tab['iteration_no_prior']=iteration
sobject_id_existing=[]
new_variables=[]
new_sigmas=[]
change_ratio=[]
convergence=[]
iteration=[]
radial_velocties=[]
for star in cross_data:
    name=str(star['sobject_id'])
    print(name)
    filename = cluster_name+'_reduction/fixed_iron_run_prior_'+str(name)
    if not (os.path.exists(filename+'_mask_finding_loop.h5') and os.path.exists(filename+'_main_loop.h5') and os.path.exists(filename+'_radial_velocities')):
          print('already done '+ name)
          continue

    # star=cross_data[0]
    radial_velocities_temp=np.load(filename+'_radial_velocities')
    radial_velocties.append(radial_velocities_temp)
    sobject_id_existing.append(star['sobject_id'])
    sampler_elem=emcee.backends.HDFBackend(filename+'_main_loop.h5')
    tags=np.load(filename+'_tags.npy')
    autocorellation=np.round(sampler_elem.get_autocorr_time(tol=0,discard=100)).astype(int)
    iteration.append(sampler_elem.iteration)
    data=sampler_elem.get_chain(flat=True,discard=100)
    value=0
    convergence_temp=[]
    new_variables_temp=[]
    new_sigmas_temp=[]
    non_zeros_elements=len(data[0])
    change_ratio_temp=[]
    for param in parameters:
        if param in tags:
            corelation=autocorellation[value]
            correlation_max=len(data)//corelation
            mean_older=np.mean(data[-int(len(data[:,value])*0.2):,value])
            mean_younger=np.mean(data[-int(len(data[:,value])*0.5):-int(len(data)*0.3),value])
            if np.isnan(mean_younger):
                print('oh no')
            change=abs(mean_older-mean_younger)
            standard_deviation=np.sqrt(np.var(data[corelation*4*non_zeros_elements:,value]))
            change_ratio_temp.append(standard_deviation/change)
            if standard_deviation/change>3:
                convergence_temp.append(0)
            else:
                convergence_temp.append(1)
            new_variables_temp.append(np.mean(data[corelation*4*non_zeros_elements:,value]))
            new_sigmas_temp.append(change)
            value+=1
        else:
            change_ratio_temp.append(np.nan)
            convergence_temp.append(2)
            new_variables_temp.append(np.nan)
            new_sigmas_temp.append(np.nan)
    change_ratio.append(change_ratio_temp)
    new_variables.append(new_variables_temp)
    new_sigmas.append(new_sigmas_temp)
    convergence.append(convergence_temp)
change_ratio=np.array(change_ratio)  
change_ratio_prior=[np.nanmean(x) for x in change_ratio.T]
fig,ax=plt.subplots(1,1)
ax.plot(np.linspace(0,len(parameters),num=len(parameters)),change_ratio_prior,c='Blue')
ax.plot(np.linspace(0,len(parameters),num=len(parameters)),change_ratio_no_prior,c='Red')
ax.set_xticks(range(0,len(parameters)))
ax.set_xticklabels(parameters)

sobject_id_existing=np.vstack(sobject_id_existing)
print(sobject_id_existing)
new_sigmas=np.vstack(new_sigmas)
new_variables=np.vstack(new_variables)
convergence=np.vstack(convergence)

parameters_prior=[x+'_prior' for x in parameters_tags]
for conv,var,sig,nam in zip(convergence.T,new_variables.T,new_sigmas.T,parameters_prior):
    tab['e_'+nam]=sig
    tab['flag_'+nam]=conv
    tab[nam]=var
tab['radial_velocities_prior']=radial_velocties
tab['iteration_prior']=iteration
tab['sobject_id']=tab['sobject_id'].reshape(len(tab))
tab=join(tab,cross_data,keys='sobject_id')
votable=from_table(tab)
writeto(votable,cluster_name+"_3_fixed_iron.xml")

    
