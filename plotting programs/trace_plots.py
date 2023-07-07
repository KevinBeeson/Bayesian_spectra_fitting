#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 10:41:47 2022

@author: kevin
"""

import arviz as az
import emcee
import numpy as np
import matplotlib.pyplot as plt
import corner
from astropy.io.votable import parse

plt.rcParams['font.family']='Times New Roman'
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

cluster_name='NGC_2682'
votable = parse(cluster_name+"_cross_galah_pdf.xml")
cross_data=votable.get_first_table().to_table(use_names_over_ids=True)

for star in cross_data:
    try:
        name=str(star['sobject_id'])
        # name='160106004101363'
        sampler=emcee.backends.HDFBackend('NGC_2682_reduction_dr61_2/'+name+'_all_elements.h5')
        parameters_no_vrad=['teff','logg','fe_h','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
        idata1 = az.from_emcee(sampler, var_names=parameters_no_vrad)
        idata1=idata1.sel(draw=slice(sampler.iteration//3, None))
        
        
        # az.plot_pair(idata1,var_names=parameters_no_vrad[:5],divergences=True)
        axes=az.plot_trace(idata1,var_names=parameters_no_vrad[:12])
        axes[0,0].title.text = str(name)
        # plt.set_size_inches(10,20)
        plt.tight_layout()
        # plt.title(str(name))
        plt.savefig('NGC_2682_reduction_dr61_2/plots/'+name+'_trace_1.pdf')
        axes=az.plot_trace(idata1,var_names=parameters_no_vrad[12:24],compact=True)
        axes[0,0].title.text = str(name)
        # plt.title(str(name))
        
        # plt.set_size_inches(10,20)
        plt.tight_layout()
        plt.savefig('NGC_2682_reduction_dr61_2/plots/'+name+'_trace_2.pdf')
        axes=az.plot_trace(idata1,var_names=parameters_no_vrad[-12:],compact=True)
        axes[0,0].title.text = str(name)
        # plt.title(str(name))
        
        # plt.set_size_inches(10,20)
        plt.tight_layout()
        plt.savefig('NGC_2682_reduction_dr61_2/plots/'+name+'_trace_3.pdf')
        
        
        
        fig=corner.corner(idata1,var_names=parameters_no_vrad[:5])
        fig.suptitle(str(name))
        
        plt.savefig('NGC_2682_reduction_dr61_2/plots/'+name+'_corner.pdf')
        
        
        
        sampler=emcee.backends.HDFBackend('NGC_2682_reduction_dr61_2/'+name+'_radial_velocities.h5')
        params_short=['teff','logg','fe_h','vmic','vsini','vrad_Blue','vrad_Green','vrad_Red','vrad_IR']
        idata1 = az.from_emcee(sampler, var_names=params_short)
        fig=corner.corner(idata1,var_names=np.hstack((params_short[:3],params_short[-4:])))
        fig.suptitle(str(name))
        
        
        plt.savefig('NGC_2682_reduction_dr61_2/plots/'+name+'_radial_velocity_corner.pdf')
        print('reduced '+name)
    except:
        print('havent reduced '+str(star['sobject_id']))

# coords = {"school": ["teff", "logg"]}
# ax = az.plot_pair(
#     idata1,
#     var_names=parameters_no_vrad[:5],
#     kind=["scatter", "kde"],
#     kde_kwargs={"fill_last": False},
#     marginals=True,
#     # coords=coords,
#     point_estimate="median",
#     figsize=(11.5, 5),
# )
