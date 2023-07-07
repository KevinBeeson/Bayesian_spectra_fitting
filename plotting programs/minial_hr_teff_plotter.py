#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 10:53:33 2023

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
import csv
from astropy.io import fits

cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)
print(cluster_details_all[:,0])
for name in cluster_details_all[:,0]:
    print(name)
    name='NGC_2682_bp_rp_run.xml'
    votable = parse(name+'_bp_rp_run.xml')
    data=votable.get_first_table().to_table(use_names_over_ids=True)

    plt.figure()
    plt.plot(data['bp_rp'],data['phot_g_mean_mag'],c=data['teff'])
    
