#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 15:53:11 2023

@author: kevin
"""
from scipy.stats import kde
from scipy import integrate

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
from astropy.io.votable import parse,from_table,writeto

import numpy as np
from multiprocessing import Pool
from astropy.io import fits
import copy
# import Galah_tool_py3 as gtools
import os
# import random
import functools
from multiprocessing.dummy import Pool as ThreadPool 
from scipy.stats import chi2
import warnings
global large_data
import csv
votable = parse("open_cluster_photometric_cross.xml")
photometric_data=votable.get_first_table().to_table(use_names_over_ids=True)
cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)
print(cluster_details_all[:,0])
names=cluster_details_all[:,0]

for name in names:
    mask=photometric_data['cluster_name']==name
    photometric_data_to_save=photometric_data[mask]
    votable=from_table(photometric_data_to_save)
    writeto(votable,name+"_photometric_cross.xml")

    print(name)
    print(len(photometric_data_to_save))

    np.savetxt(name+'_targets.txt',photometric_data_to_save['sobject_id'].astype(str),fmt='%s')