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
cluster_name='Melotte_22'
global large_data
votable = parse(cluster_name+"_photometric_cross.xml")
photometric_data=votable.get_first_table().to_table(use_names_over_ids=True)

np.savetxt('melotte_22_targets.txt',photometric_data['sobject_id'].astype(str),fmt='%s')