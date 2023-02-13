#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 11:08:03 2022

@author: kevin
"""

from astropy.io import fits
from astropy.table import Table,vstack
import emcee
import numpy as np 

data=fits.get_data('Open_cluster_kevin.fits',1)
data=Table(data[0].data)