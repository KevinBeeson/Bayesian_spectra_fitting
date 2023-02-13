#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 08:34:11 2022

@author: kevin
"""

import emcee
import csv
import numpy as np 
import matplotlib.pyplot as plt
import corner

cluster_details_all = list(csv.reader(open('my_targets.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)

for name in cluster_details_all:
    data=emcee.backends.HDFBackend(name+'try_3.h5')
    corner.corner(data.get_chain(flat=True)[300:])
    plt.title(name)
    plt.figure()
    plt.plot(data.get_chain(flat=True)[:,1])
    plt.title(name)