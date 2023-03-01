#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 15:39:26 2023

@author: kevin
"""


import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sobject_id_name", type=int)
parser.add_argument("-p", "--prior", type=str)
parser.add_argument("-n","--ncpu", required=False,type=int,default=1)
parser.add_argument("-c","--cluster", required=False,type=str,default=None)
args = parser.parse_args()
name = args.sobject_id_name
ncpu=args.ncpu
prior=args.prior
cluster=args.cluster
if prior=='False':
    prior=False
elif prior=='True':
    prior=True
else:
    print('need to give false or true for prior')
if isinstance(prior,bool):
    print('importing analysis program')
    
    import Payne_machine_solar_finder
    print('Analysing ' + str(name) +' with prior = ' + str(prior) +' with '+str(ncpu)+' cpus')
    Payne_machine_solar_finder.main_analysis(name, prior,ncpu,cluster)