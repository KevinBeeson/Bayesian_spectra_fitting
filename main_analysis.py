#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 15:39:26 2023

@author: kevin
"""


import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sobject_id_name", type=int)
parser.add_argument("-p", "--prior", type=bool)
parser.add_argument("-n","--ncpu", required=False,type=int,default=1)
args = parser.parse_args()
name = args.sobject_id_name
ncpu=args.ncpu
prior=args.prior
print('importing analysis program')

import Payne_machine_solar_finder
print('Analysing ' + str(name) +' with prior = ' + str(prior) +' with '+str(ncpu)+' cpus' )
Payne_machine_solar_finder.main_analysis(name, prior,ncpu)