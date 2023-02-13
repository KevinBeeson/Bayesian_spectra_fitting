#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 11:49:07 2021

@author: kevin
"""

from astropy.table import Table,vstack
import csv
import numpy as np
from astropy.io.votable import parse,from_table,writeto

largeData=Table.read('GALAH_DR3_main_200331.fits')
name='Berkeley_32'
votable = parse(name+"Final.xml")
data=votable.get_first_table().to_table(use_names_over_ids=True)

hermes=[x for x in largeData if x['source_id'] in data['source_id']]
hermes=vstack(hermes)