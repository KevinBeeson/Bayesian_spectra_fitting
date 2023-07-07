#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:41:44 2023

@author: kevin
"""
from astropy.table import Table,vstack
from astropy.io.votable import parse,from_table,writeto
import numpy as np
masks = Table.read('solar_spectrum_mask.fits')
masks.add_row((7745.9,7746.87,'kevin_bad_spectrums'))
masks.add_row((5720.59,5721.42,'kevin_bad_spectrums'))
masks.add_row((5743.59,5744.61,'kevin_bad_spectrums'))
# masks.add_row((5660.1,5661.1,'kevin_bad_spectrums'))

votable=from_table(masks)
writeto(votable,"spectrum_mask_kevin.fits")

