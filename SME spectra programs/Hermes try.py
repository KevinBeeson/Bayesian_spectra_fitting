#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 11:40:36 2020

@author: kevin
"""
import os.path

from pysme.gui import plot_plotly
from pysme import sme as SME
from pysme import util
from pysme.solve import solve
from pysme.synthesize import synthesize_spectrum

from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
from pysme.persistence import save_as_idl
import os
import matplotlib
if os.environ.get('DISPLAY') is None:
    # enables figure saving on clusters with ssh conection
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astropy.io import fits
from glob import glob
from scipy.interpolate import splrep, splev
from scipy.signal import savgol_filter, argrelextrema

def get_spectra_dr52(object, bands=[1,2,3,4], root='', read_sigma=False, remove_nan=False,
                     extension=4, individual=False):
    if individual:
        subfolder = 'all'
    else:
        subfolder = 'com'
    fits_path = root + object[0:6] + '/standard/'+subfolder+'/' + object
    fits_path_2 = root + object[0:6] + '/standard/'+subfolder+'/' + object
    # detersme path of spectra
    if not(os.path.isfile(fits_path + '1.fits')):
        fits_path = fits_path_2
        if not(os.path.isfile(fits_path + '1.fits')):
            print ('Spectra file not found')
            if read_sigma:
                return np.array([]), np.array([]), np.array([])
            else:
                return np.array([]), np.array([])
    # read selected spectrum bands
    spect_all = list([])
    if read_sigma:
        sigma_all = list([])
    wvls_all = list([])

    for i_b in range(len(bands)):
        band = bands[i_b]
        fits_path_band = fits_path + str(band) + '.fits'
        fits_data = fits.open(fits_path_band, memmap=False)
        if len(fits_data) < (extension+1):
            return list([]), list([])
        data_len = len(fits_data[extension].data)
        spect_all.append(fits_data[extension].data)
        if read_sigma:
            noise_data_len = len(fits_data[1].data)
            sigma_all.append(fits_data[1].data[noise_data_len-data_len:noise_data_len])
        header = fits_data[extension].header
        # calculate wavelengths of observed spectra
        wvls_all.append(header.get('CRVAL1') + np.float64(range(0, data_len)) * header.get('CDELT1'))
        # print header.get('CDELT1')
        fits_data.close()
    if read_sigma:
        return spect_all, wvls_all, sigma_all
    else:
        return spect_all, wvls_all

sme=SME.SME_Structure()
hermes=fits.open("1408240023010341.fits")
fstart= hermes[4].header.get('CRVAL1')
fend= hermes[4].header.get('CRVAL1')+len(hermes[4].data)* hermes[4].header.get('CDELT1')



length=np.array([fstart+y*hermes[4].header.get('CDELT1') for y in range(0,len(hermes[4].data))])
sme.wave= [length]
sme.spec=[ hermes[4].data]
sme.uncs=[hermes[1].data]
sme.wran =[fstart,fend]
sme.abund = Abund(0, "asplund2009")
sme.mask = [np.full(i.size, 1) for i in sme.spec]
sme.mask[sme.spec == 0] = SME.SME_Structure.mask_values["bad"]



sme.atmo.source = "marcs2012.sav"
sme.atmo.method = "grid"
sme.atmo.geom = "PP"

# Change parameters if your want
sme.vsini = 0
sme.vrad = 0.35
sme.vrad_flag = "whole"
sme.cscale_flag = "linear"
sme.cscale_type = "whole"
# sme.linelist = ValdFile(
sme.linelist=ValdFile("KevinLBeeson.4201751")

sme.nlte.set_nlte("Fe", "marcs2012_Fe2016.grd")
sme.nlte.set_nlte("Si","marcs2012_Si2016.grd")
sme.nlte.set_nlte("Ca","marcs2012p_t1.0_Ca.grd")
sme.nlte.set_nlte("O","marcs2012_O2015.grd")
# sme.nlte.set_nlte("Ba","marcs2012p_t1.0_Ba.grd ")
# sme.nlte.set_nlte("Na","marcs2012p_t1.0_Na.grd")
sme.nlte.set_nlte("Ti","marcs2012s_t2.0_Ti.grd")

fitparameters = [ "Abund Li","Abund Be","abund B","abund C","abund N","abund o","abund f","abund Ne","abund Na","abund Mg","abund Al","abund Si","abund P","abund S","abund Cl","abund Ar","abund K","abund Ca","abund Sc","abund Ti","abund V",
                 "abund Cr","abund Mn","abund Fe","abund Co","abund Ni","abund Cu","abund Zn","abund Ga","abund Ge","abund As","abund Se","abund Br","abund Kr","abund Rb","abund Sr","abund Y","abund Zr","abund Nb","abund Mo","abund Ru","abund Rh","abund Pd",
                 "abund Ag","abund Cd","abund In","abund Sn","abund Sb","abund Te","abund I","abund Xe","abund Cs","abund Ba","abund La","abund Ce","abund Pr","abund Nd","abund Sm","abund Eu","abund Gd","abund Tb","abund Dy","abund Ho","abund Er","abund Tm"
                 ,"abund Yb","abund Lu","abund Hf","abund Ta","abund W","abund Re","abund Os","abund Ir","abund Pt","abund Au","abund Hg","abund Tl","abund Pb","abund Bi","abund Th","abund U"
                 ,"monh","teff", "logg","vrad","vmic","vmac","vsini"]

# Start SME solver
# sme = synthesize_spectrum(sme)
sme = solve(sme, fitparameters)

print(sme.citation())

# Save results
sme.save("second")

# Plot results
fig = plot_plotly.FinalPlot(sme)
fig.save(filename="plot_file2.html")
