#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 11:36:25 2022

@author: kevin
"""

import numpy as np
from The_Payne import utils
from The_Payne import spectral_model
from The_Payne import fitting
from multiprocessing import Pool

def payne_sythesize(solar_values,x_min,x_max,NN_coeffs):

    scaled_labels = (solar_values-x_min)/(x_max-x_min) - 0.5

    real_spec = spectral_model.get_spectrum_from_neural_net(scaled_labels = scaled_labels, NN_coeffs = NN_coeffs)
    return real_spec
def difference(data_in):
    label,spectra=data_in
    synthetic_spectra=payne_sythesize(label, x_min, x_max, NN_coeffs)
    diffrence=synthetic_spectra-spectra
    return diffrence
        

training_data=np.load('Training_data_no_degredation_all_clean_Fe_Alpha.npy',allow_pickle=True)
print(len(training_data))

print(len(training_data))

training_labels=np.vstack(training_data[:,1])
training_labels=np.hstack((training_labels[:,0:5],training_labels[:,6:]))
training_spectra=np.vstack(np.vstack(training_data[:,0]))

validation_spectra=training_spectra[-len(training_spectra)//4:]
validation_labels=training_labels[-len(training_labels)//4:]

colours=['Blue','Green','Red','IR']
for values,x in enumerate(colours):
    tmp = np.load("NN_normalized_spectra_no_v_rad_"+x+".npz")
    w_array_0 = tmp["w_array_0"]
    w_array_1 = tmp["w_array_1"]
    w_array_2 = tmp["w_array_2"]
    b_array_0 = tmp["b_array_0"]
    b_array_1 = tmp["b_array_1"]
    b_array_2 = tmp["b_array_2"]
    x_min = tmp["x_min"]
    x_max = tmp["x_max"]
    tmp.close()
    NN_coeffs= (w_array_0, w_array_1, w_array_2, b_array_0, b_array_1, b_array_2, x_min, x_max)
    pool=Pool(6)
    diffrence_final=pool.map(difference,zip(validation_spectra,validation_labels))
    # length=len(validation_spectra)
    # for spectra,label in zip(validation_spectra,validation_labels):
    #     synthetic_spectra=payne_sythesize(label, x_min, x_max, NN_coeffs)
        