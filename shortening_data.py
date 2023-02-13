#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 12:23:34 2022

@author: kevin
"""

import numpy as np
data=np.load('Training_data_Blue_clean.npy',allow_pickle=True)
# data=data[:000]
training_labels=np.vstack(data[:,1])
training_labels=np.hstack((training_labels[:,0:5],training_labels[:,6:]))
training_spectra=np.vstack(np.vstack(data[:,0]))


validation_spectra=training_spectra[-len(data)//10:]
validation_labels=training_labels[-len(data)//10:]

training_labels=training_labels[:-len(data)//10]
training_spectra=training_spectra[:-len(data)//10]
np.save('Training_data_Blue',training_spectra)
np.save('Training_labels_Blue',training_labels)

np.save('Validation_data_Blue',validation_spectra)
np.save('Validation_labels_Blue',validation_labels)
