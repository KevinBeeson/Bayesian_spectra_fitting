#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 08:31:48 2022

@author: kevin
"""

import numpy as np
colour='IR'
data=np.load('Training_data_Final_'+colour+'.npy',allow_pickle=True)
data=data[:-len(data)//10]
data_training=data[:round(len(data)*(6.4/7.5))]
training_data=[x[0][0] for x in data_training]

data_labels=[np.hstack(x[1]) for x in data_training ]
data_labels=[np.array(x,dtype=float) for x in data_labels]
data_labels=np.vstack(data_labels)
data_validation=data[-len(data)//10:]
validation_data=[x[0][0] for x in data_validation ]
validation_data_labels=[np.hstack(x[1]) for x in data_validation ]
validation_data_labels=[np.array(x,dtype=float) for x in validation_data_labels]
validation_data_labels=np.vstack(validation_data_labels)

np.save('traning_data_all_elements_'+colour+'.npy',training_data)
np.save('traning_labels_all_elements_'+colour+'.npy',data_labels)
np.save('validation_data_all_elements_'+colour+'.npy',validation_data)
np.save('validation_labels_all_elements_'+colour+'.npy',validation_data_labels)