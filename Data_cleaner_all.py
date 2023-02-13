#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 16:11:14 2022

@author: kevin
"""

import numpy as np
data=np.load('Training_data_alpha.npy',allow_pickle=True)

labels=[np.hstack((x[1][:4],x[1][5:])) for x in data ]
labels_validation=labels[-1500:]
labels_training=labels[:-1500]
data_Blue=[x[0][0] for x in data if not( np.any(x[0][0]<0) or np.any(x[0][0]>2))]
data_Blue_validation=data_Blue[-1500:]
data_Blue_training=data_Blue[:-1500]
data_Green=[x[0][1] for x in data if not( np.any(x[0][0]<0) or np.any(x[0][0]>2))]
data_Green_validation=data_Green[-1500:]
data_Green_training=data_Green[:-1500]
data_Red=[x[0][2] for x in data if not( np.any(x[0][0]<0) or np.any(x[0][0]>2))]
data_Red_validation=data_Red[-1500:]
data_Red_training=data_Red[:-1500]
data_IR=[x[0][3] for x in data if not( np.any(x[0][0]<0) or np.any(x[0][0]>2))]
data_IR_validation=data_IR[-1500:]
data_IR_training=data_IR[:-1500]


np.save('Training_data_Blue_clean_alpha',data_Blue_training)
np.save('Validation_data_Blue_clean_alpha',data_Blue_validation)
np.save('Training_data_Green_clean_alpha',data_Green_training)
np.save('Validation_data_Green_clean_alpha',data_Green_validation)
np.save('Training_data_Red_clean_alpha',data_Red_training)
np.save('Validation_data_Red_clean_alpha',data_Red_validation)
np.save('Training_data_IR_clean_alpha',data_IR_training)
np.save('Validation_data_IR_clean_alpha',data_IR_validation)
np.save('labels_training',labels_training)
np.save('labels_validation',labels_validation)