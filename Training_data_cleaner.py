#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 10:51:35 2022

@author: kevin
"""

import numpy as np
large_data=np.load('Training_data_no_degredation_all_Fe_Alpha.npy',allow_pickle=True)
# large_data=large_data
large_data_new=[]
for value,x in enumerate(large_data):

    if value%100==0:
        print(value/len(large_data))
    a=[y for y in x[0]]
    this_max=np.max([np.nanmax(w) for w in a])
    this_min=np.min([np.nanmin(w) for w in a])

    if not (np.any([np.isnan(w).any() for w in a ]) or  np.isnan(x[1]).any() or this_max>2 or this_min<0):
        large_data_new.append(x)
    else:
        print(x,'   ',this_min,'    ',this_max)
    # this_max=np.max([np.nanmax(w) for w in a])
    # this_min=np.min([np.nanmin(w) for w in a])
    # if this_max>current_max:
    #     print('max ', this_max)
    #     current_max=this_max
    # if this_min<current_min:
    #     print('min ',this_min)
    #     current_min=this_min

np.save('Training_data_no_degredation_all_clean_Fe_Alpha.npy',large_data_new)