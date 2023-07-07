#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 13:18:51 2022

@author: kevin
"""
from scipy.optimize import curve_fit
import numpy as np 
# np.loadtxt('machine_learning_rate.txt',format=str)
import matplotlib.pyplot as plt
with open('machine_learning_rate_2.txt') as f:
    lines = f.readlines()
training_lost=[]
validation_lost=[]
for x in lines:
    for value,y in enumerate(x):
        if y=='=' and value<45:
            training_lost.append(float(x[value+1:value+8]))
        elif y=='=':
            validation_lost.append(float(x[value+1:value+8]))
            
def moving_average(a, n=20) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
def func(x, a, b, c):
     return a*np.exp(-b*x) + c
validation_lost=validation_lost[10:]
x=np.linspace(0,100,len(validation_lost))
ppt=curve_fit(func, x, validation_lost)[0]

best_fit=func(x, ppt[0], ppt[1], ppt[2])
plt.figure()
plt.plot(x,validation_lost)
plt.plot(x,best_fit)
plt.yscale('log')
training_lost=training_lost[10:]
x=np.linspace(0,100,len(training_lost))
ppt=curve_fit(func, x, training_lost)[0]

best_fit=func(x, ppt[0], ppt[1], ppt[2])
plt.figure()
plt.plot(x,training_lost)
plt.plot(x,best_fit)
# plt.figure()
# plt.plot(rate)
# plt.plot(moving_rate)
plt.yscale('log')