#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 14:11:03 2022

@author: kevin

This plots pdf plots that are made for the training data
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kde


data=np.load('traning_labels_all_elements_Blue_H.npy',allow_pickle=True)
data_val=np.load('validation_labels_all_elements_Blue_H.npy',allow_pickle=True)
size_font=10
shape=(8,5)


fig = plt.figure(figsize=shape)
# fig.set_figwidth(shape[1])
# fig.set_figheight(shape[0])
ax=[]
ax1=plt.subplot2grid(shape=shape, loc=(0, 0), colspan=3,rowspan=2)
ax=[]
for row in range(2):
    for col in range(3,shape[1]):
        ax.append(plt.subplot2grid(shape=shape,loc=(row,col)))
for row in range(2,shape[0]):
    for col in range(shape[1]):
        ax.append(plt.subplot2grid(shape=shape,loc=(row,col)))
parameters=['teff','logg','Fe/H','vmic','vbroad','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
# for x in range(5,len(parameters)):
#     data[:,x]+=data[:,2]
    
# np.save('traning_labels_all_elements_Blue_H.npy',data)
# for x in range(5,len(parameters)):
#     data_val[:,x]+=data_val[:,2]
# np.save('validation_labels_all_elements_Blue_H.npy',data_val)

ax1.scatter(data_val[:,0],data_val[:,1],s=0.3,color='black',alpha=0.3)
ax1.set_xlabel(r'$T_{\rm{eff}}$')
ax1.xaxis.set_label_position('top') 
ax1.tick_params('x', top=True, labeltop=True,bottom=False,labelbottom=False)
x=np.array(data_val[:,0],dtype=float)
y=np.array(data_val[:,1],dtype=float)
nbins=300
# k = kde.gaussian_kde([x,y])
# xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
# zi = k(np.vstack([xi.flatten(), yi.flatten()]))
# ax1.contour(xi, yi, zi.reshape(xi.shape),levels=5)

ax1.text(-0.125 ,0.5,r'log($g$)' ,horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes,rotation=90)

ax1.invert_yaxis()
ax1.invert_xaxis()

for value,(axis,col_dat) in enumerate(zip(ax[:3],data_val.T[2:5]),2):
    axis.hist(col_dat,density=False)
    # axis.set_xlabel(parameters[value])
ax[0].text(0.8,0.8,parameters[2],horizontalalignment='center',verticalalignment='center', transform=ax[0].transAxes,c='red',size=size_font,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5))
ax[1].text(0.83,0.8,parameters[3],horizontalalignment='center',verticalalignment='center', transform=ax[1].transAxes,c='red',size=size_font,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5))
ax[2].text(0.77,0.8,parameters[4],horizontalalignment='center',verticalalignment='center', transform=ax[2].transAxes,c='red',size=size_font,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5))

for value,(axis,col_dat) in enumerate(zip(ax[3:],data_val.T[5:]),5):
    axis.scatter(data_val[:,2],col_dat,s=0.1,color='black',alpha=0.3)
    nbins=300
    x=np.array(data_val[:,2],dtype=float)
    y=np.array(col_dat,dtype=float)
    axis.text(0.88,0.8,parameters[value],horizontalalignment='center',verticalalignment='center', transform=axis.transAxes,c='red',bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5),size=size_font)
    # k = kde.gaussian_kde([x,y])
    # xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    # zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    # axis.contour(xi, yi, zi.reshape(xi.shape),levels=5)
for axis in ax[-5:]:
    axis.set_xlabel('Fe/H')
for axis in ax[3:4]:
    axis.set_xlabel('Fe/H')
for value,axis in enumerate(ax[4:]):
    if value%shape[1]==0:
        axis.text(-0.45,0.5,'[x/Fe]' ,horizontalalignment='center',verticalalignment='center', transform=axis.transAxes,rotation=90)
fig.set_size_inches(6.97384806974,6.97384806974/4*8/1.6)
plt.tight_layout(w_pad=0.3,h_pad=-1.2)
# plt.savefig('/home/kevin/Documents/Paper/training data pdfs.png',dpi=200)