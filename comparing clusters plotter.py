#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 13:46:08 2022

@author: kevin
"""
from functools import  partial
from astropy.io.votable import parse

# from The_Payne import spectral_model

from pathlib import Path

from astropy.table import Table,vstack,join
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
import csv
plt.rcParams['font.family']='Times New Roman'
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)
print(cluster_details_all[:,0])
# name=input('what isochrone do you want?')
name='NGC_2682'
cluster_details=[x for x in cluster_details_all if x[0]==name][0]


name_all=['NGC_2682','NGC_2516']
# NGC_2682_summary_dr61
size_font=10
iso_type='gaia'
shape=(8,4)

fig = plt.figure(figsize=shape)
# fig.set_figwidth(shape[1])
# fig.set_figheight(shape[0])
ax=[]
for row in range(shape[0]):
    for col in range(shape[1]):
        ax.append(plt.subplot2grid(shape=shape,loc=(row,col)))
parameters_index=['fe_h','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']

for name in name_all:
    
    votable = parse(name+"_no_prior.xml")
    data=votable.get_first_table().to_table(use_names_over_ids=True)
    
    
    for value,(axis,col_dat) in enumerate(zip(ax,parameters_index)):
        axis.scatter(data['teff'],data[col_dat],s=0.5)
        axis.text(0.9,0.70,parameters_index[value],horizontalalignment='center',verticalalignment='center', transform=axis.transAxes,c='red',size=size_font,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5))

# for value,(axis,col_dat) in enumerate(zip(ax[:3],parameters_index[2:5]),2):
#     axis.hist(data_no_prior[col_dat],density=True)
# for value,(axis,col_dat) in enumerate(zip(ax[3:],data_no_prior.T[5:]),5):
#     axis.scatter(data_no_prior[:,2],col_dat,s=0.1,color='black',alpha=0.3)
#     nbins=300
#     x=np.array(data_no_prior[:,2],dtype=float)
#     y=np.array(col_dat,dtype=float)
#     k = kde.gaussian_kde([x,y])
#     xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
#     zi = k(np.vstack([xi.flatten(), yi.flatten()]))
#     axis.contour(xi, yi, zi.reshape(xi.shape),levels=5)
for axis in ax[-5:]:
    axis.set_xlabel(r'$T_{\rm{eff}}$')
for axis in ax[3:4]:
    axis.set_xlabel(r'$T_{\rm{eff}}$')
for value,axis in enumerate(ax):
    if value%shape[1]==0:
        axis.text(-0.45,0.5,r'[x/Fe]' ,horizontalalignment='center',verticalalignment='center', transform=axis.transAxes,rotation=90)

plt.tight_layout(w_pad=0.0,h_pad=-1.)
fig.set_size_inches(6.97384806974,6.97384806974/4*8/1.6)
plt.tight_layout(w_pad=0.3,h_pad=-1.2)

plt.savefig('/home/kevin/Documents/Paper/comparing different clusters.pdf', format='pdf')


# for value,axis in enumerate(ax[4:]):
#     if value%shape[1]==0:
#         axis.text(-0.45,0.5,r'frecuency' ,horizontalalignment='center',verticalalignment='center', transform=axis.transAxes,rotation=90,size=9)

# plt.tight_layout(w_pad=0.0,h_pad=-1.)
# # plt.savefig('/home/kevin/Documents/Paper/NGC_2682_dr6vsd61hist.pdf')

# shape=(8,5)

# fig = plt.figure(figsize=shape)
# # fig.set_figwidth(shape[1])
# # fig.set_figheight(shape[0])
# ax=[]
# ax1=plt.subplot2grid(shape=shape, loc=(0, 0), colspan=3,rowspan=2)
# ax=[]
# for row in range(2):
#     for col in range(3,shape[1]):
#         ax.append(plt.subplot2grid(shape=shape,loc=(row,col)))
# for row in range(2,shape[0]):
#     for col in range(shape[1]):
#         ax.append(plt.subplot2grid(shape=shape,loc=(row,col)))
# parameters=['teff','logg','Fe/H','vmic','vbroad','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
# parameters_index=['teff','logg','fe_h','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
# ax1.scatter(data_no_prior['teff'],data_no_prior['logg'],s=0.3,color='blue')
# ax1.scatter(data_prior['teff'],data_prior['logg'],s=0.3,color='red')

# ax1.set_xlabel(r'$T_{\rm{eff}}$')
# ax1.xaxis.set_label_position('top') 
# ax1.tick_params('x', top=True, labeltop=True,bottom=False,labelbottom=False)
# ax1.set_xlim(min(min(data_no_prior['teff']),min(data_prior['teff'])),max(max(data_no_prior['teff']),max(data_prior['teff'])))
# ax1.set_ylim(min(min(data_no_prior['logg']),min(data_prior['logg'])),max(max(data_no_prior['logg']),max(data_prior['logg'])))

# # x=np.array(data_no_prior['teff'],dtype=float)
# # y=np.array(data_no_prior['logg'],dtype=float)
# # nbins=300
# # k = kde.gaussian_kde([x,y])
# # xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
# # zi = k(np.vstack([xi.flatten(), yi.flatten()]))
# # ax1.contour(xi, yi, zi.reshape(xi.shape),levels=5)
# ax1.plot(10**iso[:,2],iso[:,-2],label=r'Best fit',c='black',linewidth=0.8)

# ax1.text(-0.125 ,0.5,r'log($g$)' ,horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes,rotation=90)
# circle_size=0.5
# ax1.invert_yaxis()
# ax1.invert_xaxis()
# for value,(axis,col_dat) in enumerate(zip(ax[:3],parameters_index[2:5]),2):
#     axis.scatter(data_no_prior['teff'],data_no_prior[col_dat],color='blue',s=circle_size)
#     axis.scatter(data_prior['teff'],data_prior[col_dat],color='red',s=circle_size)
# for value,(axis,col_dat) in enumerate(zip(ax[3:],parameters_index[5:]),5):
#     axis.scatter(data_no_prior['teff'],data_no_prior[col_dat],color='blue',s=circle_size)
#     axis.scatter(data_prior['teff'],data_prior[col_dat],color='red',s=circle_size)
#     axis.text(0.9,0.70,parameters[value],horizontalalignment='center',verticalalignment='center', transform=axis.transAxes,c='red',size=size_font,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5))

#     # axis.set_xlabel(parameters[value])
# ax[0].text(0.85,0.7,parameters[2],horizontalalignment='center',verticalalignment='center', transform=ax[0].transAxes,c='red',size=size_font,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5))
# ax[1].text(0.85,0.7,parameters[3],horizontalalignment='center',verticalalignment='center', transform=ax[1].transAxes,c='red',size=size_font,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5))
# ax[2].text(0.8,0.7,parameters[4],horizontalalignment='center',verticalalignment='center', transform=ax[2].transAxes,c='red',size=size_font,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5))

# # for value,(axis,col_dat) in enumerate(zip(ax[:3],parameters_index[2:5]),2):
# #     axis.hist(data_no_prior[col_dat],density=True)
# # for value,(axis,col_dat) in enumerate(zip(ax[3:],data_no_prior.T[5:]),5):
# #     axis.scatter(data_no_prior[:,2],col_dat,s=0.1,color='black',alpha=0.3)
# #     nbins=300
# #     x=np.array(data_no_prior[:,2],dtype=float)
# #     y=np.array(col_dat,dtype=float)
# #     k = kde.gaussian_kde([x,y])
# #     xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
# #     zi = k(np.vstack([xi.flatten(), yi.flatten()]))
# #     axis.contour(xi, yi, zi.reshape(xi.shape),levels=5)
# for axis in ax[-5:]:
#     axis.set_xlabel(r'$T_{\rm{eff}}$')
# for axis in ax[3:4]:
#     axis.set_xlabel(r'$T_{\rm{eff}}$')
# for value,axis in enumerate(ax[4:]):
#     if value%shape[1]==0:
#         axis.text(-0.45,0.5,r'[x/Fe]' ,horizontalalignment='center',verticalalignment='center', transform=axis.transAxes,rotation=90)

# plt.tight_layout(w_pad=0.0,h_pad=-1.)
# # plt.savefig('/home/kevin/Documents/Paper/NGC_2682_dr6vsdr61_teff.pdf')