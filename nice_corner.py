#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 09:45:27 2021

@author: kevin
"""
import emcee
import matplotlib.pyplot as plt
import corner
import numpy as np
from scipy.stats import kde
from astropy.io.votable import parse
from astropy.table import vstack
from matplotlib.lines import Line2D
from astropy.io import fits
from astropy.table import Table,vstack

from matplotlib import rc
rc('text', usetex=False)

# data=np.vstack(np.hstack((data[:,:21,:],data[:,22:24,:])))

# figure = corner.corner(reader.get_chain(flat=True)[10000:], labels=[r"teff", r"vbroad", r"log(g)", r"rv",r"[M]"],
#                        quantiles=[0.16, 0.5, 0.84],
#                        show_titles=True, title_kwargs={"fontsize": 12})

# name='prior'
# figure.suptitle(name)
name=160106001601282
cut=20000

large_data=fits.getdata('gaia_galah_dr6.fits',1)
large_data=Table(large_data)
old_abundances=[x for x in large_data if x['sobject_id']==str(name)]
old_abundances=vstack(old_abundances)

name='160106001601282_machine_learning'
reader= emcee.backends.HDFBackend(name+'_no_prior.h5')
data=reader.get_chain(flat=False)

off_value=np.max((np.mean(data,axis=0)-np.median(np.mean(data[500:,:,:],axis=0),axis=0))/np.median(np.mean(data[500:,:,:],axis=0),axis=0),axis=1)
deleted_times=0

for value,x in enumerate(off_value):
    if x>3000:
        data=np.delete(data,value-deleted_times,axis=1)
        deleted_times+=1
data=np.vstack(data)[cut:]



fig,a =plt.subplots(7,7)
fig.set_size_inches(20,20)
ticks=np.array(np.ones(7),dtype=object)


for x in range(7):
    for y in range(7):
        if x==y:
            a[x][x].hist(data[:,x],facecolor="none", 
              edgecolor='black',density=True)
            ticks[x]=np.percentile(data[:,x],[16,50,84])
            ticks[x][1]-=ticks[x][0]
            ticks[x][2]-=ticks[x][0]
        elif x<y:
            a[y][x].scatter(data[:,x],data[:,y],s=0.01,c='green')
            # Calculate the point density
            nbins=20
            data_temp=np.stack((data[:,x],data[:,y])).T
            k = kde.gaussian_kde(data_temp.T)
            xi, yi = np.mgrid[data[:,x].min():data[:,x].max():nbins*1j, data[:,y].min():data[:,y].max():nbins*1j]
            zi = k(np.vstack([xi.flatten(), yi.flatten()]))
            
            # a[y][x].pcolormesh(xi, yi, zi.reshape(xi.shape))
            a[y][x].contour(xi,yi,zi.reshape(xi.shape),levels=5,colors='black')
            # a[y][x].contour(data[:,x],data[:,y],z,levels=20)
        else:
            a[y][x].axis(False)


# reader= emcee.backends.HDFBackend('171027003801055_2limit5_priore10.h5')
# data=reader.get_chain(flat=False)

# off_value=np.max((np.mean(data,axis=0)-np.median(np.mean(data[500:,:,:],axis=0),axis=0))/np.median(np.mean(data[500:,:,:],axis=0),axis=0),axis=1)
# deleted_times=0

# for value,x in enumerate(off_value):
#     if x>1:
#         data=np.delete(data,value-deleted_times,axis=1)
#         deleted_times+=1
# data=np.vstack(data)[10000:]


# for x in range(7):
#     for y in range(7):
#         if x==y:
#             a[x][x].hist(data[:,x],facecolor="none", 
#               edgecolor='purple',density=True)
#             ticks[x]=np.percentile(data[:,x],[16,50,84])
#             ticks[x][1]-=ticks[x][0]
#             ticks[x][2]-=ticks[x][0]
#         elif x<y:
#             a[y][x].scatter(data[:,x],data[:,y],s=0.01,c='yellow')
#             # Calculate the point density
#             nbins=20
#             data_temp=np.stack((data[:,x],data[:,y])).T
#             k = kde.gaussian_kde(data_temp.T)
#             xi, yi = np.mgrid[data[:,x].min():data[:,x].max():nbins*1j, data[:,y].min():data[:,y].max():nbins*1j]
#             zi = k(np.vstack([xi.flatten(), yi.flatten()]))
            
#             # a[y][x].pcolormesh(xi, yi, zi.reshape(xi.shape))
#             a[y][x].contour(xi,yi,zi.reshape(xi.shape),levels=5,colors='purple')
#             # a[y][x].contour(data[:,x],data[:,y],z,levels=20)\
#         # else:
#         #     a[y][x].axis(False)



reader= emcee.backends.HDFBackend(name+'.h5')
data=reader.get_chain()

off_value=np.max((np.mean(data,axis=0)-np.median(np.mean(data[500:,:,:],axis=0),axis=0))/np.median(np.mean(data[500:,:,:],axis=0),axis=0),axis=1)

deleted_times=0
# for value,x in enumerate(off_value):
#     if x>5:
#         data=np.delete(data,value-deleted_times,axis=1)
#         deleted_times+=1
data=np.vstack(data)[cut:]

# data=np.vstack(np.hstack((data[:,:3,:],data[:,4:22,:],data[:,23:,:])))


ticks2=np.array(np.ones(7),dtype=object)
for x in range(7):
    for y in range(7):
        if x==y:
            a[x][x].hist(data[:,x],facecolor="none", 
              edgecolor='blue',density=True)
            ticks2[x]=np.percentile(data[:,x],[16,50,84])
            ticks2[x][1]-=ticks2[x][0]
            ticks2[x][2]-=ticks2[x][0]
            a[x][x].set_title((r'prior '+str(np.round(ticks[x][0],2))+'$^{+'+str(np.round(ticks[x][1],2))+'}_{-'+str(np.round(ticks[x][2],2))+'}$ no prior $'+str(np.round(ticks2[x][0],2))+'^{+'+str(np.round(ticks2[x][1],2))+'}_{-'+str(np.round(ticks2[x][2],2))+'}$' ))
            if x!=0:
                a[x][x].yaxis.tick_right()
            if x==0:
                teff=float(old_abundances['teff'])
                teff_e=float(old_abundances['e_teff'])
                xlim=a[x][x].get_xlim()
                x_line=np.linspace(teff-40,teff+40,50)
                y_line=1/(teff_e*np.sqrt(2*np.pi))*np.exp(-0.5*((x_line-teff)/teff_e)**2)
                a[x][x].plot(x_line,y_line,c='green')
            elif x==2:
                logg=float(old_abundances['logg'])
                logg_e=float(old_abundances['e_logg'])
                xlim=a[x][x].get_xlim()
                x_line=np.linspace(logg-0.2,logg+0.2,50)
                y_line=1/(logg_e*np.sqrt(2*np.pi))*np.exp(-0.5*((x_line-logg)/logg_e)**2)
                a[x][x].plot(x_line,y_line,c='green')
        elif x<y:
            a[y][x].scatter(data[:,x],data[:,y],s=0.01,c='red')
            # Calculate the point density
            nbins=20
            data_temp=np.stack((data[:,x],data[:,y])).T
            k = kde.gaussian_kde(data_temp.T)
            xi, yi = np.mgrid[data[:,x].min():data[:,x].max():nbins*1j, data[:,y].min():data[:,y].max():nbins*1j]
            zi = k(np.vstack([xi.flatten(), yi.flatten()]))
            
            # a[y][x].pcolormesh(xi, yi, zi.reshape(xi.shape))
            a[y][x].contour(xi,yi,zi.reshape(xi.shape),levels=5,colors='blue')
            # a[y][x].contour(data[:,x],data[:,y],z,levels=20)\

        else:
            a[y][x].axis(False)
            # print("skip")
        if y!=6 and x!=y:
            a[y][x].tick_params(labelbottom=False)

        if x!=0:
            a[y][x].tick_params( labelleft=False)
        if y==0:
            if x==1:
                a[x][y].set_ylabel('vbroad')
            elif x==2:
                a[x][y].set_ylabel('logg')
            elif x==3:
                a[x][y].set_ylabel('vrad')
            elif x==4:
                a[x][y].set_ylabel('monh')
            elif x==5:
                a[x][y].set_ylabel('vmac')

            elif x==6:
                a[x][y].set_ylabel('vmic')
        if x==6:
            if y==0:
                a[x][y].set_xlabel('teff')
            elif y==1:
                a[x][y].set_xlabel('vbroad')
            elif y==2:
                a[x][y].set_xlabel('logg')
            elif y==3:
                a[x][y].set_xlabel('vrad')
            elif y==4:
                a[x][y].set_xlabel('monh')
            elif y==5:
                a[x][y].set_xlabel('vmac')

plt.tight_layout(pad=0.5)
legend_elements = [Line2D([0], [0], color='blue', lw=2, label='no prior'),Line2D([0], [0],color='black', lw=2, label='with prior '),Line2D([0], [0], color='green', lw=2,label='prior')]


a[0][2].legend(handles=legend_elements,ncol=3,loc='center')

fig.suptitle(name)
fig.savefig(name)

# fig.show()
# p.show()

# a[0][1].scatter(data[])
