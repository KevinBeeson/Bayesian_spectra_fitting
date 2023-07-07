#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 14:10:24 2022

@author: kevin
"""

from astropy.io import fits
from astropy.table import Table,vstack,join
import matplotlib.pyplot as plt
import numpy as np
import copy
from pathlib import Path
import matplotlib.pyplot as plt
from astropy.io.votable import parse,from_table,writeto
from astropy.table import Table, vstack
from functools import partial
import  numpy as np 
import requests
import re
import csv
from os.path import exists
from pathlib import Path
from matplotlib.collections import LineCollection
from scipy.interpolate import interp1d
from multiprocessing import Pool
from scipy.optimize import curve_fit
import scipy.stats as sps
import copy
import random
from time import sleep

import seaborn
def get_isochrone(low_age,high_age,low_metalicty,high_metalicty,a_v,age_spacing=0.01,metalicty_spacing=0.01,iso_type='gaia'):
    """
    asks the padava server for an isochrone or a list of isochrones

    Parameters
    ----------
    low_age : float
        DESCRIPTION.
    high_age : float
        DESCRIPTION.
    low_metalicty : float
        DESCRIPTION.
    high_metalicty : float
        DESCRIPTION.
    a_v : float
        extinction.
    age_spacing : float, optional
        DESCRIPTION. The default is 0.01.
    metalicty_spacing : float, optional
        DESCRIPTION. The default is 0.01.

    an array of the isochrone
    -------
    None.

    """
    #if metallicity is given as a log ([M/H]) convert it here:


#    mass=[]
#    label=[]
#    mags=[]
#    imf=[]
#    teff=[]
#    logg=[]

    #parameters other than default
    d={
    'track_parsec': 'parsec_CAF09_v1.2S',
    'track_colibri':'parsec_CAF09_v1.2S_S35',
    'track_postagb':'no',
    'n_inTPC': 10,
    'eta_reimers': 0.2,
     # 'photsys_file': 'tab_mag_odfnew/tab_mag_gaiaDR2weiler.dat',
     'photsys_file':'tab_mag_odfnew/tab_mag_gaiaEDR3.dat',
   #   'photsys_file':'tab_mag_odfnew/tab_mag_gaia.dat',
  # 'photsys_file':'tab_mag_odfnew/tab_mag_2mass_spitzer_wise.dat',
    # 'photsys_file':'tab_mag_odfnew/tab_mag_panstarrs1.dat',
     #'photsys_file':'tab_mag_odfnew/tab_mag_gaiaDR2.dat',
    'photsys_version': 'OBC',
    'dust_sourceM': 'nodustM',
    'dust_sourceC': 'nodustC',
    'extinction_av': a_v,
    'extinction_coeff':'constant',
    'extinction_curve':'cardelli',   
    'imf_file': 'tab_imf/imf_kroupa_orig.dat',
    'isoc_isagelog':'1',
    'isoc_lagelow':low_age,
    'isoc_lageupp':high_age,
    'isoc_dlage':age_spacing, #steps ages
    'isoc_ismetlog':'1',
    'isoc_metlow':low_metalicty,
    'isoc_metupp':high_metalicty,
    'isoc_dmet':metalicty_spacing, #steps M/H
    'output_kind': 0,
    'submit_form': 'Submit'}
    if iso_type=='gaia':
        d['photsys_file']='tab_mag_odfnew/tab_mag_gaiaEDR3.dat'
    elif iso_type=='panstarrs':
        d['photsys_file']='tab_mag_odfnew/tab_mag_panstarrs1.dat'
    elif iso_type=='allwise':
        d['photsys_file']='tab_mag_odfnew/tab_mag_2mass_spitzer_wise.dat'
    #Check if we already downloaded this isochrone.
    #Isochrones are saved as txt files and the filename is the hash of the dictionary values.
    webserver = 'http://stev.oapd.inaf.it'
    c = requests.get(webserver + '/cgi-bin/cmd_3.7', params=d).text
#    print(c)
    aa = re.compile('output\d+')
    fname = aa.findall(c)
    if len(fname) > 0:
        url = '{0}/tmp/{1}.dat'.format(webserver, fname[0])
        #print url
        r = requests.get(url).text
        print('successfully gotten the ischrone')
    return r
def get_and_save_ishochrones(name,low_age,high_age,low_metalicty,high_metalicty,a_v,age_spacing=0.01,metalicty_spacing=0.01,special=False,iso_type='gaia'):
    """
    You can ask the function for a isochrone it will check if we already have it saved and if not it will ask from the padava servers

    Parameters
    ----------
    name : str
        name of the cluster.
    low_age : TYPE
        DESCRIPTION.
    high_age : TYPE
        DESCRIPTION.
    low_metalicty : TYPE
        DESCRIPTION.
    high_metalicty : TYPE
        DESCRIPTION.
    a_v : TYPE
        DESCRIPTION.
    age_spacing : TYPE, optional
        DESCRIPTION. The default is 0.01.
    metalicty_spacing : TYPE, optional
        DESCRIPTION. The default is 0.01.
    special : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    Path('ALL_ISO').mkdir(parents=True,exist_ok=True)
    #gets the path name of the isochrone and checks if its already been downloaded

    if (low_age==high_age or high_age==None) and (low_metalicty==high_metalicty or high_metalicty==None):
        if special:
            path_name='ALL_ISO/'+name+'_'+special+'_age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)
        else:
            path_name='ALL_ISO/age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)
    elif low_metalicty==high_metalicty or high_metalicty==None:
        path_name='ALL_ISO/age_'+str(low_age)+'_to_'+str(high_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)

    elif low_age==high_age or high_age==None:
        path_name='ALL_ISO/age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_to_'+str(high_metalicty)+'_extinction_'+str(a_v)


    else:
        
        if special:
            path_name='ALL_ISO/'+name+'_'+special+'_mixed'
        else:
            path_name='ALL_ISO/mixed'   
        if iso_type!='gaia':
            path_name+='_'+iso_type
        path_name+='.txt'
        #asks for the isochorne
        data=get_isochrone(low_age,high_age,low_metalicty,high_metalicty,a_v,age_spacing,metalicty_spacing,iso_type=iso_type)
        data=data.split('\n')
        header=data[11]
        data=data[0:len(data)-1]
        ammendedData=''
        for x in data:
            if x[0]!='#':
                numbers=np.array(x.split(' '))
                numbers=numbers[numbers!='']
                ammendedData+=' '.join(numbers)+'\n'
            else:
                numbers=np.array(x.split(' '))
                numbers=numbers[numbers!='']
                ammendedData+=' '.join(numbers)+'\n'
                
        ammendedData=ammendedData.split('\n')
        ammendedData=[i.split(' ') for i in ammendedData]
        ammendedData=ammendedData[0:len(ammendedData)-1]
        # ammendedData=np.array(ammendedData)
    
        toWrite=''
        for x in ammendedData:
            line=''
            for y in x:
                line+=str(y)+' '
            line=line[:-1]
            toWrite+=line+'\n'    
        file=open(path_name,'w')
        file.write(toWrite)
    if iso_type!='gaia':
        path_name+='_'+iso_type
    path_name+='.txt'
    if not exists(path_name):
        
        data=get_isochrone(low_age,high_age,low_metalicty,high_metalicty,a_v,iso_type=iso_type)
        data=data.split('\n')
        header=data[11]
        data=data[0:len(data)-1]
        ammendedData=''
        for x in data:
            if x[0]!='#':
                numbers=np.array(x.split(' '))
                numbers=numbers[numbers!='']
                ammendedData+=' '.join(numbers)+'\n'
            else:
                numbers=np.array(x.split(' '))
                numbers=numbers[numbers!='']
                ammendedData+=' '.join(numbers)+'\n'
                
        ammendedData=ammendedData.split('\n')
        ammendedData=[i.split(' ') for i in ammendedData]
        ammendedData=ammendedData[0:len(ammendedData)-1]
        # ammendedData=np.array(ammendedData)
    
        toWrite=''
        for x in ammendedData:
            line=''
            for y in x:
                line+=str(y)+' '
            line=line[:-1]
            toWrite+=line+'\n'    
        file=open(path_name,'w')
        file.write(toWrite)
def iso_reader(name,low_age,low_metalicty,a_v,special=False,high_age=None, high_metalicty=None,iso_type='gaia'):
    if (low_age==high_age or high_age==None) and (low_metalicty==high_metalicty or high_metalicty==None):
        if special:
            path_name='ALL_ISO/'+name+'_'+special+'_age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)
        else:
            path_name='ALL_ISO/age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)
        if iso_type!='gaia':
            path_name+='_'+iso_type
        path_name+='.txt'
        iso = list(csv.reader(open(path_name, 'rt'), delimiter=' '))
        
        iso=iso[13:-1]
        iso=np.array(iso,dtype=float)
        #crops the isochrone down to initial ini mass, int_imf, log_T, G, G_Bp_bright, G_BP_faint, G_RP ,log_g , cur mass
        if iso_type=='gaia':
            iso=np.column_stack((iso[:,3],iso[:,4],iso[:,7],iso[:,25],iso[:,26],iso[:,27],iso[:,8],iso[:,5]))
            
                
            iso=[x for x in iso if x[3]<15]
        return [low_age,low_metalicty],np.vstack(iso)
    elif low_metalicty==high_metalicty or high_metalicty==None:
        path_name='ALL_ISO/age_'+str(low_age)+'_to_'+str(high_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)+'.txt'

    elif low_age==high_age or high_age==None:
        path_name='ALL_ISO/age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_to_'+str(high_metalicty)+'_extinction_'+str(a_v)+'.txt'

    else:
        if special:
            path_name='ALL_ISO/'+name+'_'+special+'_mixed.txt'
        else:
            path_name='ALL_ISO/mixed.txt'
    iso = list(csv.reader(open(path_name, 'rt'), delimiter=' '))
    iso_full=[]
    iso_parameters=[]
    iso_temp=[]
    
    for x in iso:
        if x[0]!='#':
            iso_temp.append(x)
        else:
            if len(iso_temp):
                iso_temp=np.array(iso_temp,dtype=float)
                age_temp=iso_temp[0][2]
                metalicity_temp=iso_temp[0][1]

                iso_temp=np.column_stack((iso_temp[:,3],iso_temp[:,4],iso_temp[:,7],iso_temp[:,25],iso_temp[:,26],iso_temp[:,27],iso_temp[:,8],iso_temp[:,5]))
                iso_temp=[x for x in iso_temp if x[3]<15]

                iso_full.append(np.vstack(iso_temp))
                iso_parameters.append([age_temp,metalicity_temp])
                iso_temp=[]
    return iso_parameters,iso_full
def interpolate(oldIsochrone,scale=5):   
    mass=oldIsochrone[:,0]
    #changes log temperature to temp.
    oldIsochrone[:,2]=10**(oldIsochrone[:,2])

    
    newIso=[[None for y in range(len(oldIsochrone[0]))] for x in range((len(oldIsochrone)-1)*scale)]
    newIso=np.array(newIso)
    newIso[:,0]=np.hstack(np.transpose(np.linspace(mass[:-1],mass[1:],num=scale,endpoint=True)))
    newIso=np.array(newIso,dtype=float)
    for x in range(len(oldIsochrone[0])-1):
        f = interp1d(mass,oldIsochrone[:,x+1])
        newIso[:,x+1]=f(newIso[:,0]) 
    newIso[:,2]=np.log10(newIso[:,2])

    return np.array(newIso,dtype=float)
class isochrone:
    def __init__(self,name,age,metalicity,extinction,special=False,interpolate_scale=10,high_age=None,high_metalicity=None,limits=None,iso_type='gaia'):
        self.name=name
        self.special=special
        self.extinction=extinction
        self.iso_type=iso_type
        get_and_save_ishochrones(name,age,high_age,metalicity,high_metalicity,extinction,special=special,iso_type=iso_type)
        self.parameters,self.isochrone=iso_reader(name, age, metalicity, extinction,special=special,high_age=high_age,high_metalicty=high_metalicity,iso_type=iso_type)
        if iso_type=='allwise':
            iso_temp=isochrone(name,age,metalicity,extinction,special,interpolate_scale,high_age,high_metalicity,limits,iso_type='gaia')
            self.iso_gaia=iso_temp.isochrone
        if limits:
            if len(limits)!=2:
                raise ValueError("wrong shape for limits")
            self.limit_g=limits
        if limits:
            if isinstance(self.isochrone, np.ndarray):
                self.isochrone=np.vstack([x for x in self.isochrone if x[3]-limits[0]>0 and x[3]-limits[1]<0])  

            else:
                iso_temp=[]
                for x in self.isochrone:
                    iso_temp.append(np.vstack([x for x in isochrone if x[3]-limits[0]>0 and x[3]-limits[1]<0]))
                self.isochrone=iso_temp
        if isinstance(self.isochrone, np.ndarray) and iso_type=='gaia' and interpolate_scale:
            self.isochrone=interpolate(self.isochrone,interpolate_scale)
        elif iso_type=='gaia'and interpolate_scale:
            iso_temp=[]
            for x in self.isochrone:
                iso_temp.append(interpolate(x,interpolate_scale))
            self.isochrone=iso_temp
    def plot(self):
        plt.figure()
        if isinstance(self.isochrone,np.ndarray):
            plt.plot(self.isochrone[:,5]-self.isochrone[:,4],self.isochrone[:,3])
        else:
            for (x,label) in zip(self.isochrone,self.parameters):
                temp_label='age '+ str(label[0]) +' metalicty '+str(label[1])
            
                plt.plot(x[:,4]-x[:,5],x[:,3],label=temp_label)
                plt.legend(loc='best')
        plt.xlabel('bp_rp')
        plt.ylabel('G')
        plt.gca().invert_yaxis()
def inside_area(teff,logg,gradient,zero_point,upper=True):
    logg_line=teff*gradient+zero_point
    if logg_line>logg:
        out_side_area=True
    else:
        out_side_area=False
    if  upper:
        return out_side_area
    else:
        return not out_side_area

i_dr4=fits.getdata('galah_dr4_allspec_220713.fits',1)
i_dr4=Table(i_dr4)
mask = i_dr4['flag_sp']==0
i_dr4=i_dr4[mask]
mask =i_dr4['e_fe_h']!=np.inf
i_dr4=i_dr4[mask]

large_data_GALAH_official=fits.open('gaia_galah_cross_values.fits')
large_data_GALAH_official=Table(large_data_GALAH_official[1].data)
large_data_GALAH_official['sobject_id']=np.array(large_data_GALAH_official['sobject_id'],dtype='int64')
dr3_data=fits.getdata('GALAH_DR3_main_200331.fits',1)
dr3_data=Table(dr3_data)

cross_data=join(large_data_GALAH_official,dr3_data,keys_left='sobject_id',keys_right='sobject_id')

cross_data.sort('logg_2')

plt.figure()
plt.scatter(cross_data['teff_2'][::18],cross_data['logg_2'][::18])
parameters_names=['teff_2','logg_2','fe_h_2','vmic_2','vbroad_2']
elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
elem_mean=[]
elem_sig=[]
for elem in elements:
    elem_mean.append(np.nanmean(i_dr4[elem.lower()+'_fe']))
    elem_sig.append(np.sqrt(np.nanvar(i_dr4[elem.lower()+'_fe'])))
fe_h_mean=np.nanmean(i_dr4['fe_h'])
fe_h_sig=np.nanmean(i_dr4['e_fe_h'])
elem_max=[ 3.39219427e+00, 1.50000000e+00, 1.50000000e+00,
2.00000000e+00,0.7e+00,0.7e+00,1.2e+00,
1.0e+00,1.0e+00,1.0e+00,0.7e+00,
0.9e+0,0.7e+00,0.7e+00,1.0e+00,
0.8e+0,0.7e+00,0.7e+00,1.0e+00,
1.50000000e+00,1.0e+00, 1.00000000e+00,1.0e+00,
0.7e+0,1.0e+00, 1.20000000e+00, 1.30000000e+00,
1.30000000e+00, 1.20000000e+00,1.0e+00,1.0e+00]
elem_min=[-1.5e0,-1.0e+00,-1.5e+00,
-1.50000000e+00,-0.7e+00, -1.00000000e-0,-0.7e+00,
-5.00000000e-01,-1.0e+00,-0.7e+00,-0.7e+00,
-5.00000000e-01,-0.7e+00,-0.7e+00,-0.7e+00,
-0.7e+0,-0.7e+00,-0.7e+00,-1.0e+00,
-1.80000000e+00,-1.0e+00,-1.0e+00,-1.0e+00,
-1.0e+0, -2.00000000e+00,-1.0e+00,-0.7e+00,
-1.0e+0,-1.0e+00,-1.0e+00,-1.0e+00]
parameters=[100,0.5,0.5,5,5]
hard_limits_low=[3000,0.5,-1.5,0,0.1]
hard_limits_high=[8500,5,1.0,5,150]
np.random.seed(0)
cross_short=cross_data[::13]
param_to_compute=[]
# for star in cross_short:
#     star=cross_short[0]
#     param_temp=[star['teff_2'],star['logg_2'],np.random.normal(fe_h_mean,fe_h_sig*3,1)[0],star['vmic_2'],star['vbroad_2']]
#     param_temp=np.hstack((param_temp,np.zeros(31)))
#     if len (param_to_compute)==0:
#         param_to_compute=param_temp
#     param_to_compute=np.vstack((param_to_compute,param_temp))
#     for value,shift in enumerate(parameters):
#         to_shift=copy.copy(param_temp)
#         to_shift[value]+=shift
#         param_to_compute=np.vstack((param_to_compute,to_shift))
#         to_shift=copy.copy(param_temp)
#         to_shift[value]-=shift
#         if value==1:
#             if to_shift[value]<0.5:
#                 to_shift[value]=0.5
#         param_to_compute=np.vstack((param_to_compute,to_shift))
        
#     for value,(shift_max,shift_min) in enumerate(zip(elem_max,elem_min),len(parameters)):
#         to_shift=copy.copy(param_temp)
#         to_shift[value]=shift_max
#         param_to_compute=np.vstack((param_to_compute,to_shift))
#         to_shift=copy.copy(param_temp)
#         to_shift[value]=shift_min        
#         param_to_compute=np.vstack((param_to_compute,to_shift))
# fig, ax = plt.subplots(nrows=6, ncols=6,figsize=(16.0,9.0) )
# x=4
# for row in ax:
#     for col in row:
#         x+=1
#         if x<36:
#             # seaborn.kdeplot(param_to_compute[:,2],param_to_compute[:,x],ax=col)
#             col.scatter(param_to_compute[:,2],param_to_compute[:,x],s=0.1,c='red')
#             # col.set_xlabel('Fe/H')
#             col.set_ylabel(elements[x-5]+'/fe')
#             col.set_title(elements[x-5],y=0.7,x=0.1)
#         else:
#             col.set_axis_off()
# plt.tight_layout()
# fig.savefig('abundances distribution with limits')
logg_low=1.0
teff_low=4467
logg_high=3.2
teff_high=6000
gradient_1=(logg_high-logg_low)/(teff_high-teff_low)
zero_point_1=logg_high-teff_high*gradient_1

logg_low=1.14
logg_high=3.92
teff_low=3574
teff_high=4952
gradient_2=(logg_high-logg_low)/(teff_high-teff_low)
zero_point_2=logg_high-teff_high*gradient_2

param_to_compute=[]     
while len(param_to_compute)<2e4:
    for star in cross_data:
        print('loop_done')
        param_temp=[]
        for names,lower,upper,spread in zip(parameters_names,hard_limits_low,hard_limits_high,parameters):
            temp=-np.inf
            if 'vmic_2'==names:
                sigma=spread
            else:
                sigma=star['e_'+names]
            mean=star[names]

            if np.isnan(mean):
                mean=(upper+lower)/2
            if np.isnan(sigma):
                sigma=(upper-lower)/2
            while temp<lower or temp>upper:
                temp=np.random.normal(mean,sigma,1)
    
            param_temp.append(temp)
        for names,lower,upper,mean_zero,spread_zero in zip(elements,elem_min,elem_max,elem_mean,elem_sig):
            temp=-np.inf
            while temp<lower or temp>upper:
                if names=='N':
                    mean=0.0
                    sigma=0.5
                else:
                    mean=star[names+'_fe']
                    sigma=star['e_'+names+'_fe']
                if np.isnan(mean):
                    mean=mean_zero
                if mean<lower:
                    mean=lower
                elif mean>upper:
                    mean=upper
                if np.isnan(sigma):
                    sigma=spread_zero
                temp=-np.inf
                while temp<lower or temp>upper:
                    temp=np.random.normal(mean,sigma,1)
            param_temp.append(temp)
        while( (param_temp[0]>6000 and param_temp[1]<3.06) or
              (param_temp[0]<6000 and param_temp[1]<3.06 and inside_area(param_temp[0],param_temp[1],gradient_1,zero_point_1))or
               param_temp[0]<4952 and param_temp[1]<3.92 and inside_area(param_temp[0],param_temp[1],gradient_2,zero_point_2,upper=False)or
               param_temp[1]>5.21 or param_temp[0]>8300 or 
               (param_temp[0]>6500 and param_temp[1]>4.7)):
            # random_number=np.random.randint(0,len(cross_data))
            # star_temp=cross_data[random_number]
            mean=star['teff_2']
            sigma=star['e_teff_2']*4
            param_temp[0]=np.random.normal(mean,sigma,1)
            mean=star['logg_2']
            sigma=star['e_logg_2']*4
            param_temp[1]=np.random.normal(mean,sigma,1)
            param_temp=np.hstack(param_temp)
        # while param_temp[0]
        if len (param_to_compute)==0:

            param_to_compute=np.hstack(param_temp)
        param_to_compute=np.vstack((param_to_compute,np.hstack(param_temp)))
param_to_compute=np.array(param_to_compute)

parameters_all_names=['teff','logg','monh','vmic','vbroad','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
shift_all=[]
for x in param_to_compute:
     shift_temp={y:para for y,para in zip(parameters_all_names,x)}
     shift_all.append(shift_temp)


plt.figure()
plt.scatter(param_to_compute[:,0],param_to_compute[:,1],s=0.1,c='blue')
plt.scatter(cross_data['teff_2'],cross_data['logg_2'],s=0.3,c='red')

plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)

for name in cluster_details_all[:,0]:
    cluster_details=[x for x in cluster_details_all if x[0]==name][0]
    best_age=float(cluster_details[4])
    best_metalicty=float(cluster_details[5])
    best_extinction=float(cluster_details[6])

    isochrone_melotte=isochrone('Melotte_22',best_age,best_metalicty,best_extinction,high_age=best_age,high_metalicity=best_metalicty,interpolate_scale=1)
    plt.plot(10**isochrone_melotte.isochrone[:,2],isochrone_melotte.isochrone[:,-2],label=name)
plt.legend(loc='best')

# all_names=np.hstack((parameters_names,elements))
# fig, ax = plt.subplots(nrows=6, ncols=6,figsize=(16.0,9.0) )
# x=4
# for row in ax:
#     for col in row:
#         x+=1
#         if x<36:
#             seaborn.kdeplot(param_to_compute[:,2],param_to_compute[:,x],ax=col)
#             col.scatter(param_to_compute[:,2],param_to_compute[:,x],s=0.1,c='red')
#             # col.set_xlabel('Fe/H')
#             col.set_ylabel(elements[x-5]+'/fe')
#             col.set_title(elements[x-5],y=0.7,x=0.1)
#         else:
#             col.set_axis_off()
# plt.tight_layout()
# fig.savefig('abundances distribution')