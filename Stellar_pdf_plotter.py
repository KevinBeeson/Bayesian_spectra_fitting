#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 11:48:10 2022

@author: kevin
"""

from astropy.io.votable import parse,from_table,writeto
from matplotlib.collections import LineCollection
import scipy.stats as st
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kde
from pathlib import Path
import requests
import re
from os.path import exists
import csv
from scipy.interpolate import interp1d
from numba import jit
@jit(nopython=True,cache=True)
def normal(x,mu,sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-1/2*((x-mu)/sigma)**2)

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
        if isinstance(self.isochrone, np.ndarray) and iso_type=='gaia' and interpolate_scale!=1:
            self.isochrone=interpolate(self.isochrone,interpolate_scale)
        elif iso_type=='gaia'and interpolate_scale!=1:
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
class all_isochrones:
    def __init__(self,numbers,iso_bestfit,iso_low,iso_high,limits=None,gaussian=False):
        self.name=iso_bestfit.name
        self.iso_type=iso_bestfit.iso_type
        self.iso_bestfit=iso_bestfit
        self.iso_low=iso_low
        self.iso_high=iso_high
        self.extinction=iso_bestfit.extinction
        self.limits=limits
        if numbers!=0:
            self.difference_low=difference_iso(self.iso_low.isochrone,self.iso_bestfit.isochrone)
            self.difference_high=difference_iso(self.iso_high.isochrone,self.iso_bestfit.isochrone)
            if gaussian:
                metalicty_wanted=[]
                while len(metalicty_wanted)<numbers:
                    metalicty_temp=np.random.normal(iso_bestfit.parameters[1],abs(iso_bestfit.parameters[1]-iso_low.parameters[1]),1)
                    if metalicty_temp<iso_bestfit.parameters[1]:
                        metalicty_wanted.append(metalicty_temp)
                while len(metalicty_wanted)<numbers*2:
                    metalicty_temp=np.random.normal(iso_bestfit.parameters[1],abs(iso_bestfit.parameters[1]-iso_high.parameters[1]),1)
                    if metalicty_temp>iso_bestfit.parameters[1]:
                        metalicty_wanted.append(metalicty_temp)
                age_wanted=[]
                while len(age_wanted)<numbers:
                    age_temp=np.random.normal(iso_bestfit.parameters[0],abs(iso_bestfit.parameters[0]-iso_low.parameters[0]),1)
                    if age_temp<iso_bestfit.parameters[0]:
                        age_wanted.append(age_temp)
                while len(age_wanted)<numbers*2:
                    age_temp=np.random.normal(iso_bestfit.parameters[0],abs(iso_bestfit.parameters[0]-iso_high.parameters[0]),1)
                    if age_temp>iso_bestfit.parameters[0]:
                        age_wanted.append(age_temp)
                extinction_wanted=[]
                while len(extinction_wanted)<numbers:
                    extinction_temp=np.random.normal(iso_bestfit.extinction,abs(iso_bestfit.extinction-iso_low.extinction),1)
                    if extinction_temp<iso_bestfit.extinction:
                        extinction_wanted.append(extinction_temp)
                while len(extinction_wanted)<numbers*2:
                    extinction_temp=np.random.normal(iso_bestfit.extinction,abs(iso_bestfit.extinction-iso_high.extinction),1)
                    if extinction_temp>iso_bestfit.extinction:
                        extinction_wanted.append(extinction_temp)
                random.shuffle(extinction_wanted)
                random.shuffle(age_wanted)

                self.parameters=parameters=np.hstack((age_wanted,metalicty_wanted,extinction_wanted))
                isochrones_list=[]
                current_isochrone=0
                with Pool(processes=10) as pool:
                    results=pool.map(partial(iso_pool,self.extinction),parameters)
                self.isochrone=results
    
            else:
            # self.isochrone_list=[[iso_low.metalicty,iso_low.age],[iso_bestfit.]]
                # define the normal distribution and PDF
                difference_wanted=[]
                dist = sps.norm(loc=0, scale=self.difference_low)
                percentile_pdf=np.linspace(0.5,1,num=numbers)
                for i in percentile_pdf:
                    difference_wanted.append(dist.ppf(i))
                self.points=difference_wanted
                parameters=[[self.iso_bestfit.parameters[0],self.iso_bestfit.parameters[1],0],[self.iso_low.parameters[0],self.iso_low.parameters[1],self.difference_low]]
                isochrones_list=[]
                iso_param_have=np.copy(parameters[0])
                isochrones_list.append(iso_bestfit.isochrone)
                difference_iso_have=[0,self.difference_low,np.inf]
                while len(difference_wanted)-2:
                    print('there are '+ str(len(difference_wanted)+numbers/2)+' isocrhones left')
                    print(str(len(isochrones_list))+' isocrhones are done')
                    digitized = np.digitize(difference_wanted, difference_iso_have)
                    bin_amounts = [len(digitized[digitized == i]) for i in range(1, len(difference_iso_have))]
                    bin_max_index=bin_amounts.index(max(bin_amounts))
                    point_max_index=int(len([x for x in digitized if x<bin_max_index]))
                    point_max_index+=int(np.round(bin_amounts[bin_max_index]/2))
                    point_max_wanted=difference_wanted[point_max_index]
                    
                    if bin_max_index==len(difference_iso_have)-2:
                        bin_low=min(5,len(parameters))
                        parameters_closest=np.array(parameters[-bin_low:])
                        
                        trying_loop=True
                        while trying_loop:
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,0])
                            age_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,1])
                            metalicty_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
        
                            if abs(age_guess-parameters_closest[1][0])>0.5 or abs(metalicty_guess-parameters_closest[1][1])>0.2:
                                diff=(point_max_wanted-parameters[-1][2])/2
                                point_max_wanted-=diff
                            elif [age_guess,metalicty_guess] in np.array(parameters)[:,:2]:
                                bin_max_index=np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]
        
                                sign=np.sign(point_max_wanted-parameters[np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]][2])
                                ratio=point_max_wanted/parameters[np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]][2]
                                point_max_wanted*=ratio
                            else:
                                trying_loop=False
                    else:
                        loop=True
                        while loop:
                            bin_low=max(0,bin_max_index-2)
                            if bin_max_index+4<len(parameters):
                                bin_high=bin_max_index+4
                            else:
                                bin_high=len(parameters)
                            parameters_closest=np.array(parameters[bin_low:bin_high])
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,0])
                            age_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,1])
                            metalicty_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
                            if [age_guess,metalicty_guess] in np.array(parameters)[:,:2]:
                                bin_max_index=np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]
                                ratio=point_max_wanted/parameters[np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]][2]
                                point_max_wanted*=ratio
                            else:
                                loop=False
        
                        
                    if metalicty_guess>0.68:
                        metalicty_guess=0.68
                    elif metalicty_guess<-2.1:
                        metalicty_guess=-2.1
                    
                    if age_guess>10:
                        age_guess=10
                    if not ((metalicty_guess==0.68 or metalicty_guess==-2.1) and age_guess==10 and [age_guess,metalicty_guess] in np.array(parameters)[:,:2]):
                        iso_temp=isochrone('temp',age_guess,metalicty_guess,self.extinction,limits=limits,high_age=age_guess,high_metalicity=metalicty_guess)
                        difference_temp=difference_iso(iso_temp.isochrone,iso_bestfit.isochrone)
                        
                        difference_closest_index=np.where(abs(difference_wanted-difference_temp)==np.amin(abs(difference_wanted-difference_temp)))[0][0]
                        difference_closest=difference_wanted[difference_closest_index]
                        if abs(difference_closest-difference_temp)<self.difference_low/10:
                            del difference_wanted [difference_closest_index]
                            iso_param_have=np.vstack((iso_param_have,np.array([age_guess,metalicty_guess,difference_temp])))
                            isochrones_list.append(iso_temp.isochrone)
                            parameters.insert(bin_max_index+1,[age_guess,metalicty_guess,difference_temp])
                            difference_iso_have.insert(bin_max_index+1  ,difference_temp)
                        else:
                            parameters.insert(bin_max_index+1,[age_guess,metalicty_guess,difference_temp])
                            difference_iso_have.insert(bin_max_index+2,difference_temp)
                        parameters=sorted(parameters,key=lambda x:x[2])
                        difference_iso_have=sorted(difference_iso_have)
                    else:
                        del difference_wanted [point_max_index]
                    
                difference_wanted=[]
                dist = sps.norm(loc=0, scale=self.difference_high)
                percentile_pdf=np.linspace(0.5,1,num=numbers)
                for i in percentile_pdf:
                    difference_wanted.append(dist.ppf(i))
                self.points=difference_wanted
                
                parameters=[[self.iso_bestfit.parameters[0],self.iso_bestfit.parameters[1],0],[self.iso_high.parameters[0],self.iso_high.parameters[1],self.difference_high]]
        
                
                difference_iso_have=[0,self.difference_high,np.inf]
                while len(difference_wanted)-2:
                    print('there are '+ str(len(difference_wanted))+' isocrhones left')
                    print(str(len(isochrones_list))+' isocrhones are done')
                    digitized = np.digitize(difference_wanted, difference_iso_have)
                    bin_amounts = [len(digitized[digitized == i]) for i in range(1, len(difference_iso_have))]
                    bin_max_index=bin_amounts.index(max(bin_amounts))
                    point_max_index=int(len([x for x in digitized if x<bin_max_index]))
                    point_max_index+=int(np.round(bin_amounts[bin_max_index]/2))
                    point_max_wanted=difference_wanted[point_max_index]
                    
                    if bin_max_index==len(difference_iso_have)-2:
                        bin_low=min(5,len(parameters))
                        parameters_closest=np.array(parameters[-bin_low:])
                        
                        trying_loop=True
                        while trying_loop:
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,0])
                            age_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,1])
                            metalicty_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
        
                            if abs(age_guess-parameters_closest[1][0])>0.5 or abs(metalicty_guess-parameters_closest[1][1])>0.2:
                                diff=(point_max_wanted-parameters[-1][2])/2
                                point_max_wanted-=diff
                            elif [age_guess,metalicty_guess] in np.array(parameters)[:,:2]:
                                
                                bin_max_index=np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]
        
                                sign=np.sign(point_max_wanted-parameters[np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]][2])
                                ratio=point_max_wanted/parameters[np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]][2]
                                point_max_wanted*=ratio
                            else:
                                trying_loop=False
                    else:
                        loop=True
                        while loop:
                            bin_low=max(0,bin_max_index-2)
                            if bin_max_index+4<len(parameters):
                                bin_high=bin_max_index+4
                            else:
                                bin_high=len(parameters)
                            parameters_closest=np.array(parameters[bin_low:bin_high])
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,0])
                            age_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,1])
                            metalicty_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
                            if [age_guess,metalicty_guess] in np.array(parameters)[:,:2]:
                                bin_max_index=np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]
                                ratio=point_max_wanted/parameters[np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]][2]
                                point_max_wanted*=ratio
                            else:
                                loop=False
                    if metalicty_guess>0.68:
                        metalicty_guess=0.68
                    elif metalicty_guess<-2.1:
                        metalicty_guess=-2.1
                    
                    if age_guess>10:
                        age_guess=10
                    elif age_guess<1.5:
                        age_guess=1.5
                    if not ((metalicty_guess==0.68 or metalicty_guess==-2.1) and (age_guess==10 or age_guess==1.5) and [age_guess,metalicty_guess] in np.array(parameters)[:,:2]):
                        iso_temp=isochrone('temp',age_guess,metalicty_guess,self.extinction,limits=limits,high_age=age_guess,high_metalicity=metalicty_guess)
                        difference_temp=difference_iso(iso_temp.isochrone,iso_bestfit.isochrone)
                        
                        difference_closest_index=np.where(abs(difference_wanted-difference_temp)==np.amin(abs(difference_wanted-difference_temp)))[0][0]
                        difference_closest=difference_wanted[difference_closest_index]
                        if abs(difference_closest-difference_temp)<self.difference_low/10:
                            del difference_wanted [difference_closest_index]
                            iso_param_have=np.vstack((iso_param_have,np.array([age_guess,metalicty_guess,difference_temp])))
                            isochrones_list.append(iso_temp.isochrone)
                            parameters.insert(bin_max_index+1,[age_guess,metalicty_guess,difference_temp])
                            difference_iso_have.insert(bin_max_index+1  ,difference_temp)
                        else:
                            parameters.insert(bin_max_index+1,[age_guess,metalicty_guess,difference_temp])
                            difference_iso_have.insert(bin_max_index+2,difference_temp)
                        parameters=sorted(parameters,key=lambda x:x[2])
                        difference_iso_have=sorted(difference_iso_have)
                    else:
                        del difference_wanted [point_max_index]
                self.isochrone=isochrones_list
                self.parameters=iso_param_have
        def plot(self):
            plt.figure()
            if isinstance(self.isochrone,np.ndarray):
                plt.plot(self.isochrone[:,5]-self.isochrone[:,4],self.isochrone[:,3])
            else:
                for (x,label) in zip(self.isochrone,self.parameters):
                    temp_label='age '+ str(label[0]) +' metalicty '+str(label[1])
                
                    plt.plot(x[:,4]-x[:,5],x[:,3],label=temp_label)
                if len(self.isochrone)<5:
                    plt.legend(loc='best')
            plt.xlabel('bp_rp')
            plt.ylabel('G')
            plt.gca().invert_yaxis()
            # for x in         
            # while x>
        # def splitting(self,)
def plot_hr(isochrones,data,axis,save=False,paper=True,name=False,title=False):
    if not name:
        name=isochrones.name
    name=re.sub('_',' ', name)
    # plt.rc('font', size=15)
    iso_bestfit=isochrones.iso_bestfit.isochrone
    iso_low=isochrones.iso_low.isochrone
    iso_high=isochrones.iso_high.isochrone
    if isochrones.iso_type=='gaia':
        
        axis.plot(iso_bestfit[:,4]-iso_bestfit[:,5],iso_bestfit[:,3],label=r'Best fit',c='blue',linewidth=1.5,alpha=0.7)
        axis.plot(iso_low[:,4]-iso_low[:,5],iso_low[:,3],label=r'$\pm1 \sigma$',c='green',zorder=2,linewidth=1.5,alpha=0.7)
        axis.plot(iso_high[:,4]-iso_high[:,5],iso_high[:,3],c='green',zorder=2,linewidth=1.5,alpha=0.7)
        sigma=[2.5/np.log(10)/data['phot_g_mean_flux_over_error'],2.5/np.log(10)/data['phot_bp_mean_flux_over_error'],2.5/np.log(10)/data['phot_rp_mean_flux_over_error']]
        sigmaDCluster=5/np.log(10)*(data['r_sigma_cluster'])/data['r_est_cluster']
        # 
        errorYCluster=np.sqrt(sigmaDCluster**2+sigma[0]**2)
        errorX=np.sqrt(sigma[2]**2+sigma[1]**2)

        axis.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],s=0.3,zorder=0,c='black' )
        name_save=name
        if title:
            axis.set_title(title)
        if not paper:
            isochrone_fill=isochrones.isochrone
            for x in isochrone_fill:
                axis.plot(x[:,4]-x[:,5],x[:,3],alpha=0.1,c='blue',zorder=0)


    axis.invert_yaxis()
    axis.set_xlim((0.7,1.1))
    axis.set_ylim((5.0,1.8))
    
        
plt.rcParams['font.family']='Times New Roman'
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['text.usetex'] = True
# isochrone('temp',7.790753776549865,-0.02616414896675308,0.05,high_age=7.790753776549865,high_metalicity=-0.02616414896675308)
cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)
print(cluster_details_all[:,0])
# name=input('what isochrone do you want?')
name='NGC_2682'
cluster_details=[x for x in cluster_details_all if x[0]==name][0]
iso_type='gaia'

low_age=float(cluster_details[1])
low_metalicty=float(cluster_details[2])
low_extinction=float(cluster_details[3])
best_age=float(cluster_details[4])
best_metalicty=float(cluster_details[5])
best_extinction=float(cluster_details[6])
high_age=float(cluster_details[7])
high_metalicty=float(cluster_details[8])
high_extinction=float(cluster_details[9])


np.random.seed(0)
temperatureLimit=0.1
interpolate_scale=1
limit_g=None

iso_bestfit=isochrone(name,best_age,best_metalicty,best_extinction,special='Best_fit',limits=limit_g,high_age=best_age,high_metalicity=best_metalicty,iso_type=iso_type,interpolate_scale=interpolate_scale)


def inner_plotter(axes,star,corners,fig):

    left, bottom, width, height = axes
    ax2 = fig.add_axes([left, bottom, width, height])
    
    con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
                          xyB=(corners[0],corners[1]), coordsB=ax2.transAxes)
    fig.add_artist(con)
    con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
                          xyB=(corners[2],corners[3]), coordsB=ax2.transAxes)
    fig.add_artist(con)
    
    temperature_star=star['teff_raw']
    logg_star=star['logg_raw']
    ax2.scatter(temperature_star,logg_star,s=0.1)
    # ax2.errorbar(temperature_star,logg_star,xerr=teff_error, yerr=logg_error,fmt='o',ms=0.2,elinewidth=0.4,zorder=0)
    xx=star['teff_array']
    yy=star['logg_array']
    zz=star['probability_grid']
    ax2.contour(xx,yy,zz)
    ax2.set_xlabel(r'$T_{\rm{eff}}$')
    #ax2.set_ylabel(r'$\rm{log}(g)$') 
plt.rcParams['font.family']='Times New Roman'
    
from matplotlib.patches import ConnectionPatch
iso_bestfit=iso_bestfit.isochrone
name='NGC_2682_bp_rp_run.xml'
votable = parse(name)
data_small_run=votable.get_first_table().to_table(use_names_over_ids=True)
# plot_hr_teff(iso_all_age,data)
name=re.sub('_',' ', name)
Absmag=data_small_run['phot_g_mean_mag']-5*(np.log10(data_small_run['r_est_cluster'])-1)
data_small_run['abs_phot_g_mean_mag']=Absmag

data_small_run.sort(['abs_phot_g_mean_mag'])


plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'


plt.rcParams['ytick.labelsize']='medium'
plt.rcParams['xtick.labelsize']='medium'

fig=plt.figure(figsize=(7.0,4.0))
plt.rc('font',size=15)
plt.title(r'NGC 2682')
global ax
ax=fig.gca()
ax.tick_params(direction='in')
plt.plot(iso_bestfit[:,4]-iso_bestfit[:,5],iso_bestfit[:,3],label=r'Best fit',linewidth=1.5)
# points = np.array([iso_bestfit[:,4]-iso_bestfit[:,5],iso_bestfit[:,3]]).T.reshape(-1, 1, 2)
# segments = np.concatenate([points[:-1], points[1:]], axis=1)
# lc = LineCollection(segments, cmap='viridis', array=10**iso_bestfit[:,2],zorder=3)
# # Set the values used for colormapping
# lc.set_array(10**iso_bestfit[:-2,2])
 
# segments = np.concatenate([points[:-1], points[1:]], axis=1)
# # lc = LineCollection(segments, cmap='viridis', norm=norm)
# line = ax.add_collection(lc)

sigma=[2.5/np.log(10)/data_small_run['phot_g_mean_flux_over_error'],2.5/np.log(10)/data_small_run['phot_bp_mean_flux_over_error'],2.5/np.log(10)/data_small_run['phot_rp_mean_flux_over_error']]
sigmaDCluster=5/np.log(10)*(data_small_run['r_sigma_cluster'])/data_small_run['r_est_cluster']
# 
errorYCluster=np.sqrt(sigmaDCluster**2+sigma[0]**2)
errorX=np.sqrt(sigma[2]**2+sigma[1]**2)

limitsX=[0.53,max(data_small_run['bp_rp'])+0.1]
limitsY=[min(data_small_run['abs_phot_g_mean_mag']),max(data_small_run['abs_phot_g_mean_mag'])]


# plt.errorbar(data['bp_rp'],data['abs_phot_g_mean_mag'],xerr=errorX, yerr=errorYCluster,fmt='o',ms=1,elinewidth=1,zorder=3,c='black' )
plt.scatter(data_small_run['bp_rp'],data_small_run['abs_phot_g_mean_mag'],s=1.0,c='black')
# plt.clim(min(10**iso_bestfit[:-2,2]),max(10**iso_bestfit[:-2,2]))

plt.xlabel(r'$G_{\rm{BP}}-G_{\rm{RP}}$')
plt.ylabel(r'$M_{\rm{G}}$')
plt.xlim(limitsX)
plt.ylim(limitsY)

ax.set_xlim((0.5,2.4))

plt.gca().invert_yaxis()
plt.tight_layout()
plt.rcParams['ytick.labelsize']=8
plt.rcParams['xtick.labelsize']=8

# inner_plotter([0.15, 0.70, 0.2, 0.25],data_small_run[len(data_small_run)//30],[0,0,0,1],fig)
axes=[0.15, 0.67, 0.15, 0.20]
star=data_small_run[len(data_small_run)//30]
corners=[1,1,1,0]

left, bottom, width, height = axes
ax2 = fig.add_axes([left, bottom, width, height])

con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
                      xyB=(corners[0],corners[1]), coordsB=ax2.transAxes)
fig.add_artist(con)
con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
                      xyB=(corners[2],corners[3]), coordsB=ax2.transAxes)
fig.add_artist(con)

temperature_star=star['teff_raw']
logg_star=star['logg_raw']
ax2.scatter(temperature_star,logg_star,s=0.1,c='black')
# ax2.errorbar(temperature_star,logg_star,xerr=teff_error, yerr=logg_error,fmt='o',ms=0.2,elinewidth=0.4,zorder=0)

spread_temperature=np.sqrt(np.var(star['teff_raw']))
sig_teff=np.mean(star['e_teff_raw'])
spread_logg=np.sqrt(np.var(star['logg_raw']))
sig_logg=np.mean(star['e_logg_raw'])
spread_number=spread_temperature*spread_logg/(sig_logg*sig_teff*len(star['e_logg_raw']))
teff_raw=star['teff_raw']
logg_raw=star['logg_raw']
e_teff_raw=star['e_teff_raw']
e_logg_raw=star['e_logg_raw']
if spread_number>1:
    k=st.gaussian_kde([teff_raw,logg_raw])
    # teff_line=np.linspace(min(teff_raw)-50,max(logg_raw)+50,num=50)
    # logg_line=np.linspace(min(logg_raw)-0.25,max(logg_raw)+0.25,num=50)
    xx, yy = np.mgrid[min(teff_raw)-10:max(teff_raw)+10:100j, min(logg_raw)-0.1:max(logg_raw)+0.1   :100j]
    
    zz = k(np.vstack([xx.flatten(), yy.flatten()]))
    zz = np.reshape(zz, xx.shape)
    ax2.contour(xx,yy,zz)
else:
    xx, yy = np.mgrid[min(teff_raw)-10:max(teff_raw)+10:100j, min(logg_raw)-0.1:max(logg_raw)+0.1   :100j]

    prior_parameters=np.column_stack((teff_raw,e_teff_raw,logg_raw,e_logg_raw))
    
    zz=[]
    for x,y in zip(xx,yy):
        for temp_param in prior_parameters:           
            zz.append(normal(x,temp_param[0],temp_param[1])*normal(y,temp_param[2],temp_param[3]))
            
# ax2.set_xlabel(r'$T_{\rm{eff}}$')
#ax2.set_ylabel(r'$\rm{log}(g)$') 
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()    

axes=[0.15, 0.25, 0.15, 0.20]
star=data_small_run[140]
corners=[0,1,1,1]


left, bottom, width, height = axes
ax2 = fig.add_axes([left, bottom, width, height])

con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
                      xyB=(corners[0],corners[1]), coordsB=ax2.transAxes)
fig.add_artist(con)
con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
                      xyB=(corners[2],corners[3]), coordsB=ax2.transAxes)
fig.add_artist(con)

temperature_star=star['teff_raw']
logg_star=star['logg_raw']
ax2.scatter(temperature_star,logg_star,s=0.1,c='black')
# ax2.errorbar(temperature_star,logg_star,xerr=teff_error, yerr=logg_error,fmt='o',ms=0.2,elinewidth=0.4,zorder=0)

spread_temperature=np.sqrt(np.var(star['teff_raw']))
sig_teff=np.mean(star['e_teff_raw'])
spread_logg=np.sqrt(np.var(star['logg_raw']))
sig_logg=np.mean(star['e_logg_raw'])
spread_number=spread_temperature*spread_logg/(sig_logg*sig_teff*len(star['e_logg_raw']))
teff_raw=star['teff_raw']
logg_raw=star['logg_raw']
e_teff_raw=star['e_teff_raw']
e_logg_raw=star['e_logg_raw']
if spread_number>1:
    k=st.gaussian_kde([teff_raw,logg_raw])
    # teff_line=np.linspace(min(teff_raw)-50,max(logg_raw)+50,num=50)
    # logg_line=np.linspace(min(logg_raw)-0.25,max(logg_raw)+0.25,num=50)
    xx, yy = np.mgrid[min(teff_raw)-10:max(teff_raw)+10:100j, min(logg_raw)-0.1:max(logg_raw)+0.1   :100j]
    
    zz = k(np.vstack([xx.flatten(), yy.flatten()]))
    zz = np.reshape(zz, xx.shape)
    ax2.contour(xx,yy,zz)
else:
    xx, yy = np.mgrid[min(teff_raw)-10:max(teff_raw)+10:100j, min(logg_raw)-0.1:max(logg_raw)+0.1   :100j]

    prior_parameters=np.column_stack((teff_raw,e_teff_raw,logg_raw,e_logg_raw))
    
    zz=[]
    for x,y in zip(xx,yy):
        for temp_param in prior_parameters:           
            zz.append(normal(x,temp_param[0],temp_param[1])*normal(y,temp_param[2],temp_param[3]))



ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()    

axes=[0.50, 0.55, 0.15, 0.20]
star=data_small_run[900]
corners=[0,0,1,0]


left, bottom, width, height = axes
ax2 = fig.add_axes([left, bottom, width, height])

con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
                      xyB=(corners[0],corners[1]), coordsB=ax2.transAxes)
fig.add_artist(con)
con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
                      xyB=(corners[2],corners[3]), coordsB=ax2.transAxes)
fig.add_artist(con)

temperature_star=star['teff_raw']
logg_star=star['logg_raw']
ax2.scatter(temperature_star,logg_star,s=0.1,c='black')
# ax2.errorbar(temperature_star,logg_star,xerr=teff_error, yerr=logg_error,fmt='o',ms=0.2,elinewidth=0.4,zorder=0)

spread_temperature=np.sqrt(np.var(star['teff_raw']))
sig_teff=np.mean(star['e_teff_raw'])
spread_logg=np.sqrt(np.var(star['logg_raw']))
sig_logg=np.mean(star['e_logg_raw'])
spread_number=spread_temperature*spread_logg/(sig_logg*sig_teff*len(star['e_logg_raw']))
teff_raw=star['teff_raw']
logg_raw=star['logg_raw']
e_teff_raw=star['e_teff_raw']
e_logg_raw=star['e_logg_raw']
if spread_number>1:
    k=st.gaussian_kde([teff_raw,logg_raw])
    # teff_line=np.linspace(min(teff_raw)-50,max(logg_raw)+50,num=50)
    # logg_line=np.linspace(min(logg_raw)-0.25,max(logg_raw)+0.25,num=50)
    xx, yy = np.mgrid[min(teff_raw)-10:max(teff_raw)+10:100j, min(logg_raw)-0.1:max(logg_raw)+0.1   :100j]
    
    zz = k(np.vstack([xx.flatten(), yy.flatten()]))
    zz = np.reshape(zz, xx.shape)
    ax2.contour(xx,yy,zz)
else:
    xx, yy = np.mgrid[min(teff_raw)-10:max(teff_raw)+10:100j, min(logg_raw)-0.1:max(logg_raw)+0.1   :100j]

    prior_parameters=np.column_stack((teff_raw,e_teff_raw,logg_raw,e_logg_raw))
    
    zz=[]
    for x,y in zip(xx,yy):
        temp_probability=0
        for temp_param in prior_parameters: 
            temp_probability+=normal(x,temp_param[0],temp_param[1])*normal(y,temp_param[2],temp_param[3])        
        zz.append(temp_probability)
    zz=np.reshape(zz,xx.shape)
    ax2.contour(xx,yy,zz)
ax2.set_xlim((4150,4350))
ax2.set_ylim((4.63,4.66))

ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()    

axes=[0.74, 0.38, 0.15, 0.20]
star=data_small_run[1200]
corners=[0,0,1,0]


left, bottom, width, height = axes
ax2 = fig.add_axes([left, bottom, width, height])

con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
                      xyB=(corners[0],corners[1]), coordsB=ax2.transAxes)
fig.add_artist(con)
con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
                      xyB=(corners[2],corners[3]), coordsB=ax2.transAxes)
fig.add_artist(con)

temperature_star=star['teff_raw']
logg_star=star['logg_raw']
ax2.scatter(temperature_star,logg_star,s=0.1,c='black')
# ax2.errorbar(temperature_star,logg_star,xerr=teff_error, yerr=logg_error,fmt='o',ms=0.2,elinewidth=0.4,zorder=0)

spread_temperature=np.sqrt(np.var(star['teff_raw']))
sig_teff=np.mean(star['e_teff_raw'])
spread_logg=np.sqrt(np.var(star['logg_raw']))
sig_logg=np.mean(star['e_logg_raw'])
spread_number=spread_temperature*spread_logg/(sig_logg*sig_teff*len(star['e_logg_raw']))
teff_raw=star['teff_raw']
logg_raw=star['logg_raw']
e_teff_raw=star['e_teff_raw']
e_logg_raw=star['e_logg_raw']
if spread_number>1:
    k=st.gaussian_kde([teff_raw,logg_raw])
    # teff_line=np.linspace(min(teff_raw)-50,max(logg_raw)+50,num=50)
    # logg_line=np.linspace(min(logg_raw)-0.25,max(logg_raw)+0.25,num=50)
    xx, yy = np.mgrid[min(teff_raw)-10:max(teff_raw)+10:100j, min(logg_raw)-0.1:max(logg_raw)+0.1   :100j]
    
    zz = k(np.vstack([xx.flatten(), yy.flatten()]))
    zz = np.reshape(zz, xx.shape)
    ax2.contour(xx,yy,zz)
else:
    xx, yy = np.mgrid[min(teff_raw)-10:max(teff_raw)+10:100j, min(logg_raw)-0.1:max(logg_raw)+0.1   :100j]

    prior_parameters=np.column_stack((teff_raw,e_teff_raw,logg_raw,e_logg_raw))
    
    zz=[]
    for x,y in zip(xx,yy):
        temp_probability=0
        for temp_param in prior_parameters: 
            temp_probability+=normal(x,temp_param[0],temp_param[1])*normal(y,temp_param[2],temp_param[3])        
        zz.append(temp_probability)
    zz=np.reshape(zz,xx.shape)
    ax2.contour(xx,yy,zz)
ax2.set_ylim((4.67,4.75))

plt.tight_layout()
fig.savefig('NGC_2682_BP_RP_full_stellar_parameters.pdf')


# # ax2.set_xlabel(r'$T_{\rm{eff}}$')
# #ax2.set_ylabel(r'$\rm{log}(g)$') 
# ax2.yaxis.set_label_position("right")
# ax2.yaxis.tick_right()    
  


# axes=[0.48, 0.48, 0.15, 0.20]
# star=data_small_run[len(data_small_run)*5//10]
# corners=[0,1,0,0]


# left, bottom, width, height = axes
# ax2 = fig.add_axes([left, bottom, width, height])

# con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
#                       xyB=(corners[0],corners[1]), coordsB=ax2.transAxes)
# fig.add_artist(con)
# con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
#                       xyB=(corners[2],corners[3]), coordsB=ax2.transAxes)
# fig.add_artist(con)

# temperature_star=star['teff_raw']
# logg_star=star['logg_raw']
# ax2.scatter(temperature_star,logg_star,s=0.1)
# # ax2.errorbar(temperature_star,logg_star,xerr=teff_error, yerr=logg_error,fmt='o',ms=0.2,elinewidth=0.4,zorder=0)

# spread_temperature=np.sqrt(np.var(star['teff_raw']))
# sig_teff=np.mean(star['e_teff_raw'])
# spread_logg=np.sqrt(np.var(star['logg_raw']))
# sig_logg=np.mean(star['e_logg_raw'])
# spread_number=spread_temperature*spread_logg/(sig_logg*sig_teff*len(star['e_logg_raw']))
# teff_raw=star['teff_raw']
# logg_raw=star['logg_raw']
# e_teff_raw=star['e_teff_raw']
# e_logg_raw=star['e_logg_raw']
# if spread_number>1:
#     k=st.gaussian_kde([teff_raw,logg_raw])
#     # teff_line=np.linspace(min(teff_raw)-50,max(logg_raw)+50,num=50)
#     # logg_line=np.linspace(min(logg_raw)-0.25,max(logg_raw)+0.25,num=50)
#     xx, yy = np.mgrid[min(teff_raw)-10:max(teff_raw)+10:100j, min(logg_raw)-0.1:max(logg_raw)+0.1   :100j]
    
#     zz = k(np.vstack([xx.flatten(), yy.flatten()]))
#     zz = np.reshape(zz, xx.shape)
#     ax2.contour(xx,yy,zz)
# else:
#     xx, yy = np.mgrid[min(teff_raw)-10:max(teff_raw)+10:100j, min(logg_raw)-0.1:max(logg_raw)+0.1   :100j]

#     prior_parameters=np.column_stack((teff_raw,e_teff_raw,logg_raw,e_logg_raw))
    
#     zz=[]
#     for x,y in zip(xx,yy):
#         for temp_param in prior_parameters:           
#             zz.append(normal(x,temp_param[0],temp_param[1])*normal(y,temp_param[2],temp_param[3]))
            
# # ax2.set_xlabel(r'$T_{\rm{eff}}$')
# #ax2.set_ylabel(r'$\rm{log}(g)$') 
# ax2.yaxis.set_label_position("right")
# ax2.yaxis.tick_right()    


# axes=[0.7, 0.40, 0.2, 0.25]
# star=data_small_run[len(data_small_run)*7//10]
# corners=[0,0,1,0]


# left, bottom, width, height = axes
# ax2 = fig.add_axes([left, bottom, width, height])

# con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
#                       xyB=(corners[0],corners[1]), coordsB=ax2.transAxes)
# fig.add_artist(con)
# con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
#                       xyB=(corners[2],corners[3]), coordsB=ax2.transAxes)
# fig.add_artist(con)

# temperature_star=star['teff_raw']
# logg_star=star['logg_raw']
# ax2.scatter(temperature_star,logg_star,s=0.1)
# # ax2.errorbar(temperature_star,logg_star,xerr=teff_error, yerr=logg_error,fmt='o',ms=0.2,elinewidth=0.4,zorder=0)
# xx=star['teff_array']
# yy=star['logg_array']
# zz=star['probability_grid']
# ax2.contour(xx,yy,zz)
# ax2.set_xlabel(r'$T_{\rm{eff}}$',size=8)
# ax2.set_ylabel(r'$\rm{log}(g)$',size=8)
# ax2.xaxis.set_label_position("top")
# ax2.xaxis.tick_top()    

# # inner_plotter([0.05, 0.10, 0.2, 0.25],data_small_run[len(data_small_run)//10],[0,1,1,1],fig)

# # # inner_plotter([0.1, 0.70, 0.2, 0.25],data_small_run[len(data_small_run)//10+1],[0,0,1,0],fig)
# # # inner_plotter([0.10, 0.170, 0.17, 0.21],data_small_run[len(data_small_run)*3//10],[1,0,1,1],fig)
# # # inner_plotter([0.80, 0.270, 0.17, 0.21],data_small_run[len(data_small_run)*4//10],[1,0,0,0],fig)
# # inner_plotter([0.40, 0.400, 0.2, 0.25],data_small_run[len(data_small_run)*5//10],[0,0,1,0],fig)

# # plt.savefig(name+' pdf examples.pdf',dpi=100,format='pdf')
# # # inner_plotter( [0.25, 0.65,0.15, 0.2],data[len(data)//10],[1,0,0,0],fig)

# # # inner_plotter( [0.35, 0.55,0.15, 0.2],data_small_run[len(data_small_run)*6//10],[1,0,0,0],fig)

# # # inner_plotter( [0.15, 0.55,0.15, 0.2],data[2],[1,0,0,0])

# # # inner_plotter( [0.45, 0.20, 0.15, 0.2],data_small_run[len(data_small_run)*7//10],[0,1,1,1],fig)
# # inner_plotter( [0.65, 0.35, 0.2, 0.25],data_small_run[len(data_small_run)*8//10],[1,0,0,0],fig)
# # # inner_plotter( [0.65, 0.25, 0.15, 0.2],data[len(data)*5//10],[1,0,0,0],fig)
# # # inner_plotter( [0.75, 0.20, 0.15, 0.2],data[len(data)*6//10],[0,1,1,1],fig)
# # # inner_plotter( [0.85, 0.35, 0.15, 0.2],data[len(data)*7//10],[1,0,0,0],fig)
# # # inner_plotter( [0.95, 0.25, 0.15, 0.2],data[len(data)*8//10],[1,0,0,0],fig)

# fig.tight_layout()

# # iso_low=isochrones.iso_low.isochrone
# # iso_high=isochrones.iso_high.isochrone
# fig=plt.figure()
# ax=fig.gca()
# plt.scatter(data_small_run['bp_rp'],data_small_run['abs_phot_g_mean_mag'],c=data_small_run['teff'])
# numberDensity=np.subtract(iso_bestfit[1:,1],iso_bestfit[:-1,1])
# np.append(numberDensity,numberDensity[-1])
# numberDensity*=1e2
# points = np.array([iso_bestfit[:,4]-iso_bestfit[:,5],iso_bestfit[:,3]]).T.reshape(-1, 1, 2)
# segments = np.concatenate([points[:-1], points[1:]], axis=1)
# lc = LineCollection(segments, cmap='viridis', array=10**iso_bestfit[:,2],zorder=3)
# # Set the values used for colormapping
# lc.set_array(10**iso_bestfit[:-2,2])
 
# segments = np.concatenate([points[:-1], points[1:]], axis=1)
# # lc = LineCollection(segments, cmap='viridis', norm=norm)
# line = ax.add_collection(lc)
# plt.colorbar(line, label='teff')
# plt.clim(vmin=min(10**iso_bestfit[:-2,2]),vmax=max(10**iso_bestfit[:-2,2]))
# # plt.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],color='b',s=2)
# plt.gca().invert_yaxis()
# plt.tight_layout()
# # fig.savefig('NGC_2682_BP_RP_teff.pdf')

# fig=plt.figure()
# ax=fig.gca()
# plt.scatter(data_small_run['bp_rp'],data_small_run['abs_phot_g_mean_mag'],c=data_small_run['logg'])
# numberDensity=np.subtract(iso_bestfit[1:,1],iso_bestfit[:-1,1])
# np.append(numberDensity,numberDensity[-1])
# numberDensity*=1e2
# points = np.array([iso_bestfit[:,4]-iso_bestfit[:,5],iso_bestfit[:,3]]).T.reshape(-1, 1, 2)
# segments = np.concatenate([points[:-1], points[1:]], axis=1)
# lc = LineCollection(segments, cmap='viridis', array=iso_bestfit[:,6],zorder=3)
# # Set the values used for colormapping
# lc.set_array(iso_bestfit[:,6])
 
# segments = np.concatenate([points[:-1], points[1:]], axis=1)
# # lc = LineCollection(segments, cmap='viridis', norm=norm)
# line = ax.add_collection(lc)
# plt.colorbar(line, label='logg')
# plt.clim(vmin=min(iso_bestfit[:,6]),vmax=max(iso_bestfit[:,6]))
# # plt.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],color='b',s=2)
# plt.gca().invert_yaxis()
# plt.tight_layout()
# # fig.savefig('NGC_2682_BP_RP_logg.pdf')

