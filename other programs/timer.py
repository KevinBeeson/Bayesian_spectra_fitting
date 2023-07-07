#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 14:17:31 2022

@author: kevin
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 13:39:38 2021

@author: kevin
"""
from numba import jit
import time
from functools import  partial
from astropy.io.votable import parse
import emcee
from The_Payne import spectral_model
import scipy
from scipy import signal
from os.path import exists
import subprocess
from pathlib import Path

import os.path
import logging
from astropy.table import Table,vstack
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
from astropy.io import fits
import copy
import Galah_tool_py3 as gtools
import os
import random
import functools
from multiprocessing.dummy import Pool as ThreadPool 
from The_Payne import utils

import warnings

#ignore by message
warnings.filterwarnings("ignore", message="Polyfit may be poorly conditioned")
warnings.filterwarnings("ignore")

logger=logging.getLogger('pysme')
logger.setLevel(-100)

os.environ["OMP_NUM_THREADS"] = "1"

# def synthesize_multi(move,values,changed_values):
#     shift={values:move}
#     shift.update(changed_values)
#     return synthesize(shift)
def how_good_gaussian_fit(data):
    diff=0
    for x in range(np.shape(data)[1]):
        temp_data=data[:,x]
        mean=np.mean(temp_data)
        standard_deviation=np.std(temp_data)
        counts, bins = np.histogram(temp_data,density=True)
        bins=(bins[1:]+bins[:-1])/2
        gaussian=1/(standard_deviation*np.sqrt(2*np.pi))*np.exp(-0.5*((mean-bins)/standard_deviation)**2)
        gaussian, counts=getting_rid_zeros(gaussian,counts)
        diff+=np.sqrt(np.sum((gaussian-counts)**2/((counts+gaussian)/2)**2)/len(gaussian))
    return diff/np.shape(data)[1]
def getting_rid_zeros(data1,data2):
        data1_temp=[]
        data2_temp=[]
        for synthetic,data_bins in zip(data1,data2):
            if (abs(synthetic)+abs(data_bins)>1e-100).all():
                data1_temp.append(synthetic)
                data2_temp.append(data_bins)
        return np.array(data1_temp),np.array(data2_temp)
def shift_maker(solar_values):
    # shift={'vrad':solar_values}
    shift={'teff':solar_values[0],'logg':solar_values[1],'monh':solar_values[2],'Fe':solar_values[3],'alpha':solar_values[4],'vrad':solar_values[5],'vsini':solar_values[6],'vmac':solar_values[7],'vmic':solar_values[8]}

    return shift

def log_prior(shift):
    error=0
    for x in shift:
        if x in ('logg','teff'):
            error-=(shift[x]-float(old_abundances[x]))**2/(2*(float(old_abundances['e_'+x]))**2)
           # error-=(shift[x]-photometric_values[x])**2/(2*(photometric_values[x+'_e']*10)**2)
        # elif x=='vrad':
        #     error-=(shift[x]-float(old_abundances['rv_galah'])**2/(2*(float(old_abundances['e_rv_galah']))**2))
    if error<-1e100:
        return -np.inf
    else:
        return error
def log_posterior(solar_values,prior=True):
# pos=np.column_stack((teff_random,vsini_random,logg_random,v_rad_random,monh_random))
    if not starting_test(solar_values):
         return -np.inf
    shift=shift_maker(solar_values)
    try:
        synthetic_spectras=spectras.synthesize(shift,give_back=True)
        # print('shift' , shift)
        # print("Spectra",synthetic_spectras)
    except:
        print('couldnt finish synthesizing')
        return -np.inf
    # if not np.all([np.all(x) for x in synthetic_spectras]):
    #     print("Didnt synthesize properly")
    #     return -np.inf
    # if np.max([np.max(abs(x)) for x in synthetic_spectras])>100 or np.max([np.max(abs(x)) for x in synthetic_spectras])==np.nan:
    #     print("Didnt synthesize properly one of the spectra blew up")
    #     return -np.inf
    # if np.any(np.array_equal(synthetic_spectras,np.nan)):
    #     print("Didnt synthesize properly Nan values yay!")
    #     return -np.inf
    # print("Before normalized",synthetic_spectras, "Shift    ",shift)
    # print("max value",np.max([np.max(abs(x)) for x in synthetic_spectras]),"shift",shift)
    # normalized_spectra=spectras.normalize(data=synthetic_spectras)
    probability=spectras.log_fit(synthetic_spectra=synthetic_spectras)
    # print("This is the probability", probability, "shift prob", log_prior(shift), "shift", shift)
    # print("This is the Uncs", normalized_uncertainty,"This is the Observed",normalized_spectra)
    # print("Spectra",synthetic_spectras)
    if prior:
        probability+=log_prior(shift)
    return probability
def starting_test(solar_values):
    
    shift={'teff':solar_values[0],'logg':solar_values[1],'monh':solar_values[2],'Fe':solar_values[3],'alpha':solar_values[4],'vrad':solar_values[5],'vsini':solar_values[6],'vmac':solar_values[7],'vmic':solar_values[8]}
    if np.any([True for x in shift if shift[x]<0.0 and  x!='monh' and x!='Fe' and x!='alpha']):
        return False
    
    if 'teff' in shift:
         if abs(shift['teff']-old_abundances['teff_reduction'])>1000 and old_abundances['teff_reduction']<10000 and old_abundances['teff_reduction']>3500:
                print("temperature is too small or large")
                return False
         elif shift['teff']<3500 or shift['teff']>10000:
            return False
    if 'logg' in shift:
        if shift['logg']<1.0 or shift['logg']>5.0:
            print("logg is too small or big")
            return False  
    if 'monh' in shift:
        if abs(shift['monh']-old_abundances['fe_h'])>0.5 or abs(shift['monh'])>5.0:
            print("monh is too small or big")
            return False   
    if 'vsini' in shift:
        if shift['vsini']<0 or shift['vsini']>100:
            print("vsini is too small or big")
            return False       
    if 'vmic' in shift:
        if shift['vmic']<0 or shift['vmic']>20:
            return False 

            print("vmic is negative")
    if 'vmac' in shift:
        if shift['vmac']<0 or shift['vmac']>100:
            print("vmac is negative")
            return False 
    return True
def rsetattr(obj, attr, val):
    pre, _, post = attr.rpartition('.')
    return setattr(rgetattr(obj, pre) if pre else obj, post, val)

def rgetattr(obj, attr, *args):
    def _getattr(obj, attr):
        return getattr(obj, attr, *args)
    return functools.reduce(_getattr, [obj] + attr.split('.'))
def runGetattr(obj, attr):
    def _getattr(obj, attr):
        return getattr(obj, attr)
    return functools.reduce(_getattr, [obj] + attr.split('.'))
def payne_sythesize(solar_values,x_min,x_max,NN_coeffs):

    scaled_labels = (solar_values-x_min)/(x_max-x_min) - 0.5

    real_spec = spectral_model.get_spectrum_from_neural_net(scaled_labels = scaled_labels, NN_coeffs = NN_coeffs)
    return real_spec
def galah_kern(fwhm, b):
    """ Returns a normalized 1D kernel as is used for GALAH resolution profile """
    size=2*(fwhm/2.355)**2
    size_grid = int(size) # we limit the size of kernel, so it is as small as possible (or minimal size) for faster calculations
    if size_grid<7: size_grid=7
    x= scipy.mgrid[-size_grid:size_grid+1]
    g = scipy.exp(-0.693147*np.power(abs(2*x/fwhm), b))
    return g / np.sum(g)
def dopler(original_wavelength,spectra,v_rad,crop=False):
    c=299792.458
    delta_v=v_rad/c
    if len(original_wavelength)!=len(spectra):
        difference=abs(len(original_wavelength)-len(spectra))
        original_wavelength_strenched=np.interp(np.linspace(-difference/2,difference/2+len(original_wavelength),num=len(spectra)),range(len(original_wavelength)),original_wavelength)
        observed_wavelength=original_wavelength_strenched*(1+delta_v)
    else:
        observed_wavelength=original_wavelength*(1+delta_v)

    spectra_new=np.interp(original_wavelength,observed_wavelength,spectra)
    
    return spectra_new
@jit(nopython=True,cache=True)
def numba_syth_resolution(coefficients,l_new,sampl,min_sampl,last_frequency):
        while l_new[-1]<last_frequency+sampl:
            poly=sum([l_new[-1]**power*x for power,x in enumerate(coefficients[::-1])])
            # poly=0
            # for power,x in enumerate(coefficients[::-1]):
            #     poly+=l_new[-1]**power*x
                                     
            l_new.append(l_new[-1]+poly/sampl/min_sampl)
            # print(l_new[-1])
        
        return l_new
def synth_resolution_degradation(wave, synth,res_map,res_b,synth_res=30000000.0):
        """
        Take a synthetic spectrum with a very high  resolution and degrade its resolution to the resolution profile of the observed spectrum. The synthetic spectrum should not be undersampled, or the result of the convolution might be wrong.
        Parameters:
            synth np array or similar: an array representing the synthetic spectrum. Must have size m x 2. First column is the wavelength array, second column is the flux array. Resolution of the synthetic spectrum must be constant and higher than that of the observed spectrum.
            synth_res (float): resolving power of the synthetic spectrum
        Returns:
            Convolved syntehtic spectrum as a np array of size m x 2.
        """

        synth=np.vstack((wave,synth)).T
        
        l_original=synth[:,0]
        #check if the shape of the synthetic spectrum is correct
        if synth.shape[1]!=2: logging.error('Syntehtic spectrum must have shape m x 2.')

        #check if the resolving power is high enough
        sigma_synth=synth[:,0]/synth_res
        if max(sigma_synth)>=min(res_map)*0.95: logging.error('Resolving power of the synthetic spectrum must be higher.')

        #check if wavelength calibration of the synthetic spectrum is linear:
        if abs((synth[:,0][1]-synth[:,0][0])-(synth[:,0][-1]-synth[:,0][-2]))/abs(synth[:,0][1]-synth[:,0][0])>1e-6:
            logging.error('Synthetic spectrum must have linear (equidistant) sampling.')        

        #current sampling:
        sampl=galah_sampl=synth[:,0][1]-synth[:,0][0]

        #original sigma
        s_original=sigma_synth

        #required sigma (resample the resolution map into the wavelength range of the synthetic spectrum)
        s_out=res_map

        #the sigma of the kernel is:
        s=np.sqrt(s_out**2-s_original**2)

        #fit it with the polynomial, so we have a function instead of sampled values:
        map_fit=np.poly1d(np.polyfit(synth[:,0], s, deg=6))

        #create an array with new sampling. The first point is the same as in the spectrum:
        l_new=[synth[:,0][0]]

        #oversampling. If synthetic spectrum sampling is much finer than the size of the kernel, the code would work, but would return badly sampled spectrum. this is because from here on the needed sampling is measured in units of sigma.
        oversample=galah_sampl/sampl*5.0

        #minimal needed sampling
        min_sampl=max(s_original)/sampl/sampl*oversample

        #keep adding samples until end of the wavelength range is reached
        l_new=numba_syth_resolution(map_fit.coef,l_new,sampl,min_sampl,synth[:,0][-1])

        # while l_new[-1]<synth[:,0][-1]+sampl:
        #     l_new.append(l_new[-1]+map_fit(l_new[-1])/sampl/min_sampl)

        #interpolate the spectrum to the new sampling:
        new_f=np.interp(np.array(l_new),synth[:,0],synth[:,1])

        kernel_=galah_kern(max(s_original)/sampl*oversample, res_b)
    
        con_f=signal.fftconvolve(new_f,kernel_,mode='same')

        #inverse the warping:
        synth[:,1]=np.interp(l_original,np.array(l_new),con_f)
        return synth

class individual_spectrum:
    def __init__(self,name,interpolation,x,count):
            if x=='IR':
                starting_fraction=2048/4096
                length_fraction=2048/(4096*(1-starting_fraction))
            elif x=='Red':
                starting_fraction=0/4096
                length_fraction=4096/(4096*(1-starting_fraction))
            elif x=='Green':
                starting_fraction=0/4096
                length_fraction=4096/(4096*(1-starting_fraction))
            elif x=='Blue':
                starting_fraction=0/4096
                length_fraction=4096/(4096*(1-starting_fraction))
            else:
                starting_fraction=2085/4096
                length_fraction=2/(4096*(1-starting_fraction))
            tmp = np.load("NN_normalized_spectra_fe_alpha_shorter_"+x+".npz")
            w_array_0 = tmp["w_array_0"]
            w_array_1 = tmp["w_array_1"]
            w_array_2 = tmp["w_array_2"]
            b_array_0 = tmp["b_array_0"]
            b_array_1 = tmp["b_array_1"]
            b_array_2 = tmp["b_array_2"]
            x_min = tmp["x_min"]
            x_max = tmp["x_max"]
            tmp.close()
            NN_coeffs= (w_array_0, w_array_1, w_array_2, b_array_0, b_array_1, b_array_2, x_min, x_max)
            self.NN_coeff=NN_coeffs
            self.x_min=x_min
            self.x_max=x_max
            Path(name[0:6]+'/spectra/com/').mkdir(parents=True,exist_ok=True)
            if not exists(name[0:6]+'/spectra/com/'+name+str(count)+'.fits'):
                # source='/media/storage/HERMES_REDUCED/dr6.0/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'
                source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.0/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'
                destination=name[0:6]+'/spectra/com/'
                subprocess.run(["rsync",'-av',source,destination])


            hermes=fits.open(name[0:6]+'/spectra/com/'+name+str(count)+'.fits')
            
            if hermes[1].header.get('TEFF_R')=='None' or hermes[1].header.get('LOGG_R')=='None':
                print('Warning the reduction didnt produce an TEFF spectra might not reduce properly')
                print('enter to continue')
                input()

            x0= float(hermes[1].header.get('CRVAL1'))
            x1=float( hermes[1].header.get('CRVAL1')+len(hermes[1].data)* hermes[1].header.get('CDELT1'))
            fstart= x0+(x1-x0)*starting_fraction
            
            new_start=int(starting_fraction*len(hermes[1].data))
            new_end=new_start+int(len(hermes[1].data)*(1-starting_fraction)*length_fraction)
            
            
            length=np.linspace(x0+(x1-x0)*starting_fraction,x1,num=int(4096*interpolation*(1-starting_fraction)))
            
            fend=length[-1]
            #Increases the length of synthetic spectra so it over interpolates  
            self.wave=length
            self.spec=[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])]
            self.uncs=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[2].data[new_start:new_end])
            hermes[0].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[0].data[new_start:new_end])
            hermes[1].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])
            hermes[2].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[2].data[new_start:new_end])
            if len(hermes[7].data[new_start:new_end])!=new_end-new_start:
                hermes[7].data=np.ones(len(hermes[2].data))

            else:
                hermes[7].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[7].data[new_start:new_end])

                
            hermes[1].header['CDELT1']/=interpolation
            
            self.hermes=hermes
            self.teff=float(old_abundances['teff'])
            self.wran=[fstart,fend]
            self.vmic=float(old_abundances['vmic_reduction'])
            self.vsini=float(old_abundances['vbroad_reduction'])
            self.Fe=float(old_abundances['fe_h'])
            self.alpha=float(old_abundances['alpha_fe_r'])
            if not np.isnan(float(old_abundances['rv_'+x])):
                self.vrad=float(old_abundances['rv_'+x])
            else:
                self.vrad=0.0

                
            self.vmac=6.0
            self.logg=float(old_abundances['logg'])
            self.monh=float(old_abundances['fe_h'])
            random_number=str(np.random.randint(100000000))
            synthetic_temp=copy.deepcopy(self.hermes)
            if os.path.exists(name+x+random_number+'_synthetic.fits'):
                os.remove(name+x+random_number+'_synthetic.fits')
            synthetic_temp.writeto(name+x+random_number+'_synthetic.fits')

                        
            syn=gtools.read(name+x+random_number+'_synthetic')
            syn.equalize_resolution()
            os.remove(name+x+random_number+'_synthetic.fits')
            self.spec=syn.f
class spectrum_all:
    # bands=['IR']
    bands=['Blue','Green','Red','IR']

    def __init__(self,name,interpolation):
        name=str(name)
        bands=['Blue','Green','Red','IR']
        self.name=name
        self.interpolation=interpolation
        
        large_data=fits.getdata('gaia_galah_dr6.fits',1)
        large_data=Table(large_data)

        old_abundances=[x for x in large_data if x['sobject_id']==str(name)]
        old_abundances=vstack(old_abundances)
        for count,x in enumerate(bands,1):
            setattr(self,x,individual_spectrum(name,interpolation,x,count))

    def solar_value_maker(self,shift,colour,keys=['teff','logg','monh','Fe','alpha','vsini','vmac','vmic']):
        solar_values=[]
        for x in keys:
            if x in shift:
                solar_values.append(shift[x])
            else:
                solar_values.append(rgetattr(self,colour+'.'+x))
        return solar_values

    def synthesize(self,shift,multi=False,colours=bands,give_back=False):
        if not (multi or give_back):
            for x in colours:       
                solar_values=self.solar_value_maker(shift,x)
                spectrum = payne_sythesize(solar_values,rgetattr(self,x+'.x_min'),rgetattr(self,x+'.x_max'),rgetattr(self,x+'.NN_coeff'))
                if 'vrad' in shift:
                    dopler_shifted_spectra=dopler(rgetattr(self,x+'.wave'),spectrum,shift['vrad'])
                else:
                    dopler_shifted_spectra=dopler(rgetattr(self,x+'.wave'),spectrum,rgetattr(self,x+'.vrad'))
                dopler_shifted_spectra=synth_resolution_degradation(rgetattr(self,x+'.wave'),dopler_shifted_spectra,rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'])[:,1]

                rsetattr(self,x+'.synth',dopler_shifted_spectra)
                
        elif give_back:
            returning_spectra=np.array(np.ones(len(colours)),dtype=object)
            for number,x in enumerate(colours):       
                solar_values=self.solar_value_maker(shift,x)
                spectrum = payne_sythesize(solar_values,rgetattr(self,x+'.x_min'),rgetattr(self,x+'.x_max'),rgetattr(self,x+'.NN_coeff'))
                if 'vrad' in shift:
                    dopler_shifted_spectra=dopler(rgetattr(self,x+'.wave'),spectrum,shift['vrad'])
                else:
                    dopler_shifted_spectra=dopler(rgetattr(self,x+'.wave'),spectrum,rgetattr(self,x+'.vrad'))
                dopler_shifted_spectra=synth_resolution_degradation(rgetattr(self,x+'.wave'),dopler_shifted_spectra,rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'])[:,1]
                returning_spectra[number]=dopler_shifted_spectra
            return returning_spectra
        else:
            with Pool() as p:

                p.map(partial(getattr(self,'synthesize_multi'),shift=shift),['Blue','Green','Red','IR'])

            # pool = ThreadPool(4)
            # pool.map(partial(getattr(self,'synthesize_multi'),shift=shift),['Blue','Green','Red','IR'])
    def synthesize_multi(self,colours,shift):
        self.synthesize(shift=shift,multi=False,colours=[colours])
    def normalize(self,colours=bands,data=None):
        returning_spectra=np.array(np.ones(len(colours)),dtype=object)

        for value,x in enumerate(colours):
            if not np.array_equal(data,None):
                                
                original_line=rgetattr(self,x+'.spec')
                x_line=rgetattr(self,x+'.wave')
                synth_line=data[value]
                
                
                poly_order=4
                syth_coeff=np.nan
                while np.any(np.isnan(syth_coeff)) and poly_order>0:
                    syth_coeff=np.polyfit(x_line,synth_line,poly_order)
                    poly_order-=1
                if np.any(np.isnan(syth_coeff)) and poly_order==0:
                    np.savetxt('Not normal.csv', np.array([synth_line,x_line]), delimiter=',')
                    print("Cant be normalized")
                    break 
                ynew=np.poly1d(syth_coeff)
                poly_synth=ynew(x_line)
                
                original_coeff=np.polyfit(x_line,original_line,poly_order)
                ynew=np.poly1d(original_coeff)
                poly_original=ynew(x_line)
                # if x=='Blue':
                #     print("this is normalization",  original_line,"Line", x_line , "synthetic",  synth_line,"coeff",syth_coeff,"poly order", poly_order + 1, "spectra after",poly_synth/poly_original*original_line)                  
                returning_spectra[value]=poly_synth/poly_original*original_line
                
            elif rgetattr(self,x+'.synth').any():
                
                original_line=rgetattr(self,x+'.spec')

                x_line=rgetattr(self,x+'.wave')  
                synth_line=rgetattr(self,x+'.synth')
                
                
                poly_order=5
                syth_coeff=np.polyfit(x_line,synth_line,poly_order)
                ynew=np.poly1d(syth_coeff)
                poly_synth=ynew(x_line)
         
                original_coeff=np.polyfit(x_line,original_line,poly_order)
                ynew=np.poly1d(original_coeff)
                poly_original=ynew(x_line)
                
                synth_normal=poly_synth/poly_original*original_line
                
                # uncs_normal=poly_synth/poly_original*rgetattr(self,x+'.uncs')
                rsetattr(self,x+'.spec',synth_normal)
                # rsetattr(self,x+'.uncs',uncs_normal)
                
            else:
                print(x+' Hasnt been synthesized')
        if not np.array_equal(data,None):
            return returning_spectra

    def log_fit(self,colours=bands,synthetic_spectra=None,normal=None):
        probability=0
        return probability
        if np.array_equal(normal,None):
            for value,x in enumerate(colours):
                if not(np.array_equal(synthetic_spectra,None)):
                    probability+= -np.sum((synthetic_spectra[value]-rgetattr(self,x+'.spec'))**2/(rgetattr(self,x+'.uncs'))**2)
                else:
                    probability+= -np.sum((rgetattr(self,x+'.synth')-rgetattr(self,x+'.spec'))**2/(rgetattr(self,x+'.uncs'))**2)
        else:
            for value,x in enumerate(colours):
                if not(np.array_equal(synthetic_spectra,None)):
                    probability+= -np.sum((synthetic_spectra[value]-normal[value])**2/(rgetattr(self,x+'.uncs'))**2)
                else:
                    probability+= -np.sum((rgetattr(self,x+'.synth')-normal[value])**2/(rgetattr(self,x+'.uncs'))**2)
        return probability
    def plot(self,colours=bands):
        for x in colours:
            plt.figure()
            x_line=np.linspace(0,len(runGetattr(self,x+'.synth'))-1,num=len(runGetattr(self,x+'.synth')))
            if rgetattr(self,x+'.synth').any():
                    # labels='synthetic  chi squared= '+str(self.log_fit([x]))
                    plt.plot(x_line, runGetattr(self,x+'.synth'), label='Synthetic')
                    
            plt.plot(x_line, runGetattr(self,x+'.spec'), label='Observed')
            plt.title(x+" Band")
            plt.legend(loc='best')
            plt.tight_layout()
           

# x='Blue'
# tmp = np.load("NN_normalized_spectra_v_rad_"+x+".npz")
# w_array_0 = tmp["w_array_0"]
# w_array_1 = tmp["w_array_1"]
# w_array_2 = tmp["w_array_2"]
# b_array_0 = tmp["b_array_0"]
# b_array_1 = tmp["b_array_1"]
# b_array_2 = tmp["b_array_2"]
# wavelength=np.linspace(4718,4903,num=8192)

# x_min = tmp["x_min"]
# x_max = tmp["x_max"]
# tmp.close()
# NN_coeffs= (w_array_0, w_array_1, w_array_2, b_array_0, b_array_1, b_array_2, x_min, x_max)
# solar_values=[6700,3.2,1.1,30.0,1,3,2]
# a=payne_sythesize(solar_values,x_min,x_max,NN_coeffs)
# pos_try=np.column_stack((teff_random,vsini_random,logg_random,v_rad_random,monh_random,v_mac,v_mic))name=170517002801129
prior=True

name=160106004101254

c=299792  #km/s
# bands=['IR']

# interpolation=10
large_data=fits.getdata('gaia_galah_dr6.fits',1)
large_data=Table(large_data)

# spectras=spectrum_all(name,2)
# spectras.synthesize({})


# iron_H=float(old_abundances['fe_h'])
# names=old_abundances.colnames
# names_fe=[x for x in names if  '_fe' in x and ('e_' not in x) and ('flag' not in x) and ('nr_' not in x) and ('alpha' not in x) and ('Ti2' not in x)]
# names=[x[:-3] for x in names_fe ]

# votable = parse("cross_gaia.xml")
# photometric_values=votable.get_first_table().to_table(use_names_over_ids=True)
# photometric_values=[x for x in photometric_values if x['source_id']==old_abundances['source_id']][0]

# photometric_values.rename_column('t_eff','teff')
# photometric_values.rename_column('log_gSigma','logg_e')
# photometric_values.rename_column('t_effSigma','teff_e')
# photometric_values.rename_column('log_g','logg')






# radial_velocity
#Dopler
# fstart= float(hermes[1].header.get('CRVAL1')*(1+old_abundances['rv_galah']/c))
# fend=float( (hermes[1].header.get('CRVAL1')+len(hermes[1].data)* hermes[1].header.get('CDELT1'))*(1+old_abundances['rv_galah']/c))

# fstart= float(hermes[1].header.get('CRVAL1'))
# fend=float( (hermes[1].header.get('CRVAL1')+len(hermes[1].data)* hermes[1].header.get('CDELT1')))

# length=np.array([fstart+y*hermes[1].header.get('CDELT1')/interpolation for y in range(0,len(hermes[1].data)*interpolation)])


# #Increases the length of synthetic spectra so it over interpolates
# synthetic[0].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,len(hermes[1].data))]),hermes[0].data)
# synthetic[1].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,len(hermes[1].data))]),hermes[1].data)
# synthetic[2].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,len(hermes[1].data))]),hermes[2].data)
# synthetic[7].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,len(hermes[1].data))]),hermes[7].data)


# sme.wave= [length]
# sme.spec=[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,len(hermes[1].data))]),hermes[1].data)]
# sme.uncs=[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,len(hermes[1].data))]),hermes[2].data)]


# #Loads the spectra lines
# if str(name)[-1]=='1':
#     sme.linelist=ValdFile('Galah_Blue_2.lin')
# elif str(name)[-1]=='2':
#     sme.linelist=ValdFile('Galah_Green_2.lin')
# elif str(name)[-1]=='3':
#     sme.linelist=ValdFile('Galah_Red_2.lin')
# else:
#     sme.linelist=ValdFile('Galah_IR_2.lin')

# sme.wran =[fstart,fend]
# sme.abund = Abund(float(old_abundances['fe_h']),'asplund2009' ,type='H=12')


# sme.vmic=float(old_abundances['vmic'])
# sme.vsini=float(old_abundances['vbroad'])
# sme.vrad = float(old_abundances['rv_galah'])
# sme.vmac = 6.0
# sme.teff, sme.logg, sme.monh = old_abundances['teff'], old_abundances['logg'], old_abundances['fe_h']
# sme.h2broad = True

# sme.vrad_flag = "fix"
# sme.cscale_flag = "none"
# sme.cscale_type = "mask"


# sme.atmo.source = "marcs2012.sav"
# sme.atmo.method = "grid"
# sme.atmo.geom = "PP"

# sme.nlte.set_nlte("Fe", "marcs2012_Fe2016.grd")
# sme.nlte.set_nlte("Si", "marcs2012_Si2016.grd")
# sme.nlte.set_nlte("Ca", 'marcs2012p_t1.0_Ca.grd')
# sme.nlte.set_nlte("Na","marcs2012_Na2011.grd")
# sme.nlte.set_nlte("O",'marcs2012_O2015.grd')
# sme.nlte.set_nlte("Ti","marcs2012s_t2.0_Ti.grd")
# sme.nlte.set_nlte("Ba","marcs2012p_t1.0_Ba.grd")
# sme.nlte.set_nlte("Mg","marcs2012_Mg2016.grd")

# fitparameters = ['vrad']
# sme = solve(sme, fitparameters)


#one core
# synthesized_spectra = synthesize({})


# Multiprocessing


# steps=np.linspace(115,130,11)
# with Pool() as p:

#       synthesized_spectra=p.map(partial(synthesize_multi,values='vrad',changed_values={}),steps)





# #normalization

# if len(np.shape(synthesized_spectra))==1:
#     synth_line=synthesized_spectra
# else:
#     synth_line=synthesized_spectra[int(np.ceil(len(synthesized_spectra)/2))]
# original_line=sme.spec[0]
# x=np.linspace(fstart,fend,len(original_line))



# poly_order=1
# syth_coeff=np.polyfit(x,synth_line,poly_order)
# ynew=np.poly1d(syth_coeff)
# poly_synth=ynew(x)

# original_coeff=np.polyfit(x,original_line,poly_order)
# ynew=np.poly1d(original_coeff)
# poly_original=ynew(x)

# synth_normal=poly_synth/poly_original*original_line

# uncs_normal=poly_synth/poly_original*sme.uncs[0]

# original_coeff_new=np.polyfit(x,synth_normal,poly_order)
# ynew=np.poly1d(original_coeff_new)
# poly_original_new=ynew(x)

# # plt.figure()
# # plt.plot(x,original_line,label='original')
# # plt.plot(x,synth_normal,label='synth Normalized')
# # plt.plot(x,average_normal,label='average normalized')
# # plt.legend()
# # plt.xlabel('frequency/A')
# # plt.ylabel('Intensity')
# # plt.savefig('Normal_Blue')




# sme.spec=[synth_normal]

# sme.uncs=[uncs_normal]

# #EMCEE,
prior=False
final_results=[]
expected_results=[]
# for value,x in enumerate(large_data[8:560]):
# final_results_temp=[]
x=large_data[0]
name=str(x['sobject_id'])
global old_abundances
prior=True



old_abundances=[x for x in large_data if x['sobject_id']==name]
old_abundances=vstack(old_abundances)
filename = str(name)+'machine_learning_prior.h5'

spectras=spectrum_all(name,10)

print(filename)

# spec_janez=spectras.synthesize({},give_back=True)
spectras.synthesize({})
spectras.normalize()

nwalkers=100
pos=[]
while np.shape(pos)[0]<nwalkers:
    if old_abundances['teff'] and not np.isnan( old_abundances['teff']):
        teff_random=abs(np.random.normal(old_abundances['teff'],1000,1))
    else:
        teff_random=-100
    if old_abundances['logg'] and not np.isnan( old_abundances['logg']):
        logg_random=abs(np.random.normal(old_abundances['logg'],0.4,1))
    else:
        logg_random=-100
    if np.isnan(old_abundances['rv_av']):
        v_rad_random=np.random.normal(0,4*4,1)
    else:
        v_rad_random=np.random.normal(float(old_abundances['rv_av']),4*4,1)
    if old_abundances['fe_h']  and not np.isnan( old_abundances['fe_h']):
        monh_random=np.random.normal(old_abundances['fe_h'],0.2,1)
    else:
        monh_random=-100
    if old_abundances['vmic_reduction'] and not np.isnan( old_abundances['vmic_reduction']):
        v_mac=abs(np.random.normal(old_abundances['vmic_reduction'],8,1))+np.random.randint(-10,10)
        v_mic=abs(np.random.normal(old_abundances['vmic_reduction'],8,1))
    else:
        v_mic=np.random.normal(10,1,1)
        v_mac=np.random.normal(10,1,1)
    if old_abundances['vbroad_reduction'] and not np.isnan( old_abundances['vbroad_reduction']) and old_abundances['vbroad_reduction']<100 :
        vsini_random=abs(np.random.normal(old_abundances['vbroad_reduction'],1*4,1))
    else:
        vsini_random=np.random.normal(10,1,1)
    if old_abundances['alpha_fe_r'] and not np.isnan( old_abundances['alpha_fe_r']):
        alpha_random=np.random.normal(old_abundances['alpha_fe_r'],0.1,1)
    else:
        alpha_random=np.random.normal(0.148,0.317,1)
    fe_random=monh_random+np.random.normal(0,0.1,1)


    pos_try=np.column_stack((teff_random,logg_random,monh_random,fe_random,alpha_random,v_rad_random,vsini_random,v_mac,v_mic))
    if starting_test(pos_try[0]):
        if len(pos)==0:
            pos=pos_try
        else:
            pos=np.vstack((pos_try,pos))
print('starting')
time0=time.time()
for x in pos:
    log_posterior(x,prior=True)
time1=time.time()
print('time',(time1-time0)/len(pos))

# ndim=np.shape(pos)[1]
# backend = emcee.backends.HDFBackend(filename)
# log_posterior(pos[0])
# backend.reset(nwalkers, ndim)
# step_iteration=50
# with Pool(processes=nwalkers) as pool:
#     sampler_burn_in = emcee.EnsembleSampler(nwalkers, ndim, log_posterior,pool=pool,a=5,backend=backend,args=[prior])
#     sampler_burn_in.run_mcmc(pos, 10, progress=True)
#     pos_burn_in=sampler_burn_in.get_chain(flat=False)
#     pos=pos_burn_in[-1]
#     sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior,pool=pool,a=3,backend=backend,args=[prior])
#     # sampler.run_mcmc(pos, 2, progress=True)
    
    
#     autocorr=[]
#     index=0
#     oldTau=np.inf
#     old_mean=[]
#     old_standard_deviation=[]
    
#     for sample in sampler.sample(pos,iterations=100, progress=True):
#         if sampler.iteration % step_iteration:
#                 continue
#         tau=sampler.get_autocorr_time(tol=0)
#         print("Autocorelation Time is :",tau)
#         if not len(old_mean) and sampler.iteration==step_iteration*2:
#             old_mean=np.mean(sampler.get_chain()[-step_iteration:],axis=0)
#             old_standard_deviation=np.std(sampler.get_chain()[-step_iteration:],axis=0)
#         if sampler.iteration>np.nanmax(tau)*5 and sampler.iteration>step_iteration*3:

#             new_mean=np.mean(sampler.get_chain()[-step_iteration:],axis=0)
#             new_standard_deviation=np.std(sampler.get_chain()[-step_iteration:],axis=0)
#             temp_new_mean,temp_old_mean=getting_rid_zeros(new_mean,old_mean)
#             temp_new_standard_deviation,temp_old_standard_deviation=getting_rid_zeros(new_standard_deviation,old_standard_deviation)

#             difference=np.sqrt(np.sum((temp_new_mean-temp_old_mean)**2/((abs(temp_new_mean)+abs(temp_old_mean))/2)**2)
#                                 /(np.shape(temp_new_mean)[0]*np.shape(temp_new_mean)[1]))
#             difference+=np.sqrt(np.sum((temp_new_standard_deviation-temp_old_standard_deviation)**2/((temp_new_standard_deviation+temp_old_standard_deviation)/2)**2)
#                                 /(np.shape(temp_new_standard_deviation)[0]*np.shape(temp_new_standard_deviation)[1]))
#             difference/=2
#             temp_data=sampler.get_chain(flat=True)
#             print('gaussian and moving difference')
#             print(how_good_gaussian_fit(temp_data),print(difference))
#             if how_good_gaussian_fit(temp_data)<0.5 and difference/2<0.1:
#                 print('converged_stoped')
#                 break
                
#             else:
#                 old_mean=np.copy(new_mean)
#                 old_standard_deviation=np.copy(new_standard_deviation)
            
#         if not np.any(np.isnan(tau)):
#             autocorr.append(tau)
#             converged = np.all(tau * 100 < sampler.iteration)
#             converged &= np.all(np.abs(oldTau - tau) / tau < 0.01)
#             if converged:
#                 break
#             oldTau = tau
