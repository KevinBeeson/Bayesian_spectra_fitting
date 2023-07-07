#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 16:26:48 2022
synthesizes spectra then plots them
@author: kevin
"""

from numba import jit
import collections

from functools import  partial
from astropy.io.votable import parse
import emcee
from The_Payne import spectral_model
import scipy
from scipy import signal
from os.path import exists
import subprocess
from pathlib import Path
from datetime import date
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
from scipy.stats import chi2
import warnings
from os.path import exists
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
        """
        uses chi squared to see if the n dimentional data fits n gaussian each 
    
        Parameters
        ----------
        data : you can insert a data matrix where the np.shape[1] dimention shows the dimentions of the data 
    
        Returns
        -------
        FLOAT: probability that the data is gaussian 
    
        """
        prob=0
        for x in range(np.shape(data)[1]):
            if x==9:
                continue
            temp_data=data[:,x]
            
            #fits a gaussian
            mean=np.mean(temp_data)
            standard_deviation=np.std(temp_data)
            counts, bins = np.histogram(temp_data,density=True)
            bins=(bins[1:]+bins[:-1])/2
            gaussian=1/(standard_deviation*np.sqrt(2*np.pi))*np.exp(-0.5*((mean-bins)/standard_deviation)**2)
            gaussian, counts=getting_rid_zeros(gaussian,counts)
            diff=np.sum((gaussian-counts)**2/((counts+gaussian)/2))
    
            free_parameters=len(gaussian)
            
            #uses p value for chisquared to see if the gaussian fits
            prob_temp=chi2.cdf(diff,free_parameters)
            if prob_temp>prob:
                prob=np.copy(prob_temp)
        return prob
def getting_rid_zeros(data1,data2):
        """
        checks boths data and see if the sum of the data points are zeros and returns the non zero data points 
    
        Parameters
        ----------
        data1 , data2: arrays of equal length 1xn 
        
        Returns
        -------
        TYPE
        returns the two arrays with non zero data points
    
        """
        data1_temp=[]
        data2_temp=[]
        for synthetic,data_bins in zip(data1,data2):
            #rather than ussing 0 1e-100 is better because we are sqauring the data
            if (abs(synthetic)+abs(data_bins)>1e-100).all():
                data1_temp.append(synthetic)
                data2_temp.append(data_bins)
        return np.array(data1_temp),np.array(data2_temp)
def shift_maker(solar_values):
        """
        Create a dictionary of solar parameters from an array of values
        
    
        Parameters
        ----------
        solar_values : array of 9 length in the order of  teff,logg,monh,fe,alpha,vrad,vsini,vmac,vmic
    
        Returns
        -------
        shift : dictionary of solar parameters
    
        """
        labels=['teff','logg','monh','Fe','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']
        shift={}
        for y,x in zip(labels,solar_values):
            shift[y]=x
        return shift

def log_prior(shift):
        """
        Gives the log prior of the current star given the photometric teff and logg 
    
        Parameters
        ----------
        shift : dictionary of the solar parameters 
    
        Returns
        -------
        float of logg of the prior.
    
        """
        error=0
        for x in shift:
            if x in ('logg','teff'):
                error-=(shift[x]-float(old_abundances[x]))**2/(2*(float(old_abundances['e_'+x]))**2)
        if error<-1e100:
            return -np.inf
        else:
            return error
def log_posterior(solar_values,prior=True):
        """
        Gives how good the fit of the inserted solar parameters are compared to the observed spectra 
    
        Parameters
        ----------
        solar_values : 1x9 array of solar parameters with teff,logg,monh,fe,alpha,vrad,vsini,vmac,vmic  order
            DESCRIPTION.
        prior : BOOL, optional
            if True it takes into account the photometric prior. The default is True.
    
        Returns
        -------
        TYPE
        float of logg of how good the fit is .
    
        """
        #Sanity check to see if the solar values are in a certain section
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
        
        #Sanity checks Dont think needs to be done anymore as we're using payne
        # if not np.all([np.all(x) for x in synthetic_spectras]):
        #     print("Didnt synthesize properly")
        #     return -np.inf
        # if np.max([np.max(abs(x)) for x in synthetic_spectras])>100 or np.max([np.max(abs(x)) for x in synthetic_spectras])==np.nan:
        #     print("Didnt synthesize properly one of the spectra blew up")
        #     return -np.inf
        # if np.any(np.array_equal(synthetic_spectras,np.nan)):
        #     print("Didnt synthesize properly Nan values yay!")
        #     return -np.inf
    
    
        # normalized_spectra=spectras.normalize(data=synthetic_spectras)
        
        #gets the probability of how good the fits are and also adds the prior if wanted
        probability=spectras.log_fit(synthetic_spectra=synthetic_spectras)
    
        if prior:
            probability+=log_prior(shift)
        return probability
def starting_test(solar_values):
        """
        Sanity check to see if the solar values inserted passes a sanity check
    
        Parameters
        ----------
        solar_values : 1x9 array with the order of teff,logg,monh,fe,alpha,vrad,vsini,vmac,vmic
    
        Returns
        -------
        bool
            If it doesnt pass the sanity check it will return False and if passes will return True
    
        """
        labels=['teff','logg','monh','Fe','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']
        shift={}
        for y,x in zip(labels,solar_values):
            shift[y]=x
        labels_with_limits=['teff','logg','monh','Fe','alpha','vsini','vmac','vmic']
        for value,x in enumerate(labels_with_limits):
            if shift[x]<x_min[value]-abs(x_min[value])*0.2 or shift[x]>x_max[value]*1.2:
                print('outside of the models limits ',x,' value ', shift[x])
                return False
        vrad=['vrad_Blue','vrad_Green','vrad_Red','vrad_IR']
        vrad_r=['rv_Blue_r','rv_Green_r','rv_Red_r','rv_IR_r']
        radial_velocities=[shift[x] for x in vrad]
        radial_velocities_r=[old_abundances[x] for x in vrad_r]
        for rad,rad_r in zip(vrad,vrad_r):
            if not np.isnan(old_abundances[rad_r]):
                mean=float(old_abundances[rad_r])
            elif not np.isnan(old_abundances['rv_r']):
                mean=float(old_abundances['rv_r'])
            else:
                break
            if not np.isnan(old_abundances['e_'+rad_r]):
                sig=float(old_abundances['e_'+rad_r])*3
            elif not np.isnan(old_abundances['e_rv_r']):
                sig=float(old_abundances['e_rv_r'])*3
            else:
                sig=5
            if abs(shift[rad]-mean)>sig+1:
                print(rad,' is wrong by ',str(float(abs(shift[rad]-old_abundances[rad_r]))))
                return False
        if np.count_nonzero(~np.isnan(radial_velocities))>1.0:
            if max(radial_velocities)-min(radial_velocities)>np.nanmax(radial_velocities_r)-np.nanmin(radial_velocities_r)+2:
                print('radial velocities too different')
                return False
        else:
            if max(radial_velocities)-min(radial_velocities)>5+2:
                print('radial velocities too different')
                return False
        return True
#Three function so one can set and get attributed in classes in mutiple levels i.e. spectras.Blue.wave
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

def equalize_resolution(wavelength,res_map,observed_spectra,observed_error,hermes, r='lowest', precision=0.005, kernel='galah',synthetic=False):
    """
    Convolves a spectrum with a kernel with a variable width. Works by warping the data, performing the convolution and unwarping the data, so it is vectorized (mostly) and fast
    Parameters:
        r (str or float): resolution of resulting spectra. If r=='lowest', the resulting spectra will have the lowest resolution given in the resolution map
        precision (float): precision of intermitten sampling. Lower number means better precision.
        kernel (str): 'gauss' or 'galah' (default). Will either use a gaussian kernel or a kernel derived for GALAH as the konvolution kernel.
    """


    #current sampling:
    sampl=wavelength[1]-wavelength[0]

    #target sigma coresponds to the R=22000. We want the constant sigma, not constant R, so take the sigma that corresponds to the average R=22000
    #s_target=np.ones(len(self.map_l))*np.average(self.map_l)/18400.
    if r=='lowest':
        s_target=np.ones(len(wavelength))*max(res_map)*(1+precision)
    elif isinstance(r, (int,  float)):
        s_target=res_map*(1+precision)
        s_target[r>s_target]=r
    else:
        logging.error('Parameter r must be either \'lowest\' or a number.')
        
        
    if type(synthetic)==int or type(synthetic)==float:
        theoretical_resolution=synthetic
    else:
        theoretical_resolution=0.0001
    #original sigma
    if synthetic==False:
        s_original=res_map
    else:
        s_original=np.ones(len(res_map))*theoretical_resolution
    #the sigma of the kernel is:
    s=s_original


    #fit it with the polynomial, so we have a function instead of sampled values:
    map_fit=np.poly1d(np.polyfit(wavelength, s, deg=6))

    #create an array with new sampling. The first point is the same as in the spectrum:
    l_new=[wavelength[0]]

    #minimal needed sampling
    min_sampl=max(s_original)/sampl/sampl

    #keep adding samples until end of the wavelength range is reached
    l_new=numba_syth_resolution(map_fit.coef,l_new,sampl,min_sampl,wavelength[-1])

    # while l_new[-1]<wavelength[-1]+sampl:
    #     l_new.append(l_new[-1]+map_fit(l_new[-1])/sampl/min_sampl)

    #interpolate the spectrum to the new sampling:
    new_f=np.interp(np.array(l_new),wavelength,observed_spectra)
    new_fe=np.interp(np.array(l_new),wavelength,observed_error)

    #plot(l_new, l_new-np.roll(l_new, 1), 'r,')
    #show()

    #convolve
    res_b=float(hermes[7].header['b'])
    if  kernel=='galah':
        kernel_=galah_kern(max(s_original)/sampl, res_b)
    else:
        logging.error('Kernel %s is not available. Use gauss or galah.' % kernel)
        return False
    con_f=signal.fftconvolve(new_f,kernel_,mode='same')
    con_fe=signal.fftconvolve(new_fe**2,kernel_,mode='same')

    #inverse the warping:
    new_observed_spectra=np.interp(wavelength,np.array(l_new),con_f)
    new_error=np.sqrt(np.interp(wavelength,np.array(l_new),con_fe))
    
    return new_observed_spectra,new_error
def payne_sythesize(solar_values,x_min,x_max,NN_coeffs):
        """
        Synthesizes the spectra using Payne
    
        Parameters
        ----------
        solar_values : a 1x8 array( the solar value arrays without vrad )  using teff,logg,monh,fe,alpha,vrad,vsini,vmac,vmic order
        x_min : min value for payne
        x_max : max values for payne
        NN_coeffs :Matrix coefficients gotten from Payne 
    
        Returns
        -------
        real_spec : 1xn array which is the payne sythesized spectra
    
        """
    
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
def dopler(original_wavelength,spectra,v_rad):
        """
        Shifts and crops your spectra to be the same length as the observed spectra
    
        Parameters
        ----------
        original_wavelength : 1xn array of the observed wavelength range (i.e. observed wavelengths)
        spectra : 1xm synthesized spectra that you want to shift
        v_rad : float
            the radial velocity you want to shift your spectra by.
        Returns
        -------
        spectra_new : TYPE
            returns a radially shifted spectra with the length of the observed spectra .
    
        """
        c=299792.458
        delta_v=v_rad/c
        #my spectra are syntheszied to be larger than the galah spectra so now it crops see by how much 
        if len(original_wavelength)!=len(spectra):
            difference=abs(len(original_wavelength)-len(spectra))
            original_wavelength_strenched=np.interp(np.linspace(-difference/2,difference/2+len(original_wavelength),num=len(spectra)),range(len(original_wavelength)),original_wavelength)
            observed_wavelength=original_wavelength_strenched*(1+delta_v)
        else:
            observed_wavelength=original_wavelength*(1+delta_v)
            
        #crops and interpolates 
        spectra_new=np.interp(original_wavelength,observed_wavelength,spectra)
        
        return spectra_new
@jit(nopython=True,parallel=True,cache=True)
def numba_syth_resolution(coefficients,l_new,sampl,min_sampl,last_frequency):
        """
        Slow part of resolution degradation is done in JIT to make it faster
    
        Parameters
        ----------
        coefficients : TYPE
            DESCRIPTION.
        l_new : TYPE
            DESCRIPTION.
        sampl : TYPE
            DESCRIPTION.
        min_sampl : TYPE
            DESCRIPTION.
        last_frequency : TYPE
            DESCRIPTION.
    
        Returns
        -------
        l_new : TYPE
            DESCRIPTION.
    
        """
        while l_new[-1]<last_frequency+sampl:
            poly=sum([l_new[-1]**power*x for power,x in enumerate(coefficients[::-1])])
                                     
            l_new.append(l_new[-1]+poly/sampl/min_sampl)
        
        return l_new
def synth_resolution_degradation(wave, synth,res_map,res_b,synth_res=300000.0):
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
            tmp = np.load("NN_normalized_spectra_fe_alpha_shortest_"+x+".npz")
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
                #source='/media/storage/HERMES_REDUCED/dr6.0/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'
                source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.0/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'

                destination=name[0:6]+'/spectra/com/'
                subprocess.run(["rsync",'-av',source,destination])


            hermes=fits.open(name[0:6]+'/spectra/com/'+name+str(count)+'.fits')
            
            if hermes[1].header.get('TEFF_R')=='None' or hermes[1].header.get('LOGG_R')=='None':
                print('Warning the reduction didnt produce an TEFF spectra might not reduce properly')
                print('enter to continue')
                self.bad_reduction=True
            else:
                self.bad_reduction=False
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
            self.spec_original=[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])][0]

            self.uncs=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[2].data[new_start:new_end])
            hermes[0].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[0].data[new_start:new_end])
            hermes[1].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])
            hermes[2].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[2].data[new_start:new_end])
            # not_resolved=[np.int64(name)]
            if len(hermes[7].data[new_start:new_end])==new_end-new_start:
            #     hermes[7].data=False
            #     print("Resolution map is wrong ",name)
            #     temp_data=[x for x in all_reduced_data if x['sobject_id']==np.int64(name)]
            #     temp_data=np.vstack(temp_data)
            #     fibre_number=temp_data['fibre']
            #     plate=temp_data['plate']
            #     time=temp_data['utdate'][0][0]
            #     fibre_sisters=[x for x in all_reduced_data if x['fibre']==fibre_number and x['plate']==plate and not x['sobject_id'] in not_resolved]
            #     fibre_sisters=vstack(fibre_sisters)
            #     difference=difference_UT_vector(time, fibre_sisters['utdate'])
            #     fibre_sisters['difference']=difference
            #     fibre_sisters.sort('difference')
            #     name_target=str(fibre_sisters[0]['sobject_id'])
            #     print('copying from '+name_target)
            #     not_resolved.append(name_target)
            #     Path(name_target[0:6]+'/spectra/com/').mkdir(parents=True,exist_ok=True)
            #     if not exists(name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'):
            #         # source='/media/storage/HERMES_REDUCED/dr6.0/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'
            #         source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.0/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'

            #         destination=name_target[0:6]+'/spectra/com/'
            #         subprocess.run(["rsync",'-av',source,destination])


            #     hermes_temp=fits.open(name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits')
            #     hermes[7].data=hermes_temp[7].data
            # else:
                hermes[7].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[7].data[new_start:new_end])

                
            hermes[1].header['CDELT1']/=interpolation
            
            self.hermes=hermes
            self.teff=float(old_abundances['teff'])
            self.wran=[fstart,fend]
            self.vmic=float(old_abundances['vmic_reduction'])
            self.vsini=float(old_abundances['vbroad_reduction'])
            self.Fe=float(old_abundances['fe_h'])
            self.alpha=float(old_abundances['alpha_fe_r'])
            if not np.isnan(float(old_abundances['rv_'+x+'_r'])):
                self.vrad=float(old_abundances['rv_'+x+'_r'])
            elif not np.isnan(float(old_abundances['rv_r'])):
                self.vrad=float(old_abundances['rv_r'])
            else:
                self.vrad=0.0

                
            self.vmac=6.0
            self.logg=float(old_abundances['logg'])
            self.monh=float(old_abundances['fe_h'])
@jit(nopython=True,parallel=True,cache=True)
def array_limiter(mother_array,baby_array_1,baby_array_2,limit=1.05):
    """
    returns the both arrays where the value of the mother_array is bellow 1.05

    Parameters
    ----------
    mother_array : 1xn array
    baby_array : 1xn array
        DESCRIPTION.
    limit : float, optional
        whats the cut of limit. The default is 1.05.

    Returns
    -------
    2 arrays 
        both arrays where the mother array's value is bellow 1.05.

    """
    all_arrays=np.column_stack((mother_array,baby_array_1,baby_array_2))
    all_temp=[]
    [all_temp.append(x) for x in all_arrays if x[0]<1.05]
    all_temp=np.array(all_temp)
    return all_temp[:,0],all_temp[:,1],all_temp[:,2]
def difference_UT_vector(date_str1,date_vector):
    # date_1=datetime(int(date_str1[:4]),int(date_str1[5:7]),int(date_str1[8:10]),int(date_str1[11:13]),int(date_str1[14:16]),int(date_str1[17:19]))
    date_1=date(int(date_str1[:4]),int(date_str1[5:7]),int(date_str1[8:10]))

    year=[int(x[:4]) for x in date_vector]
    month=[int(x[5:7]) for x in date_vector]
    day=[int(x[8:10]) for x in date_vector]
    # hour=[int(x[11:13]) for x in date_vector]
    # minute=[int(x[14:16]) for x in date_vector]
    # second=[int(x[17:19]) for x in date_vector]
    date_2=[]
    for x in range(len(year)):
        date_2.append(date(year[x],month[x],day[x]))
        # date_2.append(datetime(year[x],month[x],day[x],hour[x],minute[x],second[x]))

    difference=[]
    for x in date_2:
        diffrence_temp=x-date_1
        difference.append(abs(diffrence_temp.days))
    return difference
class spectrum_all:
    # bands=['IR']
    bands=['Blue','Green','Red','IR']

    def __init__(self,name,interpolation):
        name=str(name)
        bands=['Blue','Green','Red','IR']
        self.bands=bands
        self.name=name
        self.interpolation=interpolation
        self.sister_stars=None
        large_data=fits.getdata('photometric_reduction_dr6_cross.fits',1)
        large_data=Table(large_data)

        old_abundances=[x for x in large_data if x['sobject_id']==str(name)]
        old_abundances=vstack(old_abundances)
        for count,x in enumerate(bands,1):
            setattr(self,x,individual_spectrum(name,interpolation,x,count))
        self.correct_resolution_map()
        self.equilize_spectras()
        self.limit_array()
    def equilize_spectras(self,colours=bands):
        for x in colours:
            wavelength=rgetattr(self, x+'.wave')
            hermes=rgetattr(self, x+'.hermes')
            resmap=hermes[7].data
            observed_spectra=rgetattr(self, x+'.spec')[0]
            observed_error=rgetattr(self, x+'.uncs')
            new_observed,new_error=equalize_resolution(wavelength,resmap,observed_spectra,observed_error,hermes)
            rsetattr(self, x+'.spec',new_observed)
            rsetattr(self, x+'.uncs',new_error)
    def dip_priority(self,colours=bands,sigfig=0.5):
        for x in colours:
            dip_array=[np.exp((y-1.0)**2/2*sigfig**2) for y in rgetattr(self,x+'.spec')]
            rsetattr(self, x+'.dip_array',dip_array)
    def limit_array(self,colours=bands,limit=1.05):
        for x in colours:
            limit_array=[]
            for y in rgetattr(self,x+'.spec'):
                if y>limit:
                    limit_array.append(0)
                else:
                    limit_array.append(1)
            rsetattr(self,x+'.limit',limit_array)
    def solar_value_maker(self,shift,colour,keys=['teff','logg','monh','Fe','alpha','vsini','vmac','vmic']):
        solar_values=[]
        for x in keys:
            if x in shift:
                solar_values.append(shift[x])
            else:
                solar_values.append(rgetattr(self,colour+'.'+x))
        return solar_values
    def hermes_checker(self,limit=20,colours=bands):
        for x in colours:
            spectra=rgetattr(self, x+'.spec')
            if max(spectra)>limit:
                print (max(spectra),' ',x)
                return True
        return False

    def synthesize(self,shift={},colours=bands,multi=False,give_back=False,full=False):
        if full and not give_back:
            print('you probably want the spectrum back run again with give_back=True')
            return False
        if not give_back:
            if not multi:
                for x in colours:       
                    solar_values=self.solar_value_maker(shift,x)
                    spectrum = payne_sythesize(solar_values,rgetattr(self,x+'.x_min'),rgetattr(self,x+'.x_max'),rgetattr(self,x+'.NN_coeff'))
                    if 'vrad_'+x in shift:
                        dopler_shifted_spectra=dopler(rgetattr(self,x+'.wave'),spectrum,shift['vrad_'+x])
                    else:
                        dopler_shifted_spectra=dopler(rgetattr(self,x+'.wave'),spectrum,rgetattr(self,x+'.vrad'))
                    dopler_shifted_spectra=synth_resolution_degradation(rgetattr(self,x+'.wave'),dopler_shifted_spectra,rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'])[:,1]
    
                    rsetattr(self,x+'.synth',dopler_shifted_spectra)
            else:
                with ThreadPool(4) as pool:
                    inputs=[(shift,x,False,False) for x in colours]
                    pool.map(getattr(self,'synthesize'),inputs)
                
        else:
            if not multi:
                returning_spectra=np.array(np.ones(len(colours)),dtype=object)
                for number,x in enumerate(colours):       
                    solar_values=self.solar_value_maker(shift,x)

                    spectrum = payne_sythesize(solar_values,rgetattr(self,x+'.x_min'),rgetattr(self,x+'.x_max'),rgetattr(self,x+'.NN_coeff'))
                    if full:
                        returning_spectra[number]=spectrum
                        continue                   
                    if 'vrad_'+x in shift:
                        dopler_shifted_spectra=dopler(rgetattr(self,x+'.wave'),spectrum,shift['vrad_'+x])
                    else:
                        dopler_shifted_spectra=dopler(rgetattr(self,x+'.wave'),spectrum,rgetattr(self,x+'.vrad'))
                    dopler_shifted_spectra=synth_resolution_degradation(rgetattr(self,x+'.wave'),dopler_shifted_spectra,rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'])[:,1]
                    returning_spectra[number]=dopler_shifted_spectra
                return returning_spectra
            else :
                with ThreadPool() as pool:
                    inputs=[(shift,x,False,True) for x in colours]
                    return pool.map(partial(self.synthesize_multi,shift=shift),colours)
                    # pool.map(partial(getattr(self,'synthesize_multi'),give_back=True,shift=shift),colours=colours)
    def synthesize_multi(self,colours,shift):
        return self.synthesize(shift,multi=False,give_back=True,colours=[colours])
    def correct_resolution_map(self,colours=bands):
        name=self.name
        not_resolved=[np.int64(name)]
        for count,x in enumerate(colours,1):
            hermes=rgetattr(self,x+'.hermes')
            if len(hermes[7].data)!=len(hermes[1].data):
                print('Resolution map of ',name +str(count),' is wrong')
                if self.sister_stars==None:
                    while len(hermes[7].data)!=len(hermes[1].data)/self.interpolation:
                        temp_data=[x for x in all_reduced_data if x['sobject_id']==np.int64(name)]
                        temp_data=np.vstack(temp_data)
                        fibre_number=temp_data['fibre']
                        plate=temp_data['plate']
                        time=temp_data['utdate'][0][0]
                        fibre_sisters=[x for x in all_reduced_data if x['fibre']==fibre_number and x['plate']==plate and not x['sobject_id'] in not_resolved]
                        fibre_sisters=vstack(fibre_sisters)
                        difference=difference_UT_vector(time, fibre_sisters['utdate'])
                        fibre_sisters['difference']=difference
                        fibre_sisters.sort('difference')
                        name_target=str(fibre_sisters[0]['sobject_id'])
                        print('copying from '+name_target)
                        not_resolved.append(name_target)
                        Path(name_target[0:6]+'/spectra/com/').mkdir(parents=True,exist_ok=True)
                        if not exists(name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'):
                            # source='/media/storage/HERMES_REDUCED/dr6.0/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'
                            source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.0/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'
        
                            destination=name_target[0:6]+'/spectra/com/'
                            subprocess.run(["rsync",'-av',source,destination])
        
        
                        hermes_temp=fits.open(name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits')
                        hermes[7].data=hermes_temp[7].data
                    self.sister_stars=name_target
                else:
                    name_target=self.sister_stars
                    print('copying from '+name_target)

                    hermes_temp=fits.open(name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits')
                    hermes[7].data=hermes_temp[7].data
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
                x0= float(hermes[3].header.get('CRVAL1'))
                x1=float( hermes[3].header.get('CRVAL1')+len(hermes[3].data)* hermes[3].header.get('CDELT1'))
                fstart= x0+(x1-x0)*starting_fraction
                
                new_start=int(starting_fraction*len(hermes[3].data))
                new_end=new_start+int(len(hermes[3].data)*(1-starting_fraction)*length_fraction)
                
                
                length=np.linspace(x0+(x1-x0)*starting_fraction,x1,num=int(4096*self.interpolation*(1-starting_fraction)))
                hermes[7].data=np.interp(length,np.array([fstart+y*hermes[3].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[7].data[new_start:new_end])

    def normalize(self,colours=bands,data=None):
        returning_spectra=np.array(np.ones(len(colours)),dtype=object)

        for value,x in enumerate(colours):
            if not np.array_equal(data,None):
                                
                original_line=rgetattr(self,x+'.spec_original')
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
                
                original_line=rgetattr(self,x+'.spec_original')

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

    def log_fit(self,colours=bands,dip_array=False,synthetic_spectra=None,normal=None):
        probability=0
        if not dip_array:
            if np.array_equal(normal,None):
                for value,x in enumerate(colours):
                    if not(np.array_equal(synthetic_spectra,None)):
                        probability+= -np.sum((synthetic_spectra[value]-rgetattr(self,x+'.spec'))**2/(rgetattr(self,x+'.uncs'))**2*rgetattr(self,x+'.limit'))
                    else:
                        probability+= -np.sum((rgetattr(self,x+'.synth')-rgetattr(self,x+'.spec'))**2/(rgetattr(self,x+'.uncs'))**2*rgetattr(self,x+'.limit'))
            else:
                for value,x in enumerate(colours):
                    if not(np.array_equal(synthetic_spectra,None)):
                        probability+= -np.sum((synthetic_spectra[value]-normal[value])**2/(rgetattr(self,x+'.uncs'))**2*rgetattr(self,x+'.limit'))
                    else:
                        probability+= -np.sum((rgetattr(self,x+'.synth')-normal[value])**2/(rgetattr(self,x+'.uncs'))**2*rgetattr(self,x+'.limit'))
        else:
            if np.array_equal(normal,None):
                for value,x in enumerate(colours):
                    if not(np.array_equal(synthetic_spectra,None)):
                        probability+= -np.sum((synthetic_spectra[value]-rgetattr(self,x+'.spec'))**2/(rgetattr(self,x+'.uncs'))**2*rgetattr(self,x+'.limit')*rgetattr(self,x+'.dip_array'))
                    else:
                        probability+= -np.sum((rgetattr(self,x+'.synth')-rgetattr(self,x+'.spec'))**2/(rgetattr(self,x+'.uncs'))**2*rgetattr(self,x+'.limit')*rgetattr(self,x+'.dip_array'))
            else:
                for value,x in enumerate(colours):
                    if not(np.array_equal(synthetic_spectra,None)):
                        probability+= -np.sum((synthetic_spectra[value]-normal[value])**2/(rgetattr(self,x+'.uncs'))**2*rgetattr(self,x+'.limit')*rgetattr(self,x+'.dip_array'))
                    else:
                        probability+= -np.sum((rgetattr(self,x+'.synth')-normal[value])**2/(rgetattr(self,x+'.uncs'))**2*rgetattr(self,x+'.limit')*rgetattr(self,x+'.dip_array'))
        return probability
    def observed_spectra_giver(self,colours=bands):
        returning_spectra=np.array(np.ones(len(colours)),dtype=object)
        wavelengths=np.array(np.ones(len(colours)),dtype=object)
        uncs=np.array(np.ones(len(colours)),dtype=object)
        for value,x in enumerate(colours):
            returning_spectra[value]=rgetattr(self,x+'.spec')
            wavelengths[value]=rgetattr(self,x+'.wave')
            uncs[value]=rgetattr(self,x+'.uncs')
        return returning_spectra,wavelengths,uncs
    def plot(self,colours=bands):
        for x in colours:
            plt.figure()
            x_line=np.linspace(0,len(runGetattr(self,x+'.synth'))-1,num=len(runGetattr(self,x+'.synth')))
            if rgetattr(self,x+'.synth').any():
                    # labels='synthetic  chi squared= '+str(self.log_fit([x]))
                    plt.plot(x_line, runGetattr(self,x+'.synth'), label='Synthetic',c='r')
                    
            plt.fill_between(x_line, runGetattr(self,x+'.spec')- runGetattr(self,x+'.uncs'),runGetattr(self,x+'.spec')+ runGetattr(self,x+'.uncs'),alpha=0.5)
            plt.plot(x_line, runGetattr(self,x+'.spec'), label='observing',c='b')
            plt.title(x+" Band")
            plt.legend(loc='best')
            plt.tight_layout()
            
def starter_walkers_maker(nwalkers):
    """
    Creates an 9 x nwalkers dimentional array thats is a good starter posistion for the walkers 

    Parameters
    ----------
    nwalkers : float 

    Returns
    -------
     the 9xn dimentional array 
    """
    pos=[]
    # labels=['teff','logg','monh','Fe','alpha','vsini','vmac','vmic']
    while np.shape(pos)[0]<nwalkers:
        if old_abundances['teff'] and not np.isnan( old_abundances['teff']) and old_abundances['teff']>x_min[0] and old_abundances['teff']<x_max[0]:
            teff_random=abs(np.random.normal(old_abundances['teff'],1000,1))
        else:
            teff_random=np.random.normal(5000,100,1)
        if old_abundances['logg'] and not np.isnan( old_abundances['logg'])and old_abundances['logg']>x_min[1] and old_abundances['logg']<x_max[1]:
            logg_random=abs(np.random.normal(old_abundances['logg'],0.4,1))
        else:
            logg_random=np.random.normal(4.5,1,1)
            
        any_radial_velocity=np.isnan(old_abundances['rv_r'])
        
        if not np.isnan(old_abundances['rv_Blue_r']):
            v_rad_random_Blue=np.random.normal(old_abundances['rv_Blue_r'],1,1)
        elif not any_radial_velocity:
            v_rad_random_Blue=np.random.normal(old_abundances['rv_r'],1,1)
        else:
            v_rad_random_Blue=np.random.normal(0,4,1)
            
        if not np.isnan(old_abundances['rv_Green_r']):
            v_rad_random_Green=np.random.normal(old_abundances['rv_Green_r'],1,1)
        elif not any_radial_velocity:
            v_rad_random_Green=np.random.normal(old_abundances['rv_r'],1,1)
        else:
            v_rad_random_Green=np.random.normal(0,4,1)
            
        if not np.isnan(old_abundances['rv_Red_r']):
             v_rad_random_Red=np.random.normal(old_abundances['rv_Red_r'],1,1)
        elif not any_radial_velocity:
             v_rad_random_Red=np.random.normal(old_abundances['rv_r'],1,1)
        else:
             v_rad_random_Red=np.random.normal(0,4,1)           

        if not np.isnan(old_abundances['rv_IR_r']):
             v_rad_random_IR=np.random.normal(old_abundances['rv_IR_r'],1,1)
        elif not any_radial_velocity:
             v_rad_random_IR=np.random.normal(old_abundances['rv_r'],1,1)
        else:
             v_rad_random_IR=np.random.normal(0,4,1) 

        if old_abundances['fe_h']  and not np.isnan( old_abundances['fe_h']) and abs(old_abundances['fe_h'])<2.0 and old_abundances['fe_h']>x_min[3] and old_abundances['fe_h']<x_max[3]:
            monh_random=np.random.normal(old_abundances['fe_h'],0.2,1)
        else:
            monh_random=np.random.normal(0,0.1,1)
        if old_abundances['vmic_reduction'] and not np.isnan( old_abundances['vmic_reduction']) and old_abundances['vmic_reduction']>x_min[7] and old_abundances['vmic_reduction']<x_max[7]:
            v_mac=abs(abs(np.random.normal(old_abundances['vmic_reduction'],8,1))+np.random.randint(-1,1))
            v_mic=abs(np.random.normal(old_abundances['vmic_reduction'],8,1))
        else:
            v_mic=np.random.normal(10,1,1)
            v_mac=np.random.normal(10,1,1)
        if old_abundances['vbroad_reduction'] and not np.isnan( old_abundances['vbroad_reduction']) and old_abundances['vbroad_reduction']<100 and old_abundances['vbroad_reduction']>x_min[5] and old_abundances['vbroad_reduction']<x_max[5]:
            vsini_random=abs(np.random.normal(old_abundances['vbroad_reduction'],1*4,1))
        else:
            vsini_random=np.random.normal(10,1,1)
        if old_abundances['alpha_fe_r'] and not np.isnan( old_abundances['alpha_fe_r'])  and old_abundances['alpha_fe_r']>x_min[4] and old_abundances['alpha_fe_r']<x_max[4]:
            alpha_random=np.random.normal(old_abundances['alpha_fe_r'],0.1,1)
        else:
            alpha_random=np.random.normal(0.148,0.317,1)
        fe_random=monh_random+np.random.normal(0,0.1,1)
    
    
        pos_try=np.hstack((teff_random,logg_random,monh_random,fe_random,alpha_random,v_rad_random_Blue,v_rad_random_Green,v_rad_random_Red,v_rad_random_IR,vsini_random,v_mac,v_mic))
        if starting_test(pos_try):
            if len(pos)==0:
                pos=pos_try
            else:
                pos=np.vstack((pos_try,pos))
    return pos
global old_abundances
colours=['Blue','Green','Red','IR']
large_data=fits.getdata('photometric_reduction_dr6_cross.fits',1)
large_data=Table(large_data)
global x_min,x_max
tmp = np.load("NN_normalized_spectra_fe_alpha_shortest_Blue.npz")
x_min=tmp['x_min']
x_max=tmp['x_max']
path=r'machine_learning_long/'

directories=os.listdir(path)

names=[x[:15] for x in directories]
names=np.sort(names)
double_names=[]
[double_names.append(x) for x in names if collections.Counter(names)[x]>1.0 and x not in double_names]

# bands=['IR']
#Needs to be read because some spectra doesnt have 
global all_reduced_data
all_reduced_data=fits.getdata('dr6.0.fits',1)
all_reduced_data=Table(all_reduced_data)
# interpolation=10
svens_mine=fits.getdata('svens_fe_alpha_with_mine.fits',1)
svens_mine=Table(svens_mine)
min_fit=0.01
def spectra_makers(shift):
    return spectras.synthesize(shift,give_back=True)
for x in double_names:
    sobject_id=x
    file_name='machine_learning_long/'+sobject_id+'_old_moves_machine_learning_prior.h5'
    file_name2='machine_learning_long/'+sobject_id+'_old_moves_machine_learning_no_prior.h5'
    # if not (exists(file_name) and exists(file_name2)):
    #     print('File doesnt exists Skipping ', sobject_id)
    #     continue
    #if exists('machine_learning_results_alpha_fe_svens_new/spectras/'+sobject_id+'4.fits'):
    #	print('already synthesized', sobject_id)
    #	continue
        
        
    old_abundances=[y for y in large_data if y['sobject_id']==str(sobject_id)]
    # if len(old_abundances)==0:
    #  	print('spectra doesnt exists skipping'+sobject_id)
    #  	continue
    old_abundances=vstack(old_abundances)

    print('synthesizing '+sobject_id)
    # spectras=spectrum_all(name,2)
    # spectras.synthesize({})
    spectras=spectrum_all(sobject_id,10)
    
    
    # spec_janez=spectras.synthesize({},give_back=True)
    before_spectras=spectras.synthesize({},give_back=True)
    spectras.synthesize({})
    spectras.normalize()
    
    # training_data=np.load('training_data_30.npy',allow_pickle=True)
    # which=10
    # validation=np.load('Validation_data_all_no_v_rad_fixed.npy',allow_pickle=True)
    # validation_labels=np.vstack(validation[:,1])
    # shift=shift_maker(validation_labels[which])
    # validation_spectra=np.vstack(np.vstack(validation[:,0]))
    
    # synthesized_spectra=spectras.synthesize(shift,give_back=True)
    
    cut=20000
    reader= emcee.backends.HDFBackend('machine_learning_long/'+sobject_id+'_old_moves_machine_learning_prior.h5')
    fitted_data=reader.get_chain(flat=True)[72*20:]
    cut=1000
    total_cut=0
    fit=how_good_gaussian_fit(fitted_data)
    print(fit)
    max_cut=0.5*len(fitted_data)
    while fit>min_fit and total_cut<max_cut:
        
        fitted_data=fitted_data[cut:]
        total_cut+=cut
        fit=how_good_gaussian_fit(fitted_data)
    fitted_data_mean=np.mean(fitted_data,axis=0)
    shift=shift_maker(fitted_data_mean)
    # data_to_synthesize=fitted_data[::len(fitted_data)//100]
    # with Pool(12) as pool:
    #     prior_data_full=pool.map(spectra_makers,data_to_synthesize)
    synthetic_prior=spectras.synthesize(shift,give_back=True)
    observed_spectra,wavelengths,uncs=spectras.observed_spectra_giver()
    
    reader= emcee.backends.HDFBackend('machine_learning_long/'+sobject_id+'_old_moves_machine_learning_no_prior.h5')
    fitted_data=reader.get_chain(flat=True)[72*20:]
    cut=1000
    total_cut=0
    fit=how_good_gaussian_fit(fitted_data)
    print(fit)
    max_cut=0.5*len(fitted_data)
    while fit>min_fit and total_cut<max_cut:
        
        fitted_data=fitted_data[cut:]
        total_cut+=cut
        fit=how_good_gaussian_fit(fitted_data)
    fitted_data_mean=np.mean(fitted_data,axis=0)
    shift=shift_maker(fitted_data_mean)
    # data_to_synthesize=fitted_data[::len(fitted_data)//100]
    # with Pool(12) as pool:
    #     no_prior_data_full=pool.map(spectra_makers,data_to_synthesize)
    no_prior_data=spectras.synthesize(shift,give_back=True)
    to_save=[]
    keys=['wavelength','observed_spectra','synthetic_spectra_with_prior','synthetic_spectra_without_prior','uncs','original_synthetic_spectra']
    for value,y in enumerate(colours,1):
        table_to_save=Table({keys[0]:wavelengths[value-1]})
        table_to_save[keys[1]]=observed_spectra[value-1]
        table_to_save[keys[2]]=synthetic_prior[value-1]
        table_to_save[keys[3]]=no_prior_data[value-1]
        table_to_save[keys[4]]=uncs[value-1]
        table_to_save[keys[5]]=before_spectras[value-1]
        # table_to_save.write('machine_learning_long/spectras/'+sobject_id+str(y+1)+'.fits',overwrite=True)
        
        plt.plot(table_to_save['wavelength'], table_to_save['observed_spectra'],color='green',label='observed_spectra')
        plt.plot(table_to_save['wavelength'],table_to_save['synthetic_spectra_with_prior'],color='blue',label='synthetic_spectra_prior')
        plt.plot(table_to_save['wavelength'],table_to_save['synthetic_spectra_without_prior'],color='red',label='synthetic_spectra_without_prior')
        # plt.plot(table_to_save['wavelength'],np.vstack(np.array(prior_data_full)[:,value-1]).T,alpha=0.1,c='blue')
        plt.legend(loc='best')
        plt.title(y +' '+str(x))
        plt.xlabel('wavelength')
        plt.show()
        plt.draw()   
        # input(('drawed ' ,y))
        # plt.figure()        
