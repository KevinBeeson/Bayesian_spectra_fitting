#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 16:55:30 2022

@author: kevin
"""

import pymc as pm 
import aesara.tensor as at
import numpy as np
import arviz as az
from numba import jit

from functools import  partial
from astropy.io.votable import parse
# from The_Payne import spectral_model
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
def shift_maker(solar_value,given_labels=['teff','logg','monh','Fe','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']):
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
        skip=0
        for value,x in enumerate(labels):
            if x in given_labels:
                shift[x]=solar_value[value-skip]
            else:
                skip+=1
                bands=spectras.bands
                colour=[y for y in bands if y in x]

                if not 'vrad' in x:
                    shift[x]=rgetattr(spectras,bands[0]+'.'+x)
                elif len(colour): 
                    shift[x]=rgetattr(spectras, colour[0]+'.vrad')
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

def log_posterior(solar_values,parameters,prior=False):
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
        full_parameters=['teff','logg','monh','Fe','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']
        if len(solar_values)!= len(parameters):
            print('youre missing parameters')
        
        # if not starting_test(solar_values):
        #       return -np.inf
        # print('passed')
        shift=shift_maker(solar_values,parameters)
        if not 'Fe' in parameters: 
            shift['Fe']=shift['monh']
        # shift['teff']=(shift['teff']+6000)*10
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
        probability=spectras.log_fit(synthetic_spectra=synthetic_spectras,solar_shift=shift)
    
        if prior:
            probability+=log_prior(shift)
        if probability>0:
            print('oh no prob ',probability)
            print(solar_values)
            # return False
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
def leaky_relu(z,grad=False):
    '''
    This is the activation function used by default in all our neural networks.
    '''
    limits=(z > 0)
    leaky=z*limits + 0.01*z*np.invert(limits)
    if grad:
        return leaky,limits
    return leaky
def payne_sythesize(solar_values,x_min,x_max,NN_coeffs,grad=False):
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
        w_array_0, w_array_1, w_array_2, b_array_0, b_array_1, b_array_2, x_min, x_max = NN_coeffs
        if grad is False:
            
            inside =np.einsum('ij,j->i', w_array_0, scaled_labels)+ b_array_0
            inside1=np.einsum('ij,j->i', w_array_1, leaky_relu(inside))+ b_array_1
            real_spec=np.einsum('ij,j->i', w_array_2, leaky_relu(inside1))+ b_array_2
            return real_spec

        else:
            inside =np.einsum('ij,j->i', w_array_0, scaled_labels)+ b_array_0
            limit_1,leaked=leaky_relu(inside,True)
            # leaked=np.ones(len(w_array_0))
            # grad_return=w_array_0.T*leaked 
            # grad_return=w_array_0.T
            
            grad_return=w_array_0.T*leaked + 0.01*w_array_0.T*np.invert(leaked)
            
            inside =np.einsum('ij,j->i', w_array_1, limit_1)+ b_array_1
            limit_2,leaked=leaky_relu(inside,True)
            grad_return=np.dot(grad_return,w_array_1.T)
            # leaked=np.ones(len(grad_return.T))
            # grad_return=grad_return*leaked 
            grad_return=grad_return*leaked + 0.01*grad_return*np.invert(leaked)
            
            real_spec =np.einsum('ij,j->i', w_array_2,limit_2)+ b_array_2
            grad_return=np.dot(grad_return,w_array_2.T)
            return real_spec,grad_return
        # real_spec = spectral_model.get_spectrum_from_neural_net(scaled_labels = scaled_labels, NN_coeffs = NN_coeffs)
def galah_kern(fwhm, b):
        """ Returns a normalized 1D kernel as is used for GALAH resolution profile """
        size=2*(fwhm/2.355)**2
        size_grid = int(size) # we limit the size of kernel, so it is as small as possible (or minimal size) for faster calculations
        if size_grid<7: size_grid=7
        x= scipy.mgrid[-size_grid:size_grid+1]
        g = scipy.exp(-0.693147*np.power(abs(2*x/fwhm), b))
        return g / np.sum(g)
def dopler(original_wavelength,spectra,v_rad,synthetic_wavelength=None,grad=None):
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
        if len(original_wavelength)!=len(spectra) and synthetic_wavelength is None:
            difference=abs(len(original_wavelength)-len(spectra))
            original_wavelength_strenched=np.interp(np.linspace(-difference/2,difference/2+len(original_wavelength),num=len(spectra)),range(len(original_wavelength)),original_wavelength)
            
            observed_wavelength=original_wavelength_strenched*(1+delta_v)
        elif not synthetic_wavelength is None:
            observed_wavelength=synthetic_wavelength*(1+delta_v)
        else:
            observed_wavelength=original_wavelength*(1+delta_v)
            
        #crops and interpolates 
        spectra_new=np.interp(original_wavelength,observed_wavelength,spectra)
        if not grad is None:
            grad_new=[np.interp(original_wavelength,observed_wavelength,x) for x in grad]
            return spectra_new,grad_new
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
def synth_resolution_degradation(wave_synth, synth,res_map,res_b,wave_original,l_new_premade=None,kernel_=None,synth_res=300000.0,grad=None):
        """
        Take a synthetic spectrum with a very high  resolution and degrade its resolution to the resolution profile of the observed spectrum. The synthetic spectrum should not be undersampled, or the result of the convolution might be wrong.
        Parameters:
            synth np array or similar: an array representing the synthetic spectrum. Must have size m x 2. First column is the wavelength array, second column is the flux array. Resolution of the synthetic spectrum must be constant and higher than that of the observed spectrum.
            synth_res (float): resolving power of the synthetic spectrum
        Returns:
            Convolved syntehtic spectrum as a np array of size m x 2.
        """
        synth=np.vstack((wave_synth,synth)).T
        
        l_original=wave_synth
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
        galah_sampl=wave_original[1]-wave_original[0]


        #original sigma
        s_original=sigma_synth





        #oversampling. If synthetic spectrum sampling is much finer than the size of the kernel, the code would work, but would return badly sampled spectrum. this is because from here on the needed sampling is measured in units of sigma.
        oversample=galah_sampl/sampl*5.0

        #minimal needed sampling

        #keep adding samples until end of the wavelength range is reached
        if l_new_premade is None:
            #required sigma (resample the resolution map into the wavelength range of the synthetic spectrum)
            s_out=np.interp(synth[:,0], wave_original, res_map)


            #the sigma of the kernel is:
            s=np.sqrt(s_out**2-s_original**2)

            #fit it with the polynomial, so we have a function instead of sampled values:
            map_fit=np.poly1d(np.polyfit(synth[:,0], s, deg=6))

            #create an array with new sampling. The first point is the same as in the spectrum:
            l_new=[synth[:,0][0]]

            min_sampl=max(s_original)/sampl/sampl*oversample

            l_new=numba_syth_resolution(map_fit.coef,l_new,sampl,min_sampl,synth[:,0][-1])
            kernel_=galah_kern(max(s_original)/sampl*oversample, res_b)

        else:
            l_new=np.copy(l_new_premade)
        # while l_new[-1]<synth[:,0][-1]+sampl:
        #     l_new.append(l_new[-1]+map_fit(l_new[-1])/sampl/min_sampl)

        #interpolate the spectrum to the new sampling:
        new_f=np.interp(np.array(l_new),synth[:,0],synth[:,1])

        
        con_f=signal.fftconvolve(new_f,kernel_,mode='same')

        #inverse the warping:
        synth[:,1]=np.interp(l_original,np.array(l_new),con_f)
        if l_new_premade is None:
            return synth[:,1],l_new,kernel_
        if not grad is None:
            new_grad=[np.interp(np.array(l_new),synth[:,0],x) for x in grad]
            con_grad=[signal.fftconvolve(x,kernel_,mode='same') for x in new_grad]
            grad=[np.interp(l_original,np.array(l_new),x) for x in con_grad]
            return synth[:,1],grad

        return synth[:,1]

class individual_spectrum:
    
    def __init__(self,name,interpolation,x,count,old_abundances,cluster=True):
            limits={'Blue':[4713,4908],'Green':[5643,5879],'Red':[6474.0,6746.0],'IR':[7730.0,7900.0]}

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
            tmp = np.load("NN_normalized_spectra_no_fe_"+x+".npz")
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
            self.grad=(w_array_0, w_array_1, w_array_2)
            self.NN_coeff=NN_coeffs
            self.x_min=x_min
            self.x_max=x_max
            Path(name[0:6]+'/spectra/com/').mkdir(parents=True,exist_ok=True)
            if not exists(name[0:6]+'/spectra/com/'+name+str(count)+'.fits'):
                #source='/media/storage/HERMES_REDUCED/dr6.0/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'
                source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.0/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'

                destination=name[0:6]+'/spectra/com/'
                subprocess.run(["rsync",'-av',source,destination])

            try:
                hermes=fits.open(name[0:6]+'/spectra/com/'+name+str(count)+'.fits')
                x0= float(hermes[1].header.get('CRVAL1'))
                x1=float( hermes[1].header.get('CRVAL1')+len(hermes[1].data)* hermes[1].header.get('CDELT1'))
                fstart= x0+(x1-x0)*starting_fraction
                
                new_start=int(starting_fraction*len(hermes[1].data))
                new_end=new_start+int(len(hermes[1].data)*(1-starting_fraction)*length_fraction)
                
                
                length=np.linspace(x0+(x1-x0)*starting_fraction,x1,num=int(4096*(1-starting_fraction)))
                length_synthetic=np.linspace(limits[x][0], limits[x][1],num=int((4096*(1-starting_fraction)+160)*interpolation))
                fend=length[-1]
                #Increases the length of synthetic spectra so it over interpolates  
                self.wave_synth=length_synthetic
                self.wave=length
                self.spec=[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])]
                self.spec_original=[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])][0]

                self.uncs=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[2].data[new_start:new_end])
                hermes[0].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[0].data[new_start:new_end])
                hermes[1].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])
                hermes[2].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[2].data[new_start:new_end])
                # not_resolved=[np.int64(name)]
                try:
                    if len(hermes[7].data[new_start:new_end])==new_end-new_start and min(hermes[7].data)>0:
                        hermes[7].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[7].data[new_start:new_end])
                except TypeError:
                    hermes[7].data=[0]
                    
                # hermes[1].header['CDELT1']/=interpolation
                self.wran=[fstart,fend]
                self.hermes=hermes
                if cluster:
                    self.teff=float(old_abundances['teff'])
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
                else:
                    self.teff=float(old_abundances['teff_r'])
                    self.vmic=float(old_abundances['vmic_r'])
                    self.vsini=float(old_abundances['vbroad_r'])
                    self.Fe=float(old_abundances['fe_h_r'])
                    self.alpha=float(old_abundances['alpha_fe_r'])
                    if not np.isnan(float(old_abundances['rv'][count-1])):
                        self.vrad=float(old_abundances['rv'][count-1])
                    elif not np.isnan(float(old_abundances['rv_com'])):
                        self.vrad=float(old_abundances['rv_com'])
                    else:
                        self.vrad=0.0
                    self.vmac=6.0
                    self.logg=float(old_abundances['logg_r'])
                    self.monh=float(old_abundances['fe_h_r'])
                if hermes[1].header.get('TEFF_R')=='None' or hermes[1].header.get('LOGG_R')=='None':
                    print('Warning the reduction didnt produce an TEFF spectra might not reduce properly')
                    print('enter to continue')
                    self.bad_reduction=True
                else:
                    self.bad_reduction=False

            except FileNotFoundError:
                print('No hermes found')
                self.hermes=None
                self.bad_reduction=True
            


            self.normal_value=None

            self.l_new=None
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

    def __init__(self,name,interpolation=10,cluster=True):
        name=str(name)
        bands=['Blue','Green','Red','IR']
        self.rv_shift=1e-10
        self.bands=bands
        self.name=name
        self.interpolation=interpolation
        self.sister_stars=None
        if cluster:
            old_abundances=[x for x in large_data if x['sobject_id']==str(name)]
            old_abundances=vstack(old_abundances)
        else:

            old_abundances=[x for x in all_reduced_data if x['sobject_id']==np.int64(name)]
            old_abundances=vstack(old_abundances)
            old_abundances=old_abundances[0]
        for count,x in enumerate(bands,1):
            setattr(self,x,individual_spectrum(name,interpolation,x,count,old_abundances,cluster))
        self.correct_resolution_map()
        self.equilize_spectras()
        self.limit_array()
        self.old_abundances=old_abundances 
    def equilize_spectras(self,colours=bands):
        for x in colours:
            if rgetattr(self, x+'.hermes')!=None:
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
            if rgetattr(self, x+'.hermes')!=None:
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

    def synthesize(self,shift={},colours=bands,multi=False,give_back=False,full=False,grad=False):
        if full and not give_back:
            print('you probably want the spectrum back run again with give_back=True')
            return False
        if not give_back:
            if not multi:
                for x in colours:       
                    solar_values=self.solar_value_maker(shift,x)
                    spectrum = payne_sythesize(solar_values,rgetattr(self,x+'.x_min'),rgetattr(self,x+'.x_max'),rgetattr(self,x+'.NN_coeff'))
                    if 'vrad_'+x in shift:
                        dopler_shifted_spectra=dopler(rgetattr(self,x+'.wave'),spectrum,shift['vrad_'+x],synthetic_wavelength=rgetattr(self,x+'.wave_synth'))
                    else:
                        dopler_shifted_spectra=dopler(rgetattr(self,x+'.wave'),spectrum,rgetattr(self,x+'.vrad'),synthetic_wavelength=rgetattr(self,x+'.wave_synth'))
                    if rgetattr(self,x+'.l_new') is None:
                        dopler_shifted_spectra,l_new,kernel=synth_resolution_degradation(rgetattr(self,x+'.wave'),dopler_shifted_spectra,rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'],rgetattr(self,x+'.wave'))
                        rsetattr(self,x+'.l_new',l_new)
                        rsetattr(self,x+'.kernel',kernel)
                    else:
                        dopler_shifted_spectra=synth_resolution_degradation(rgetattr(self,x+'.wave'),dopler_shifted_spectra,rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'],rgetattr(self,x+'.wave'),rgetattr(self,x+'.l_new'),rgetattr(self,x+'.kernel'))
                    rsetattr(self,x+'.synth',dopler_shifted_spectra)
            else:
                with Pool(4) as pool:
                    inputs=[(shift,x,False,False) for x in colours]
                    pool.map(getattr(self,'synthesize'),inputs)
                # with ThreadPool(4) as pool:
                #     inputs=[(shift,x,False,False) for x in colours]
                #     pool.map(getattr(self,'synthesize'),inputs)
                
        else:
            if not multi:
                if grad:
                    
                    returning_spectra=np.array(np.ones(len(colours)),dtype=object)
                    returning_grad=np.array(np.ones(len(colours)),dtype=object)
                    for number,x in enumerate(colours):       
                        solar_values=self.solar_value_maker(shift,x)

                        spectrum,grad_spec = payne_sythesize(solar_values,rgetattr(self,x+'.x_min'),rgetattr(self,x+'.x_max'),rgetattr(self,x+'.NN_coeff'),grad=True)
                        if full:
                            returning_spectra[number]=spectrum
                            continue   
                        if 'vrad_'+x in shift:
                            dopler_shifted_spectra,grad_spec=dopler(rgetattr(self,x+'.wave'),spectrum,shift['vrad_'+x],synthetic_wavelength=rgetattr(self,x+'.wave_synth'),grad=grad_spec)
                        else:
                            dopler_shifted_spectra,grad_spec=dopler(rgetattr(self,x+'.wave'),spectrum,rgetattr(self,x+'.vrad'),synthetic_wavelength=rgetattr(self,x+'.wave_synth'),grad=grad_spec)

                        if rgetattr(self,x+'.l_new') is None:
                            
                            dopler_shifted_spectra,l_new,kernel,grad_dopler=synth_resolution_degradation(rgetattr(self,x+'.wave'),dopler_shifted_spectra,rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'],rgetattr(self,x+'.wave'),grad=grad_spec)
                            rsetattr(self,x+'.l_new',l_new)
                            rsetattr(self,x+'.kernel',kernel)

                        else:
                            dopler_shifted_spectra,grad_dopler=synth_resolution_degradation(rgetattr(self,x+'.wave'),dopler_shifted_spectra,rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'],rgetattr(self,x+'.wave'),rgetattr(self,x+'.l_new'),rgetattr(self,x+'.kernel'),grad=grad_spec)
                        returning_grad[number]=grad_dopler
                        returning_spectra[number]=dopler_shifted_spectra
                    return returning_spectra,returning_grad
                else:
                    returning_spectra=np.array(np.ones(len(colours)),dtype=object)
                    for number,x in enumerate(colours):       
                        solar_values=self.solar_value_maker(shift,x)
    
                        spectrum = payne_sythesize(solar_values,rgetattr(self,x+'.x_min'),rgetattr(self,x+'.x_max'),rgetattr(self,x+'.NN_coeff'))
                        if full:
                            returning_spectra[number]=spectrum
                            continue                   
                        if 'vrad_'+x in shift:
                            dopler_shifted_spectra=dopler(rgetattr(self,x+'.wave'),spectrum,shift['vrad_'+x],synthetic_wavelength=rgetattr(self,x+'.wave_synth'))
                        else:
                            dopler_shifted_spectra=dopler(rgetattr(self,x+'.wave'),spectrum,rgetattr(self,x+'.vrad'),synthetic_wavelength=rgetattr(self,x+'.wave_synth'))
                        if rgetattr(self,x+'.l_new') is None:
                            dopler_shifted_spectra,l_new,kernel=synth_resolution_degradation(rgetattr(self,x+'.wave'),dopler_shifted_spectra,rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'],rgetattr(self,x+'.wave'))
                            rsetattr(self,x+'.l_new',l_new)
                            rsetattr(self,x+'.kernel',kernel)
    
                        else:
                            dopler_shifted_spectra=synth_resolution_degradation(rgetattr(self,x+'.wave'),dopler_shifted_spectra,rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'],rgetattr(self,x+'.wave'),rgetattr(self,x+'.l_new'),rgetattr(self,x+'.kernel'))
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
       correct_spectra=[]
       # print(not_resolved)

       for count,x in enumerate(colours,1):
           if rgetattr(self,x+'.hermes')!=None:
               hermes=rgetattr(self,x+'.hermes')
               if len(hermes[7].data)!=len(hermes[1].data) or min(hermes[7].data)<0:
                   correct_spectra.append(0)
                   tried=True
                   print('Resolution map of ',name +str(count),' is wrong')
                   pivot_sisters=[]
                   temp_data=[x for x in all_reduced_data if x['sobject_id']==np.int64(name)]
                   temp_data=temp_data[0]
                   plate=temp_data['plate']
                   time=temp_data['utdate']
                   epoch=temp_data['epoch']
                   pivot_number=temp_data['pivot']
                   min_difference=0
                   while len(hermes[7].data)!=len(hermes[1].data) or min(hermes[7].data)<=0:
                       if self.sister_stars==None or tried==False:
   
                           if pivot_number!=0 and min_difference*365<15 and len(not_resolved)<15:
                               if len(pivot_sisters)==0:
                                   pivot_sisters=[x for x in all_reduced_data if x['pivot']==pivot_number and x['plate']==plate and not x['sobject_id'] in not_resolved and abs(x['epoch']-epoch)<0.1 and x['res'][count-1]>0]
                               else:
                                   pivot_sisters=[x for x in pivot_sisters if x['pivot']==pivot_number and x['plate']==plate and not x['sobject_id'] in not_resolved and abs(x['epoch']-epoch)<0.1 and x['res'][count-1]>0]
                               if len(pivot_sisters)>0:
                                   pivot_sisters=vstack(pivot_sisters)
                                   # difference=difference_UT_vector(time, pivot_sisters['utdate'])
                                   pivot_sisters['difference']=abs(pivot_sisters['epoch']-epoch)
                                   pivot_sisters.sort('difference')
                                   name_target=str(pivot_sisters[0]['sobject_id'])
                                   min_difference=pivot_sisters[0]['difference']
                                   # print('here')
                                   print('copying from '+name_target+' for '+ name)
                                   not_resolved.append(np.int64(name_target))
                                   # print(not_resolved)
                                   Path(name_target[0:6]+'/spectra/com/').mkdir(parents=True,exist_ok=True)
                                   if not exists(name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'):
                                       # source='/media/storage/HERMES_REDUCED/dr6.0/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'
                                       source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.0/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'
                   
                                       destination=name_target[0:6]+'/spectra/com/'
                                       subprocess.run(["rsync",'-av',source,destination])
                                   try:
                                       if os.path.getsize(name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits')>380000:
                                           hermes_temp=fits.open(name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits')
                                           if x=='IR':
                                               try:
                                                   hermes[7].data=hermes_temp[7].data[-2048:]
                                               except TypeError:
                                                   print('cant be openned')
                                           else:
                                               try:
                                                   hermes[7].data=hermes_temp[7].data
                                               except TypeError:
                                                   print('cant be openned')
                                       else:
                                           print('file is too small')
                                   except FileNotFoundError:
                                       print('hermes_error')
                               else:
                                   print('no more sister stars ', name)
                                   min_difference=30
                           else:
                               print('pivot is ',pivot_number,' setting res_map to 0.4 A and the difference is ',min_difference*365)
                               hermes[7].data=0.40*(np.ones(len(hermes[1].data)))
                               name_target=None
   
   
   
                       elif tried:
                           tried=False
                           name_target=self.sister_stars
                           try:
                               Path(name_target[0:6]+'/spectra/com/').mkdir(parents=True,exist_ok=True)
                               if not exists(name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'):
                                   # source='/media/storage/HERMES_REDUCED/dr6.0/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'
                                   source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.0/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'
               
                                   destination=name_target[0:6]+'/spectra/com/'
                                   subprocess.run(["rsync",'-av',source,destination])
       
                               if os.path.getsize(name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits')>380000:
                                   hermes_temp=fits.open(name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits')
                                   if x=='IR':
                                       try:
                                           hermes[7].data=hermes_temp[7].data[-2048:]
                                       except TypeError:
                                           print('cant be openned')
                                   else:
                                       try:
                                           hermes[7].data=hermes_temp[7].data
                                       except TypeError:
                                           print('cant be openned')
                               else:
                                   print('file is too small')
                           except FileNotFoundError:
                               print('file not found')
                   self.sister_stars=name_target
               else:
                   correct_spectra.append(1)
           else:
               correct_spectra.append(0)
       self.resolution_spectra=correct_spectra
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
    def gradient(self,shift={},colours=bands,self_normal=False,labels=['teff','logg','monh','Fe','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']):
        labels_payne=['teff','logg','monh','Fe','alpha','vsini','vmac','vmic']
        labels_payne_front=['teff','logg','monh','Fe','alpha']
        vrad_labels=['vrad_Blue','vrad_Green','vrad_Red','vrad_IR']
        labels_payne_front_wanted=[x for x in labels_payne_front if x in labels]
        labels_payne_wanted=[x for x in labels_payne if x in labels]

        rv_wanted=[x for x in vrad_labels if x in labels]
        synthetic_spectra,grad_all=self.synthesize(shift,give_back=True,grad=True)
        grad_out=np.ones(len(colours),dtype=object)
        if len(labels)-len(rv_wanted)!=0:
            for count,(syn,grad,x,current_label) in enumerate(zip(synthetic_spectra,grad_all,colours,labels_payne)):
                if self_normal:
                    grad_out[count]=-np.sum(np.array(grad)*rgetattr(self,x+'.limit')*(syn-rgetattr(self,x+'.spec'))/(rgetattr(self,x+'.uncs')**2),axis=1)/(rgetattr(self,x+'.x_max')-rgetattr(self,x+'.x_min'))
                else:
                    grad_out[count]=-np.sum(np.array(grad)*rgetattr(self,x+'.limit')*(syn-rgetattr(self,x+'.spec'))/(rgetattr(self,x+'.uncs')**2),axis=1)/(rgetattr(self,x+'.x_max')-rgetattr(self,x+'.x_min'))
            grad_out=np.sum(np.array(grad_out),axis=0)

        if len(rv_wanted)!=0:
            count=0
            grad_out_rv=np.ones(len(rv_wanted),dtype=object)
            for syn,x,current_labels in zip(synthetic_spectra,colours,vrad_labels):
                if current_labels in labels:
                    shift_temp=copy.copy(shift)
                    if 'vrad_'+x in shift_temp:
                        shift_temp['vrad_'+x]-=self.rv_shift
                    else:
                        shift_temp['vrad_'+x]=rgetattr(self,x+'.vrad')-self.rv_shift
                    grad_low=self.log_fit([x],shift_temp)
                    if 'vrad_'+x in shift_temp:
                        shift_temp['vrad_'+x]+=2*self.rv_shift
                    else:
                        shift_temp['vrad_'+x]=rgetattr(self,x+'.vrad')+self.rv_shift
                    grad_high=self.log_fit([x],shift_temp)
                    grad_out_rv[count]=(grad_high-grad_low)/(2*self.rv_shift)
                    count+=1
        if len(labels)-len(rv_wanted)!=0:
            grad_out_new=[x for x,y in zip(grad_out,labels_payne) if y in labels_payne_wanted]

            if len(rv_wanted)!=0:
                grad_out_new=np.array(np.hstack((grad_out_new[:len(labels_payne_front_wanted)],grad_out_rv,grad_out_new[len(labels_payne_front_wanted):])),dtype=np.float64)
        else:
            grad_out_new=grad_out_rv
        # print(grad_out_new)
        return np.array(grad_out_new,dtype=np.float64)
    def log_fit(self,colours=bands,shift=None,dip_array=False,synthetic_spectra=None,normal=None,self_normal=False,solar_shift=None):
        if not shift is None:
            synthetic_spectra=self.synthesize(shift,colours=colours,give_back=True)
        probability=np.ones(len(colours))
        if not dip_array:
            if np.array_equal(normal,None):
                for value,x in enumerate(colours):
                    if not(np.array_equal(synthetic_spectra,None)):
                        probability[value]= -np.sum((synthetic_spectra[value]-rgetattr(self,x+'.spec'))**2/(2*(rgetattr(self,x+'.uncs'))**2)*rgetattr(self,x+'.limit'))
                    else:
                        probability[value]= -np.sum((rgetattr(self,x+'.synth')-rgetattr(self,x+'.spec'))**2/(2*(rgetattr(self,x+'.uncs'))**2)*rgetattr(self,x+'.limit'))
            else:
                for value,x in enumerate(colours):
                    if not(np.array_equal(synthetic_spectra,None)):
                        probability[value]= -np.sum((synthetic_spectra[value]-normal[value])**2/(2*(rgetattr(self,x+'.uncs'))**2)*rgetattr(self,x+'.limit'))
                    else:
                        probability[value]= -np.sum((rgetattr(self,x+'.synth')-normal[value])**2/(2*(rgetattr(self,x+'.uncs'))**2)*rgetattr(self,x+'.limit'))
        else:
            if np.array_equal(normal,None):
                for value,x in enumerate(colours):
                    if not(np.array_equal(synthetic_spectra,None)):
                        probability[value]= -np.sum((synthetic_spectra[value]-rgetattr(self,x+'.spec'))**2/(2*(rgetattr(self,x+'.uncs'))**2)*rgetattr(self,x+'.limit')*rgetattr(self,x+'.dip_array'))
                    else:
                        probability[value]= -np.sum((rgetattr(self,x+'.synth')-rgetattr(self,x+'.spec'))**2/(2*(rgetattr(self,x+'.uncs'))**2)*rgetattr(self,x+'.limit')*rgetattr(self,x+'.dip_array'))
            else:
                for value,x in enumerate(colours):
                    if not(np.array_equal(synthetic_spectra,None)):
                        probability[value]= -np.sum((synthetic_spectra[value]-normal[value])**2/(2*(rgetattr(self,x+'.uncs'))**2)*rgetattr(self,x+'.limit')*rgetattr(self,x+'.dip_array'))
                    else:
                        probability[value]= -np.sum((rgetattr(self,x+'.synth')-normal[value])**2/(2*(rgetattr(self,x+'.uncs'))**2)*rgetattr(self,x+'.limit')*rgetattr(self,x+'.dip_array'))
        for value,x in enumerate(colours):
            if rgetattr(self,x+'.normal_value') is None:
                rsetattr(self, x+'.normal_value',probability[value]+1)
        if self_normal:
            normal_shift=np.ones(len(colours))
            for value,x in enumerate(colours):
                normal_shift[value]=rgetattr(self, x+'.normal_value')

            probability=probability-normal_shift
            if np.any(probability>0):
                for prob,x in zip(probability,colours):
                    if prob>0:
                        print(prob,x)
                        print('old',rgetattr(self, x+'.normal_value'))
                        rsetattr(self, x+'.normal_value',rgetattr(self, x+'.normal_value')+prob+1)
                        print('new',rgetattr(self, x+'.normal_value'))
                        print(solar_shift)
                return False      
        return sum(probability)
    def observed_spectra_giver(self,colours=bands):
        returning_spectra=np.array(np.ones(len(colours)),dtype=object)
        wavelengths=np.array(np.ones(len(colours)),dtype=object)
        uncs=np.array(np.ones(len(colours)),dtype=object)
        for value,x in enumerate(colours):
            returning_spectra[value]=rgetattr(self,x+'.spec')
            wavelengths[value]=rgetattr(self,x+'.wave')
            uncs[value]=rgetattr(self,x+'.uncs')
        return returning_spectra,wavelengths,uncs
    def plot(self,colours=bands,dopler_shift=True):
        for x in colours:
            plt.figure(self.name+x) 
            # x_line=np.linspace(0,len(runGetattr(self,x+'.synth'))-1,num=len(runGetattr(self,x+'.synth')))
            if rgetattr(self,x+'.synth').any():
                    # labels='synthetic  chi squared= '+str(self.log_fit([x]))
                    if dopler_shift:
                        synth=dopler(rgetattr(self,x+'.wave'),rgetattr(self,x+'.synth'),-rgetattr(self,x+'.vrad'))
                    else:
                        synth=rgetattr(self,x+'.synth')
                    plt.plot(runGetattr(self,x+'.wave'), synth, label='Synthetic',c='r')
            if dopler_shift:
                observed=dopler(rgetattr(self,x+'.wave'),rgetattr(self,x+'.spec'),-rgetattr(self,x+'.vrad'))
                uncs=dopler(rgetattr(self,x+'.wave'),rgetattr(self,x+'.uncs'),-rgetattr(self,x+'.vrad'))
            else:
                observed=rgetattr(self,x+'.spec')
                uncs=rgetattr(self,x+'.uncs')
            plt.fill_between(runGetattr(self,x+'.wave'), observed- uncs,observed+ uncs,alpha=0.5)
            plt.plot(runGetattr(self,x+'.wave'), observed, label='observing',c='b')
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
def rv(rv,full_synth):
    bands=['Blue','Green','Red','IR']
    returning_spectra_rv=np.array(np.ones(len(bands)),dtype=object)
    error_bands=np.ones(4)
    for value,(x,colour) in enumerate(zip(full_synth,bands)):
        dopler_shifted_spectra=dopler(rgetattr(spectras,colour+'.wave'),x,rv)
        returning_spectra_rv[value]=synth_resolution_degradation(rgetattr(spectras,colour+'.wave'),dopler_shifted_spectra,rgetattr(spectras,colour+'.hermes')[7].data,rgetattr(spectras,colour+'.hermes')[7].header['b'])[:,1]
        error_bands[value]=spectras.log_fit(colours=[colour],synthetic_spectra=[returning_spectra_rv[value]])
    # error=spectras.log_fit(synthetic_spectra=returning_spectra_rv)
    return error_bands


def my_model(theta, x):
    m, c = theta
    return m * x + c


def my_loglike(theta, x, data, sigma):
    model = my_model(theta, x)
    return -(0.5 / sigma**2) * np.sum((data - model) ** 2)
class LogLike(at.Op):

    """
    Specify what type of object will be passed and returned to the Op when it is
    called. In our case we will be passing it a vector of values (the parameters
    that define our model) and returning a single "scalar" value (the
    log-likelihood)
    """

    itypes = [at.dvector]  # expects a vector of parameter values when called
    otypes = [at.dscalar]  # outputs a single scalar value (the log likelihood)

    def __init__(self, loglike):
        """
        Initialise the Op with various things that our log-likelihood function
        requires. Below are the things that are needed in this particular
        example.

        Parameters
        ----------
        loglike:
            The log-likelihood (or whatever) function we've defined
        data:
            The "observed" data that our log-likelihood function takes in
        x:
            The dependent variable (aka 'x') that our model requires
        sigma:
            The noise standard deviation that our function requires.
        """

        # add inputs as class attributes
        self.likelihood = loglike
    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        (theta,) = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(theta)

        outputs[0][0] = np.array(logl)  # output the log-likelihood
class LogLikeWithGrad(at.Op):

    itypes = [at.dvector]  # expects a vector of parameter values when called
    otypes = [at.dscalar]  # outputs a single scalar value (the log likelihood)

    def __init__(self, loglike,param):
        """
        Initialise with various things that the function requires. Below
        are the things that are needed in this particular example.

        Parameters
        ----------
        loglike:
            The log-likelihood (or whatever) function we've defined
        data:
            The "observed" data that our log-likelihood function takes in
        x:
            The dependent variable (aka 'x') that our model requires
        sigma:
            The noise standard deviation that out function requires.
        """

        # add inputs as class attributes
        self.param=param
        self.likelihood = loglike
        # initialise the gradient Op (below)
        self.logpgrad = LogLikeGrad(param)

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        (theta,) = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(theta,self.param)

        outputs[0][0] = np.array(logl) 

    def grad(self, inputs, g):
        # the method that calculates the gradients - it actually returns the
        # vector-Jacobian product - g[0] is a vector of parameter values
        (theta,) = inputs  # our parameters
        return [g[0] * self.logpgrad(theta)]
class LogLikeGrad(at.Op):

    """
    This Op will be called with a vector of values and also return a vector of
    values - the gradients in each dimension.
    """

    itypes = [at.dvector]
    otypes = [at.dvector]

    def __init__(self,param):
        """
        Initialise with various things that the function requires. Below
        are the things that are needed in this particular example.

        Parameters
        ----------
        data:
            The "observed" data that our log-likelihood function takes in
        x:
            The dependent variable (aka 'x') that our model requires
        sigma:
            The noise standard deviation that out function requires.
        """
        self.param=param
    def perform(self, node, inputs, outputs):
        (theta,) = inputs

        # calculate gradients
        shift_temp=shift_maker(theta,self.param)
        grads = spectras.gradient(shift_temp,labels=self.param)

        outputs[0][0] = grads  

def rv_mean_sig():
    vrad=['vrad_Blue','vrad_Green','vrad_Red','vrad_IR']
    vrad_r=['rv_Blue_r','rv_Green_r','rv_Red_r','rv_IR_r']
    mean_all=[]
    sig_all=[]
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
        mean_all.append(mean)
        sig_all.append(sig)
    return mean_all,sig_all
labels=['teff','logg','monh','Fe','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']   
global large_data         
large_data=fits.getdata('photometric_reduction_dr6_cross.fits',1)
large_data=Table(large_data)

#Needs to be read because some spectra doesnt have 
global all_reduced_data
all_reduced_data=fits.getdata('dr6.0.fits',1)
all_reduced_data=Table(all_reduced_data)
global x_min,x_max
tmp = np.load("NN_normalized_spectra_all_Blue.npz")
x_min=tmp['x_min']
x_max=tmp['x_max']
w_array_0 = tmp["w_array_0"]
w_array_1 = tmp["w_array_1"]
w_array_2 = tmp["w_array_2"]
b_array_0 = tmp["b_array_0"]
b_array_1 = tmp["b_array_1"]
b_array_2 = tmp["b_array_2"]
x_min = tmp["x_min"]
x_max = tmp["x_max"]
NN_coeffs= (w_array_0, w_array_1, w_array_2, b_array_0, b_array_1, b_array_2, x_min, x_max)


# #EMCEE,
# prior=False
np.random.seed(589403)
final_results=[]
expected_results=[]
large_data_GALAH_official=fits.open('gaia_galah_cross_values.fits')
large_data_GALAH_official_unchanged=large_data_GALAH_official[1].data
large_data_GALAH_official=Table(large_data_GALAH_official[1].data)
nwalkers=10



# name=210115002201239
# old_abundances=[y for y in all_reduced_data if y['sobject_id']==name]
# old_abundances=old_abundances[0]

# log_fit_all=[]
# spectras_fit=[]
# names=[]
# uncs_all=[]
# observed_spectra=[]
# to_be_reduced=all_reduced_data[::1000]
# for value,x in enumerate(to_be_reduced[1300:]):
#     name=x['sobject_id']
#     print('now_doing',name)
#     spectras=spectrum_all(name,cluster=False)
#     spectras.synthesize()
#     bands=spectras.bands
#     bad_reduction=np.any([rgetattr(spectras,x+'.bad_reduction') for x in bands ])
#     hermes_flag=np.any([rgetattr(spectras, x+'.bad_reduction') is None for x in spectras.bands ])

#     if bad_reduction or hermes_flag:
#         fit=None
#     else:
#         spectras.normalize()
#         fit=spectras.log_fit()
#     log_fit_all.append(fit)
#     spectra=spectras.synthesize(give_back=(True))
#     spectra_moved=[]
#     spectra_moved=[]
#     uncs_moved=[]
#     original_moved=[]
#     for colour,unshifted_spectra in zip(bands,spectra):
#         uncs=rgetattr(spectras,colour+'.uncs')
#         wave=rgetattr(spectras,colour+'.wave')
#         vrad=rgetattr(spectras,colour+'.vrad')
#         spec=rgetattr(spectras,colour+'.spec')
#         spectra_moved.append(dopler(wave,unshifted_spectra,-vrad))
#         uncs_moved.append(dopler(wave,uncs,-vrad))
#         original_moved.append(dopler(wave,spec,-vrad))

#     names.append(spectras.name)
#     spectras_fit.append(spectra_moved)
#     observed_spectra.append(original_moved)
#     uncs_all.append(uncs_moved)
                                                                  
#     if value%20==0:
#         np.save('testing_fit_5',np.array((names,spectras_fit,log_fit_all,spectras_fit,observed_spectra,uncs_all)))
# def making_spectra_and_fit(name):
#     print('doing',name)
#     spectras=spectrum_all(name,cluster=False)

#     reduction_flags=np.any([rgetattr(spectras, x+'.bad_reduction') for x in spectras.bands ])
#     hermes_flag=np.any([rgetattr(spectras, x+'.bad_reduction') is None for x in spectras.bands ])
#     if reduction_flags or hermes_flag:
#         print('bad',name)
#         fit=None
#         spectra_fit=np.ones(4)
#         spectra_moved=np.ones(4)
#         original_moved=np.ones(4)
#         uncs_moved=np.ones(4)
#     else:    
#         spectras.synthesize()

#         spectras.normalize()
#         fit=spectras.log_fit()
#         spectra_fit=spectras.synthesize(give_back=(True))
#         spectra_moved=[]
#         spectra_moved=[]
#         uncs_moved=[]
#         original_moved=[]
#         bands=spectras.bands
#         for colour,unshifted_spectra in zip(bands,spectra_fit):
#             uncs=rgetattr(spectras,colour+'.uncs')
#             wave=rgetattr(spectras,colour+'.wave')
#             vrad=rgetattr(spectras,colour+'.vrad')
#             spec=rgetattr(spectras,colour+'.spec')
#             spectra_moved.append(dopler(wave,unshifted_spectra,-vrad))
#             uncs_moved.append(dopler(wave,uncs,-vrad))
#             original_moved.append(dopler(wave,spec,-vrad))

        
        
#     return fit,spectra_moved,original_moved,uncs_moved,name
# ncpu=12
# batch=10
# large_results=[]
# for x in range(len(to_be_reduced)//(ncpu*batch)):
#     pool=Pool(processes=ncpu)
#     results=pool.map(making_spectra_and_fit,to_be_reduced[x*ncpu*batch:(x+1)*ncpu*batch]['sobject_id'])
#     pool.close()
#     if len(large_results)==0:
#         large_results=results
#     else:
#         large_results=np.vstack((large_results,results))
#     np.save('testing_fit_3',large_results)
#     print('Saving at '+str((x+1)*ncpu)+' iteration. It is '+ str((x+1)*ncpu*batch/len(to_be_reduced)))
# for x in all_reduced_data:
    
# spectras.synthesize()
# spectras.gradient()
# spectras.normalize()
# spectras.log_fit()
# spectras.gradient({},labels=['teff'])
# spectras.gradient()
# spectras.gradient({'teff':500},labels=['teff'])
# spectras.gradient({'teff':5000},labels=['teff','vrad_Blue'])
# a=spectras.synthesize(give_back=True,grad=True)
# spectras.gradient()
# spectras.normalize()
# spectras=spectrum_all(name, 10)
# spectras.synthesize({})
# spectras.normalize()
# define a aesara Op for our likelihood function

        
# set up our data



# labels=['teff','logg','monh','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']
# plt.figure()
# how_many_starting=30
# starting_radial=[]
# labels=['vrad_Blue','vrad_Green','vrad_Red','vrad_IR']
# for lab in labels:
#     log_temp=np.ones(how_many_starting)
#     grad_temp=np.ones(how_many_starting)
#     for value,rad in enumerate(np.linspace(20, 50,num=how_many_starting)):
        
#         log_temp[value]=spectras.log_fit(shift={lab:rad},self_normal=(True))
#         grad_temp[value]=spectras.gradient(shift={lab:rad},labels=[lab])
#     plt.plot(np.linspace(20, 50,num=how_many_starting),log_temp,label=lab)
#     plt.plot(np.linspace(20, 50,num=how_many_starting),grad_temp,label=lab+'grad')

#     starting_radial.append(np.linspace(20, 50,num=how_many_starting)[np.where(log_temp==max(log_temp))][0])
# plt.legend(loc='best')
# plt.figure()
# log_temp=np.ones(how_many_starting)
# grad_temp=np.ones(how_many_starting)
# for value,rad in enumerate( np.linspace(5000, 7000,num=how_many_starting)):
#     log_temp[value]=spectras.log_fit(shift={'teff':rad})
#     grad_temp[value]=spectras.gradient(shift={'teff':rad},labels=['teff'])
   
# plt.plot(np.linspace(5000, 7000,num=how_many_starting),log_temp/max(abs(log_temp)),label='log')
# plt.plot(np.linspace(5000, 7000,num=how_many_starting),grad_temp/max(abs(grad_temp)),label='grad')
# plt.legend(loc='best')
# # shift_temp={y:x for y,x in zip(labels,starting_radial)}
# # spectras.synthesize()
# # # spectras.log_fit(shift=shift_temp)
# # # spectras.gradient(shift={'teff': 6244.1780897915805, 'logg': 4.246294816412452, 'monh': 0.053370714, 'Fe': 0.053370714, 'alpha': 0.0064120926, 'vrad_Blue': 73.23569035759559, 'vrad_Green': 83.3159958109764, 'vrad_Red': 80.64146650027877, 'vrad_IR': 81.24084720075764, 'vsini': 6.3641806, 'vmac': 6.0, 'vmic': 1.2957745},labels=['vrad_Red'])

# # # # create our Op
# # # logl = LogLike(log_posterior)
labels=['teff','logg','monh','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']
# labels=['teff','logg','monh','alpha','vsini','vmac','vmic']
logl=LogLikeWithGrad(log_posterior,labels)
# keys=['teff','logg','monh','Fe','alpha','vsini','vmac','vmic']
# # use PyMC to sampler from log-likelihood
# idata_mh=[]


# teff=np.normal(old_abundances['teff'],old_abundances['e_teff'])

spectras=spectrum_all(150107004201330,cluster=False)
spectras.synthesize()
spectras.normalize()
# old_abundances=spectras.old_abundances
# with pm.Model():
    
    
#     # # solar_values=pm.Uniform('teff',lower=x_min[0],upper=x_max[0])
#     # # teff=pm.Uniform('teff',lower=old_abundances['teff']-200,upper=old_abundances['teff']+200,testval=old_abundances['teff'])
#     #uniform priors
#     # Normal priors
#     # teff=pm.Normal('teff',mu=old_abundances['teff'],sigma=old_abundances['teff'],testval=old_abundances['teff'])
#     # logg=pm.Normal('logg',mu=old_abundances['logg'],sigma=old_abundances['e_logg'],testval=old_abundances['logg'])   
#     # #uniform
#     teff=pm.Uniform('teff',lower=old_abundances['teff_r']-200,upper=old_abundances['teff_r']+200,testval=old_abundances['teff_r'])
#     logg=pm.Uniform('logg',lower=old_abundances['logg_r']-1,upper=old_abundances['logg_r']+1,testval=old_abundances['logg_r'])
#     monh=pm.Uniform('monh',lower=old_abundances['fe_h_r']-1,upper=old_abundances['fe_h_r']+1,testval=old_abundances['fe_h_r'])
#     alpha=pm.Uniform('alpha',lower=old_abundances['alpha_fe_r']-1,upper=old_abundances['alpha_fe_r']+1,testval=old_abundances['alpha_fe_r'])
#     # rv_mean,rv_sig=rv_mean_sig() 
#     vrad_Blue=pm.Uniform('vrad_Blue',lower=old_abundances['rv_com']-5,upper=old_abundances['rv_com']+5,testval=old_abundances['rv_com'])
#     vrad_Red=pm.Uniform('vrad_Red',lower=old_abundances['rv_com']-5,upper=old_abundances['rv_com']+5,testval=old_abundances['rv_com'])
#     vrad_Green=pm.Uniform('vrad_Green',lower=old_abundances['rv_com']-5,upper=old_abundances['rv_com']+5,testval=old_abundances['rv_com'])
#     vrad_IR=pm.Uniform('vrad_IR',lower=old_abundances['rv_com']-5,upper=old_abundances['rv_com']+5,testval=old_abundances['rv_com'])
#     vsini=pm.Uniform('vsini',lower=x_min[5],upper=old_abundances['vbroad_r']+20,testval=old_abundances['vbroad_r'])
#     vmac=pm.Uniform('vmac',lower=x_min[6],upper=x_max[6],testval=old_abundances['vmic_r'])
#     vmic=pm.Uniform('vmic',lower=x_min[7],upper=x_max[7],testval=old_abundances['vmic_r'])
    
    
#     #stellar
#     # teff=pm.Uniform('teff',lower=old_abundances['teff_r']-200,upper=old_abundances['teff_r']+200,testval=old_abundances['teff_r'])
#     # monh=pm.Uniform('monh',lower=old_abundances['fe_h_r']-1,upper=old_abundances['fe_h_r']+1,testval=old_abundances['fe_h_r'])
#     # fe=pm.Uniform('fe',lower=old_abundances['fe_h_r']-1,upper=old_abundances['fe_h_r']+1,testval=old_abundances['fe_h_r'])
#     # logg=pm.Uniform('logg',lower=old_abundances['logg_r']-1,upper=old_abundances['logg_r']+1,testval=old_abundances['logg_r'])
#     # alpha=pm.Uniform('alpha',lower=old_abundances['alpha_fe_r']-1,upper=old_abundances['alpha_fe_r']+1,testval=old_abundances['alpha_fe_r'])   
#     # vrad_Blue=pm.Uniform('vrad_Blue',lower=old_abundances['rv'][0]-5,upper=old_abundances['rv'][0]+5,testval=old_abundances['rv'][0])
#     # vrad_Red=pm.Uniform('vrad_Red',lower=old_abundances['rv'][1]-5,upper=old_abundances['rv'][1]+5,testval=old_abundances['rv'][1])
#     # vrad_Green=pm.Uniform('vrad_Green',lower=old_abundances['rv'][2]-5,upper=old_abundances['rv'][2]+5,testval=old_abundances['rv'][2])
#     # vrad_IR=pm.Uniform('vrad_IR',lower=old_abundances['rv'][3]-5,upper=old_abundances['rv'][3]+5,testval=old_abundances['rv'][3])
#     # vsini=pm.Uniform('vsini',lower=x_min[5],upper=x_max[5],testval=old_abundances['vbroad_r'])
#     # vmac=pm.Uniform('vmac',lower=x_min[6],upper=x_max[6],testval=old_abundances['vmic_r'])
#     # vmic=pm.Uniform('vmic',lower=x_min[7],upper=x_max[7],testval=old_abundances['vmic_r'])    
#     # solar_values=at.as_tensor_variable([teff,logg,monh,fe,alpha,vrad_Blue,vrad_Green,vrad_Red,vrad_IR,vsini,vmac,vmic])
    
#     #weak prior
    
#     # teff=pm.Normal('teff',mu=old_abundances['teff_r'],sigma=50,testval=old_abundances['teff_r'])
#     # monh=pm.Normal('monh',mu=old_abundances['fe_h_r'],sigma=1,testval=old_abundances['fe_h_r'])
#     # # fe=pm.Normal('fe',mu=old_abundances['fe_h_r'],sigma=1,testval=old_abundances['fe_h_r'])
#     # logg=pm.Normal('logg',mu=old_abundances['logg_r'],sigma=1,testval=old_abundances['logg_r'])
#     # alpha=pm.Normal('alpha',mu=old_abundances['alpha_fe_r'],sigma=1,testval=old_abundances['alpha_fe_r'])   
#     # # vrad_Blue=pm.Normal('vrad_Blue',mu=old_abundances['rv'][0],sigma=old_abundances['e_rv'][0]*4,testval=old_abundances['rv'][0])
#     # # vrad_Red=pm.Normal('vrad_Red',mu=old_abundances['rv'][1],sigma=old_abundances['e_rv'][1]*4,testval=old_abundances['rv'][1])
#     # # vrad_Green=pm.Normal('vrad_Green',mu=old_abundances['rv'][2],sigma=old_abundances['e_rv'][2]*4,testval=old_abundances['rv'][2])
#     # # vrad_IR=pm.Normal('vrad_IR',mu=old_abundances['rv'][3],sigma=old_abundances['e_rv'][3]*4,testval=old_abundances['rv'][3])
#     # vsini=pm.Normal('vsini',mu=old_abundances['vbroad_r'],sigma=10,testval=old_abundances['vbroad_r'])
#     # vmac=pm.Normal('vmac',mu=old_abundances['vmic_r'],sigma=10,testval=old_abundances['vmic_r'])
#     # vmic=pm.Normal('vmic',mu=old_abundances['vmic_r'],sigma=10,testval=old_abundances['vmic_r'])    
#     # teff=pm.Normal('teff',mu=old_abundances['teff_r'],sigma=50)
#     # monh=pm.Normal('monh',mu=old_abundances['fe_h_r'],sigma=1)
#     # # fe=pm.Normal('fe',mu=old_abundances['fe_h_r'],sigma=1,testval=old_abundances['fe_h_r'])
#     # logg=pm.Normal('logg',mu=old_abundances['logg_r'],sigma=1)
#     # alpha=pm.Normal('alpha',mu=old_abundances['alpha_fe_r'],sigma=1)
#     # # vrad_Blue=pm.Normal('vrad_Blue',mu=old_abundances['rv'][0],sigma=old_abundances['e_rv'][0]*4,testval=old_abundances['rv'][0])
#     # # vrad_Red=pm.Normal('vrad_Red',mu=old_abundances['rv'][1],sigma=old_abundances['e_rv'][1]*4,testval=old_abundances['rv'][1])
#     # # vrad_Green=pm.Normal('vrad_Green',mu=old_abundances['rv'][2],sigma=old_abundances['e_rv'][2]*4,testval=old_abundances['rv'][2])
#     # # vrad_IR=pm.Normal('vrad_IR',mu=old_abundances['rv'][3],sigma=old_abundances['e_rv'][3]*4,testval=old_abundances['rv'][3])
#     # vsini=pm.Normal('vsini',mu=old_abundances['vbroad_r'],sigma=10)
#     # vmac=pm.Normal('vmac',mu=old_abundances['vmic_r'],sigma=10)
#     # vmic=pm.Normal('vmic',mu=old_abundances['vmic_r'],sigma=10)   

#     solar_values=at.as_tensor_variable([teff,logg,monh,alpha,vrad_Blue,vrad_Green,vrad_Red,vrad_IR,vsini,vmac,vmic])
#     # solar_values=at.as_tensor_variable([teff,logg,monh,alpha,vsini,vmac,vmic])
    
#     # use a Potential to "call" the Op and include it in the logp computation
#     pm.Potential("likelihood", logl(solar_values))
#     # func = pm.Model.logp_dlogp_function(profile=True)
#     # func.set_extra_values({})
#     # x = np.random.randn(func.size)
#     # # %timeit func(x)
    
#     # func.profile.print_summary()  
#     # Use custom number of draws to replace the HMC based defaults
#     # func = model.logp_dlogp_function(profile=True)
#     # a=pm.find_MAP()
#     idata_mh = pm.sample(200, tune=200,return_inferencedata=True,cores=4)

# # # # # plot the traces
# # # # idata_mh.to_netcdf
# # # # idata_mh.posterior
# # # # az.from_netcdf
# az.plot_trace(idata_mh)
# # fig, axes = plt.subplots(3,2, figsize=(12,6))
# # for i, kind in enumerate(("trace", "rank_bars", "rank_vlines")):
# #     az.plot_trace(idata_mh, var_names="teff", kind=kind,compact=True);
# # fig.tight_layout()