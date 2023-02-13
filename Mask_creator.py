
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 13:39:38 2021

@author: kevin
"""
from scipy.stats import kde
from scipy import integrate

from numba import jit
import time
from functools import  partial
from astropy.io.votable import parse
import emcee
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
# import Galah_tool_py3 as gtools
import os
# import random
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
def shift_maker(solar_value,given_labels=['teff','logg','monh','fe_h','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic'],fe_hold=False,all_labels=['teff','logg','monh','fe_h','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']):
        """
        Create a dictionary of solar parameters from an array of values
        
    
        Parameters
        ----------
        solar_values : array of 9 length in the order of  teff,logg,monh,fe,alpha,vrad,vsini,vmac,vmic
    
        Returns
        -------
        shift : dictionary of solar parameters
    
        """
        shift={}
        skip=0
        for value,x in enumerate(all_labels):
            if x in given_labels:
                if x=='fe_h' and fe_hold:
                    shift[x]=shift['monh']
                else:
                    shift[x]=solar_value[value-skip]
            else:
                skip+=1
                bands=spectras.bands
                colour=[y for y in bands if y in x]
                if x=='fe_h' and fe_hold:
                    shift[x]=shift['monh']
                elif not 'vrad' in x:

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

def log_posterior(solar_values,parameters,prior=False,first_try=False,full_parameters=['teff']):
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
        if len(solar_values)!= len(parameters):
            print('youre missing parameters')
        
        if not starting_test(solar_values,spectras.old_abundances,parameters,cluster=prior):
              return -np.inf
        # print('passed')
        shift=shift_maker(solar_values,given_labels=parameters,all_labels=full_parameters)
        # shift['teff']=(shift['teff']+6000)*10
        # try:
        synthetic_spectras=spectras.synthesize(shift,give_back=True)
            # print('shift' , shift)
            # print("Spectra",synthetic_spectras)
        # except:
        #     print('couldnt finish synthesizing')
        #     return -np.inf
        
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
    
    
        normalized_spectra,normalized_uncs=spectras.normalize(data=synthetic_spectras)
        if first_try==True:
            normalized_limit_array=spectras.create_masks(clean=True,shift=shift)
        else:
            normalized_limit_array=first_try
        probability=spectras.log_fit(synthetic_spectra=synthetic_spectras,solar_shift=shift,normal=normalized_spectra,uncertainty=normalized_uncs,limit_array=normalized_limit_array,combine_masks=False)
        # print(prior_2d(shift))
        if prior:
            probability+=prior_2d(shift)
        if probability>0:
            print('oh no prob ',probability)
            print(repr(solar_values))
            # return False
        # print('probability', probability)
        return probability
@jit(nopython=True,cache=True)
def normal(x,mu,sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-1/2*((x-mu)/sigma)**2)
def prior_2d(shift):
    """
    Creates a prior from a 2d probabilty distribution gotten from photometric temperatures and logg

    Parameters
    ----------
    shift :DICT
        a dictionary that has temperature and logg which we would like to calculate the probability of it using the prior.

    Returns
    -------
    FLOAT
        log of the probability.

    """
    #Spread number defines if the sampling of the prior has been fine enough. if it has we use a sumation method(spread number is bellow 1) 
    #and if it hasnt we use a kde method (spread number above 1)
    if spectras.spread_number>1:
        #calculates the kde if hasnt been done
        if spectras.kde==None:
            spectras.kde=kde.gaussian_kde([old_abundances['teff_raw'],old_abundances['logg_raw']])
        photometric_probability_density=spectras.kde
        e_teff=spectras.e_teff_photometric
        e_logg=spectras.e_logg_photometric
        # f= lambda teff_var,logg_var,teff_0,logg_0:photometric_probability_density([teff_var,logg_var])*normal(teff_0,teff_var,e_teff)*normal(logg_0,logg_var,e_logg)
        total_int=0
        teff_prior=shift['teff']
        #applies the shift from spectroscopic temperature to photometric
        polynomial_coeff=old_abundances['coeff']
        spectroscopic_shift=np.poly1d(polynomial_coeff)
        teff_prior+=spectroscopic_shift(teff_prior)
        
        logg_prior=shift['logg']
        # teff_line=np.linspace(teff_prior-e_teff*10, teff_prior+e_teff*10,2)
        # logg_line=np.linspace(logg_prior-e_logg*10, logg_prior+e_logg*10,2)
        
        total_int=photometric_probability_density([teff_prior,logg_prior])[0]
        # for x in range(len(teff_line)-1):
        #     for y in range(len(logg_line)-1):
        #         total_int+=integrate.dblquad(f, logg_line[y], logg_line[y+1], lambda teff_var: teff_line[x], lambda teff_var:teff_line[x+1],args=(teff_prior,logg_prior))[0]
    else:
        # teff_raw=fast_abundances['teff_raw']
        # logg_raw=old_abundances['logg_raw']
        # e_teff_raw=old_abundances['e_teff_raw']
        # e_logg_raw=old_abundances['e_logg_raw']
        
        teff_raw=fast_abundances[0]
        logg_raw=fast_abundances[1]
        e_teff_raw=fast_abundances[2]
        e_logg_raw=fast_abundances[3]

        prior_parameters=np.column_stack((teff_raw,e_teff_raw,logg_raw,e_logg_raw))
        teff_prior=shift['teff']
        logg_prior=shift['logg']
        #applied the shift from spectroscopic temperature to photometric

        # polynomial_coeff=fast_abundances[5]
        spectroscopic_shift=np.poly1d(fast_coefficients)

        teff_prior+=spectroscopic_shift(teff_prior)

        total_int=0
        for temp_param in prior_parameters:           
            total_int+=normal(teff_prior,temp_param[0],temp_param[1])*normal(logg_prior,temp_param[2],temp_param[3])
        total_int/=len(old_abundances['teff_raw'])
    #if the probability is zero it means that the asked parameters are too far away from the photometric temperatures. 
    #Thus we have to use a simple normal approximation
    if total_int==0:
        e_teff=spectras.e_teff_photometric
        e_logg=spectras.e_logg_photometric
        return np.log(1/(e_teff*e_logg*2*np.pi))-(1/2*((teff_prior-old_abundances['teff_photometric'])/e_teff)**2)-\
            (1/2*((logg_prior-old_abundances['logg_photometric'])/e_logg)**2)
    
    return np.log(total_int)
            
            
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
def starting_test(solar_values,old_abundances,parameters=['teff','logg','monh','fe_h','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic'],cluster=False):
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
        shift={}
        for y,x in zip(parameters,solar_values):
            shift[y]=x
        labels_with_limits=['teff','logg','fe_h','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
        labels_with_limits=[x for x in labels_with_limits if x in parameters]
        for x in labels_with_limits:
            if shift[x]<x_min[x]-abs(x_min[x])*0.2 or shift[x]>x_max[x]*1.2:
                print('outside of the models limits ',x,' value ', shift[x])
                return False
        vrad=['vrad_Blue','vrad_Green','vrad_Red','vrad_IR']
        vrad=[x for x in vrad if x in parameters]
        vrad_r=['rv'+x[4:]+'_r' for x in vrad]
        radial_velocities=[shift[x] for x in vrad]
        if len (radial_velocities):
            if cluster:
                for rad,rad_r in zip(vrad,vrad_r):
                    if not np.isnan(old_abundances['red_rv_ccd'][vrad.index(rad)]):
                        mean=float(old_abundances['red_rv_ccd'][vrad.index(rad)])
                    elif not np.isnan(old_abundances['red_rv_com']):
                        mean=float(old_abundances['red_rv_com'])
                    else:
                        break
                    if not np.isnan(old_abundances['red_e_rv_ccd'][vrad.index(rad)]):
                        sig=float(old_abundances['red_e_rv_ccd'][vrad.index(rad)])*3
        
                    elif not np.isnan(old_abundances['red_e_rv_com']):
                        sig=float(old_abundances['red_e_rv_com'])*3
                    else:
                        sig=5
                    if abs(shift[rad]-mean)>sig+2:
                        print(rad,' is wrong by ',str(float(abs(shift[rad]-old_abundances['red_rv_com']))))
                        return False
                if np.count_nonzero(~np.isnan(radial_velocities))>1.0:
                    if max(radial_velocities)-min(radial_velocities)>(np.nanmax(radial_velocities)-np.nanmin(radial_velocities))*4:
                        print('radial velocities too different')
                        return False
                else:
                    if max(radial_velocities)-min(radial_velocities)>5.0:
                        print('radial velocities too different')
                        return False
            else:
                for rad,rad_r in zip(vrad,vrad_r):
                    if not np.isnan(old_abundances['rv'][vrad.index(rad)]):
                        mean=float(old_abundances['rv'][vrad.index(rad)])
                    elif not np.isnan(old_abundances['rv_com']):
                        mean=float(old_abundances['rv_com'])
                    else:
                        break
                    if not np.isnan(old_abundances['e_rv'][vrad.index(rad)]):
                        sig=float(old_abundances['e_rv'][vrad.index(rad)])*3
        
                    elif not np.isnan(old_abundances['e_rv_com']):
                        sig=float(old_abundances['e_rv_com'])*3
                    else:
                        sig=5
                    if abs(shift[rad]-mean)>sig+2:
                        print(rad,' is wrong by ',str(float(abs(shift[rad]-old_abundances['rv_com']))))
                        return False
                if np.count_nonzero(~np.isnan(radial_velocities))>1.0:
                    if max(radial_velocities)-min(radial_velocities)>np.nanmax(radial_velocities)-np.nanmin(radial_velocities)+np.nanmax(old_abundances['e_rv']):
                        print('radial velocities too different')
                        return False
                else:
                    if max(radial_velocities)-min(radial_velocities)>5.0:
                        print('radial velocities too different')
                        return False
        return True
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

# @jit(nopython=True,parallel=False,cache=True)
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
            limits={'Blue':[4705,4908],'Green':[5643,5879],'Red':[6470,6743],'IR':[7577.0,7894.0]}

            if x=='IR':
                starting_fraction=1000/4096
                length_fraction=3096/(4096*(1-starting_fraction))
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
            tmp = np.load("NN_normalized_spectra_all_elements_"+x+".npz")
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
            Path('dr6.1/'+name[0:6]+'/spectra/com/').mkdir(parents=True,exist_ok=True)
            if not exists('dr6.1/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'):
                #source='/media/storage/HERMES_REDUCED/dr6.1/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'
                source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.1/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'

                destination='dr6.1/'+name[0:6]+'/spectra/com/'
                subprocess.run(["rsync",'-av',source,destination])

            try:
                hermes=fits.open('dr6.1/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits')
                x0= float(hermes[1].header.get('CRVAL1'))
                x1=float( hermes[1].header.get('CRVAL1')+len(hermes[1].data)* hermes[1].header.get('CDELT1'))
                fstart= x0+(x1-x0)*starting_fraction
                
                new_start=int(starting_fraction*len(hermes[1].data))
                new_end=new_start+int(len(hermes[1].data)*(1-starting_fraction)*length_fraction)
                
                
                length=np.linspace(x0+(x1-x0)*starting_fraction,x1,num=int(4096*(1-starting_fraction)))
                length_synthetic=np.linspace(limits[x][0], limits[x][1],num=len(b_array_2))
                fend=length[-1]
                #Increases the length of synthetic spectra so it over interpolates  
                self.wave_synth=length_synthetic
                self.wave=length
                self.spec=hermes[0].data[new_start:new_end]
                self.spec_original=hermes[0].data[new_start:new_end]

                self.uncs=hermes[2].data[new_start:new_end]
                self.uncs*=self.spec[0]
                self.uncs_original=copy.deepcopy(self.uncs)
                hermes[0].data=hermes[0].data[new_start:new_end]
                hermes[1].data=hermes[1].data[new_start:new_end]
                hermes[2].data=hermes[2].data[new_start:new_end]
                hermes[7].data=hermes[7].data[new_start:new_end]
                # not_resolved=[np.int64(name)]
                # try:
                #     if len(hermes[7].data[new_start:new_end])==new_end-new_start and min(hermes[7].data)>0:
                #         hermes[7].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[7].data[new_start:new_end])
                # except TypeError:
                #     hermes[7].data=[0]
                    
                # hermes[1].header['CDELT1']/=interpolation
                self.wran=[fstart,fend]
                self.hermes=hermes
                if cluster:
                    self.teff=float(old_abundances['teff_spectroscopic'])
                    self.vmic=float(old_abundances['vmic'])
                    self.vsini=float(old_abundances['vsini'])
                    self.Fe=float(old_abundances['fe_h'])
                    self.fe_h=float(old_abundances['fe_h'])
                    if not np.isnan(old_abundances['red_rv_ccd'][count-1]):
                        self.vrad=old_abundances['red_rv_ccd'][count-1]
                    elif not np.isnan(float(old_abundances['red_rv_com'])):
                        self.vrad=float(old_abundances['red_rv_com'])
                    else:
                        self.vrad=0.0
                    self.vmac=6.0
                    self.logg=float(old_abundances['logg_photometric'])
                    self.monh=float(old_abundances['fe_h'])
                    elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
                    for individual_element in elements:
                        if not isinstance(old_abundances[individual_element.lower()+'_fe'],np.float32):
                            setattr(self,individual_element,0.0)
                        else:
                            setattr(self,individual_element,old_abundances[individual_element.lower()+'_fe'])

                    
                else:
                    self.teff=float(old_abundances['teff_r'])
                    self.vmic=float(old_abundances['vmic_r'])
                    self.vsini=float(old_abundances['vbroad_r'])
                    self.Fe=float(old_abundances['fe_h_r'])
                    self.fe_h=float(old_abundances['fe_h_r'])

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
                    elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
                    for individual_element in elements:
                        setattr(self,individual_element,0.0)

                if hermes[1].header.get('TEFF_R')=='None' or hermes[1].header.get('LOGG_R')=='None':
                    print('Warning the reduction didnt produce an TEFF spectra might not reduce properly')
                    print('enter to continue')
                    self.bad_reduction=True
                else:
                    self.bad_reduction=False
                line, wavelength = np.loadtxt('important_lines',usecols=(0,1),unpack=True,dtype=str, comments=';')
                np.append(wavelength,'4861.3230')
                np.append(wavelength,'6562.7970')
                np.append(line,r'H$_\beta$')
                np.append(line,r'H$_\alpha$')
                wavelength=wavelength.astype(float)
                important_lines_temp=np.vstack([[elem_temp,wave_temp] for elem_temp,wave_temp in zip(line,wavelength) if wave_temp>length_synthetic[0] and wave_temp<length_synthetic[-1]])
                self.important_lines=important_lines_temp
                rsetattr(self,'.important_lines',important_lines_temp)


            except FileNotFoundError:
                print('No hermes found')
                self.hermes=None
                self.bad_reduction=True
            


            self.normal_value=None

            self.l_new=None
            # Load spectrum masks
            masks = Table.read('solar_spectrum_mask.fits')
            masks_temp=vstack([x for x in masks if x['mask_begin']>self.wran[0]-200 and x['mask_end']<self.wran[1]+200])
            self.masks=masks_temp
            vital_lines = Table.read('vital_lines.fits')
            vital_lines_temp=vstack([x for x in vital_lines if x['line_begin']>self.wran[0]-200 and x['line_end']<self.wran[1]+200])
            self.vital_lines=vital_lines_temp
def sclip(p,fit,n,ye=[],sl=99999,su=99999,min=0,max=0,min_data=1,grow=0,global_mask=None,verbose=True):
    """
    p: array of coordinate vectors. Last line in the array must be values that are fitted. The rest are coordinates.
    fit: name of the fitting function. It must have arguments x,y,ye,and mask and return an array of values of the fitted function at coordinates x
    n: number of iterations
    ye: array of errors for each point
    sl: lower limit in sigma units
    su: upper limit in sigma units
    min: number or fraction of rejected points below the fitted curve
    max: number or fraction of rejected points above the fitted curve
    min_data: minimal number of points that can still be used to make a constrained fit
    global_mask: if initial mask is given it will be used throughout the whole fitting process, but the final fit will be evaluated also in the masked points
    grow: number of points to reject around the rejected point.
    verbose: print the results or not
    """

    nv,dim=np.shape(p)

    #if error vector is not given, assume errors are equal to 0:
    if ye==[]: ye=np.zeros(dim)
    #if a single number is given for y errors, assume it means the same error is for all points:
    if isinstance(ye, (int, float)): ye=np.ones(dim)*ye

    if global_mask==None: global_mask=np.ones(dim, dtype=bool)
    else: pass

    f_initial=fit(p,ye,global_mask)
    s_initial=np.std(p[-1]-f_initial)

    f=f_initial
    s=s_initial

    tmp_results=[]

    b_old=np.ones(dim, dtype=bool)

    for step in range(n):
        #check that only sigmas or only min/max are given:
        if (sl!=99999 or su!=99999) and (min!=0 or max!=0):
            raise RuntimeError('Sigmas and min/max are given. Only one can be used.')

        #if sigmas are given:
        if sl!=99999 or su!=99999:
            b=np.zeros(dim, dtype=bool)
            if sl>=99999 and su!=sl: sl=su#check if only one is given. In this case set the other to the same value
            if su>=99999 and sl!=su: su=sl

            good_values=np.where(((f-p[-1])<(sl*(s+ye))) & ((f-p[-1])>-(su*(s+ye))))#find points that pass the sigma test
            b[good_values]=True

        #if min/max are given
        if min!=0 or max!=0:
            b=np.ones(dim, dtype=bool)
            if min<1: min=dim*min#detect if min is in number of points or percentage
            if max<1: max=dim*max#detect if max is in number of points or percentage

            bad_values=np.concatenate(((p[-1]-f).argsort()[-int(max):], (p[-1]-f).argsort()[:int(min)]))
            b[bad_values]=False

        #check the grow parameter:
        if grow>=1 and nv==2:
            b_grown=np.ones(dim, dtype=bool)
            for ind,val in enumerate(b):
                if val==False:
                    ind_l=ind-int(grow)
                    ind_u=ind+int(grow)+1
                    if ind_l<0: ind_l=0
                    b_grown[ind_l:ind_u]=False

            b=b_grown

        tmp_results.append(f)

        #check that the minimal number of good points is not too low:
        if len(b[b])<min_data:
            step=step-1
            b=b_old
            break

        #check if the new b is the same as old one and break if yes:
        if np.array_equal(b,b_old):
            step=step-1
            break

        #fit again
        f=fit(p,ye,b&global_mask)
        s=np.std(p[-1][b]-f[b])
        b_old=b

    if verbose:
        print('')
        print('FITTING RESULTS:')
        print('Number of iterations requested:    ',n)
        print('Number of iterations performed:    ', step+1)
        print('Initial standard deviation:        ', s_initial)
        print('Final standard deviation:          ', s)
        print('Number of rejected points:         ',len(np.invert(b[np.invert(b)])))
        print('')

    return f,tmp_results,b
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

    def __init__(self,name,interpolation=10,cluster=True,bands=None):
        name=str(name)
        if bands==None:
            bands=['Blue','Green','Red','IR']
        self.rv_shift=1e-10
        self.bands=bands
        self.name=name
        self.interpolation=interpolation
        self.sister_stars=None
        if cluster:
            mask=photometric_data['sobject_id']==np.int64(name)
            old_abundances=photometric_data[mask]
            old_abundances=old_abundances[0]
            spread_temperature=np.sqrt(np.var(old_abundances['teff_raw']))
            sig_teff=np.mean(old_abundances['e_teff_raw'])
            spread_logg=np.sqrt(np.var(old_abundances['logg_raw']))
            sig_logg=np.mean(old_abundances['e_logg_raw'])
            self.spread_number=spread_temperature*spread_logg/(sig_logg*sig_teff*len(old_abundances['e_logg_raw']))
            self.e_teff_photometric=sig_teff
            self.e_logg_photometric=sig_logg
            self.kde=None
        else:
            mask=all_reduced_data['sobject_id']==np.int64(name)
            old_abundances=all_reduced_data[mask]
            old_abundances=old_abundances[0]
        self.old_abundances=old_abundances
        
        for count,x in enumerate(bands,1):
            setattr(self,x,individual_spectrum(name,interpolation,x,count,old_abundances,cluster))
        self.correct_resolution_map()
        self.equilize_spectras()
        self.limit_array()
    def equilize_spectras(self,colours=None):
        if colours==None:
            colours=self.bands
        for x in colours:
            if rgetattr(self, x+'.hermes')!=None:
                wavelength=rgetattr(self, x+'.wave')
                hermes=rgetattr(self, x+'.hermes')
                resmap=hermes[7].data
                observed_spectra=rgetattr(self, x+'.spec')
                observed_error=rgetattr(self, x+'.uncs')
                new_observed,new_error=equalize_resolution(wavelength,resmap,observed_spectra,observed_error,hermes)
                rsetattr(self, x+'.spec',new_observed)
                rsetattr(self, x+'.uncs',new_error)
                rsetattr(self, x+'.spec_equalized',new_observed)
                rsetattr(self, x+'.uncs_equalized',new_error)
                
    def dip_priority(self,colours=None,sigfig=0.5):
        if colours==None:
            colours=self.bands

        for x in colours:
            dip_array=[np.exp((y-1.0)**2/2*sigfig**2) for y in rgetattr(self,x+'.spec')]
            rsetattr(self, x+'.dip_array',dip_array)
    def limit_array(self,colours=None,limit=1.05,observed_spectra=None):
        if colours==None:
            colours=self.bands
        if not np.array_equal(observed_spectra,None):
            returning_limit=np.array(np.ones(len(colours)),dtype=object)
            for value,spec in enumerate(observed_spectra):
                limit_array_temp=[]
                for y in spec:
                    if y>limit:
                        limit_array_temp.append(0)
                    else:
                        limit_array_temp.append(1)
                returning_limit[value]=limit_array_temp
            return returning_limit
        for x in colours:
            if rgetattr(self, x+'.hermes')!=None:
                limit_array=[]
                for y in rgetattr(self,x+'.spec'):
                    if y>limit:
                        limit_array.append(0)
                    else:
                        limit_array.append(1)
                rsetattr(self,x+'.limit',limit_array)
            
    def solar_value_maker(self,shift,colours=None,keys=['teff','logg','monh','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']):
        if colours==None:
            colours=self.bands

        solar_values=[]
        for x in keys:
            if x in shift:
                solar_values.append(shift[x])
            else:
                solar_values.append(rgetattr(self,colours+'.'+x))
        return solar_values
    def hermes_checker(self,limit=20,colours=None):
        if colours==None:
            colours=self.bands
        for x in colours:
            spectra=rgetattr(self, x+'.spec')
            if min(spectra)<0:
                print (min(spectra),' ',x)
                return True
        return False
    def create_masks(self,colours=None,synthetic_spectra_insert=None,uncs_insert=None,normalized_observed_spectra_insert=None,shift=None,clean=False,limits=[5,0.4]):
        if colours==None:
            colours=self.bands

        masks_all=np.array(np.ones(len(colours)),dtype=object)
        #Clean is a mask without the difference between synthetik spectra and observed
        if clean:
            masks_all=np.array(np.ones(len(colours)),dtype=object)

            for value,x in enumerate(colours):
                
                masks=rgetattr(self,x+'.masks')
                vital_lines=rgetattr(self,x+'.vital_lines')
                original_wavelength_strenched=rgetattr(self,x+'.wave')
                v_rad=shift['vrad_'+x]
                c=299792.458
                delta_v=v_rad/c
                wave_opt=original_wavelength_strenched*(1-delta_v)
                masks_temp2=[]
                masks_temp=(
                    (
                        # Not in unreliable synthesis region
                        (~np.any(np.array([((wave_opt >= mask_beginning) & (wave_opt <= mask_end)) for (mask_beginning, mask_end) in zip(masks['mask_begin'],masks['mask_end'])]),axis=0))
                    ) |
                    # or is in vital line wavelengths
                    np.any(np.array([((wave_opt >= line_beginning) & (wave_opt <= line_end)) for (line_beginning, line_end) in zip(vital_lines['line_begin'],vital_lines['line_end'])]),axis=0)
                )
                for y in masks_temp:
                    if y:
                        masks_temp2.append(1)
                    else:
                        masks_temp2.append(0)
                masks_all[value]=masks_temp2
            return masks_all

        if  np.array_equal(synthetic_spectra_insert,None):
            for value,x in enumerate(colours):
                synthetic_spectra=rgetattr(self,x+'.synth')
                observed_spectra=rgetattr(self,x+'.spec')
                uncs=rgetattr(self,x+'.uncs')
                masks=rgetattr(self,x+'.masks')
                vital_lines=rgetattr(self,x+'.vital_lines')
                original_wavelength_strenched=rgetattr(self,x+'.wave')
                v_rad=rgetattr(self,x+'.vrad')
                c=299792.458
                delta_v=v_rad/c
                wave_opt=original_wavelength_strenched*(1-delta_v)
                masks_temp2=[]
                masks_temp=(
                    (
                        # Not too large difference between obs and synthesis
                        (~((np.abs(synthetic_spectra-observed_spectra)/np.sqrt(uncs) > 5) & (np.abs(synthetic_spectra-observed_spectra) > 0.4))) & 
                        # Not in unreliable synthesis region
                        (~np.any(np.array([((wave_opt >= mask_beginning) & (wave_opt <= mask_end)) for (mask_beginning, mask_end) in zip(masks['mask_begin'],masks['mask_end'])]),axis=0))
                    ) |
                    # or is in vital line wavelengths
                    np.any(np.array([((wave_opt >= line_beginning) & (wave_opt <= line_end)) for (line_beginning, line_end) in zip(vital_lines['line_begin'],vital_lines['line_end'])]),axis=0)
                )
                for y in masks_temp:
                    if y:
                        masks_temp2.append(1)
                    else:
                        masks_temp2.append(0)
                rsetattr(self, x+'.masked_area',masks_temp2)
                limits_temp=[]
                limits_first_loop=[]
                first=True
                for wave_temp,mask_temp in zip(wave_opt,masks_temp2):
                    if not mask_temp:
                        if first:
                            limits_first_loop.append(wave_temp)
                            first=False
                    else:
                        if not first:
                            limits_first_loop.append(wave_temp)
                            limits_temp.append(limits_first_loop)
                            limits_first_loop=[]
                            first=True
                rsetattr(self, x+'.masked_limits',limits_temp)
                    
        else:
            masks_all=np.array(np.ones(len(colours)),dtype=object)

            for value,x in enumerate(colours):
        
                synthetic_spectra=synthetic_spectra_insert[value]
                observed_spectra=normalized_observed_spectra_insert[value]
                uncs=uncs_insert[value]
                masks=rgetattr(self,x+'.masks')
                vital_lines=rgetattr(self,x+'.vital_lines')
                original_wavelength_strenched=rgetattr(self,x+'.wave')
                v_rad=shift['vrad_'+x]
                c=299792.458
                delta_v=v_rad/c
                wave_opt=original_wavelength_strenched*(1-delta_v)
                masks_temp2=[]
                masks_temp=(
                    (
                        # Not too large difference between obs and synthesis
                        (~((np.abs(synthetic_spectra-observed_spectra)/np.sqrt(uncs) > limits[0]) & (np.abs(synthetic_spectra-observed_spectra) > limits[1]))) & 
                        # Not in unreliable synthesis region
                        (~np.any(np.array([((wave_opt >= mask_beginning) & (wave_opt <= mask_end)) for (mask_beginning, mask_end) in zip(masks['mask_begin'],masks['mask_end'])]),axis=0))
                    ) |
                    # or is in vital line wavelengths
                    np.any(np.array([((wave_opt >= line_beginning) & (wave_opt <= line_end)) for (line_beginning, line_end) in zip(vital_lines['line_begin'],vital_lines['line_end'])]),axis=0)
                )
                for y in masks_temp:
                    if y:
                        masks_temp2.append(1)
                    else:
                        masks_temp2.append(0)
                masks_all[value]=masks_temp2
            return masks_all
    def synthesize(self,shift={},colours=None,multi=False,give_back=False,full=False,grad=False):
        if colours==None:
            colours=self.bands

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
    def correct_resolution_map(self,colours=None):
        if colours==None:
            colours=self.bands
        name=self.name
        not_resolved=[np.int64(name)]
        correct_spectra=[]
        # print(not_resolved)

        for count,x in enumerate(colours,1):
            if rgetattr(self,x+'.hermes')!=None:
                hermes=rgetattr(self,x+'.hermes')
                if len(hermes[7].data)!=len(hermes[1].data) or min(hermes[7].data)<0.02:
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
                    while len(hermes[7].data)!=len(hermes[1].data) or min(hermes[7].data)<=0.02:
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
                                    Path('dr6.1/'+name_target[0:6]+'/spectra/com/').mkdir(parents=True,exist_ok=True)
                                    if not exists('dr6.1/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'):
                                        #source='/media/storage/HERMES_REDUCED/dr6.1/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'
                                        source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.1/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'
                    
                                        destination='dr6.1/'+name_target[0:6]+'/spectra/com/'
                                        subprocess.run(["rsync",'-av',source,destination])
                                    try:
                                        if os.path.getsize('dr6.1/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits')>380000:
                                            hermes_temp=fits.open('dr6.1/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits')
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
                                Path('dr6.1/'+name_target[0:6]+'/spectra/com/').mkdir(parents=True,exist_ok=True)
                                if not exists('dr6.1/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'):
                                    #source='/media/storage/HERMES_REDUCED/dr6.1/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'
                                    source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.1/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits'
                
                                    destination='dr6.1/'+name_target[0:6]+'/spectra/com/'
                                    subprocess.run(["rsync",'-av',source,destination])
        
                                if os.path.getsize('dr6.1/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits')>380000:
                                    hermes_temp=fits.open('dr6.1/'+name_target[0:6]+'/spectra/com/'+name_target+str(count)+'.fits')
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
    def normalize(self,colours=None,data=None):
        if colours==None:
            colours=self.bands

        returning_spectra=np.array(np.ones(len(colours)),dtype=object)
        returning_uncs=np.array(np.ones(len(colours)),dtype=object)
        for value,x in enumerate(colours):
            if not np.array_equal(data,None):
                
                
                
                original_line=copy.deepcopy(rgetattr(self,x+'.spec_equalized'))
                x_line=rgetattr(self,x+'.wave')
                synth_line=data[value]
                uncertainty=copy.deepcopy(rgetattr(self,x+'.uncs_equalized'))
                
                renormalisation_fit = sclip((x_line,original_line/synth_line),chebyshev,int(3),ye=uncertainty,su=5,sl=5,min_data=100,verbose=False)
                
                returning_spectra[value]=original_line/renormalisation_fit[0]
                returning_uncs[value]=uncertainty/renormalisation_fit[0]
            elif rgetattr(self,x+'.synth').any():
                
                original_line=copy.deepcopy(rgetattr(self,x+'.spec_equalized'))
                x_line=rgetattr(self,x+'.wave')
                synth_line=rgetattr(self,x+'.synth')
                uncertainty=copy.deepcopy(rgetattr(self,x+'.uncs_equalized'))
                
                renormalisation_fit = sclip((x_line,original_line/synth_line),chebyshev,int(3),ye=uncertainty,su=5,sl=5,min_data=100,verbose=False)
                
                # uncs_normal=poly_synth/poly_original*rgetattr(self,x+'.uncs')
                rsetattr(self,x+'.spec',original_line/renormalisation_fit[0])
                rsetattr(self,x+'.uncs',uncertainty/renormalisation_fit[0])
                self.limit_array()
            else:
                print(x+' Hasnt been synthesized')
            
        if not np.array_equal(data,None):
            return returning_spectra,returning_uncs
    def gradient(self,shift={},colours=None,self_normal=True,labels=['teff','logg','monh','Fe','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']):
        if colours==None:
            colours=self.bands
        labels_payne=['teff','logg','monh','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
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
    def log_fit(self,colours=None,shift=None,dip_array=False,synthetic_spectra=None,normal=None,self_normal=False,solar_shift=None,uncertainty=None,limit_array=None,combine_masks=False):
        if colours==None:
            colours=self.bands
 
        limit_combine_array=np.array(np.ones(len(colours)),dtype=object)
        for value,x in enumerate(colours):
            if combine_masks:
                   if np.array_equal(limit_array,None):
                       limit_combine_array[value]=rgetattr(self,x+'.limit')*rgetattr(self,x+'.masked_area')
                   else:
                       limit_combine_array[value]=limit_array[value]*rgetattr(self,x+'.masked_area')
            else:
               if np.array_equal(limit_array,None):
                   limit_combine_array[value]=rgetattr(self,x+'.limit')
               else:
                   limit_combine_array[value]=limit_array[value]
        if not shift is None:
           synthetic_spectra=self.synthesize(shift,colours=colours,give_back=True)
        probability=np.ones(len(colours))
        if not dip_array:
           if np.array_equal(normal,None):
               for value,x in enumerate(colours):
                   if not(np.array_equal(synthetic_spectra,None)):
                       probability[value]= -np.sum((synthetic_spectra[value]-rgetattr(self,x+'.spec'))**2/(2*(rgetattr(self,x+'.uncs'))**2)*limit_combine_array[value])
                   else:
                       probability[value]= -np.sum((rgetattr(self,x+'.synth')-rgetattr(self,x+'.spec'))**2/(2*(rgetattr(self,x+'.uncs'))**2)*limit_combine_array[value])
           else:
               for value,x in enumerate(colours):
                   if not(np.array_equal(synthetic_spectra,None)):
                       probability[value]= -np.sum((synthetic_spectra[value]-normal[value])**2/(2*(uncertainty[value])**2)*limit_combine_array[value])
                   else:
                       probability[value]= -np.sum((rgetattr(self,x+'.synth')-normal[value])**2/(2*(uncertainty[value])**2)*limit_combine_array[value])
        else:
           if np.array_equal(normal,None):
               for value,x in enumerate(colours):
                   if not(np.array_equal(synthetic_spectra,None)):
                       probability[value]= -np.sum((synthetic_spectra[value]-rgetattr(self,x+'.spec'))**2/(2*(rgetattr(self,x+'.uncs'))**2)*limit_combine_array[value]*rgetattr(self,x+'.dip_array'))
                   else:
                       probability[value]= -np.sum((rgetattr(self,x+'.synth')-rgetattr(self,x+'.spec'))**2/(2*(rgetattr(self,x+'.uncs'))**2)*limit_combine_array[value]*rgetattr(self,x+'.dip_array'))
           else:
               for value,x in enumerate(colours):
                   if not(np.array_equal(synthetic_spectra,None)):
                       probability[value]= -np.sum((synthetic_spectra[value]-normal[value])**2/(2*(rgetattr(self,x+'.uncs'))**2)*limit_combine_array[value]*rgetattr(self,x+'.dip_array'))
                   else:
                       probability[value]= -np.sum((rgetattr(self,x+'.synth')-normal[value])**2/(2*(rgetattr(self,x+'.uncs'))**2)*limit_combine_array[value]*rgetattr(self,x+'.dip_array'))
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
    def mass_setter(self,colours=None,shift={}):
        if colours==None:
            colours=self.bands

        for x in colours:
            for param in shift:
                if 'vrad' in param:
                    if x in param:
                        rsetattr(self,x+'.vrad',shift[param])
                else:
                    rsetattr(self,x+'.'+param,shift[param])

    def observed_spectra_giver(self,colours=None):
        if colours==None:
            colours=self.bands

        returning_spectra=np.array(np.ones(len(colours)),dtype=object)
        wavelengths=np.array(np.ones(len(colours)),dtype=object)
        uncs=np.array(np.ones(len(colours)),dtype=object)
        for value,x in enumerate(colours):
            returning_spectra[value]=rgetattr(self,x+'.spec')
            wavelengths[value]=rgetattr(self,x+'.wave')
            uncs[value]=rgetattr(self,x+'.uncs')
        return returning_spectra,wavelengths,uncs
    def plot(self,colours=None,lines='all',masks=False):
        if colours==None:
            colours=self.bands

        c=299792.458

        for x in colours:
            plt.figure()
            rv=rgetattr(self,x+'.vrad')
            x_line=rgetattr(self,x+'.wave')
            x_shifted=(1-rv/c)*x_line
            observed= runGetattr(self,x+'.spec')
            plt.plot(x_shifted, observed, label='Observed')

            # x_line=np.linspace(0,len(runGetattr(self,x+'.synth')[0])-1,num=len(runGetattr(self,x+'.synth')[0]))
            if rgetattr(self,x+'.synth').any():
                    # labels='synthetic  chi squared= '+str(self.log_fit([x]))
                    synthetic=runGetattr(self,x+'.synth')
                    plt.plot(x_shifted, runGetattr(self,x+'.synth'), label='Synthetic')
                    min_synth=min( runGetattr(self,x+'.synth'))
            if lines=='all':
                important_lines=rgetattr(self,x+'.important_lines')
                for individual_line in important_lines:
                    individual_line_temp=float(individual_line[1])
                    minimum=min(abs(x_line-individual_line_temp))
                    min_synth=min(observed[np.where(abs(x_line-individual_line_temp)==minimum)][0],synthetic[np.where(abs(x_line-individual_line_temp)==minimum)][0])
                    plt.axvline(x=individual_line_temp,c='pink')
                    plt.text(float(individual_line_temp),min_synth-0.05,individual_line[0][:2],fontsize=20,ha='center',color='pink')
            else:
                important_lines=rgetattr(self,x+'.important_lines')
                for individual_line in important_lines:
                    if individual_line[0][:2] in lines:
                        
                        plt.axvline(x=float(individual_line[1]))
                        plt.text(float(individual_line[1]),0.5,individual_line[0][:2],fontsize=20,ha='center',color='pink')
            if masks:
                vital_lines=rgetattr(self,x+'.vital_lines')
                for line in vital_lines:
                    plt.axvspan(line['line_begin'],line['line_end'],alpha=0.7,color='blue',label='vital lines')
                masks_line=rgetattr(self,x+'.masks')
                for line in masks_line:
                    plt.axvspan(line['mask_begin'],line['mask_end'],alpha=0.2,color='red',label='masks for bad spectra overall')
                synthetic_spectra=rgetattr(self,x+'.synth')
                observed_spectra=rgetattr(self,x+'.spec')
                uncs=rgetattr(self,x+'.uncs')
                masks_bad_synth=(
                        # Not too large difference between obs and synthesis
                        (~((np.abs(synthetic_spectra-observed_spectra)/np.sqrt(uncs) > 5) & (np.abs(synthetic_spectra-observed_spectra) > 0.4))))                 
                limits_temp=[]
                limits_first_loop=[]
                first=True
                for wave_temp,mask_temp in zip(x_shifted,masks_bad_synth):
                    if not mask_temp:
                        if first:
                            limits_first_loop.append(wave_temp)
                            first=False
                    else:
                        if not first:
                            limits_first_loop.append(wave_temp)
                            limits_temp.append(limits_first_loop)
                            limits_first_loop=[]
                            first=True
                for area in limits_temp:
                    plt.axvspan(area[0],area[1],color='orange',alpha=0.2,label='bad  current synthetic spectra')
                limit_array=rgetattr(self,x+'.limit')
                limits_temp=[]
                limits_first_loop=[]
                first=True
                for wave_temp,mask_temp in zip(x_shifted,limit_array):
                    if not mask_temp:
                        if first:
                            limits_first_loop.append(wave_temp)
                            first=False
                    else:
                        if not first:
                            limits_first_loop.append(wave_temp)
                            limits_temp.append(limits_first_loop)
                            limits_first_loop=[]
                            first=True
                for area in limits_temp:
                    plt.axvspan(area[0],area[1],color='orange',alpha=0.2,label='limits')
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            plt.legend(by_label.values(), by_label.keys())
            plt.title(x+" Band")
            plt.xlim([x_shifted[0],x_shifted[-1]])
            plt.tight_layout()
def chebyshev(p,ye,mask):
    coef=np.polynomial.chebyshev.chebfit(p[0][mask], p[1][mask], 4)
    cont=np.polynomial.chebyshev.chebval(p[0],coef)
    return cont            
def starter_walkers_maker(nwalkers,old_abundances,parameters=['teff','logg','monh','fe_h','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic'],cluster=False):
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
    rv_labels={'vrad_Blue':0,'vrad_Green':1,'vrad_Red':2,'vrad_IR':3}
    parameters_Payne=['teff','logg','fe_h','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
    elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
    if cluster:
        label_change={'teff':'teff_spectroscopic','logg':'logg_spectroscopic','vrad_Green':'red_rv_ccd'}
    else:
        label_change={'alpha':'alpha_fe_r','fe_h':'fe_h_r','vsini':'vbroad_r','vmic':'vmic_r','vmac':'vmic_r','vrad_Blue':'rv',
                           'vrad_Green':'rv','vrad_Red':'rv','vrad_IR':'rv','monh':'fe_h','teff':'teff_r','logg':'logg_r'}
    limits_Payne={x:[(x_max[x]+x_min[x])/2,(x_max[x]-x_min[x])/10] for x  in parameters if x in parameters_Payne}
    limits_rv={x:[0,0] for x  in parameters if not  x in parameters_Payne}
    # parameters_Payne=[x for x in parameters_Payne if x in parameters]
    # labels=['teff','logg','monh','fe_h','alpha','vsini','vmac','vmic']
    while np.shape(pos)[0]<nwalkers or len(np.shape(pos))==1:
        dictionary_parameters={}
        for value,x in enumerate(limits_Payne):
            if x in label_change:
                    labels=label_change[x]
            else:
                labels=x
            if x in elements:
                if cluster and isinstance(old_abundances[x.lower()+'_fe'],np.float32) and old_abundances[x.lower()+'_fe'] >x_min[x] and old_abundances[x.lower()+'_fe']<x_max[x]:
                    dictionary_parameters[x]=np.random.normal(old_abundances[x.lower()+'_fe'],old_abundances['e_'+x.lower()+'_fe'],1)
                else:
                    dictionary_parameters[x]=np.random.normal(0,0.1,1)
            elif old_abundances[labels] and not np.isnan( old_abundances[labels]) and old_abundances[labels]>x_min[x] and old_abundances[labels]<x_max[x]:
                dictionary_parameters[x]=abs(np.random.normal(old_abundances[labels],limits_Payne[x][1],1))
            else:
                dictionary_parameters[x]=np.random.normal(limits_Payne[x][0],limits_Payne[x][1],1)
        for x in limits_rv:
            if cluster:
                if not np.isnan( old_abundances['red_rv_ccd'][rv_labels[x]]):
                    dictionary_parameters[x]=np.random.normal(old_abundances['red_rv_ccd'][rv_labels[x]],1.0,1)
                elif not np.isnan(old_abundances['red_rv_ccd']):
                    dictionary_parameters[x]=np.random.normal(old_abundances['red_rv_ccd'],1.0,1)
    
                else:
                    dictionary_parameters[x]=np.random.normal(0,20.0,1)            
            else:
                if not np.isnan( old_abundances['rv'][rv_labels[x]]):
                    dictionary_parameters[x]=np.random.normal(old_abundances['rv'][rv_labels[x]],1.0,1)
                elif not np.isnan(old_abundances['rv_com']):
                    dictionary_parameters[x]=np.random.normal(old_abundances['rv_com'],1.0,1)
    
                else:
                    dictionary_parameters[x]=np.random.normal(0,20.0,1)            
        pos_try_number=np.hstack([dictionary_parameters[x] for x in parameters])
        if starting_test(pos_try_number,spectras.old_abundances,parameters=parameters,cluster=cluster):
            if len(pos)==0:
                pos=pos_try_number
            else:
                pos=np.vstack((pos_try_number,pos))
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
def plot_spectrum(wave,flux,flux_uncertainty,unmasked_region,title_text,comp1_text,comp2_text,important_lines,neglect_ir_beginning=False):
    """
    Let's plot a spectrum, that is, flux over wavelenth
    
    We will plot 12 different subplot ranges (3 for each CCD) to allow better assessment of the results
    
    INPUT:
    wave : 1D-array with N pixels
    flux : 1D-array with N pixels or (M,N)-array with N pixels for M spectra (e.g. M = 2 for observed and synthetic spectrum)
    """
    
    # Let's define the wavelength beginnings and ends for each suplot
    if neglect_ir_beginning:
        subplot_wavelengths = np.array([
            [4710,4775],
            [4770,4850],
            [4840,4905],
            [5645,5730],
            [5720,5805],
            [5795,5878],
            [6470,6600],
            [6590,6670],
            [6660,6739],
            [7677,7720],
            [7710,7820],
            [7810,7890]
        ])
    else:
        subplot_wavelengths = np.array([
            [4710,4775],
            [4770,4850],
            [4840,4905],
            [5645,5730],
            [5720,5805],
            [5795,5878],
            [6470,6600],
            [6590,6670],
            [6660,6739],
            [7577,7697],
            [7677,7720],
            [7710,7820],
            [7810,7890]
        ])
    
    # How many subplots will we need?
    nr_subplots = np.shape(subplot_wavelengths)[0]
    
    f, gs = plt.subplots(nr_subplots,1,figsize=(8.3,11.7),sharey=True)
    
    try:
        # test if several spectra fed into flux
        dummy = np.shape(flux)[1] == len(wave)
        flux_array_indices = np.shape(flux)[0]
        flux = np.array(flux)
    except:
        flux_array_indices = 1

    # Let's loop over the subplots
    for subplot in range(nr_subplots):
        
        # Which part of the observed/model spectrum is in our subplot wavelength range?
        in_subplot_wavelength_range = (wave > subplot_wavelengths[subplot,0]) & (wave < subplot_wavelengths[subplot,1])

        ax = gs[subplot]
        ax.set_xlim(subplot_wavelengths[subplot,0],subplot_wavelengths[subplot,1])
        
        if len(wave[in_subplot_wavelength_range]) > 0:
            # if only 1 spectrum
            if flux_array_indices == 1:
                ax.plot(wave[in_subplot_wavelength_range],flux[in_subplot_wavelength_range],lw=0.5);
            else:
                for index in range(flux_array_indices):
                    if index == 0:
                        ax.plot(wave[in_subplot_wavelength_range],flux[0,in_subplot_wavelength_range],lw=0.5,c='k',label='data');
                        ax.plot(wave[in_subplot_wavelength_range],1.05 + flux_uncertainty[in_subplot_wavelength_range],lw=0.5,c='C3',label='scatter');
                    if index == 1:
                        ax.plot(wave[in_subplot_wavelength_range],flux[index,in_subplot_wavelength_range],lw=0.5,c='r',label='model (optimised)');
                        ax.plot(wave[in_subplot_wavelength_range],1.05 + np.abs(flux[0,in_subplot_wavelength_range]-flux[index,in_subplot_wavelength_range]),lw=0.5,c='C4',label='residuals');
                if subplot == nr_subplots-1:
                    ax.legend(ncol=2,loc='lower right',fontsize=6)

            maski = 0
            for maski, pixel in enumerate(wave[in_subplot_wavelength_range & unmasked_region]):
                if maski == 0:
                    ax.axvline(pixel,color='C0',alpha=0.1,label='Mask')
                    maski += 1
                else:
                    ax.axvline(pixel,color='C0',alpha=0.1)
            each_index = 0 
            for each_element in important_lines:
                if (each_element[0] > subplot_wavelengths[subplot,0]) & (each_element[0] < subplot_wavelengths[subplot,1]):

                    offset = -0.05+0.1*(each_index%3)
                    each_index+=1
                    ax.axvline(each_element[0],lw=0.2,ls='dashed',c='r')
                    if each_element[1] in ['Li','C','O']:
                        ax.text(each_element[0],offset,each_element[1],fontsize=10,ha='center',color='pink')
                    elif each_element[1] in ['Mg','Si','Ca','Ti','Ti2']:
                        ax.text(each_element[0],offset,each_element[1],fontsize=10,ha='center',color='b')
                    elif each_element[1] in ['Na','Al','K']:
                        ax.text(each_element[0],offset,each_element[1],fontsize=10,ha='center',color='orange')
                    elif each_element[1] in ['Sc','V', 'Cr','Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']:
                        ax.text(each_element[0],offset,each_element[1],fontsize=10,ha='center',color='brown')
                    elif each_element[1] in ['Rb', 'Sr', 'Y', 'Zr', 'Ba', 'La', 'Ce','Mo','Ru', 'Nd', 'Sm','Eu']:
                        ax.text(each_element[0],offset,each_element[1],fontsize=10,ha='center',color='purple')
        ax.set_ylim(-0.1,1.2)
        if subplot == nr_subplots-1:
            ax.set_xlabel(r'Wavelength / $\mathrm{\AA}$')
        ax.set_ylabel('Flux / norm.')
    f.suptitle(title_text+' \n '+comp1_text+' \n '+comp2_text)
    plt.tight_layout(h_pad=0)
    
    return f
def load_dr3_lines(mode_dr3_path = 'important_lines'):
    """
    
    """
    important_lines = [
        [4861.3230,r'H$_\beta$',r'H$_\beta$'],
        [6562.7970,r'H$_\alpha$',r'H$_\alpha$']
    ]
    
    important_molecules = [
        [4710,4740,'Mol. C2','Mol. C2'],
        [7594,7695,'Mol. O2 (tell.)','Mol. O2 (tell.)']
        ]

    line, wave = np.loadtxt(mode_dr3_path,usecols=(0,1),unpack=True,dtype=str, comments=';')

    for each_index in range(len(line)):
        if line[each_index] != 'Sp':
            if len(line[each_index]) < 5:
                important_lines.append([float(wave[each_index]), line[each_index], line[each_index]])
            else:
                important_lines.append([float(wave[each_index]), line[each_index][:-4], line[each_index]])
        
    return(important_lines,important_molecules)



labels=['teff','logg','monh','fe_h','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']  

open_cluster_sobject_id=np.loadtxt('target_sobject_id.txt',delimiter=',',dtype=str)

cluster_name='Melotte_22'
global large_data
votable = parse(cluster_name+"_photometric_cross.xml")
photometric_data=votable.get_first_table().to_table(use_names_over_ids=True)
          

# large_data=[x for x in large_data if x['sobject_id'] in open_cluster_sobject_id]
# large_data=vstack(large_data)
#Needs to be read because some spectra doesnt have 
global all_reduced_data
all_reduced_data=fits.getdata('dr6.1.fits',1)
all_reduced_data=Table(all_reduced_data)
global x_min,x_max
tmp = np.load("NN_normalized_spectra_all_elements_Blue.npz")   
labels_with_limits=['teff','logg','fe_h','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']

# votable = parse(cluster_name+"_cross_galah_light.xml")
file_directory = cluster_name+'_reduction_priors/plots/'
Path(file_directory).mkdir(parents=True, exist_ok=True)
# cross_data=votable.get_first_table().to_table(use_names_over_ids=True)

# cross_data=cross_data

x_min=tmp['x_min']
x_max=tmp['x_max']
x_min={x:y for x,y in zip(labels_with_limits,x_min)}
x_max={x:y for x,y in zip(labels_with_limits,x_max)}

# #EMCEE,
# prior=False
np.random.seed(589404)
final_results=[]
expected_results=[]
large_data_GALAH_official=fits.open('gaia_galah_cross_values.fits')
large_data_GALAH_official_unchanged=large_data_GALAH_official[1].data
large_data_GALAH_official=Table(large_data_GALAH_official[1].data)
nwalkers=22

spectras=spectrum_all(photometric_data['sobject_id'][0],cluster=True)
old_abundances=spectras.old_abundances
 
colours=spectras.bands
reduction_status=np.any([rgetattr(spectras,x+'.bad_reduction') for x in colours ])or spectras.hermes_checker()

# if reduction_status:
#       print('reduction failed will skip'+str(name)+ 'for now')
#       continue


parameters=['teff','logg','fe_h','vmic','vsini','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']

spectras.spread_number=10

full_parameters=['teff','logg','monh','fe_h','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']
parameters=['teff','logg','monh','alpha','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','vsini','vmac','vmic']

# log_posterior(a[0],parameters,fe_hold=True)
# name=cross_data['sobject_id'][0]
# spectras=spectrum_all(name,10,cluster=False)
# spectras.synthesize()
# spectras.limit_array()
# spectras.normalize()
# spectras.create_masks()
# spectras.plot(masks=True)




parameters=['teff','logg','fe_h','vmic','vsini','vrad_Blue','vrad_Green','vrad_Red','vrad_IR','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
params_short=['teff','logg','fe_h','vmic','vsini','vrad_Blue','vrad_Green','vrad_Red','vrad_IR']
parameters_no_vrad=['teff','logg','fe_h','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
parameters_no_vrad=['teff','logg','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']


nwalkers=len(parameters)*2
prior=True
# for star in photometric_data:
#     name=star['sobject_id']
# name=photometric_data[0]['sobject_id']
name=160106001601280
filename = cluster_name+'_reduction/prior_'+str(name)

name='Melotte_22'

# NGC_2682_summary_dr61
votable = parse(name+"_all4.xml")
data=votable.get_first_table().to_table(use_names_over_ids=True)
difference=[np.zeros(4096) for x in range(3)]
difference.append(np.zeros(3096))
masks_original=[np.ones(4096) for x in range(3)]
masks_original.append(np.zeros(3096))

for star in data:
    spectras=spectrum_all(star['sobject_id'],cluster=True)
    shift={}
    for param in parameters:
        param_index=copy.copy(param)
        if param in elements:
            param_index+='_Fe'
        shift[param]=star[param_index+'_prior']
    synthetic_spectras=spectras.synthesize(shift,give_back=True)
    normalized_spectra,normalized_uncs=spectras.normalize(data=synthetic_spectras)
    normalized_limit_array=spectras.limit_array(observed_spectra=normalized_spectra)
    normalized_mask=spectras.create_masks(clean=True,shift=shift)

    for value,(syn,norm,unc) in enumerate(zip(synthetic_spectras,normalized_spectra,normalized_uncs)):
        difference[value]+=(syn-norm)/unc
    masks_original=[ x*y*z for (x,y,z) in zip(normalized_mask,masks_original,normalized_limit_array)]
difference=[x/len(data) for x in difference]
c=299792.458
new_masks=[np.zeros(4096) for x in range(3)]
new_masks.append(np.zeros(3096))
colours=spectras.bands
for value, col in enumerate(colours):
    plt.figure()
    plt.title(col)
    vital_lines=rgetattr(spectras,col+'.vital_lines')
    x_line=rgetattr(spectras,col+'.wave')
    rv=shift['vrad_'+col]
    x_shifted=(1-rv/c)*x_line
    plt.plot(x_shifted,difference[value])
    
    for line in vital_lines:
        plt.axvspan(line['line_begin'],line['line_end'],alpha=0.7,color='blue',label='vital lines')
    masks_line=rgetattr(spectras,col+'.masks')
    for line in masks_line:
        plt.axvspan(line['mask_begin'],line['mask_end'],alpha=0.2,color='red',label='masks for bad spectra overall')
    for value2,x in enumerate(difference[value]):
        if abs(x)>3:
            new_masks.append
