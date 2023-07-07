#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 16:59:16 2022

@author: kevin
"""
from numba import jit
from functools import  partial
from astropy.io.votable import parse
import subprocess
import emcee
import os.path
from astropy.table import Table,vstack
from pysme import sme as SME
from pysme.synthesize import synthesize_spectrum
from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
from astropy.io import fits
import copy
import Galah_tool_py3 as gtools
import os
import random
import functools
import logging
from pathlib import Path
from os.path import exists
import scipy
from scipy import signal
from datetime import date
logger=logging.getLogger('pysme')
logger.setLevel(logging.DEBUG   )

os.environ["OMP_NUM_THREADS"] = "1"

# def synthesize_multi(move,values,changed_values):
#     shift={values:move}
#     shift.update(changed_values)
#     return synthesize(shift)
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
def galah_kern(fwhm, b):
        """ Returns a normalized 1D kernel as is used for GALAH resolution profile """
        size=2*(fwhm/2.355)**2
        size_grid = int(size) # we limit the size of kernel, so it is as small as possible (or minimal size) for faster calculations
        if size_grid<7: size_grid=7
        x= scipy.mgrid[-size_grid:size_grid+1]
        g = np.exp(-0.693147*np.power(abs(2*x/fwhm), b))
        return g / np.sum(g)
def shift_maker(solar_values):
    # shift={'vrad':solar_values}
    shift={'teff':solar_values[0],'logg':solar_values[1],'monh':solar_values[2],'Fe':solar_values[3],'alpha':solar_values[4],'vrad':solar_values[5],'vsini':solar_values[6],'vmac':solar_values[7],'vmic':solar_values[8]}

    return shift
def synthesize(shift):
    spectrum=copy.copy(sme)
    for key in shift:
        setattr(spectrum,key,shift[key])
    spectrum = synthesize_spectrum(spectrum)
    synthetic_temp=copy.copy(synthetic)

    synthetic_temp[1].data=spectrum.synth[0]
    if os.path.exists(name+str(shift)+'_synthetic.fits'):
        os.remove(name+str(shift)+'_synthetic.fits')
    synthetic_temp.writeto(name+str(shift)+'_synthetic.fits')

    syn=gtools.read(name+str(shift)+'_synthetic')
    syn.equalize_resolution()
    os.remove(name+str(shift)+'_synthetic.fits')
    return syn.f

def log_fit(synthetic_spectrum):
    return -0.5*np.sum((synthetic_spectrum-sme.spec[0])**2/(sme.uncs[0])**2)

def log_prior(shift):
    error=0
    for x in shift:
        if x in ('teff','logg'):
            error-=(shift[x]-photometric_values[x])**2/(2*(photometric_values[x+'_e'])**2)
    if error<-1e100:
        return -np.inf
    else:
        return error
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
def log_posterior(solar_values):
# pos=np.column_stack((teff_random,vsini_random,logg_random,v_rad_random,monh_random))
    limit=10
    shift={'teff':solar_values[0],'logg':solar_values[1],'monh':solar_values[2],'Fe':solar_values[3],'alpha':solar_values[4],'vrad':solar_values[5],'vsini':solar_values[6],'vmac':solar_values[7],'vmic':solar_values[8]}
    # print("this is SHIFT", shift)
    if np.any([True for x in shift if shift[x]<0.0 and x!='monh']):
        return -np.inf

    if 'teff' in shift:
         if abs(shift['teff']-old_abundances['teff'])/old_abundances['e_teff']>limit:
                print("temperature is too small or large")
                return -np.inf
    if 'logg' in shift:
        if shift['logg']<1.0 or shift['logg']>5.0:
            print("logg is too small or big")
            return -np.inf  
        if abs(shift['logg']-old_abundances['logg'])/old_abundances['e_logg']>limit:
            print("logg is too small or big")
            return -np.inf  

    if 'monh' in shift:
        if abs(shift['monh']-old_abundances['fe_h'])/old_abundances['e_fe_h']>limit:
            print("monh is too small or big")
            return -np.inf   
    if 'vsini' in shift:
        if abs(shift['vsini']-old_abundances['vbroad'])/old_abundances['e_vbroad']>limit:
            print("monh is too small or big")
            return -np.inf       
    if 'vmic' in shift:
        if shift['vmic']<0 or shift['vmic']>20:
            return -np.inf 

            print("vmic is negative")
        if abs(shift['vmic']-old_abundances['vmic'])/old_abundances['e_vmic']>limit:
            
            return -np.inf 
    # if 'vrad' in shift:
    #     if abs(shift['vrad']-old_abundances['rv_galah'])/old_abundances['e_rv_galah']>limit:
    #         return -np.inf
    if 'vmac' in shift:
        if shift['vmac']<0 or shift['vmac']>100:
            print("vmac is negative")
            return -np.inf 
    try:
        synthetic_spectras=spectras.synthesize(shift,give_back=True)
        # print('shift' , shift)
        # print("Spectra",synthetic_spectras)
    except:
        print('couldnt finish synthesizing')
        return -np.inf
    if not np.all([np.all(x) for x in synthetic_spectras]):
        print("Didnt synthesize properly")
        return -np.inf
    if np.max([np.max(abs(x)) for x in synthetic_spectras])>100 or np.max([np.max(abs(x)) for x in synthetic_spectras])==np.nan:
        print("Didnt synthesize properly one of the spectra blew up")
        return -np.inf
    if np.any(np.array_equal(synthetic_spectras,np.nan)):
        print("Didnt synthesize properly Nan values yay!")
        return -np.inf
    # print("Before normalized",synthetic_spectras, "Shift    ",shift)
    # print("max value",np.max([np.max(abs(x)) for x in synthetic_spectras]),"shift",shift)
    # normalized_spectra, normalized_uncertainty=spectras.normalize(data=synthetic_spectras)
    probability=spectras.log_fit(synthetic_spectra=synthetic_spectras)
    # print("This is the probability", probability, "shift prob", log_prior(shift), "shift", shift)
    # print("This is the Uncs", normalized_uncertainty,"This is the Observed",normalized_spectra)
    # print("Spectra",synthetic_spectras)
    return probability 
def starting_test(solar_values,abundances):
    
    shift={'teff':solar_values[0],'logg':solar_values[1],'monh':solar_values[2],'Fe':solar_values[3],'alpha':solar_values[4],'vrad':solar_values[5],'vsini':solar_values[6],'vmac':solar_values[7],'vmic':solar_values[8]}
    if np.any([True for x in shift if shift[x]<0.0 and x!='monh' and x!='alpha' and x!='Fe']):
        return False
    if 'Fe' in shift:
        if abs(shift['Fe'])>10:
            return False
    if 'teff' in shift:
         if abs(shift['teff']-abundances['teff_r'])>1000:
                print("temperature is too small or large")
                return False
         elif shift['teff']<2000 or shift['teff']>10000:
            return False
    if 'logg' in shift:
        if shift['logg']<1.0 or shift['logg']>5.0:
            print("logg is too small or big")
            return False  
    if 'monh' in shift:
        if abs(shift['monh']-abundances['fe_h_r'])>0.5:
            print("monh is too small or big")
            return False   
    if 'vsini' in shift:
        if shift['vsini']<0 or shift['vsini']>100:
            print("monh is too small or big")
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
def synth_resolution_degradation(wave, synth,res_map,res_b,l_new_premade=None,kernel_=None,synth_res=300000.0,grad=None):
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





        #oversampling. If synthetic spectrum sampling is much finer than the size of the kernel, the code would work, but would return badly sampled spectrum. this is because from here on the needed sampling is measured in units of sigma.
        oversample=galah_sampl/sampl*5.0

        #minimal needed sampling

        #keep adding samples until end of the wavelength range is reached
        if l_new_premade is None:
            #required sigma (resample the resolution map into the wavelength range of the synthetic spectrum)
            s_out=res_map

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
@jit(nopython=True,parallel=False,cache=True)
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
class spectrum_all:
    # bands=['Blue','Green','Red','IR']
    bands=['Blue','Green','Red','IR']
    alpha_elements_abundaces={'O':8.69,'Ne':7.93,'Mg':7.6,'Si':7.510,'S':7.120,'Ar':6.400,'Ca':6.340,'Ti':4.950,'Fe':7.450}
    fe_abund=7.450

    # bands=['Red']
    def __init__(self,name,interpolation,cluster=True):
        name=str(name)
        bands=['Blue','Green','Red','IR']

        # bands=['Blue','Green','Red','IR']

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

        file_names={'Blue':1,'Green':2,'Red':3,'IR':4}
        for x in bands:
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
            setattr(self,x,SME.SME_Structure())
            count=file_names[x]
            Path(name[0:6]+'/spectra/com/').mkdir(parents=True,exist_ok=True)
            if not exists(name[0:6]+'/spectra/com/'+name+str(count)+'.fits'):
                #source='/media/storage/HERMES_REDUCED/dr6.0/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'
                source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.0/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'

                destination=name[0:6]+'/spectra/com/'
                subprocess.run(["rsync",'-av',source,destination])


            hermes=fits.open(name[0:6]+'/spectra/com/'+name+str(count)+'.fits')
            x0= float(hermes[1].header.get('CRVAL1'))
            x1=float( hermes[1].header.get('CRVAL1')+len(hermes[1].data)* hermes[1].header.get('CDELT1'))
            fstart= x0+(x1-x0)*starting_fraction
            
            new_start=int(starting_fraction*len(hermes[1].data))
            new_end=new_start+int(len(hermes[1].data)*(1-starting_fraction)*length_fraction)
            n_points=(new_end-new_start)*interpolation

            # length=np.linspace(limits[x][0],limits[x][1],num=(1600+n_points))
            length=np.array([fstart+y*hermes[1].header.get('CDELT1')/interpolation for y in range(0,n_points)])
            fstart=length[0]
            fend=length[-1]
            #Increases the length of synthetic spectra so it over interpolates  
            rsetattr(self,x+'.wave',length)
            rsetattr(self,x+'.spec',[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])])
            rsetattr(self,x+'.uncs',[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[2].data[new_start:new_end])])
            hermes[0].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[0].data[new_start:new_end])
            hermes[1].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])
            hermes[2].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[2].data[new_start:new_end])
            if len(hermes[7].data[new_start:new_end])==new_end-new_start:
                hermes[7].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[7].data[new_start:new_end])
            hermes[1].header['CDELT1']/=interpolation
            
            rsetattr(self,x+'.hermes',hermes)
            rsetattr(self,x+'.hermes',hermes)
            rsetattr(self,x+'.l_new',None)

            rsetattr(self,x+'.wran',[length[0],length[-1]])
            rsetattr(self,x+'.linelist',ValdFile('Galah_'+x+'_6.lin'))

            if cluster:
                rsetattr(self,x+'.abund',Abund(float(old_abundances['fe_h']),'asplund2009' ))
                print(x+'_HFS.lin')
            # rsetattr(self,x+'.linelist',ValdFile(x+'_no_molecule.lin'))
                rsetattr(self,x+'.vmic',float(old_abundances['vmic_reduction']))
                rsetattr(self,x+'.vsini',float(old_abundances['vbroad_reduction']))
                rsetattr(self,x+'.vrad',float(old_abundances['rv_r']))
                rsetattr(self,x+'.vmac',float(old_abundances['vmic_reduction']))
                rsetattr(self,x+'.teff',float(old_abundances['teff_reduction']))
                rsetattr(self,x+'.logg',float(old_abundances['logg_reduction']))
                rsetattr(self,x+'.monh',float(old_abundances['fe_h']))
            else:
                rsetattr(self,x+'.abund',Abund(float(old_abundances['fe_h_r']),'asplund2009' ))
                print(x+'_HFS.lin')
            # rsetattr(self,x+'.linelist',ValdFile(x+'_no_molecule.lin'))
                rsetattr(self,x+'.vmic',float(old_abundances['vmic_r']))
                rsetattr(self,x+'.vsini',float(old_abundances['vbroad_r']))
                if not np.isnan(float(old_abundances['rv'][count-1])):
                    rsetattr(self,x+'.vrad',float(old_abundances['rv'][count-1]))
                elif not np.isnan(float(old_abundances['rv_com'])):
                    rsetattr(self,x+'.vrad',float(old_abundances['rv_com']))
                else:
                    rsetattr(self,x+'.vrad',0.0)
                rsetattr(self,x+'.vmac',float(old_abundances['vmic_r']))
                rsetattr(self,x+'.teff',float(old_abundances['teff_r']))
                rsetattr(self,x+'.logg',float(old_abundances['logg_r']))
                rsetattr(self,x+'.monh',float(old_abundances['fe_h_r']))
            rsetattr(self,x+'.h2broad',True)
            rsetattr(self,x+'.vrad_flag','whole')
            rsetattr(self,x+'.cscale_flag','none')
            rsetattr(self,x+'.cscale_type','mask')
            rsetattr(self,x+'.atmo.source',"marcs2012.sav")
            rsetattr(self,x+'.atmo.method','grid')
            rsetattr(self,x+'.atmo.geom','pp')
            runGetattr(self,x+'.nlte.set_nlte')("Fe", "nlte_Fe_ama51_Feb2022_pysme.grd")
            runGetattr(self,x+'.nlte.set_nlte')("Si", "marcs2012_Si2016.grd")
            runGetattr(self,x+'.nlte.set_nlte')("Ca", 'marcs2012p_t1.0_Ca.grd')
            runGetattr(self,x+'.nlte.set_nlte')("Na","marcs2012_Na2011.grd")
            runGetattr(self,x+'.nlte.set_nlte')("O",'marcs2012_O2015.grd')
            runGetattr(self,x+'.nlte.set_nlte')("Ti","marcs2012s_t2.0_Ti.grd")
            runGetattr(self,x+'.nlte.set_nlte')("Ba","marcs2012p_t1.0_Ba.grd")
            runGetattr(self,x+'.nlte.set_nlte')("Li","marcs2012_Li.grd")
            runGetattr(self,x+'.nlte.set_nlte')("Mg","marcs2012_Mg2016.grd")
        self.correct_resolution_map()
        self.equilize_spectras()        

    def equilize_spectras(self,colours=bands):
        for x in colours:
            wavelength=rgetattr(self, x+'.wave')[0]
            hermes=rgetattr(self, x+'.hermes')
            resmap=hermes[7].data
            observed_spectra=rgetattr(self, x+'.spec')[0]
            observed_error=rgetattr(self, x+'.uncs')[0]
            new_observed,new_error=equalize_resolution(wavelength,resmap,observed_spectra,observed_error,hermes)
            rsetattr(self, x+'.spec',new_observed)
            rsetattr(self, x+'.uncs',new_error)            
    def solar_value_maker(self,shift,colour,keys=['teff','logg','monh','Fe','alpha','vrad','vsini','vmac','vmic']):
        solar_values=[]
        for x in keys:
            if x in shift:
                solar_values.append(shift[x])
            else:
                solar_values.append(rgetattr(self,colour+'.'+x))
        return solar_values        
    def correct_resolution_map(self,colours=bands):
        name=self.name
        not_resolved=[np.int64(name)]
        # print(not_resolved)

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
                        # print('here')
                        print('copying from '+name_target)
                        not_resolved.append(np.int64(name_target))
                        # print(not_resolved)
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

                    hermes_temp=fits.open(name_target[0:6]+'/spectra/com/'+name_target+str(1)+'.fits')
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

    def synthesize(self,shift,multi=False,colours=bands,fe_abund=fe_abund,alpha_elements_abundaces=alpha_elements_abundaces,give_back=False):
        if not (multi or give_back):
            for x in colours:       
                spectrum=copy.copy(getattr(self,x))
                for key in shift:
                    if key=='Fe':
                        spectrum.abund.update_pattern({'Fe':fe_abund+shift[key]-spectrum.monh})
                    elif key=='alpha':
                        alpha_shift=spectrum.abund.get_element('Fe')-fe_abund+shift[key]
                        alpha_shift_dict={x:alpha_elements_abundaces[x]+alpha_shift-spectrum.monh for x in alpha_elements_abundaces}
                        spectrum.abund.update_pattern(alpha_shift_dict)
                    else:
                        setattr(spectrum,key,shift[key])
                spectrum = synthesize_spectrum(spectrum)
                if rgetattr(self,x+'.l_new') is None:
                    dopler_shifted_spectra,l_new,kernel=synth_resolution_degradation(rgetattr(self,x+'.wave')[0],spectrum.synth[0],rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'])
                    rsetattr(self,x+'.l_new',l_new)
                    rsetattr(self,x+'.kernel',kernel)
                else:
                    dopler_shifted_spectra=synth_resolution_degradation(rgetattr(self,x+'.wave')[0],spectrum.synth[0],rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'],rgetattr(self,x+'.l_new'),rgetattr(self,x+'.kernel'))

                rsetattr(self,x+'.synth',dopler_shifted_spectra)
        elif give_back:
            returning_spectra=np.array(np.ones(len(colours)),dtype=object)
            for number,x in enumerate(colours):
                spectrum=copy.copy(getattr(self,x))
                for key in shift:
                    if key=='Fe':
                        spectrum.abund.update_pattern({'Fe':fe_abund+shift[key]-spectrum.monh})
                    elif key=='alpha':
                        alpha_shift=spectrum.abund.get_element('Fe')-fe_abund+shift[key]
                        alpha_shift_dict={x:alpha_elements_abundaces[x]+alpha_shift-spectrum.monh for x in alpha_elements_abundaces}
                        spectrum.abund.update_pattern(alpha_shift_dict)
                    else:
                        setattr(spectrum,key,shift[key])
                spectrum = synthesize_spectrum(spectrum)
                if rgetattr(self,x+'.l_new') is None:
                    dopler_shifted_spectra,l_new,kernel=synth_resolution_degradation(rgetattr(self,x+'.wave')[0],spectrum.synth[0],rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'])
                    rsetattr(self,x+'.l_new',l_new)
                    rsetattr(self,x+'.kernel',kernel)
                else:
                    dopler_shifted_spectra=synth_resolution_degradation(rgetattr(self,x+'.wave')[0],spectrum.synth[0],rgetattr(self,x+'.hermes')[7].data,rgetattr(self,x+'.hermes')[7].header['b'],rgetattr(self,x+'.l_new'),rgetattr(self,x+'.kernel'))

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
        returning_uncs=np.array(np.ones(len(colours)),dtype=object)

        for value,x in enumerate(colours):
            if not np.array_equal(data,None):
                                
                original_line=rgetattr(self,x+'.spec')[0]
                x_line=rgetattr(self,x+'.wave')[0]  
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
                
                returning_uncs[value]=rgetattr(self,x+'.uncs')[0]
                
                
            elif rgetattr(self,x+'.synth'):
                
                original_line=rgetattr(self,x+'.spec')[0]

                x_line=rgetattr(self,x+'.wave')[0]  
                synth_line=rgetattr(self,x+'.synth')[0]
                
                
                poly_order=5
                syth_coeff=np.polyfit(x_line,synth_line,poly_order)
                ynew=np.poly1d(syth_coeff)
                poly_synth=ynew(x_line)
         
                original_coeff=np.polyfit(x_line,original_line,poly_order)
                ynew=np.poly1d(original_coeff)
                poly_original=ynew(x_line)
                
                synth_normal=poly_synth/poly_original*original_line
                
                # uncs_normal=poly_synth/poly_original*rgetattr(self,x+'.uncs')[0]
                rsetattr(self,x+'.spec',synth_normal)
                # rsetattr(self,x+'.uncs',uncs_normal)
                
            else:
                print(x+' Hasnt been synthesized')
        if not np.array_equal(data,None):
            return returning_spectra,returning_uncs

    def log_fit(self,colours=bands,synthetic_spectra=None):
        probability=0
        for value,x in enumerate(colours):
            if not(np.array_equal(synthetic_spectra,None)):
                probability+= -np.sum((synthetic_spectra[value]-rgetattr(self,x+'.spec')[0])**2/(rgetattr(self,x+'.uncs')[0])**2)
            else:
                    probability+= -np.sum((rgetattr(self,x+'.synth')[0]-rgetattr(self,x+'.spec')[0])**2/(rgetattr(self,x+'.uncs')[0])**2)
        return probability
    def plot(self,colours=bands,dopler_shift=True):
        for x in colours:
            plt.figure()
            if rgetattr(self,x+'.synth'):
                    # labels='synthetic  chi squared= '+str(self.log_fit([x]))
                    plt.plot(rgetattr(self,x+'.wave')[0], runGetattr(self,x+'.synth')[0], label='Synthetic')
            plt.fill_between(runGetattr(self,x+'.wave')[0], runGetattr(self,x+'.spec')[0]- runGetattr(self,x+'.uncs'),runGetattr(self,x+'.spec')[0]+ runGetattr(self,x+'.uncs'),alpha=0.5)
            plt.plot(rgetattr(self,x+'.wave')[0], runGetattr(self,x+'.spec')[0], label='Observed')
            plt.title(x+" Band")
            plt.legend(loc='best')
            plt.tight_layout()
            
            
            

name=160108002601369
filename = "160108002601369"


c=299792  #km/s

#Needs to be read because some spectra doesnt have 
global all_reduced_data
all_reduced_data=fits.getdata('dr6.0.fits',1)
all_reduced_data=Table(all_reduced_data)

interpolation=1
sme=SME.SME_Structure()
global large_data
large_data=fits.getdata('photometric_reduction_dr6_cross.fits',1)
large_data=Table(large_data)

name=160108002601369

# name=''
spectras=spectrum_all(name,1,cluster=False)
# spectras.synthesize({})

# svens=fits.getdata('svens_fe_alpha_with_mine.fits',1)
# svens=Table(svens)
# svens['e_teff_sven']=np.array(svens['e_teff_sven'],dtype=float)
# # spectras.synthesize({'teff':6244.17,'logg':4.24,'monh':0.053,'Fe':0.053,'alpha':0.0064,'vmic':1.29,'vmac':6.0,'vsini':6.3})









