#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 17:11:34 2021

@author: kevin
"""

from functools import  partial
from astropy.io.votable import parse
import emcee
import os.path
import logging
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
from multiprocessing.dummy import Pool as ThreadPool 

logger=logging.getLogger('pysme')
logger.setLevel(0)

os.environ["OMP_NUM_THREADS"] = "1"


def log_prior(shift,photometric_values):
    '''
    # Returns the probability of the solar parameters given the prior
    needs to give e.g shift={'teff':5700,'logg':4.45} 
    and an astropy table with the prior photometric values

    '''
    error=0
    for x in shift:
        if x in ('teff','logg'):
            #both priors asuumes that it is gaussian. The first one is a prior from photometry and the second is a backgrouhnd prior
            error-=(shift[x]-photometric_values[x+'_gaia'])**2/(2*(photometric_values[x+'_e_gaia'])**2)
           # error-=(shift[x]-photometric_values[x])**2/(2*(photometric_values[x+'_e']*10)**2)

    # if the prior is smaller than exp(-1e100) will just return -np.inf
    if error<-1e100:
        return -np.inf
    else:
        return error
    
def log_posterior(solar_values,galah_spectra):
    '''
    #Main section of the MCMC code, You enter your solar values and you will recive the probability that this spectra fits the observed spectrum

    '''
    
    #checks if the entered solar values are in the parameter range
    if not starting_test(solar_values,galah_spectra.old_abundances):
        print('Your solar parameters is out side of the range')
        return -np.inf    
    shift={'teff':solar_values[0],'vsini':solar_values[1],'logg':solar_values[2],'vrad':solar_values[3],'monh':solar_values[4],'vmac':solar_values[5],'vmic':solar_values[6]}
    
    #Synthesizes the spectra for all 4 bands
    try:
        synthetic_spectras=galah_spectra.synthesize(shift,give_back=True)
    except:
        print('couldnt finish synthesizing')
        return -np.inf
    #performs basic checks on the spectra to make sure if the spectra synthesizes properly
    if not np.all([np.all(x) for x in synthetic_spectras]):
        print("Didnt synthesize properly")
        return -np.inf
    if np.max([np.max(abs(x)) for x in synthetic_spectras])>100 or np.max([np.max(abs(x)) for x in synthetic_spectras])==np.nan:
        print("Didnt synthesize properly one of the spectra blew up")
        return -np.inf
    if np.any(np.array_equal(synthetic_spectras,np.nan)):
        print("Didnt synthesize properly Nan values yay!")
        return -np.inf
    #gets the probability that the synthetic spectra fits the observed spectra
    probability=galah_spectra.log_fit(synthetic_spectra=synthetic_spectras)
    
    
    return probability + log_prior(shift,galah_spectra.old_abundances)


def starting_test(solar_values,old_abundances):
    limit=5
    shift={'teff':solar_values[0],'vsini':solar_values[1],'logg':solar_values[2],'vrad':solar_values[3],'monh':solar_values[4],'vmac':solar_values[5],'vmic':solar_values[6]}
    if np.any([True for x in shift if shift[x]<0.0 and x!='monh']):
        return False
    
    if 'teff' in shift:
         if abs(shift['teff']-old_abundances['teff'])/old_abundances['e_teff']>limit:
                print("temperature is too small or large")
                return False
    if 'logg' in shift:
        if shift['logg']<1.0 or shift['logg']>5.0:
            print("logg is too small or big")
            return False  
        if abs(shift['logg']-old_abundances['logg'])/old_abundances['e_logg']>limit:
            print("logg is too small or big")
            return False  

    if 'monh' in shift:
        if abs(shift['monh']-old_abundances['fe_h'])/old_abundances['e_fe_h']>limit:
            print("monh is too small or big")
            return False   
    if 'vsini' in shift:
        if abs(shift['vsini']-old_abundances['vbroad'])/old_abundances['e_vbroad']>limit:
            print("monh is too small or big")
            return False       
    if 'vmic' in shift:
        if shift['vmic']<0 or shift['vmic']>20:
            return False 

            print("vmic is negative")
        if abs(shift['vmic']-old_abundances['vmic'])/old_abundances['e_vmic']>limit:
            
            return False 
    if 'vrad' in shift:
        if abs(shift['vrad']-old_abundances['rv_galah'])/old_abundances['e_rv_galah']>limit:
            return False

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
class spectrum_all:
    # bands=['IR']
    bands=['IR']

    def __init__(self,name,interpolation):
        name=str(name)
        bands=['IR']
        self.name=name
        self.interpolation=interpolation
        #Finds the prior solar parameters to enter for the star to set for initial values
        votable=parse('cross_hermes.xml')
        largeData=votable.get_first_table().to_table(use_names_over_ids=True)
        old_abundances=[x for x in largeData if x['sobject_id']==int(name)]
        old_abundances=vstack(old_abundances)
        self.old_abundances=old_abundances
        

        for count,x in enumerate(bands,1):
            if x=='IR':
                starting_fraction=2048/4096
                length_fraction=20/(4096*(1-starting_fraction))
            elif x=='Red':
                starting_fraction=0/4096
                length_fraction=20/(4096*(1-starting_fraction))
            elif x=='Green':
                starting_fraction=0/4096
                length_fraction=20/(4096*(1-starting_fraction))
            elif x=='Blue':
                starting_fraction=0/4096
                length_fraction=20/(4096*(1-starting_fraction))
            else:
                starting_fraction=2085/4096
                length_fraction=2/(4096*(1-starting_fraction))
            setattr(self,x,SME.SME_Structure())
            hermes=fits.open(name[0:6]+'/spectra/com/'+name+str(count)+'.fits')
            x0= float(hermes[1].header.get('CRVAL1'))
            x1=float( hermes[1].header.get('CRVAL1')+len(hermes[1].data)* hermes[1].header.get('CDELT1'))
            fstart= x0+(x1-x0)*starting_fraction
            
            new_start=int(starting_fraction*len(hermes[1].data))
            new_end=new_start+int(len(hermes[1].data)*(1-starting_fraction)*length_fraction)
            n_points=(new_end-new_start)*interpolation

            
            length=np.array([fstart+y*hermes[1].header.get('CDELT1')/interpolation for y in range(0,n_points)])
            
            fend=length[-1]
            #Increases the length of synthetic spectra so it over interpolates  
            rsetattr(self,x+'.wave',length)
            rsetattr(self,x+'.spec',[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])])
            rsetattr(self,x+'.uncs',[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[2].data[new_start:new_end])])
            hermes[0].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[0].data[new_start:new_end])
            hermes[1].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])
            hermes[2].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[2].data[new_start:new_end])
            hermes[7].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[7].data[new_start:new_end])
            hermes[1].header['CDELT1']/=interpolation
            
            rsetattr(self,x+'.hermes',hermes)
            
            rsetattr(self,x+'.wran',[fstart,fend])
            rsetattr(self,x+'.abund',Abund(float(old_abundances['fe_h']),'asplund2009' ,type='H=12'))
            print(x+'_HFS.lin')
            rsetattr(self,x+'.linelist',ValdFile('Galah_'+x+'_4.lin'))
            rsetattr(self,x+'.vmic',float(old_abundances['vmic']))
            rsetattr(self,x+'.vsini',float(old_abundances['vbroad']))
            rsetattr(self,x+'.vrad',float(old_abundances['rv_galah']))
            rsetattr(self,x+'.vmac',6.0)
            rsetattr(self,x+'.teff',float(old_abundances['teff']))
            rsetattr(self,x+'.logg',float(old_abundances['logg']))
            rsetattr(self,x+'.monh',float(old_abundances['fe_h']))
            rsetattr(self,x+'.h2broad',True)
            rsetattr(self,x+'.vrad_flag','whole')
            rsetattr(self,x+'.cscale_flag','none')
            rsetattr(self,x+'.cscale_type','mask')
            rsetattr(self,x+'.atmo.source',"marcs2012.sav")
            rsetattr(self,x+'.atmo.method','grid')
            rsetattr(self,x+'.atmo.geom','pp')
            runGetattr(self,x+'.nlte.set_nlte')("Fe", "marcs2012_Fe2016.grd")
            runGetattr(self,x+'.nlte.set_nlte')("Si", "marcs2012_Si2016.grd")
            runGetattr(self,x+'.nlte.set_nlte')("Ca", 'marcs2012p_t1.0_Ca.grd')
            runGetattr(self,x+'.nlte.set_nlte')("Na","marcs2012_Na2011.grd")
            runGetattr(self,x+'.nlte.set_nlte')("O",'marcs2012_O2015.grd')
            runGetattr(self,x+'.nlte.set_nlte')("Ti","marcs2012s_t2.0_Ti.grd")
            runGetattr(self,x+'.nlte.set_nlte')("Ba","marcs2012p_t1.0_Ba.grd")
            runGetattr(self,x+'.nlte.set_nlte')("Li","marcs2012_Li.grd")
            runGetattr(self,x+'.nlte.set_nlte')("Mg","marcs2012_Mg2016.grd")


            random_number=str(np.random.randint(100000000))
            synthetic_temp=copy.deepcopy(rgetattr(self,x+'.hermes'))
            if os.path.exists(self.name+x+random_number+'_synthetic.fits'):
                os.remove(self.name+x+random_number+'_synthetic.fits')
            synthetic_temp.writeto(self.name+x+random_number+'_synthetic.fits')

                    	
            syn=gtools.read(self.name+x+random_number+'_synthetic')
            syn.equalize_resolution()
            os.remove(self.name+x+random_number+'_synthetic.fits')
            rsetattr(self,x+'.spec',syn.f)
            
            
    def synthesize(self,shift,multi=False,colours=bands,give_back=False):
        if not (multi or give_back):
            for x in colours:       
                spectrum=copy.deepcopy(getattr(self,x))
                random_number=str(np.random.randint(100000000))
                for key in shift:
                    setattr(spectrum,key,shift[key])
                spectrum = synthesize_spectrum(spectrum)
                synthetic_temp=copy.deepcopy(rgetattr(self,x+'.hermes'))
                synthetic_temp[1].data=spectrum.synth[0]
                if os.path.exists(self.name+x+random_number+str(shift)+'_synthetic.fits'):
                    os.remove(self.name+x+random_number+str(shift)+'_synthetic.fits')
                synthetic_temp.writeto(self.name+x+random_number+str(shift)+'_synthetic.fits')      
                syn=gtools.read(self.name+x+random_number+str(shift)+'_synthetic')
                rsetattr(self,x+'.synth',syn.synth_resolution_degradation(np.transpose((rgetattr(self,x+'.wave')[0],syn.f)))[:,1])
                os.remove(self.name+x+random_number+str(shift)+'_synthetic.fits')

        elif give_back:
            returning_spectra=np.array(np.ones(len(colours)),dtype=object)
            for number,x in enumerate(colours):       
                spectrum=copy.deepcopy(getattr(self,x))
                random_number=str(np.random.randint(100000000))
                for key in shift:
                    setattr(spectrum,key,shift[key])
                spectrum = synthesize_spectrum(spectrum)
                synthetic_temp=copy.deepcopy(rgetattr(self,x+'.hermes'))
                synthetic_temp[1].data=spectrum.synth[0]
                if os.path.exists(self.name+x+random_number+str(shift)+'_synthetic.fits'):
                    os.remove(self.name+x+random_number+str(shift)+'_synthetic.fits')
                synthetic_temp.writeto(self.name+x+random_number+str(shift)+'_synthetic.fits')      
                syn=gtools.read(self.name+x+random_number+str(shift)+'_synthetic')
                returning_spectra[number]=syn.synth_resolution_degradation(np.transpose((rgetattr(self,x+'.wave')[0],syn.f)))[:,1]
                os.remove(self.name+x+random_number+str(shift)+'_synthetic.fits')
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
                # 	print("this is normalization",  original_line,"Line", x_line , "synthetic",  synth_line,"coeff",syth_coeff,"poly order", poly_order + 1, "spectra after",poly_synth/poly_original*original_line)                  
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
    def plot(self,colours=bands):
        for x in colours:
            plt.figure()
            x_line=np.linspace(0,len(runGetattr(self,x+'.synth')[0])-1,num=len(runGetattr(self,x+'.synth')[0]))
            if rgetattr(self,x+'.synth'):
                    # labels='synthetic  chi squared= '+str(self.log_fit([x]))
                    plt.plot(x_line, runGetattr(self,x+'.synth')[0], label='Synthetic')
                    
            plt.plot(x_line, runGetattr(self,x+'.spec')[0], label='Observed')
            plt.title(x+" Band")
            plt.legend(loc='best')
            plt.tight_layout()