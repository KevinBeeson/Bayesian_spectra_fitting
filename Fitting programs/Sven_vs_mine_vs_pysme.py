#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 14:16:24 2023

@author: kevin
"""
from functools import  partial
from astropy.io.votable import parse
import os.path
from astropy.table import Table,vstack,join

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

logger=logging.getLogger('pysme')
logger.setLevel(logging.ERROR)

os.environ["OMP_NUM_THREADS"] = "1"

# def synthesize_multi(move,values,changed_values):
#     shift={values:move}
#     shift.update(changed_values)
#     return synthesize(shift)
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
def shift_maker(solar_values):
    # shift={'vrad':solar_values}
    shift={'teff':solar_values[0],'logg':solar_values[1],'monh':solar_values[2],'alpha':solar_values[3],'vrad':solar_values[4],'vsini':solar_values[5],'vmac':solar_values[6],'vmic':solar_values[7]}

    return shift
def synthesize(shift):
    spectrum=copy.deepcopy(sme)
    for key in shift:
        setattr(spectrum,key,shift[key])
    spectrum = synthesize_spectrum(spectrum)
    synthetic_temp=copy.deepcopy(synthetic)

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
def log_posterior(solar_values):
# pos=np.column_stack((teff_random,vsini_random,logg_random,v_rad_random,monh_random))
    limit=10
    shift={'teff':solar_values[0],'logg':solar_values[1],'monh':solar_values[2],'alpha':solar_values[3],'vrad':solar_values[4],'vsini':solar_values[5],'vmac':solar_values[6],'vmic':solar_values[7]}
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
    
    shift={'teff':solar_values[0],'logg':solar_values[1],'monh':solar_values[2],'alpha':solar_values[3],'vrad':solar_values[4],'vsini':solar_values[5],'vmac':solar_values[6],'vmic':solar_values[7]}
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

class spectrum_all:
    bands=['Blue']

    alpha_elements_abundaces={'O':8.69,'Ne':7.93,'Mg':7.6,'Si':7.510,'S':7.120,'Ar':6.400,'Ca':6.340,'Ti':4.950}
    fe_abund=7.500

    # bands=['Blue','Green']
    def __init__(self,name,interpolation):
        limits={'Blue':[4705,4908],'Green':[5643,5879],'Red':[6474.0,6746.0],'IR':[7577.0,7894.0]}
        name=str(name)

        bands=['Blue']

        self.name=name
        self.interpolation=interpolation
        
        votable=parse('cross_hermes.xml')
        largeData=votable.get_first_table().to_table(use_names_over_ids=True)

        old_abundances=[x for x in largeData if x['sobject_id']==int(name)]
        old_abundances=vstack(old_abundances)
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
            hermes=fits.open(name[0:6]+'/spectra/com/'+name+str(file_names[x])+'.fits')
            x0= float(hermes[1].header.get('CRVAL1'))
            x1=float( hermes[1].header.get('CRVAL1')+len(hermes[1].data)* hermes[1].header.get('CDELT1'))
            fstart= x0+(x1-x0)*starting_fraction
            
            new_start=int(starting_fraction*len(hermes[1].data))
            new_end=new_start+int(len(hermes[1].data)*(1-starting_fraction)*length_fraction)
            n_points=(new_end-new_start)*interpolation

            length=np.linspace(limits[x][0],limits[x][1],num=int((limits[x][1]-limits[x][0])/0.004))
            # length=np.array([fstart+y*hermes[1].header.get('CDELT1')/interpolation for y in range(-80*interpolation,n_points+80*interpolation)])
            fstart=length[-1]
            fend=length[-1]
            #Increases the length of synthetic spectra so it over interpolates  
            rsetattr(self,x+'.wave',length)
            # rsetattr(self,x+'.spec',[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])])
            # rsetattr(self,x+'.uncs',[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[2].data[new_start:new_end])])
            # hermes[0].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[0].data[new_start:new_end])
            # hermes[1].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[1].data[new_start:new_end])
            # hermes[2].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[2].data[new_start:new_end])
            # hermes[7].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,new_end-new_start)]),hermes[7].data[new_start:new_end])
            # hermes[1].header['CDELT1']/=interpolation
            
            # rsetattr(self,x+'.hermes',hermes)
            rsetattr(self,x+'.hermes',hermes)
            rsetattr(self,x+'.wran',[length[0],length[-1]])
            rsetattr(self,x+'.abund',Abund(float(old_abundances['fe_h']),'asplund2009' ,type='H=12'))
            print(x+'_HFS.lin')
            rsetattr(self,x+'.linelist',ValdFile('line_list/'+x+'/combined_line_list_'+x+'.txt'))
            rsetattr(self,x+'.vmic',float(old_abundances['vmic']))
            rsetattr(self,x+'.vsini',float(old_abundances['vbroad']))
            rsetattr(self,x+'.vrad',0)
            rsetattr(self,x+'.vmac',0)
            rsetattr(self,x+'.teff',float(old_abundances['teff']))
            rsetattr(self,x+'.logg',float(old_abundances['logg']))
            rsetattr(self,x+'.monh',float(old_abundances['fe_h']))
            rsetattr(self,x+'.h2broad',True)
            rsetattr(self,x+'.vrad_flag','whole')
            rsetattr(self,x+'.cscale_flag','none')
            rsetattr(self,x+'.cscale_type','mask')
            rsetattr(self,x+'.atmo.source',"marcs2012p_t1.0.sav")
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
            line, wavelength = np.loadtxt('important_lines',usecols=(0,1),unpack=True,dtype=str, comments=';')
            wavelength=wavelength.astype(float)
            important_lines_temp=np.vstack([[elem_temp,wave_temp] for elem_temp,wave_temp in zip(line,wavelength) if wave_temp>runGetattr(self,x+'.wran')[0][0]-200 and wave_temp<runGetattr(self,x+'.wran')[0][1]+200])
            self.important_lines=important_lines_temp
            rsetattr(self,x+'.important_lines',important_lines_temp)
            
            random_number=str(np.random.randint(100000000))
                        

            
    def solar_value_maker(self,shift,colour,keys=['teff','logg','monh','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']):
        solar_values=[]
        for x in keys:
            if x in shift:
                solar_values.append(shift[x])
            else:
                solar_values.append(rgetattr(self,colour+'.'+x))
        return solar_values        
    def synthesize(self,shift,multi=False,colours=bands,fe_abund=fe_abund,alpha_elements_abundaces=alpha_elements_abundaces,give_back=False):
        elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']

        if not (multi or give_back):
            for x in colours:       
                spectrum=copy.deepcopy(getattr(self,x))
                for key in shift:
                    if key=='Fe':
                        spectrum.monh=shift[key]
                    elif key in elements:
                        
                        element_solar=Abund(0.0,'asplund2009')[key]
                        element_ratio=shift[key]
                        shift_dict={key:(element_ratio+element_solar)}
                        spectrum.abund.update_pattern(shift_dict)
                    else:
                        setattr(spectrum,key,shift[key])
                spectrum = synthesize_spectrum(spectrum)
                rsetattr(self,x+'.synth',spectrum.synth[0])
        elif give_back:
            returning_spectra=np.array(np.ones(len(colours)),dtype=object)
            for number,x in enumerate(colours):
                spectrum=copy.deepcopy(getattr(self,x))
                for key in shift:
                    if key=='Fe':
                        spectrum.monh=shift[key]
                    elif key in elements:
                        
                        element_solar=Abund(0.0,'asplund2009')[key]
                        element_ratio=shift[key]
                        shift_dict={key:(element_ratio+element_solar)}
                        spectrum.abund.update_pattern(shift_dict)
                    else:
                        setattr(spectrum,key,shift[key])
                spectrum = synthesize_spectrum(spectrum)
                returning_spectra[number]=spectrum.synth[0]
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
    def plot(self,colours=bands,lines='all'):
        for x in colours:
            plt.figure()
            # x_line=np.linspace(0,len(runGetattr(self,x+'.synth')[0])-1,num=len(runGetattr(self,x+'.synth')[0]))
            if rgetattr(self,x+'.synth'):
                    # labels='synthetic  chi squared= '+str(self.log_fit([x]))
                    plt.plot(rgetattr(self,x+'.wave')[0], runGetattr(self,x+'.synth')[0], label='Synthetic')
            if lines=='all':
                important_lines=rgetattr(self,x+'.important_lines')
                for individual_line in important_lines:
                    
                    plt.axvline(x=float(individual_line[1]))
                    plt.text(float(individual_line[1]),0.5,individual_line[0][:2],fontsize=20,ha='center',color='pink')
            else:
                important_lines=rgetattr(self,x+'.important_lines')
                for individual_line in important_lines:
                    if individual_line[0][:2] in lines:
                        
                        plt.axvline(x=float(individual_line[1]))
                        plt.text(float(individual_line[1]),0.5,individual_line[0][:2],fontsize=20,ha='center',color='pink')

            # plt.plot(x_line, runGetattr(self,x+'.spec')[0], label='Observed')
            plt.title(x+" Band")
            plt.legend(loc='best')
            plt.tight_layout()
            
            
            

name=170517002801129
filename = "170517002801129"


c=299792  #km/s



interpolation=1
sme=SME.SME_Structure()
votable=parse('cross_hermes.xml')
largeData=votable.get_first_table().to_table(use_names_over_ids=True)
global old_abundances
old_abundances=[x for x in largeData if x['sobject_id']==int(name)]
old_abundances=vstack(old_abundances)


all_galah=fits.open('dr6.0.fits')
all_galah=all_galah[1].data
# print(all_galah[100]['teff_r'])

# spectras=spectrum_all(name,10)
# spectras.synthesize({})
large_shift=[]
number_shift=[]
spectras=spectrum_all(170517002801129,10)

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
def leaky_relu(z,grad=False):
    '''
    This is the activation function used by default in all our neural networks.
    '''
    limits=(z > 0)
    leaky=z*limits + 0.01*z*np.invert(limits)
    if grad:
        return leaky,limits
    return leaky
tmp = np.load("NN_normalized_spectra_all_elements_2_Blue.npz")
w_array_0 = tmp["w_array_0"]
w_array_1 = tmp["w_array_1"]
w_array_2 = tmp["w_array_2"]
b_array_0 = tmp["b_array_0"]
b_array_1 = tmp["b_array_1"]
b_array_2 = tmp["b_array_2"]
x_min_kevin = tmp["x_min"]
x_max_kevin = tmp["x_max"]
tmp.close()
NN_coeffs_kevin= (w_array_0, w_array_1, w_array_2, b_array_0, b_array_1, b_array_2, x_min_kevin, x_max_kevin)
labels_with_limits=['teff','logg','fe_h','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
solar_values=[4500,4.5,-0.25,0.1,5]

solar_values=np.hstack((solar_values,np.zeros(len(elements))))
limits_kevin=[4705,4908]
wave_kevin=np.linspace(4705,4908,num=len(b_array_2))


tmp = np.load("test_sven_Blue.npz")
w_array_0 = tmp["w_array_0"]
w_array_1 = tmp["w_array_1"]
w_array_2 = tmp["w_array_2"]
b_array_0 = tmp["b_array_0"]
b_array_1 = tmp["b_array_1"]
b_array_2 = tmp["b_array_2"]
x_min_sven = tmp["x_min"]
x_max_sven = tmp["x_max"]
tmp.close()
NN_coeffs_sven = (w_array_0, w_array_1, w_array_2, b_array_0, b_array_1, b_array_2, x_min_sven, x_max_sven)
wave_sven=np.load('Blue_frequency_sven.npy')
solar_values=(x_max_sven+x_min_sven)/2
spectra_sven=payne_sythesize(solar_values,x_min_sven,x_max_sven,NN_coeffs_sven)
spectra_kevin=payne_sythesize(solar_values,x_min_kevin,x_max_kevin,NN_coeffs_kevin)

shift={y:x for (x,y) in zip(solar_values,labels_with_limits)}
spectras.synthesize(shift)
spectra_pysme=spectras.Blue.synth[0]
wave_pysme=spectras.Blue.wave[0]

spectra_kevin=np.interp(wave_sven,wave_kevin,spectra_kevin)
spectra_pysme=np.interp(wave_sven,wave_kevin,spectras_pysme)

plt.figure()
plt.plot(wave_sven,spectra_sven)
plt.plot(wave_kevin,spectra_kevin)
spectra_kevin=np.interp(wave_sven,wave_kevin,spectra_kevin)
plt.plot(wave_pysme,spectra_pysme)

total=np.column_stack((wave_sven,spectra_sven,spectra_kevin,spectra_pysme))
np.savetxt('comparison of spectras.txt',total,delimiter="\t",header='\t'.join(['wavelength','svens_spectra','kevins_payne_spectra','kevins_pysme_spectra']))

# i_dr4=fits.getdata('galah_dr4_allspec_220713.fits',1)
# i_dr4=Table(i_dr4)
# mask = i_dr4['flag_sp']==0
# i_dr4=i_dr4[mask]
# mask =i_dr4['e_fe_h']!=np.inf
# i_dr4=i_dr4[mask]

# large_data_GALAH_official=fits.open('gaia_galah_cross_values.fits')
# large_data_GALAH_official=Table(large_data_GALAH_official[1].data)
# large_data_GALAH_official['sobject_id']=np.array(large_data_GALAH_official['sobject_id'],dtype='int64')
# dr3_data=fits.getdata('GALAH_DR3_main_200331.fits',1)
# dr3_data=Table(dr3_data)

# cross_data=join(large_data_GALAH_official,dr3_data,keys_left='sobject_id',keys_right='sobject_id')

# cross_data.sort('logg_2')

# plt.figure()
# plt.scatter(cross_data['teff_2'][::18],cross_data['logg_2'][::18])
# parameters_names=['teff_2','logg_2','fe_h_2','vmic_2','vbroad_2']
# elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
# elem_mean=[]
# elem_sig=[]
# for elem in elements:
#     elem_mean.append(np.nanmean(i_dr4[elem.lower()+'_fe']))
#     elem_sig.append(np.sqrt(np.nanvar(i_dr4[elem.lower()+'_fe'])))
# fe_h_mean=np.nanmean(i_dr4['fe_h'])
# fe_h_sig=np.nanmean(i_dr4['e_fe_h'])
# elem_max=[ 3.39219427e+00, 1.50000000e+00, 1.50000000e+00,
# 2.00000000e+00,0.7e+00,0.7e+00,1.2e+00,
# 1.0e+00,1.0e+00,1.0e+00,0.7e+00,
# 0.9e+0,0.7e+00,0.7e+00,1.0e+00,
# 0.8e+0,0.7e+00,0.7e+00,1.0e+00,
# 1.50000000e+00,1.0e+00, 1.00000000e+00,1.0e+00,
# 0.7e+0,1.0e+00, 1.20000000e+00, 1.30000000e+00,
# 1.30000000e+00, 1.20000000e+00,1.0e+00,1.0e+00]
# elem_min=[-1.5e0,-1.0e+00,-1.5e+00,
# -1.50000000e+00,-0.7e+00, -1.00000000e-0,-0.7e+00,
# -5.00000000e-01,-1.0e+00,-0.7e+00,-0.7e+00,
# -5.00000000e-01,-0.7e+00,-0.7e+00,-0.7e+00,
# -0.7e+0,-0.7e+00,-0.7e+00,-1.0e+00,
# -1.80000000e+00,-1.0e+00,-1.0e+00,-1.0e+00,
# -1.0e+0, -2.00000000e+00,-1.0e+00,-0.7e+00,
# -1.0e+0,-1.0e+00,-1.0e+00,-1.0e+00]
# parameters=[100,0.5,0.5,5,5]
# hard_limits_low=[3000,0.5,-1.5,0,0.1]
# hard_limits_high=[8500,5,1.0,5,150]
# np.random.seed(0)
# cross_short=cross_data[::13]
# param_to_compute=[]
# # for star in cross_short:
# #     star=cross_short[0]
# #     param_temp=[star['teff_2'],star['logg_2'],np.random.normal(fe_h_mean,fe_h_sig*3,1)[0],star['vmic_2'],star['vbroad_2']]
# #     param_temp=np.hstack((param_temp,np.zeros(31)))
# #     if len (param_to_compute)==0:
# #         param_to_compute=param_temp
# #     param_to_compute=np.vstack((param_to_compute,param_temp))
# #     for value,shift in enumerate(parameters):
# #         to_shift=copy.copy(param_temp)
# #         to_shift[value]+=shift
# #         param_to_compute=np.vstack((param_to_compute,to_shift))
# #         to_shift=copy.copy(param_temp)
# #         to_shift[value]-=shift
# #         if value==1:
# #             if to_shift[value]<0.5:
# #                 to_shift[value]=0.5
# #         param_to_compute=np.vstack((param_to_compute,to_shift))
        
# #     for value,(shift_max,shift_min) in enumerate(zip(elem_max,elem_min),len(parameters)):
# #         to_shift=copy.copy(param_temp)
# #         to_shift[value]=shift_max
# #         param_to_compute=np.vstack((param_to_compute,to_shift))
# #         to_shift=copy.copy(param_temp)
# #         to_shift[value]=shift_min        
# #         param_to_compute=np.vstack((param_to_compute,to_shift))
# # fig, ax = plt.subplots(nrows=6, ncols=6,figsize=(16.0,9.0) )
# # x=4
# # for row in ax:
# #     for col in row:
# #         x+=1
# #         if x<36:
# #             # seaborn.kdeplot(param_to_compute[:,2],param_to_compute[:,x],ax=col)
# #             col.scatter(param_to_compute[:,2],param_to_compute[:,x],s=0.1,c='red')
# #             # col.set_xlabel('Fe/H')
# #             col.set_ylabel(elements[x-5]+'/fe')
# #             col.set_title(elements[x-5],y=0.7,x=0.1)
# #         else:
# #             col.set_axis_off()
# # plt.tight_layout()
# # fig.savefig('abundances distribution with limits')
# logg_low=1.0
# teff_low=4467
# logg_high=3.2
# teff_high=6000
# gradient_1=(logg_high-logg_low)/(teff_high-teff_low)
# zero_point_1=logg_high-teff_high*gradient_1

# logg_low=1.14
# logg_high=3.92
# teff_low=3574
# teff_high=4952
# gradient_2=(logg_high-logg_low)/(teff_high-teff_low)
# zero_point_2=logg_high-teff_high*gradient_2
# print('compiling_targets')
# param_to_compute=[]     
# while len(param_to_compute)<1.5e4:
#     for star in cross_data:
#         param_temp=[]
#         for names,lower,upper,spread in zip(parameters_names,hard_limits_low,hard_limits_high,parameters):
#             temp=-np.inf
#             if 'vmic_2'==names:
#                 sigma=spread
#             else:
#                 sigma=star['e_'+names]
#             mean=star[names]

#             if np.isnan(mean):
#                 mean=(upper+lower)/2
#             if np.isnan(sigma):
#                 sigma=(upper-lower)/2
#             while temp<lower or temp>upper:
#                 temp=np.random.normal(mean,sigma,1)
    
#             param_temp.append(temp)
#         for names,lower,upper,mean_zero,spread_zero in zip(elements,elem_min,elem_max,elem_mean,elem_sig):
#             temp=-np.inf
#             while temp<lower or temp>upper:
#                 if names=='N':
#                     mean=0.0
#                     sigma=0.5
#                 else:
#                     mean=star[names+'_fe']
#                     sigma=star['e_'+names+'_fe']
#                 if np.isnan(mean):
#                     mean=mean_zero
#                 if mean<lower:
#                     mean=lower
#                 elif mean>upper:
#                     mean=upper
#                 if np.isnan(sigma):
#                     sigma=spread_zero
#                 temp=-np.inf
#                 while temp<lower or temp>upper:
#                     temp=np.random.normal(mean,sigma,1)
#             param_temp.append(temp)
#         while( (param_temp[0]>6000 and param_temp[1]<3.06) or
#               (param_temp[0]<6000 and param_temp[1]<3.06 and inside_area(param_temp[0],param_temp[1],gradient_1,zero_point_1))or
#                 param_temp[0]<4952 and param_temp[1]<3.92 and inside_area(param_temp[0],param_temp[1],gradient_2,zero_point_2,upper=False)or
#                 param_temp[1]>5.21 or param_temp[0]>8300 or 
#                 (param_temp[0]>6500 and param_temp[1]>4.7) or
#                 param_temp[0]>6000 or param_temp[0]<5000 or 
#                 param_temp[1]<4.0 or param_temp[1]>5.0):
#             random_number=np.random.randint(0,len(cross_data))
#             star_temp=cross_data[random_number]
#             mean=star_temp['teff_2']
#             sigma=star_temp['e_teff_2']*4
            
#             if np.isnan(mean):
#                 mean=-1000
#             if np.isnan(sigma):
#                 sigma=0

#             param_temp[0]=np.random.normal(mean,sigma,1)
#             mean=star_temp['logg_2']
#             sigma=star_temp['e_logg_2']*4
#             if np.isnan(mean):
#                 mean=-1000
#             if np.isnan(sigma):
#                 sigma=0
#             param_temp[1]=np.random.normal(mean,sigma,1)
#             param_temp=np.hstack(param_temp)
#         # while param_temp[0]
#         if len (param_to_compute)==0:

#             param_to_compute=np.hstack(param_temp)
#         param_to_compute=np.vstack((param_to_compute,np.hstack(param_temp)))
# param_to_compute=np.array(param_to_compute)

# parameters_all_names=['teff','logg','monh','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
# shift_all=[]
# for x in param_to_compute:
#       shift_temp={y:para for y,para in zip(parameters_all_names,x)}
#       shift_all.append(shift_temp)
     

# def multi_synthetic(shift):
#     spectra=[]
#     spectras=spectrum_all(170517002801129,10)
#     while len(spectra)==0:
    
#         try:
#             pos_try=spectras.solar_value_maker(shift,colour='Blue')
#             spectra=spectras.synthesize(shift,give_back=True)
#         except:
#             print('failed',shift)
#             shift={}
#             random_number=np.random.randint(0,len(cross_data))
#             star=cross_data[random_number]
#             param_temp=[]
#             for names,lower,upper,spread in zip(parameters_names,hard_limits_low,hard_limits_high,parameters):
#                 temp=-np.inf
#                 if 'vmic_2'==names:
#                     sigma=spread
#                 else:
#                     sigma=star['e_'+names]
#                 mean=star[names]

#                 if np.isnan(mean):
#                     mean=(upper+lower)/2
#                 if np.isnan(sigma):
#                     sigma=(upper-lower)/2
#                 while temp<lower or temp>upper:
#                     temp=np.random.normal(mean,sigma,1)
        
#                 param_temp.append(temp)
#             for names,lower,upper,mean_zero,spread_zero in zip(elements,elem_min,elem_max,elem_mean,elem_sig):
#                 temp=-np.inf
#                 while temp<lower or temp>upper:
#                     if names=='N':
#                         mean=0.0
#                         sigma=0.5
#                     else:
#                         mean=star[names+'_fe']
#                         sigma=star['e_'+names+'_fe']
#                     if np.isnan(mean):
#                         mean=mean_zero
#                     if mean<lower:
#                         mean=lower
#                     elif mean>upper:
#                         mean=upper
#                     if np.isnan(sigma):
#                         sigma=spread_zero
#                     temp=-np.inf
#                     while temp<lower or temp>upper:
#                         temp=np.random.normal(mean,sigma,1)
#                 param_temp.append(temp)
#             while( (param_temp[0]>6000 and param_temp[1]<3.06) or
#                   (param_temp[0]<6000 and param_temp[1]<3.06 and inside_area(param_temp[0],param_temp[1],gradient_1,zero_point_1))or
#                    param_temp[0]<4952 and param_temp[1]<3.92 and inside_area(param_temp[0],param_temp[1],gradient_2,zero_point_2,upper=False)or
#                    param_temp[1]>5.21 or param_temp[0]>8300 or 
#                    (param_temp[0]>6500 and param_temp[1]>4.7) ):
#                 random_number=np.random.randint(0,len(cross_data))
#                 star_temp=cross_data[random_number]
#                 mean=star_temp['teff_2']
#                 sigma=star_temp['e_teff_2']*4
                
#                 if np.isnan(mean):
#                     mean=-1000
#                 if np.isnan(sigma):
#                     sigma=0

#                 param_temp[0]=np.random.normal(mean,sigma,1)
#                 mean=star_temp['logg_2']
#                 sigma=star_temp['e_logg_2']*4
#                 if np.isnan(mean):
#                     mean=-1000
#                 if np.isnan(sigma):
#                     sigma=0
#                 param_temp[1]=np.random.normal(mean,sigma,1)
#                 param_temp=np.hstack(param_temp)
#             shift={y:para for y,para in zip(parameters_all_names,param_temp)}
            
#             pos_try=spectras.solar_value_maker(shift,colour='Blue')
#             spectra=spectras.synthesize(shift,give_back=True)


#     return spectra ,pos_try
# # multi_synthetic(7678157831)
# # # multi_synthetic({})
# # spectras.synthesize({})
# large_results=np.array([],dtype=object)
# steps=1
# # spectras.synthesize(large_shift[0])
# # spectras.synthesize(large_shift[0],give_back=True)
# # spec,pos1=multi_synthetic(large_shift[0])
# ncpu=2
# sythesizing_length=len(shift_all)

# spectras.synthesize(shift_all[0])
# for x in range(int(sythesizing_length/ncpu/steps)):
#     pool=Pool(processes=ncpu)
#     results=pool.map(multi_synthetic,shift_all[x*ncpu*steps:(x+1)*ncpu*steps])
#     pool.close()
#     if len(large_results)==0:
#         large_results=np.array(results,dtype=object)
#     else:
#         large_results=np.vstack((large_results,results))
#     np.save('Training_data_small_Blue_2',large_results)
#     print('Saving at '+str((x+1)*ncpu*steps)+' iteration. It is '+ str((x+1)*ncpu*steps/sythesizing_length))


