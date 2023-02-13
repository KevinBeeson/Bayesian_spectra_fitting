#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 11:36:28 2022

@author: kevin
"""
from os.path import exists
import subprocess
from pathlib import Path
import csv
from functools import  partial
from astropy.io.votable import parse
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
logger=logging.getLogger('pysme')
logger.setLevel(0)

os.environ["OMP_NUM_THREADS"] = "1"

# def synthesize_multi(move,values,changed_values):
#     shift={values:move}
#     shift.update(changed_values)
#     return synthesize(shift)
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
def log_posterior(solar_values):
# pos=np.column_stack((teff_random,vsini_random,logg_random,v_rad_random,monh_random))
    limit=10
    shift={'teff':solar_values[0],'logg':solar_values[1],'monh':solar_values[2],'Fe':solar_values[3],'alpha':solar_values[4],'vrad':solar_values[5],'vsini':solar_values[6],'vmac':solar_values[7],'vmic':solar_values[8]}
    # print("this is SHIFT", shift)
    if np.any([True for x in shift if shift[x]<0.0 and x!='monh' and x!='Fe' and x!='alpha']):
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
    if np.any([True for x in shift if shift[x]<0.0 and  x!='monh' and x!='Fe' and x!='alpha']):
        return False
    
    if 'teff' in shift:
         if abs(shift['teff']-abundances['teff'])>1000:
                print("temperature is too small or large")
                return False
         elif shift['teff']<3500 or shift['teff']>10000:
            return False
    if 'logg' in shift:
        if shift['logg']<1.0 or shift['logg']>5.0:
            print("logg is too small or big")
            return False  
    if 'monh' in shift:
        if abs(shift['monh']-abundances['fe_h'])>0.5 or abs(shift['monh'])>3.0:
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
    solar_abundances={'H': 12.0, 'He': 10.93, 'Li': 1.05, 'Be': 1.38, 'B': 2.7, 'C': 8.43, 'N': 7.83, 'O': 8.69, 
                      'F': 4.56, 'Ne': 7.93, 'Na': 6.24, 'Mg': 7.6, 'Al': 6.45, 'Si': 7.51, 'P': 5.41, 'S': 7.12, 'Cl': 5.5, 'Ar': 6.4, 
                      'K': 5.03, 'Ca': 6.34, 'Sc': 3.15, 'Ti': 4.95, 'V': 3.93, 'Cr': 5.64, 'Mn': 5.43, 'Fe': 7.5, 'Co': 4.99, 'Ni': 6.22, 
                      'Cu': 4.19, 'Zn': 4.56, 'Ga': 3.04, 'Ge': 3.65, 'As': 2.3, 'Se': 3.34, 'Br': 2.54, 'Kr': 3.25, 'Rb': 2.52, 'Sr': 2.87, 
                      'Y': 2.21, 'Zr': 2.58, 'Nb': 1.46, 'Mo': 1.88, 'Tc': np.nan, 'Ru': 1.75, 'Rh': 0.91, 'Pd': 1.57, 'Ag': 0.94, 'Cd': 1.71, 
                      'In': 0.8, 'Sn': 2.04, 'Sb': 1.01, 'Te': 2.18, 'I': 1.55, 'Xe': 2.24, 'Cs': 1.08, 'Ba': 2.18, 'La': 1.1, 'Ce': 1.58, 
                      'Pr': 0.72, 'Nd': 1.42, 'Pm': np.nan, 'Sm': 0.96, 'Eu': 0.52, 'Gd': 1.07, 'Tb': 0.3, 'Dy': 1.1, 'Ho': 0.48, 'Er': 0.92, 
                      'Tm': 0.1, 'Yb': 0.84, 'Lu': 0.1, 'Hf': 0.85, 'Ta': -0.12, 'W': 0.85, 'Re': 0.26, 'Os': 1.4, 'Ir': 1.38, 'Pt': 1.62, 
                      'Au': 0.92, 'Hg': 1.17, 'Tl': 0.9, 'Pb': 1.75, 'Bi': 0.65, 'Po': np.nan, 'At': np.nan, 'Rn': np.nan, 'Fr': np.nan,
                      'Ra': np.nan, 'Ac': np.nan, 'Th': 0.02, 'Pa': np.nan, 'U': -0.54, 'Np': np.nan, 'Pu': np.nan, 'Am': np.nan, 
                      'Cm': np.nan, 'Bk': np.nan, 'Cf': np.nan, 'Es': np.nan}
    alpha_elements=['O', 'Ne', 'Mg', 'Si', 'S', 'Ar', 'Ca', 'Ti']
    all_elements=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 
                  'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 
                  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 
                  'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 
                  'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es']
    fe_abund=7.500
    priority=['monh','fe_h','H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 
                  'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 
                  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 
                  'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 
                  'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
                  'Bk', 'Cf', 'Es','alpha']
    # bands=['IR']
    def __init__(self,name,bands=['Blue','Green','Red','IR'],interpolation=10,limits=None):
        name=str(name)

        self.bands=bands
        self.name=name
        self.interpolation=interpolation
        
        votable=parse('cross_hermes.xml')
        largeData=votable.get_first_table().to_table(use_names_over_ids=True)

        old_abundances=[x for x in largeData if x['sobject_id']==int(name)]
        old_abundances=vstack(old_abundances)
        file_names={'Blue':1,'Green':2,'Red':3,'IR':4}
        for x,current_limit in zip(bands,limits):
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
            Path(name[0:6]+'/spectra/com/').mkdir(parents=True,exist_ok=True)
            if not exists(name[0:6]+'/spectra/com/'+name+str(file_names[x])+'.fits'):
                source='/media/storage/HERMES_REDUCED/dr6.0/'+name[0:6]+'/spectra/com/'+name+str(file_names[x])+'.fits'
                #source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.0/'+name[0:6]+'/spectra/com/'+name+str(file_names[x])+'.fits'

                destination=name[0:6]+'/spectra/com/'
                subprocess.run(["rsync",'-av',source,destination])


            hermes=fits.open(name[0:6]+'/spectra/com/'+name+str(file_names[x])+'.fits')
            if isinstance(current_limit,list):
                fstart=current_limit[0]
                fend=current_limit[1]
                length=np.array([fstart+y*hermes[1].header.get('CDELT1')/interpolation for y in range(int(np.ceil((fend-fstart)/hermes[1].header.get('CDELT1')*interpolation)))])
                wran=[fstart,fend]
            else:
                wran=[]
                length=[]
                for y in current_limit:
                    fstart_temp=y[0]
                    fend_temp=y[1]
                    length_temp=np.array([fstart_temp+y*hermes[1].header.get('CDELT1')/interpolation for y in range(int(np.ceil((fend_temp-fstart_temp)/hermes[1].header.get('CDELT1')*interpolation)))])
                    wran.append([fstart_temp,fend_temp])
                    length.append(length_temp)
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
            
            rsetattr(self,x+'.wran',wran)
            rsetattr(self,x+'.abund',Abund(float(old_abundances['fe_h']),'asplund2009' ,type='H=12'))
            print(x+'_HFS.lin')
            rsetattr(self,x+'.linelist',ValdFile('Galah_'+x+'_4.lin'))
            rsetattr(self,x+'.vmic',float(old_abundances['vmic']))
            rsetattr(self,x+'.vsini',float(old_abundances['vbroad']))
            rsetattr(self,x+'.vrad',0)
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
            
            
                        

            
    def solar_value_maker(self,shift,colour,keys=['teff','logg','monh','Fe','alpha','vrad','vsini','vmac','vmic']):
        solar_values=[]
        for x in keys:
            if x in shift:
                solar_values.append(shift[x])
            else:
                solar_values.append(rgetattr(self,colour+'.'+x))
        return solar_values        
    def synthesize(self,shift={},multi=False,colours=None,fe_abund=fe_abund,priority=priority,solar_abundaces=solar_abundances,alpha_elements=alpha_elements,give_back=False):
        if colours is None:
            colours=self.bands
        keys=shift.keys()
        [priority.append(x) for x in keys if not x in priority]
        shift_sorted={x:shift[x] for x in priority if x in shift}

        if not (multi or give_back):
            for x in colours:
                spectrum=copy.deepcopy(getattr(self,x))
                for key in shift_sorted:
                    if key in solar_abundaces:
                        spectrum.abund.update_pattern({key:solar_abundaces[key]+shift_sorted[key]-spectrum.monh})
                    elif key=='alpha':
                        alpha_fixed=[x for x in shift_sorted if x in alpha_elements ]
                        current_alpha=[x for x in alpha_elements if not x in alpha_fixed]
                        alpha_shift=(spectrum.abund.get_element('Fe')-fe_abund+shift_sorted[key])*len(alpha_elements)/len(current_alpha)
                        alpha_shift_dict={x:solar_abundaces[x]+alpha_shift-spectrum.monh for x in current_alpha}
                        spectrum.abund.update_pattern(alpha_shift_dict)
                    else:
                        setattr(spectrum,key,shift[key])
                spectrum = synthesize_spectrum(spectrum)
                rsetattr(self,x+'.synth',spectrum.synth)
        elif give_back:
            returning_spectra=np.array(np.ones(len(colours)),dtype=object)
            for number,x in enumerate(colours):
                spectrum=copy.deepcopy(getattr(self,x))
                for key in shift_sorted:
                    if key in solar_abundaces:
                        spectrum.abund.update_pattern({key:solar_abundaces[key]+shift_sorted[key]-spectrum.monh})
                    elif key=='alpha':
                        alpha_fixed=[x for x in shift_sorted if x in alpha_elements ]
                        current_alpha=[x for x in alpha_elements if not x in alpha_fixed]
                        alpha_shift=(spectrum.abund.get_element('Fe')-fe_abund+shift_sorted[key])*len(alpha_elements)/len(current_alpha)
                        alpha_shift_dict={x:solar_abundaces[x]+alpha_shift-spectrum.monh for x in current_alpha}
                        spectrum.abund.update_pattern(alpha_shift_dict)
                    else:
                        setattr(spectrum,key,shift[key])
                spectrum = synthesize_spectrum(spectrum)
                returning_spectra[number]=spectrum.synth
            return returning_spectra
        else:
            with Pool() as p:

                p.map(partial(getattr(self,'synthesize_multi'),shift=shift),['Blue','Green','Red','IR'])

            # pool = ThreadPool(4)
            # pool.map(partial(getattr(self,'synthesize_multi'),shift=shift),['Blue','Green','Red','IR'])
    def synthesize_multi(self,colours,shift):
        self.synthesize(shift=shift,multi=False,colours=[colours])
    def normalize(self,colours=None,data=None):
        if colours is None:
            colours=self.bands
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

    def log_fit(self,colours=None,synthetic_spectra=None):
        if colours is None:
            colours=self.bands
        probability=0
        for value,x in enumerate(colours):
            if not(np.array_equal(synthetic_spectra,None)):
                probability+= -np.sum((synthetic_spectra[value]-rgetattr(self,x+'.spec')[0])**2/(rgetattr(self,x+'.uncs')[0])**2)
            else:
                    probability+= -np.sum((rgetattr(self,x+'.synth')[0]-rgetattr(self,x+'.spec')[0])**2/(rgetattr(self,x+'.uncs')[0])**2)
        return probability
    def plot(self,colours=None):
        if colours is None:
            colours=self.bands
        for x in colours:
            plt.figure()
            x_line=np.linspace(0,len(runGetattr(self,x+'.synth')[0])-1,num=len(runGetattr(self,x+'.synth')[0]))
            if rgetattr(self,x+'.synth'):
                    # labels='synthetic  chi squared= '+str(self.log_fit([x]))
                    plt.plot(x_line, runGetattr(self,x+'.synth')[0], label='Synthetic')
                    
            # plt.plot(x_line, runGetattr(self,x+'.spec')[0], label='Observed')
            plt.title(x+" Band")
            plt.legend(loc='best')
            plt.tight_layout()
            
            
            

name='160401002101192'
filename = "160401002101192"

line_list = list(csv.reader(open('GALAH_DR3_line_info_table.csv', 'rt'), delimiter=','))

wanted_element='Li'

specific_list=[x for x in line_list if x[0]==wanted_element]
def out_of_limits(limit_array,test_limit):
    truth=[(y[0]>float(test_limit[8]) and y[1]>float(test_limit[9])) or (y[0]<float(test_limit[8]) and y[1]<float(test_limit[9])) for y in limit_array]
    return np.all(truth)
limits=[]
for x in specific_list:
    if not len(limits):
        limits.append([float(x[8]),float(x[9])])
    else:
        for y in limits:
            if y[0]>float(x[8]) and y[1]<float(x[9]):
                y[0]=float(x[8])
                break
            elif y[0]<float(x[8]) and y[1]>float(x[9]):
                y[1]=float(x[9])
                break
        if out_of_limits(limits,x):
            limits.append([float(x[8]),float(x[9])])

band_limits=[ [4718, 4903],[5649, 5873], [6481, 6739], [7590, 7890]]

bands=['Blue','Green','Red','IR']
bands_in=[]
bands_all=[]
for y in limits:
    for band,band_name in zip(band_limits,bands):
        if y[0]>band[0] and y[1]< band[1]:
            if not band_name in bands_in:
                bands_in.append(band_name)
            bands_all.append(band_name)
                
line_list_bands=np.array(np.ones(len(bands_in)),dtype=object)
bands=[x for x in bands if x in bands_in]
for band,limit in zip(bands_all,limits):
    if isinstance(line_list_bands[bands.index(band)],float):  
        line_list_bands[bands.index(band)]=limit
    else:
        line_list_bands[bands.index(band)]=np.vstack([line_list_bands[bands.index(band)],limit])
c=299792  #km/s

spectras=spectrum_all(name,bands=bands,limits=line_list_bands)
# spectras.synthesize({'logg':3.2,'alpha':1.2,'Si':0.1,'teff':5000})
global old_abundances


all_galah=fits.open('GALAH_DR3_main_200331.fits')
all_galah=all_galah[1].data
# print(all_galah[100]['teff_r'])

# spectras.synthesize({})
large_shift=[]
number_shift=[]
ncpu=10
sythesizing_length=ncpu*100
# spectras.synthesize({})
while len(large_shift)<sythesizing_length:
    random_num=np.random.randint(0,len(all_galah))
    old_abundances=all_galah[random_num]
    if old_abundances['teff'] and not np.isnan( old_abundances['teff']):
        teff_random=abs(np.random.normal(old_abundances['teff'],1000,1))
    else:
        teff_random=-100
    if old_abundances['logg'] and not np.isnan( old_abundances['logg']):
        logg_random=abs(np.random.normal(old_abundances['logg'],0.4,1))
    else:
        logg_random=-100
    v_rad_random= 0  
    if old_abundances['fe_h']  and not np.isnan( old_abundances['fe_h']):
        monh_random=np.random.normal(old_abundances['fe_h'],0.2,1)
    else:
        monh_random=-100
    if old_abundances['vmic'] and not np.isnan( old_abundances['vmic']):
        v_mac=abs(np.random.normal(old_abundances['vmic'],8,1))+np.random.randint(-10,10)
        v_mic=abs(np.random.normal(old_abundances['vmic'],8,1))
    else:
        v_mic=-100
    if old_abundances['vbroad'] and not np.isnan( old_abundances['vbroad']):
        vsini_random=abs(np.random.normal(old_abundances['vbroad'],1*4,1))
    else:
        v_mac=-100
        vsini_random=-100
    if old_abundances['alpha_fe'] and not np.isnan( old_abundances['alpha_fe']):
        alpha_random=np.random.normal(old_abundances['alpha_fe'],0.1,1)
    else:
        alpha_random=np.random.normal(0.148,0.317,1)
    if not(np.isnan(old_abundances[wanted_element+'_fe']) or np.isnan(old_abundances['e_'+wanted_element+'_fe']) or np.isnan(old_abundances['fe_h'])): 
        element_random=np.random.normal(old_abundances[wanted_element+'_fe']-old_abundances['fe_h'],old_abundances['e_'+wanted_element+'_fe'],1)
    else:
        element_random=np.random.normal(0,3,1)
    fe_random=monh_random+np.random.normal(0,0.1,1)


    pos_try=np.column_stack((teff_random,logg_random,monh_random,fe_random,alpha_random,v_rad_random,vsini_random,v_mac,v_mic,element_random))
    if starting_test(pos_try[0],old_abundances):
        pos_try=pos_try[0]
        shift={'teff':pos_try[0],'logg':pos_try[1],'monh':pos_try[2],'Fe':pos_try[3],'alpha':pos_try[4],'vrad':pos_try[5],'vsini':pos_try[6],'vmac':pos_try[7],'vmic':pos_try[8],wanted_element:pos_try[9]}
        number_shift.append(pos_try)
        large_shift.append(shift)
def multi_synthetic(shift):
    spectra=[]

    while len(spectra)==0:
        
        try:
            pos_try=[shift[x] for x in shift]
            spectra=spectras.synthesize(shift,give_back=True)
            
        except RuntimeError:
            
            print('probably an athmosphere errror ', shift)
            shift={}
            while len(shift)==0:
                random_num=np.random.randint(0,len(all_galah))
                old_abundances=all_galah[random_num]
                if old_abundances['teff'] and not np.isnan( old_abundances['teff']):
                    teff_random=abs(np.random.normal(old_abundances['teff'],1000,1))
                else:
                    teff_random=-100
                if old_abundances['logg'] and not np.isnan( old_abundances['logg']):
                    logg_random=abs(np.random.normal(old_abundances['logg'],0.4,1))
                else:
                    logg_random=-100
                v_rad_random= 0  
                if old_abundances['fe_h']  and not np.isnan( old_abundances['fe_h']):
                    monh_random=np.random.normal(old_abundances['fe_h'],0.2,1)
                else:
                    monh_random=-100
                if old_abundances['vmic'] and not np.isnan( old_abundances['vmic']):
                    v_mac=abs(np.random.normal(old_abundances['vmic'],8,1))+np.random.randint(-10,10)
                    v_mic=abs(np.random.normal(old_abundances['vmic'],8,1))
                else:
                    v_mic=-100
                if old_abundances['vbroad'] and not np.isnan( old_abundances['vbroad']):
                    vsini_random=abs(np.random.normal(old_abundances['vbroad'],1*4,1))
                else:
                    v_mac=-100
                    vsini_random=-100
                if old_abundances['alpha_fe'] and not np.isnan( old_abundances['alpha_fe']):
                    alpha_random=np.random.normal(old_abundances['alpha_fe'],0.1,1)
                else:
                    alpha_random=np.random.normal(0.148,0.317,1)
                fe_random=monh_random+np.random.normal(0,0.1,1)
                if not(np.isnan(old_abundances[wanted_element+'_fe']) or np.isnan(old_abundances['e_'+wanted_element+'_fe']) or np.isnan(old_abundances['fe_h'])): 
                    element_random=np.random.normal(old_abundances[wanted_element+'_fe']-old_abundances['fe_h'],old_abundances['e_'+wanted_element+'_fe'],1)
                else:
                    element_random=np.random.normal(0,3,1)
            
                pos_try=np.column_stack((teff_random,logg_random,monh_random,fe_random,alpha_random,v_rad_random,vsini_random,v_mac,v_mic,element_random))
                if starting_test(pos_try[0],old_abundances):
                    pos_try=pos_try[0]
                    shift={'teff':pos_try[0],'logg':pos_try[1],'monh':pos_try[2],'Fe':pos_try[3],'alpha':pos_try[4],'vrad':pos_try[5],'vsini':pos_try[6],'vmac':pos_try[7],'vmic':pos_try[8],wanted_element:pos_try[9]}
                    spectra=[]
    return spectra ,pos_try
# multi_synthetic({})

large_results=[]
steps=20
for x in range(int(sythesizing_length/ncpu/steps)):
    with Pool(processes=ncpu) as pool:
        results=pool.map(multi_synthetic,large_shift[x*ncpu*steps:(x+1)*ncpu*steps])

    if len(large_results)==0:
        large_results=results
    else:
        large_results=np.vstack((large_results,results))
    np.save('Training_data_no_degradation_fe_alpha'+wanted_element,large_results)
    print('Saving at '+str((x+1)*ncpu*5)+' iteration. It is '+ str((x+1)*ncpu*5/sythesizing_length))


        