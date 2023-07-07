#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 20:32:44 2021

@author: kevin
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 13:39:38 2021

@author: kevin
"""

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
from multiprocessing.dummy import Pool as ThreadPool 


os.environ["OMP_NUM_THREADS"] = "1"

def synthesize_multi(solar_values):
    shift={'teff':solar_values[0],'vsini':solar_values[1],'logg':solar_values[2],'vrad':solar_values[3],'monh':solar_values[4]}
    spectra=spectras.synthesize(shift,colours=['Blue'],give_back=True)
    normalized_spectra, normalized_uncertainty=spectras.normalize(data=spectra,colours=['Blue'])

    return spectra ,normalized_spectra,normalized_uncertainty
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
            error-=(shift[x]-photometric_values[x])**2/(photometric_values[x+'_e'])**2

    return error
def log_posterior(solar_values):
# pos=np.column_stack((teff_random,vsini_random,logg_random,v_rad_random,monh_random))

    shift={'teff':solar_values[0],'vsini':solar_values[1],'logg':solar_values[2],'vrad':solar_values[3],'monh':solar_values[4]}
    # print("this is SHIFT", shift)
    if np.any([True for x in shift if shift[x]<0.0 and x!='monh']):
        return -np.inf
    if 'teff' in shift:
         if abs(shift['teff']-old_abundances['teff'])>1000:
                # print("temperature is too small or large")
                return -np.inf
    if 'logg' in shift:
        if shift['logg']<1.0 or shift['logg']>5.0:
            # print("logg is too small or big")
            return -np.inf  
        if abs(shift['logg']-old_abundances['logg'])>1:
            # print("logg is too small or big")
            return -np.inf  

    if 'monh' in shift:
        if abs(shift['monh']-0.0)>5:
            print("monh is too small or big")
            return -np.inf   
    if 'vsini' in shift:
        if abs(shift['vsini']-old_abundances['vbroad'])>20:
            print("monh is too small or big")
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
    # print("Before normalized",synthetic_spectras, "Shift	",shift)
    # print("max value",np.max([np.max(abs(x)) for x in synthetic_spectras]),"shift",shift)
    normalized_spectra, normalized_uncertainty=spectras.normalize(data=synthetic_spectras)
    probability=spectras.log_fit(synthetic_spectra=synthetic_spectras,observed_spectra=normalized_spectra,uncs=normalized_uncertainty)
    # print("This is the probability", probability, "shift prob", log_prior(shift), "shift", shift)
    # print("This is the Uncs", normalized_uncertainty,"This is the Observed",normalized_spectra)
    # print("Spectra",synthetic_spectras)
    return probability
def starting_test(solar_values):
    shift={'teff':solar_values[0],'vsini':solar_values[1],'logg':solar_values[2],'vrad':solar_values[3],'monh':solar_values[4]}
    if np.any([True for x in shift if shift[x]<0.0 and x!='monh']):
        return False
    if 'teff' in shift:
         if abs(shift['teff']-old_abundances['teff'])>1000:
                print("temperature is too small or large")
                return False
    if 'logg' in shift:
        if shift['logg']<1.0 or shift['logg']>5.0:
            print("logg is too small or big")
            return False  
        if abs(shift['logg']-old_abundances['logg'])>1:
            print("logg is too small or big")
            return False  

    if 'monh' in shift:
        if abs(shift['monh']-0.0)>5:
            print("monh is too small or big")
            return False   
    if 'vsini' in shift:
        if abs(shift['vsini']-old_abundances['vbroad'])>20:
            print("monh is too small or big")
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
    bands=['Blue','Green','Red','IR']
    # bands=['Green','Red','IR']

    def __init__(self,name,interpolation):
        name=str(name)
        bands=['Blue','Green','Red','IR']
        self.name=name
        self.interpolation=interpolation
        
        votable=parse('cross_hermes.xml')
        largeData=votable.get_first_table().to_table(use_names_over_ids=True)

        old_abundances=[x for x in largeData if x['sobject_id']==int(name)]
        old_abundances=vstack(old_abundances)
        for count,x in enumerate(bands,1):
            if x=='IR':
                fraction=0.5
            else:
                fraction=0.0
            setattr(self,x,SME.SME_Structure())
            hermes=fits.open(name[0:6]+'/spectra/com/'+name+str(count)+'.fits')

            fstart= float(hermes[1].header.get('CRVAL1'))
            fend=float( (hermes[1].header.get('CRVAL1')+len(hermes[1].data)* hermes[1].header.get('CDELT1')))
            fstart+=fraction*(fend-fstart)
            
            original_length=len(hermes[1].data)
            new_start=int(fraction*len(hermes[1].data))
            
            length=np.array([fstart+y*hermes[1].header.get('CDELT1')/interpolation for y in range(0,len(hermes[1].data)-new_start)])
            
            #Increases the length of synthetic spectra so it over interpolates  
            rsetattr(self,x+'.wave',length)
            rsetattr(self,x+'.spec',[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,original_length-new_start)]),hermes[1].data[new_start:])])
            rsetattr(self,x+'.uncs',[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,original_length-new_start)]),hermes[2].data[new_start:])])
            hermes[0].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,original_length-new_start)]),hermes[0].data[new_start:])
            hermes[1].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,original_length-new_start)]),hermes[1].data[new_start:])
            hermes[2].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,original_length-new_start)]),hermes[2].data[new_start:])
            hermes[7].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,original_length-new_start)]),hermes[7].data[new_start:])
            hermes[1].header['CDELT1']/=interpolation
            
            rsetattr(self,x+'.hermes',hermes)
            
            rsetattr(self,x+'.wran',[fstart,fend])
            rsetattr(self,x+'.abund',Abund(float(old_abundances['fe_h']),'asplund2009' ,type='H=12'))
            rsetattr(self,x+'.linelist',ValdFile('Galah_'+x+'_2.lin'))
            rsetattr(self,x+'.vmic',float(old_abundances['vmic']))
            rsetattr(self,x+'.vsini',float(old_abundances['vbroad']))
            rsetattr(self,x+'.vrad',float(old_abundances['rv_galah']))
            rsetattr(self,x+'.vmac',6.0)
            rsetattr(self,x+'.teff',float(old_abundances['teff']))
            rsetattr(self,x+'.logg',float(old_abundances['logg']))
            rsetattr(self,x+'.monh',float(old_abundances['fe_h']))
            rsetattr(self,x+'.h2broad',True)
            rsetattr(self,x+'.vrad_flag','fix')
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
    def synthesize(self,shift,multi=False,colours=bands,give_back=False):
        if not (multi or give_back):
            for x in colours:       
                spectrum=copy.copy(getattr(self,x))
                random_number=str(np.random.randint(100000000))
                for key in shift:
                    setattr(spectrum,key,shift[key])
                spectrum = synthesize_spectrum(spectrum)
                synthetic_temp=copy.copy(rgetattr(self,x+'.hermes'))
            
                synthetic_temp[1].data=spectrum.synth[0]
                if os.path.exists(self.name+x+random_number+str(shift)+'_synthetic.fits'):
                    os.remove(self.name+x+random_number+str(shift)+'_synthetic.fits')
                synthetic_temp.writeto(self.name+x+random_number+str(shift)+'_synthetic.fits')

                        	
                syn=gtools.read(self.name+x+random_number+str(shift)+'_synthetic')
                syn.equalize_resolution()
                os.remove(self.name+x+random_number+str(shift)+'_synthetic.fits')
                rsetattr(self,x+'.synth',syn.f)
        elif give_back:
            returning_spectra=np.array(np.ones(len(colours)),dtype=object)
            for number,x in enumerate(colours):       
                spectrum=copy.copy(getattr(self,x))
                random_number=str(np.random.randint(100000000))
                for key in shift:
                    setattr(spectrum,key,shift[key])
                spectrum = synthesize_spectrum(spectrum)
               
                
                synthetic_temp=copy.copy(rgetattr(self,x+'.hermes'))
                # if x=='Blue':
                #     # print("in spectra",spectrum.synth[0], "shift", shift)             
                synthetic_temp[1].data=spectrum.synth[0]
                if os.path.exists(self.name+x+random_number+str(shift)+'_synthetic.fits'):
                    os.remove(self.name+x+random_number+str(shift)+'_synthetic.fits')
                synthetic_temp.writeto(self.name+x+random_number+str(shift)+'_synthetic.fits')
            
                syn=gtools.read(self.name+x+random_number+str(shift)+'_synthetic')
                # if x=='Blue':
                #     print("just before in spectra equalized",syn.f, "shift", shift)                 
                syn.equalize_resolution()

                # if x=='Blue':
                #     print("in spectra equalized",syn.f, "shift", shift)  
                returning_spectra[number]=syn.f
                os.remove(self.name+x+random_number+str(shift)+'_synthetic.fits')
            # print("just before giving the spectra back", returning_spectra,shift)
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
                    # print("Cant be normalized")
                    return np.nan 
                ynew=np.poly1d(syth_coeff)
                poly_synth=ynew(x_line)
                
                original_coeff=np.polyfit(x_line,original_line,poly_order)
                ynew=np.poly1d(original_coeff)
                poly_original=ynew(x_line)
                # if x=='Blue':
                # 	print("this is normalization",  original_line,"Line", x_line , "synthetic",  synth_line,"coeff",syth_coeff,"poly order", poly_order + 1, "spectra after",poly_synth/poly_original*original_line)                  
                returning_spectra[value]=poly_synth/poly_original*original_line
                
                returning_uncs[value]=poly_synth/poly_original*rgetattr(self,x+'.uncs')[0]
                
                
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
                
                uncs_normal=poly_synth/poly_original*rgetattr(self,x+'.uncs')[0]
                rsetattr(self,x+'.spec',synth_normal)
                rsetattr(self,x+'.uncs',uncs_normal)
                
            else:
                print(x+' Hasnt been synthesized')
        if not np.array_equal(data,None):
            return returning_spectra,returning_uncs

    def log_fit(self,colours=bands,synthetic_spectra=None,observed_spectra=None,uncs=None):
        probability=0
        for value,x in enumerate(colours):
            if not(np.array_equal(synthetic_spectra,None) and np.array_equal(observed_spectra,None) and np.array_equal(uncs,None)):
                probability+= -0.5*np.sum((observed_spectra[value]-synthetic_spectra[value])**2/(uncs[value]**2))
            else:
                    probability+= -0.5*np.sum((rgetattr(self,x+'.synth')[0]-rgetattr(self,x+'.spec')[0])**2/(rgetattr(self,x+'.uncs')[0])**2)
        return probability
    def plot(self,colours=bands):
        for x in colours:
            plt.figure()
            
            if rgetattr(self,x+'.synth'):
                    # labels='synthetic  chi squared= '+str(self.log_fit([x]))
                    plt.plot(runGetattr(self,x+'.wave')[0], runGetattr(self,x+'.synth')[0], label='Synthetic')
            plt.plot(runGetattr(self,x+'.wave')[0], runGetattr(self,x+'.spec')[0], label='Observed')
            plt.title(x+" Band")
            plt.legend(loc='best')
            plt.tight_layout()
            

name=180126001601066
filename = "100_teff_difference_Blue_gaia_start_shifted_no_prior.h5"


c=299792  #km/s
global old_abundances
# bands=['Blue','Green','Red','IR']

interpolation=1
sme=SME.SME_Structure()
votable=parse('cross_hermes.xml')
largeData=votable.get_first_table().to_table(use_names_over_ids=True)
old_abundances=[x for x in largeData if x['sobject_id']==int(name)]
old_abundances=vstack(old_abundances)
iron_H=float(old_abundances['fe_h'])
elem=sme.abund.elem_dict
names=old_abundances.colnames
names_fe=[x for x in names if  '_fe' in x and ('e_' not in x) and ('flag' not in x) and ('nr_' not in x) and ('alpha' not in x) and ('Ti2' not in x)]
names=[x[:-3] for x in names_fe ]
star_abundances=dict(zip(names, np.ones(len(elem))*0.01))

votable = parse("cross_gaia.xml")
photometric_values=votable.get_first_table().to_table(use_names_over_ids=True)
# photometric_values.rename_column('t_eff','teff')
# photometric_values.rename_column('log_gSigma','logg_e')
# photometric_values.rename_column('t_effSigma','teff_e')
# photometric_values.rename_column('log_g','logg')



photometric_values=[x for x in photometric_values if x['source_id']==old_abundances['source_id']][0]

for i in range(len(names)):
    if not np.isnan(old_abundances[names_fe[i]]):
        star_abundances[names[i]]=float(old_abundances[names_fe[i]])+12

hermes=fits.open(str(name)[0:6]+'/spectra/com/'+str(name)+'1.fits')
synthetic=fits.open(str(name)[0:6]+'/spectra/com/'+str(name)+'1.fits')




#Dopler
# fstart= float(hermes[1].header.get('CRVAL1')*(1+old_abundances['rv_galah']/c))
# fend=float( (hermes[1].header.get('CRVAL1')+len(hermes[1].data)* hermes[1].header.get('CDELT1'))*(1+old_abundances['rv_galah']/c))

fstart= float(hermes[1].header.get('CRVAL1'))
fend=float( (hermes[1].header.get('CRVAL1')+len(hermes[1].data)* hermes[1].header.get('CDELT1')))

length=np.array([fstart+y*hermes[1].header.get('CDELT1')/interpolation for y in range(0,len(hermes[1].data)*interpolation)])


#Increases the length of synthetic spectra so it over interpolates
synthetic[0].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,len(hermes[1].data))]),hermes[0].data)
synthetic[1].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,len(hermes[1].data))]),hermes[1].data)
synthetic[2].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,len(hermes[1].data))]),hermes[2].data)
synthetic[7].data=np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,len(hermes[1].data))]),hermes[7].data)


sme.wave= [length]
sme.spec=[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,len(hermes[1].data))]),hermes[1].data)]
sme.uncs=[np.interp(length,np.array([fstart+y*hermes[1].header.get('CDELT1') for y in range(0,len(hermes[1].data))]),hermes[2].data)]


#Loads the spectra lines
if str(name)[-1]=='1':
    sme.linelist=ValdFile('Galah_Blue_2.lin')
elif str(name)[-1]=='2':
    sme.linelist=ValdFile('Galah_Green_2.lin')
elif str(name)[-1]=='3':
    sme.linelist=ValdFile('Galah_Red_2.lin')
else:
    sme.linelist=ValdFile('Galah_IR_2.lin')

sme.wran =[fstart,fend]
sme.abund = Abund(float(old_abundances['fe_h']),'asplund2009' ,type='H=12')


sme.vmic=float(old_abundances['vmic'])
sme.vsini=float(old_abundances['vbroad'])
sme.vrad = float(old_abundances['rv_galah'])
sme.vmac = 6.0
sme.teff, sme.logg, sme.monh = old_abundances['teff'], old_abundances['logg'], old_abundances['fe_h']
sme.h2broad = True

sme.vrad_flag = "fix"
sme.cscale_flag = "none"
sme.cscale_type = "mask"


sme.atmo.source = "marcs2012.sav"
sme.atmo.method = "grid"
sme.atmo.geom = "PP"

sme.nlte.set_nlte("Fe", "marcs2012_Fe2016.grd")
sme.nlte.set_nlte("Si", "marcs2012_Si2016.grd")
sme.nlte.set_nlte("Ca", 'marcs2012p_t1.0_Ca.grd')
sme.nlte.set_nlte("Na","marcs2012_Na2011.grd")
sme.nlte.set_nlte("O",'marcs2012_O2015.grd')
sme.nlte.set_nlte("Ti","marcs2012s_t2.0_Ti.grd")
sme.nlte.set_nlte("Ba","marcs2012p_t1.0_Ba.grd")
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

n_dim=5
spectras=spectrum_all(name,2)

while True:
    teff_random=abs(np.random.normal(old_abundances['teff_gaia'],old_abundances['teff_e_gaia']*4,n_dim))
    vsini_random=abs(np.random.normal(old_abundances['vbroad'],old_abundances['e_vbroad']*4,n_dim))
    logg_random=abs(np.random.normal(old_abundances['logg_gaia'],old_abundances['logg_e_gaia']*4,n_dim))
    v_rad_random=abs(np.random.normal(old_abundances['rv_galah'],old_abundances['e_rv_galah']*4,n_dim))    
    monh_random=np.random.normal(old_abundances['fe_h'],old_abundances['e_fe_h']*4,n_dim)
    pos=np.column_stack((teff_random,vsini_random,logg_random,v_rad_random,monh_random))
    
    spectras=spectrum_all(name,2)
    
    with Pool() as p:
        n_dim=5
        teff_random=abs(np.random.normal(old_abundances['teff_gaia'],old_abundances['teff_e_gaia']*4,n_dim))
        vsini_random=abs(np.random.normal(old_abundances['vbroad'],old_abundances['e_vbroad']*4,n_dim))
        logg_random=abs(np.random.normal(old_abundances['logg_gaia'],old_abundances['logg_e_gaia']*4,n_dim))
        v_rad_random=abs(np.random.normal(old_abundances['rv_galah'],old_abundances['e_rv_galah']*4,n_dim))    
        monh_random=np.random.normal(-0.255,0.1*4,n_dim)
        pos=np.column_stack((teff_random,vsini_random,logg_random,v_rad_random,monh_random))
        
    
        synthesized_spectra=p.map(synthesize_multi,pos)
        # synthesized_spectra=[x[0] for x in synthesized_spectra]
    print("THIS IS THE MAX VALUE Spectra",np.max(synthesized_spectra), "normal spec",np.max(synthesized_spectra))






# probability_results=np.ones(nwalkers)
# teff_random=abs(np.random.normal(old_abundances['teff'],old_abundances['e_teff']*4,nwalkers))
# vsini_random=abs(np.random.normal(old_abundances['vbroad'],old_abundances['e_vbroad']*4,nwalkers))
# logg_random=abs(np.random.normal(old_abundances['logg'],old_abundances['e_logg']*4,nwalkers))
# v_rad_random=abs(np.random.normal(old_abundances['rv_galah'],old_abundances['e_rv_galah']*4,nwalkers))
# monh_random=np.random.normal(old_abundances['fe_h'],old_abundances['e_fe_h']*4,nwalkers)
# pos=np.column_stack((teff_random,vsini_random,logg_random,v_rad_random,monh_random))
# old_pos=np.column_stack((teff_random,vsini_random,logg_random,v_rad_random,monh_random))
# old_probability_results=np.ones(nwalkers)
# while np.all(probability_results):
#     print("No mess up yet")
#     teff_random=abs(np.random.normal(old_abundances['teff'],old_abundances['e_teff']*4,nwalkers))
#     vsini_random=abs(np.random.normal(old_abundances['vbroad'],old_abundances['e_vbroad']*4,nwalkers))
#     logg_random=abs(np.random.normal(old_abundances['logg'],old_abundances['e_logg']*4,nwalkers))
#     v_rad_random=abs(np.random.normal(old_abundances['rv_galah'],old_abundances['e_rv_galah']*4,nwalkers))
#     monh_random=np.random.normal(old_abundances['fe_h'],old_abundances['e_fe_h']*4,nwalkers)
#     pos=np.column_stack((teff_random,vsini_random,logg_random,v_rad_random,monh_random))
#     with Pool() as p:
#         probability_results=p.map(log_posterior, pos)
#     old_pos=np.vstack((old_pos,pos))
#     old_probability_results=np.hstack((old_probability_results,probability_results))
#     np.savetxt('problem_solving.csv',np.transpose(np.vstack((np.transpose(old_pos),old_probability_results))), delimiter=',')
# np.savetxt('Final_problem_solving.csv',np.transpose(np.vstack((old_pos,old_probability_results))), delimiter=',')

        
#Noise levels compared
# sigma_original=np.sqrt(np.sum((poly_original-original_line)**2)/len(original_line))
# noise=abs(np.random.normal(0,sigma_original,len(synth_line)))*-1

# x=np.linspace(0,len(noise),len(noise))
# gaussian=np.exp(-((x-len(noise)/2)/20)**2/2)/(20*np.sqrt(np.pi*2))
# noise=np.convolve(noise,gaussian,mode='same')
# noise=np.convolve(noise, np.ones(10)/10, mode='same')



# plt.figure()
# plt.plot(x,poly_original_new+noise,c='black',label='Observed')

# noise=abs(np.random.normal(0,sigma_original,len(synth_line)))*-1
# gaussian=np.exp(-((x-len(noise)/2)/20)**2/2)/(20*np.sqrt(np.pi*2))
# noise=np.convolve(noise,gaussian,mode='same')
# noise=np.convolve(noise, np.ones(10)/10, mode='same')


# plt.plot(x,poly_synth+noise,label='Observed')



# plt.figure()

# if len(np.shape(synthesized_spectra))==1:
#     chi_squared=np.sum((synth_normal-np.array(synthesized_spectra))**2/synth_normal)
#     labels='synthetic  chi squared= '+str(chi_squared)
#     fig=plt.figure()
#     plt.plot(x, synthesized_spectra, label=labels)

# else:
#     chi_squared=np.sum((synth_normal-np.array(synthesized_spectra))**2/synth_normal,axis=1)
#     labels=[str(x) for x in steps]
#     labels=[(x+'  chi squared= '+str(y)) for x,y in zip(labels,chi_squared)]
#     fig=plt.figure()
#     for y_arr, label in zip(synthesized_spectra, labels):
#         plt.plot(x, y_arr, label=label)
#         plt.legend(loc='best')
# plt.plot(x,synth_normal,c='black',label='Observed')

# plt.xlabel("wavelength/Angstrom")
# plt.ylabel("intensity")
# fig.set_size_inches(40,20)


# # plt.plot(x,np.transpose(synthesized_spectra),label=labels)
# # plt.plot(x,synth_normal,c='black',label='Observed')
# # plt.legend(loc='best')

# # sme.spec=synth_normal
# # fitparameters = ['teff']
# # sme = solve(sme, fitparameters)

# # plt.figure()
# # for y_arr, label in zip(synthesized_spectra, labels):
# #     plt.plot(x, y_arr, label=label)
# plt.tight_layout()
# plt.savefig(name+'_plot.pdf',dpi=10000,type='pdf')

# fig = plot_plotly.FinalPlot(sme)
# fig.save(filename=name+"_normalized.html")