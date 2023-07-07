#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates a cross file for dr6 fits and all_gaia
Created on Thu Jan 20 16:54:25 2022

@author: kevin
"""
from astropy.io import fits
import astropy
from astropy.io.votable import parse
from astropy.table import vstack
import numpy as np
from astropy.table import Table
import copy 
#reduction files
large_data_GALAH_reduction=fits.open('dr6.0.fits')
large_data_GALAH_reduction=Table(large_data_GALAH_reduction[1].data)
large_data_GALAH_reduction=large_data_GALAH_reduction
#Galah fits file
# large_data_GALAH_official=fits.open('GALAH_DR3_main_allstar_v2.fits')
# large_data_GALAH_official_unchanged=large_data_GALAH_official[1].data
# large_data_GALAH_official=Table(large_data_GALAH_official[1].data)



# keys=['logg_r','fe_h_r','vmic_r','vbroad_r']

#All gaia stars
votable=parse('all_gaia.xml')
gaia=votable.get_first_table().to_table(use_names_over_ids=True)

gaia_all=votable.get_first_table().to_table(use_names_over_ids=True)


#gaia dr2 to gaia dr3
votable = parse("dr2_to_dr3_rv.vot")
dr2_dr3=Table(votable.get_first_table().to_table(use_names_over_ids=True))
dr2_dr3.rename_column('col1', 'dr2_source_id')
cross=np.array([],dtype=str)
cross_rv=[]
cross_rv_error=[]
for count,x in enumerate(gaia):
    where=np.where(dr2_dr3['dr2_source_id']==x['source_id'])
    cross=np.append(cross,str(dr2_dr3[where[0][0]]['dr3_source_id']))
    cross_rv.append(dr2_dr3[where[0][0]]['dr2_radial_velocity'])
    cross_rv_error.append(dr2_dr3[where[0][0]]['dr2_radial_velocity_error'])
gaia['dr3_source_id']=cross
gaia_all['dr3_source_id']=cross


#only having galah members with gaia clusters
cross_galah_reduction=[x for x in large_data_GALAH_reduction if x['gaia_id'] in cross]
cross_gaia=[x for x in gaia if x['dr3_source_id'] in large_data_GALAH_reduction['gaia_id']]
cross_galah_reduction=vstack(cross_galah_reduction)
cross_gaia=vstack(cross_gaia)
cross_gaia_all=np.copy(cross_gaia)

#only having galah members with gaia clusters
# cross_galah_reduction=[x for x in large_data_GALAH_reduction if x['gaia_id'] in cross]
# dr3_source_id_str=[str(x) for x in large_data_GALAH_official['dr3_source_id']]
# cross_gaia_official=[x for x in gaia if x['dr3_source_id'] in dr3_source_id_str]
# cross_gaia_official=np.vstack(cross_gaia_official)


# delete=0
# for value,x in enumerate(cross_values_galah):
    
# cross_values_galah=[x for x in large_data_GALAH_official if str(x['dr3_source_id']) in cross]
# cross_values_galah=vstack(cross_values_galah)

# cross tables galah
galah_reduction_to_gaia_values=[]
added_gaia=[]
galah_official_keys=['gaia_id','sobject_id','teff_sven','e_teff_sven','logg_sven','e_logg_sven','fe_h','e_fe_h','vmic','vbroad','e_vbroad',
                     'alpha_fe','e_alpha_fe','rv_galah','e_rv_galah']

reduction_keys=['gaia_id','sobject_id','teff_reduction','logg_reduction','fe_h','vmic_reduction','vbroad_reduction',
            'reduction_flags_reduction','rv_Blue_r','rv_Green_r','rv_Red_r','rv_IR_r','e_rv_Blue_r','e_rv_Green_r','e_rv_Red_r',
            'e_rv_IR_r','rv_r','e_rv_r','alpha_fe_r']
# gets galah_reduction 
for count,x in enumerate(cross_galah_reduction):
    galah_reduction_to_gaia_values.append([x['gaia_id'],x['sobject_id'],x['teff_r'],x['logg_r'],x['fe_h_r'],x['vmic_r'],x['vbroad_r'],
                                           x['reduction_flags'],x['rv'][0],x['rv'][1],x['rv'][2],x['rv'][3],x['e_rv'][0],x['e_rv'][1],x['e_rv'][2],
                                           x['e_rv'][3],x['rv_com'],x['e_rv_com'],x['alpha_fe_r']])
# galah_values_to_gaia_values=[]
# for x in (cross_values_galah):
#     galah_values_to_gaia_values.append([x['dr3_source_id'],x['sobject_id'],x['teff'],x['e_teff'],x['logg'],x['e_logg'],x['fe_h'],x['e_fe_h'],x['vmic'],x['vbroad'],x['e_vbroad'],x['alpha_fe'],x['e_alpha_fe'],x['rv_galah'],x['e_rv_galah']])

#make cross file with reduction data
all_data=[]
gaia_keys=gaia.colnames
gaia_keys.remove('dr3_source_id')
gaia_keys.remove('phot_variable_flag')
for x in galah_reduction_to_gaia_values:
    where=np.where(cross_gaia['dr3_source_id']==x[0])
    temp_gala=cross_gaia[where]
    a=[temp_gala[y] for y in gaia_keys]
    temp_galah=copy.copy(x)
    for y in a:
        temp_galah.append(y[0])
    if len(temp_galah)!=len(gaia_keys)+len(reduction_keys):
        print('oh_no')
    all_data.append(temp_galah)
all_data=np.array(all_data)
final=Table({'eDR3_source_id':all_data[:,0].astype(str)})
bands=['Blue','Green','Red','IR']
for count,x in enumerate(gaia_keys):
    if  x!='sobject_id' and x!='gaia_id' and x!='dr2_source_id' and x!='dr3_source_id' and x!='source_id':
        final[x]=all_data[:,count+19].astype(np.float64)
    elif  x=='sobject_id' or x=='gaia_id' or x=='dr2_source_id'or x=='dr3_source_id':
        final[x]=all_data[:,count+19].astype(str)
    else:
        final[x]=all_data[:,count+19]
for count,x in enumerate(reduction_keys):
    if x!='sobject_id' and x!='gaia_id' and x!='dr2_source_id' and x!='dr3_source_id':
        final[x]=all_data[:,count].astype(np.float64)
    elif x=='sobject_id' or x=='gaia_id' or x=='dr2_source_id'or x=='dr3_source_id' :
        final[x]=all_data[:,count].astype(str)
    else:
        final[x]=all_data[:,count]
for x in final.colnames:
    print(final[x].dtype  )  
    print(x)
final.rename_column('t_eff','teff')
final.rename_column('t_effSigma','e_teff')

final.rename_column('log_g','logg')
final.rename_column('log_gSigma','e_logg')


final.write('photometric_reduction_dr6_cross.fits',overwrite=True)

#make cross file with values data
# all_data=[]
# gaia_keys=gaia.colnames
# # gaia_keys.remove('dr3_source_id')
# gaia_keys.remove('phot_variable_flag')
# for x in galah_values_to_gaia_values:
#     where=np.where(cross_gaia_official['dr3_source_id']==str(x[0]))
#     temp_gala=cross_gaia_official[where]
#     a=[temp_gala[y] for y in gaia_keys]
    
#     a=np.vstack(a)
#     temp_galah=copy.copy(x)
#     for y in a:
#         temp_galah.append(y[0])
#     if len(temp_galah)!=len(gaia_keys)+len(galah_official_keys):
#         print('oh_no')
#     all_data.append(temp_galah)
# all_data=np.array(all_data)
# final=Table({'eDR3_source_id':all_data[:,0].astype(str)})
# bands=['Blue','Green','Red','IR']
# for count,x in enumerate(gaia_keys):
#     if all_data[:,count+15].dtype==object and x!='sobject_id' and x!='gaia_id' and x!='dr2_source_id' and x!='dr3_source_id':
#         if len([q for q in all_data[:,count+15] if q.size==0])!=0:
#             temp_data=[]
#             for y in all_data[:,count+15]:
#                 if len(y)==0:
#                     temp_data.append(np.nan)
#                 else:
#                     temp_data.append(y)
#             temp_data=np.vstack(temp_data)
#             final[x]=temp_data
  
#         else:
#             final[x]=all_data[:,count+15].astype(np.float64)
#     elif  x!='sobject_id' or x!='gaia_id' or x!='dr2_source_id'or x!='dr3_source_id':
#         final[x]=all_data[:,count+15].astype(str)
#     else:
#         final[x]=all_data[:,count+15]
# for count,x in enumerate(galah_official_keys):
#     if all_data[:,count].dtype==object and x!='sobject_id' and x!='gaia_id' and x!='dr2_source_id' and x!='dr3_source_id':
#         if len([q for q in all_data[:,count] if q.size==0])!=0:
#             temp_data=[]
#             for y in all_data[:,count]:
#                 if len(y)==0:
#                     temp_data.append(np.nan)
#                 else:
#                     temp_data.append(y)
#             temp_data=np.vstack(temp_data)
#             final[x]=temp_data
#         else:
#             final[x]=all_data[:,count].astype(np.float64)
        
#     elif x=='sobject_id' or x=='gaia_id' or x=='dr2_source_id'or x=='dr3_source_id' :
#         final[x]=all_data[:,count].astype(str)
#     else:
#         final[x]=all_data[:,count]
# for x in final.colnames:
#     print(final[x].dtype  )  
#     print(x)
# final.rename_column('t_eff','teff')
# final.rename_column('t_effSigma','e_teff')

# final.rename_column('log_g','logg')
# final.rename_column('log_gSigma','e_logg')


# final.write('gaia_galah_cross_values.fits',overwrite=True)