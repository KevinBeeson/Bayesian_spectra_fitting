#this will load NGC_2682_no_normalization_.xml get calculate its spread factor and for all stars with a spread factor above 0.5 it will save its sobject_id
from astropy.io.votable import parse,from_table,writeto
from astropy.table import Table,vstack,join
import numpy as np
name='NGC_2682'
votable = parse(name+"_no_normalization_.xml")
photometric_data=votable.get_first_table().to_table(use_names_over_ids=True)
sig_teff=np.mean(photometric_data['e_teff_raw'],axis=1)
sig_logg=np.mean(photometric_data['e_logg_raw'],axis=1)
spread_teff=np.std(photometric_data['teff_raw'],axis=1)
spread_logg=np.std(photometric_data['logg_raw'],axis=1)
spread_factor=spread_teff*spread_logg/(sig_teff*sig_logg*len(photometric_data[0]['teff_raw']))
mask=spread_factor>0.5
photometric_data_to_save=photometric_data[mask]
sobject_id=photometric_data_to_save['sobject_id']
np.savetxt(name+'_targets_large_spread_factor.txt',sobject_id.astype(str),fmt='%s')

#check how many stars have a spread factor above 0.5 in new data
votable = parse(name+"_bp_rp_run_2.xml")
photometric_data=votable.get_first_table().to_table(use_names_over_ids=True)
sig_teff=np.mean(photometric_data['e_teff_raw'],axis=1)
sig_logg=np.mean(photometric_data['e_logg_raw'],axis=1)
spread_teff=np.std(photometric_data['teff_raw'],axis=1)
spread_logg=np.std(photometric_data['logg_raw'],axis=1)
spread_factor=spread_teff*spread_logg/(sig_teff*sig_logg*len(photometric_data[0]['teff_raw']))
mask=spread_factor>0.5
photometric_data_to_save=photometric_data[mask]
print(len(photometric_data_to_save))
