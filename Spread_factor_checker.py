
#This program will load data from no_normalization_.xml files and plot the e_fe_h vs spread factor
#load the data
import numpy as np
import matplotlib.pyplot as plt
import csv
from astropy.io.votable import parse,from_table,writeto
from astropy.table import Table,vstack,join


plt.rcParams['font.family']='Times New Roman'
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['text.usetex'] = True

name='NGC_2682'
votable = parse(name+"_no_normalization_.xml")
data=votable.get_first_table().to_table(use_names_over_ids=True)

teff_raw=data['teff_raw']
logg_raw=data['logg_raw']
e_teff_raw=data['e_teff_raw']
e_logg_raw=data['e_logg_raw']
spread_temperature=np.std(teff_raw,axis=1)
spread_logg=np.std(logg_raw,axis=1)
sig_teff=np.mean(e_teff_raw,axis=1)
sig_logg=np.mean(e_logg_raw,axis=1)
numbers=[len(x) for x in teff_raw]

spread_factor=spread_temperature*spread_logg/(sig_logg*sig_teff*numbers)

plt.figure()
plt.scatter(spread_factor,data['e_logg_prior']/data['e_logg_no_prior'],s=0.7,c='black')

plt.xlim((0,2))
plt.ylim((0,2))
plt.xlabel('spread factor')
plt.ylabel(r'$\sigma_{\rm{ [Fe]\, prior}}$/$\sigma_{\rm{ [Fe]\, no\, prior}}$')
#set figure size
plt.gcf().set_size_inches(6.97384806974/2,6.97384806974/4)
plt.tight_layout()
plt.savefig('spread_factor_vs_error_ratio.pdf')