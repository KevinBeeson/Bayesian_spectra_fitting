
#in this program we will stack all the _normalization_.xml together to make a single file

#get the names of the clusters from the dr3_clusters.txt file
import numpy as np
import matplotlib.pyplot as plt
import csv
from astropy.io.votable import parse,from_table,writeto
from astropy.table import Table,vstack,join

plt.rcParams['font.family']='Times New Roman'
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['text.usetex'] = True


cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)

#load the data from the xml files
votable = parse(cluster_details_all[0,0]+"_no_normalization_.xml")
data_all   = votable.get_first_table().to_table(use_names_over_ids=True)
#delete radial velocties column
del data_all['radial_velocities_no_prior']
del data_all['radial_velocities_prior']
for name in cluster_details_all[1:,0]:
    print(name)
    votable = parse(name+"_no_normalization_.xml")
    data=votable.get_first_table().to_table(use_names_over_ids=True)
    #delete radial velocties column
    del data['radial_velocities_no_prior']
    del data['radial_velocities_prior']
    data_all=vstack([data_all,data])
elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
other=['teff','logg','fe_h','vmic','vsini']
#do a plot of elements np.ma.median(prior / no prior) 
plt.figure()
ratio_elements=[np.ma.median(data_all['e_'+x+'_Fe_prior']/data_all['e_'+x+'_Fe_no_prior']) for x in elements]
#create a for loop and append the other ratios also append at the beginning the elements ratios
ratio_other=[np.ma.median(data_all['e_'+x+'_prior']/data_all['e_'+x+'_no_prior']) for x in other]
ratio=ratio_other+ratio_elements

plt.bar(np.arange(len(elements)+len(other)),ratio,color='black')
plt.ylabel(r'$\sigma_{\rm{prior}}  / \sigma_{\rm{no prior}}$')
#for the x ticks put element and other names 
other[0]=r'$T_{\rm{eff}}$'
other[2]=r'$\small{\rm{[Fe/H]}}$'
other[1]=r'log\tiny{(}$g$\tiny{)}'
plt.xticks(np.arange(len(elements)+len(other)),other+elements,rotation=90)
#put horizontal lines at every 0.1  tick
plt.yticks(np.arange(0,1.1,0.1))
plt.grid(axis='y')

plt.gcf().set_size_inches(6.97384806974/2,6.97384806974/4)
#adjust the limits of x axis so it stops just after the last tick and stars just at the first tick

plt.xlim(-0.5,len(elements)+len(other)-0.5)
#make x ticks smaller
plt.tick_params(axis='x', which='major', labelsize=6)
plt.tight_layout(pad=0.1)
plt.savefig('e_prior_vs_e_no_prior_ratio.pdf')


#do another plot of only the median value of the errors in the elements prior vs non prior the plot should be a scatter plot with the elements on the x axis 
median_prior=[np.ma.median(data_all['e_'+x+'_Fe_prior']) for x in elements]
median_no_prior=[np.ma.median(data_all['e_'+x+'_Fe_no_prior']) for x in elements]
plt.figure()
#on the x axis put elements
plt.scatter(np.arange(len(elements)),median_prior,color='black',label='prior')
plt.scatter(np.arange(len(elements)),median_no_prior,color='red',label='no prior')
plt.ylabel(r'$\sigma$')
plt.xticks(np.arange(len(elements)),elements,rotation=90)
plt.legend()
plt.grid(axis='y')
plt.gcf().set_size_inches(6.97384806974/2,6.97384806974/4)
plt.xlim(-0.5,len(elements)-0.5)
plt.tick_params(axis='x', which='major', labelsize=7)

#save data_all as a xml file
votable=from_table(data_all)
writeto(votable,'all_clusters_no_normalization_.xml')