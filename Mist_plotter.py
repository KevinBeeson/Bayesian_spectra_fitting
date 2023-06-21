#This file will plot the isochrones from padava and mist and compare them both to each other 
import numpy as np
import matplotlib.pyplot as plt
import csv
from astropy.io.votable import parse

plt.rcParams['font.family']='Times New Roman'
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['text.usetex'] = True

# import names of clusters from dr3_clusters.txt first row
cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)
names=cluster_details_all[:,0]

name=names[0]
for name in ['NGC_2516']:  
    #load the data from xml files
    votable = parse(name+"_bp_rp_run.xml")
    data=votable.get_first_table().to_table(use_names_over_ids=True)

    #load the data from the isochrones
    isochrone_data_mist=np.genfromtxt('Mist_isochrones/'+name+'/'+name+'_main')
    isochrone_data_padava=np.genfromtxt('ALL_ISO/'+name+'_Best_fit.txt')

    fig=plt.figure()
    plt.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],c='black',s=0.2)
    plt.plot(isochrone_data_mist[:,-3]-isochrone_data_mist[:,-2],isochrone_data_mist[:,-4],c='red',linewidth=0.9)
    plt.plot(isochrone_data_padava[:,-2]-isochrone_data_padava[:,-1],isochrone_data_padava[:,-3],linewidth=0.9)
    #flip y axis
    plt.gca().invert_yaxis()
    plt.ylabel('$M_G$')
    plt.xlabel('$G_{BP}-G_{RP}$')
    plt.xlim((-0.1,data['bp_rp'].max()+0.1))
    # plt.xlim((data['bp_rp'].min()-0.1,data['bp_rp'].max()+0.1))
    plt.ylim((data['abs_phot_g_mean_mag'].max()-0.1,data['abs_phot_g_mean_mag'].min()+0.1))
    fig.set_size_inches(3.32088003321,3.32088003321/1.61)
    plt.tight_layout()
    plt.savefig('Isochrone_comparison/'+name+'_isochrone_comparison.pdf')

