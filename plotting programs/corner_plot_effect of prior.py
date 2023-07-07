#this program will do corner plots of emcee output of spectra and compare prior vs non prior
import numpy as np
import matplotlib.pyplot as plt
import csv
from astropy.io.votable import parse,from_table,writeto
from astropy.table import Table,vstack,join
import emcee
import matplotlib.patches as mpatches
import subprocess
import scipy.stats as st
from scipy.stats import kde
from os.path import exists
from numba import jit
@jit(nopython=True,cache=True)
def normal(x,mu,sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-1/2*((x-mu)/sigma)**2)

plt.rcParams['font.family']='Times New Roman'
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['text.usetex'] = True

parameters_index=['teff','logg','fe_h','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']

cluster_name='NGC_2682'
data = parse("NGC_2682_no_normalization_fixed_errors.xml",)

data=data.get_first_table().to_table(use_names_over_ids=True)
# mask=data['sobject_id']==131217003901059
# data=data[mask]
computing_cluster='olimp'
labels=['teff','logg','fe_h','vmic','vsini']
# for name in data['sobject_id'][10:20]:
# def corner_plot(name):
name=160106004101321
#ending = labels that has join with '_'
ending="_".join(labels)+'_'
ending='all_elements_long_run_'
data_star=data[data['sobject_id']==name][0]
name=str(name)
#check if the data exists and if not gets it from the database
if not exists(f"{cluster_name}_reduction_fixed_photometric/_{ending}_prior_{name}_main_loop.h5"):
    if computing_cluster=='gigli':
        source=f'kevin@gigli.fmf.uni-lj.si:/media/hdd/home2/kevin/Pysme/{cluster_name}_reduction_fixed_photometric/_{ending}_prior_{name}_main_loop.h5'
    elif computing_cluster=='olimp':
        source=f'beeson@home.fmf.uni-lj.si:/home/beeson/Working_Folder/Pysme/{cluster_name}_reduction_fixed_photometric/_{ending}_prior_{name}_main_loop.h5'

    #source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.1/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'

    destination=cluster_name+'_reduction_fixed_photometric/'
    subprocess.run(["rsync",'-av',source,destination])
if not exists(f"{cluster_name}_reduction_fixed_photometric/_{ending}_prior_{name}_tags.npy"):
    if computing_cluster=='gigli':
        source=f"kevin@gigli.fmf.uni-lj.si:/media/hdd/home2/kevin/Pysme/{cluster_name}_reduction_fixed_photometric/_{ending}_prior_{name}_tags.npy"
    elif computing_cluster=='olimp':
        source=f'beeson@home.fmf.uni-lj.si:/home/beeson/Working_Folder/Pysme/{cluster_name}_reduction_fixed_photometric/_{ending}_prior_{name}_tags.npy'
    #source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.1/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'

    destination=cluster_name+'_reduction_fixed_photometric/'
    subprocess.run(["rsync",'-av',source,destination])
if not exists(f"{cluster_name}_reduction_fixed_photometric/_{ending}_no_prior_{name}_tags.npy"):
    if computing_cluster=='gigli':
        source=f'kevin@gigli.fmf.uni-lj.si:/media/hdd/home2/kevin/Pysme/{cluster_name}_reduction_fixed_photometric/_{ending}_no_prior_{name}_tags.npy'
    elif computing_cluster=='olimp':
        source=f'beeson@home.fmf.uni-lj.si:/home/beeson/Working_Folder/Pysme/{cluster_name}_reduction_fixed_photometric/_{ending}_no_prior_{name}_tags.npy'
    #source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.1/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'

    destination=cluster_name+'_reduction_fixed_photometric/'
    subprocess.run(["rsync",'-av',source,destination])
if not exists(f"{cluster_name}_reduction_fixed_photometric/_{ending}_no_prior_{name}_main_loop.h5"):
    if computing_cluster=='gigli':
        source=f'kevin@gigli.fmf.uni-lj.si:/media/hdd/home2/kevin/Pysme/{cluster_name}_reduction_fixed_photometric/_{ending}_no_prior_{name}_main_loop.h5'
    elif computing_cluster=='olimp':
        source=f'beeson@home.fmf.uni-lj.si:/home/beeson/Working_Folder/Pysme/{cluster_name}_reduction_fixed_photometric/_{ending}_no_prior_{name}_main_loop.h5'
    #source='kevin@gigli.fmf.uni-lj.si:/media/storage/HERMES_REDUCED/dr6.1/'+name[0:6]+'/spectra/com/'+name+str(count)+'.fits'

    destination=cluster_name+'_reduction_fixed_photometric/'
    subprocess.run(["rsync",'-av',source,destination])

#gets the data from the h5 file
reader_prior = emcee.backends.HDFBackend(f'{cluster_name}_reduction_fixed_photometric/_{ending}_prior_{name}_main_loop.h5', read_only=True)
reader_no_prior = emcee.backends.HDFBackend(f'{cluster_name}_reduction_fixed_photometric/_{ending}_no_prior_{name}_main_loop.h5', read_only=True)


kde_prior = reader_prior.get_chain(flat=True,discard=500)
tags_prior=np.load(f"{cluster_name}_reduction_fixed_photometric/_{ending}_prior_{name}_tags.npy")

kde_no_prior = reader_no_prior.get_chain(flat=True,discard=500)
tags_no_prior=np.load(f"{cluster_name}_reduction_fixed_photometric/_{ending}_no_prior_{name}_tags.npy")
# get the autocorrelation time
print('calculating autocorrelation time')
tau_prior = reader_prior.get_autocorr_time(tol=0)
tau_no_prior = reader_no_prior.get_autocorr_time(tol=0)
#get the number of steps
nsteps_prior = len(reader_prior.get_chain(discard=500))
nsteps_no_prior = len(reader_no_prior.get_chain(discard=500))


#make an astropy table of the data 
data_prior=Table(kde_prior,names=tags_prior)
data_no_prior=Table(kde_no_prior,names=tags_no_prior)

#create a new table for the autocorrelation time
tau_prior_table=Table(tau_prior,names=tags_prior)
tau_no_prior_table=Table(tau_no_prior,names=tags_no_prior)

indipendent_samples=100

#2x2 plot
fig, axes = plt.subplots(len(labels), len(labels), figsize=(10, 7))
#create a an array of the limits of each coloumn
limits_array=np.zeros((len(labels),2))
for coloumn,row in np.ndindex(axes.shape):
    print(coloumn,row)
    if coloumn==row:
        prior_to_plot=data_prior[labels[coloumn]][::int(tau_prior_table[f'{labels[coloumn]}'])][:indipendent_samples]
        
        axes[coloumn,row].hist(prior_to_plot,bins=50,alpha=0.5,color='red',density=True)
        no_prior_to_plot=data_no_prior[labels[coloumn]][::int(tau_no_prior_table[f'{labels[coloumn]}'])][:indipendent_samples]
        axes[coloumn,row].hist(no_prior_to_plot,bins=50,alpha=0.5,color='blue',density=True)
        if labels[row]=='teff':
            
            temperature_star=data_star['teff_raw']
            spread_temperature=np.sqrt(np.var(data_star['teff_raw']))
            sig_teff=np.mean(data_star['e_teff_raw'])
            spread_logg=np.sqrt(np.var(data_star['logg_raw']))
            sig_logg=np.mean(data_star['e_logg_raw'])
            teff_raw=data_star['teff_raw']
            logg_raw=data_star['logg_raw']
            e_teff_raw=data_star['e_teff_raw']
            e_logg_raw=data_star['e_logg_raw']
            k=st.gaussian_kde(teff_raw,weights=1/np.sqrt((e_teff_raw*e_logg_raw)))
            #get the limits from axes[coloumn,row]
            limits=axes[coloumn,row].get_xlim()
            x=np.linspace(limits[0],limits[1],100)
            axes[coloumn,row].plot(x,k(x),color='green',linewidth=2)
            #set the limits back
            axes[coloumn,row].set_xlim(limits)
        if labels[row]=='logg':
            temperature_star=data_star['teff_raw']
            spread_temperature=np.sqrt(np.var(data_star['teff_raw']))
            sig_teff=np.mean(data_star['e_teff_raw'])
            spread_logg=np.sqrt(np.var(data_star['logg_raw']))
            sig_logg=np.mean(data_star['e_logg_raw'])
            teff_raw=data_star['teff_raw']
            logg_raw=data_star['logg_raw']
            e_teff_raw=data_star['e_teff_raw']
            e_logg_raw=data_star['e_logg_raw']
            k=st.gaussian_kde(logg_raw,weights=1/np.sqrt((e_teff_raw*e_logg_raw)))
            #get the limits from axes[coloumn,row]
            limits=axes[coloumn,row].get_xlim()
            x=np.linspace(limits[0],limits[1],100)
            axes[coloumn,row].plot(x,k(x),color='green',linewidth=2)
            #set the limits back
            axes[coloumn,row].set_xlim(limits)
    elif coloumn>row:
        for data_emcee,tau_table,colours in zip([data_prior,data_no_prior],[tau_prior_table,tau_no_prior_table],['red','blue']):
            tau=np.max(int(tau_table[f'{labels[coloumn]}']),int(tau_table[f'{labels[row]}']))
            x=data_emcee[labels[row]][::tau][:indipendent_samples]
            y=data_emcee[labels[coloumn]][::tau][:indipendent_samples]
            axes[coloumn,row].scatter(x,y,s=0.1,color=colours,linewidth=0,alpha=0.1)  
            k = kde.gaussian_kde([x,y])
            nbins=100
            xi, yi = np.mgrid[min(x):max(x):nbins*1j, min(y):max(y):nbins*1j]
            zi = k(np.vstack([xi.flatten(), yi.flatten()]))
            con=axes[coloumn,row].contour(xi, yi, zi.reshape(xi.shape), colors=colours,alpha=0.5,levels=5)
            #store the limits of the x axis in the array
            limits_x=axes[coloumn,row].get_xlim()
            limits_array[row]=limits_x

            #axes[coloumn,row].clabel(con, inline=1, fontsize=3)
        #if both of  labels are teff or logg
        if (labels[coloumn]=='teff' and labels[row]=='logg') or (labels[coloumn]=='logg' and labels[row]=='teff'):
            #get the limits from axes[coloumn,row] for x and y 
            limits_x=axes[coloumn,row].get_xlim()
            limits_y=axes[coloumn,row].get_ylim()
            temperature_star=data_star['teff_raw']
            logg_star=data_star['logg_raw']

            spread_temperature=np.sqrt(np.var(data_star['teff_raw']))
            sig_teff=np.mean(data_star['e_teff_raw'])
            spread_logg=np.sqrt(np.var(data_star['logg_raw']))
            sig_logg=np.mean(data_star['e_logg_raw'])
            spread_number=spread_temperature*spread_logg/(sig_logg*sig_teff*len(data_star['e_logg_raw']))
            teff_raw=data_star['teff_raw']
            logg_raw=data_star['logg_raw']
            e_teff_raw=data_star['e_teff_raw']
            e_logg_raw=data_star['e_logg_raw']
            if spread_number>0.1:
                k=st.gaussian_kde([teff_raw,logg_raw],weights=1/np.sqrt((e_teff_raw*e_logg_raw)))
                # teff_line=np.linspace(min(teff_raw)-50,max(logg_raw)+50,num=50)
                # logg_line=np.linspace(min(logg_raw)-0.25,max(logg_raw)+0.25,num=50)
                #xx, yy = np.mgrid[min(teff_raw)-10:max(teff_raw)+10:100j, min(logg_raw)-0.3:max(logg_raw)+0.3   :100j]
                #use the limits to make xx and yy using np.mgrid
                xx, yy = np.mgrid[limits_x[0]:limits_x[1]:100j, limits_y[0]:limits_y[1]:100j]
                zz = k(np.vstack([xx.flatten(), yy.flatten()]))
                zz = np.reshape(zz, xx.shape)
                # polynomial_coeff=data_star['coeff']
                # spectroscopic_shift=np.poly1d(polynomial_coeff)
                # xx-=spectroscopic_shift(xx)
                con=axes[coloumn,row].contour(xx,yy,zz,color='green',alpha=0.5,levels=50)
                #set the limits back
                axes[coloumn,row].set_xlim(limits_x)
                axes[coloumn,row].set_ylim(limits_y)
                #extract the contour levels
            else:
                xx, yy = np.mgrid[min(teff_raw)-10:max(teff_raw)+10:100j, min(logg_raw)-0.1:max(logg_raw)+0.1   :100j]

                prior_parameters=np.column_stack((teff_raw,e_teff_raw,logg_raw,e_logg_raw))
                
                zz=[]
                for x,y in zip(xx,yy):
                    temp_probability=0
                    for temp_param in prior_parameters: 
                        temp_probability+=normal(x,temp_param[0],temp_param[1])*normal(y,temp_param[2],temp_param[3])        
                    zz.append(temp_probability)
                zz=np.reshape(zz,xx.shape)
                # polynomial_coeff=data_star['coeff']
                # spectroscopic_shift=np.poly1d(polynomial_coeff)
                # xx-=spectroscopic_shift(xx)
                
                

    else:
        axes[coloumn,row].axis('off')
    if row==0:
        axes[coloumn,row].set_ylabel(labels[coloumn])
    if coloumn==len(labels)-1:
        axes[coloumn,row].set_xlabel(labels[row])
#set the limits of column == row to the limits of the other plots
for coloumn in range(len(labels)-1):

    axes[coloumn,coloumn].set_xlim(limits_array[coloumn])
red_patch = mpatches.Patch(color='red', label=('prior').replace("_", " "))
blue_patch = mpatches.Patch(color='blue', label=('no_prior').replace("_", " "))
green_patch = mpatches.Patch(color='green', label=('phtometric_prior').replace("_", " "))
fig.legend(handles=[red_patch,blue_patch,green_patch],loc=(0.46,0.755),borderaxespad=0)


fig.savefig(f'example_figs/{name}_{indipendent_samples}_indipendent_samples.png',dpi=300)

# from multiprocessing import Pool
# #use multiprocessing to plot the data
# if __name__ == '__main__':
#     pool = Pool(10)
#     pool.map(corner_plot, [160107004101259])
#     pool.close()
#     pool.join()