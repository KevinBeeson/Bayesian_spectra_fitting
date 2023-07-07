#This program will do a running average of the data and plot it but only for individual elements

import numpy as np
import matplotlib.pyplot as plt
import csv
from astropy.io.votable import parse,from_table,writeto
from astropy.table import Table,vstack,join
import matplotlib.patches as mpatches


plt.rcParams['font.family']='Times New Roman'
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['text.usetex'] = True
cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)
print(cluster_details_all[:,0])
# name=input('what isochrone do you want?')

convergence=1
parameters=['logg','logg','Fe/H','vmic','vbroad','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
parameters_index=['logg','logg','fe_h','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
parameters_new=[]
for x in parameters_index:
    if x in elements:
        parameters_new.append(x+'_Fe')
    else:
        parameters_new.append(x)
parameters_index=parameters_new
cluster_1='NGC_2682'
cluster_2='NGC_2682'
file_name_1=cluster_1+'_no_normalization_.xml'
file_name_2=cluster_2+'_no_normalization_.xml'
circle_size=1
colours=['red','blue']
endings=['_prior','_no_prior']
element_to_do=['Co','Mn','Cr','K','O','Si']
original_element_to_do=element_to_do.copy()
to_remove_temp=[]
for x in element_to_do:
    if x in elements:
        element_to_do.append(x+'_Fe')
        to_remove_temp.append(x)
for x in to_remove_temp:
    element_to_do.remove(x)
converged=True
size_font=10
for name in ['NGC_2682']:
    spread=0

    # fig.set_figwidth(shape[1])
    # fig.set_figheight(shape[0])
    
    cluster_1=cluster_2=name
    file_name_1=cluster_1+'_no_normalization_.xml'
    file_name_2=cluster_2+'_no_normalization_.xml'
    votable = parse(file_name_1)
    data_1=votable.get_first_table().to_table(use_names_over_ids=True)
    if file_name_1==file_name_2:
        data_2=data_1.copy()
    else:
        votable = parse(file_name_2)
        data_2=votable.get_first_table().to_table(use_names_over_ids=True)
    mask=data_1['vsini'+endings[0]]<20
    data_1=data_1[mask]
    mask=data_2['vsini'+endings[0]]<20
    data_2=data_2[mask]
    # #do another mask for e_teff larger than 25
    # mask=data_1['e_teff'+endings[0]]<25
    # data_1=data_1[mask]
    # mask=data_2['e_teff'+endings[1]]<25
    # data_2=data_2[mask]
    # #do another mask for logg larger than 0.06
    # mask=data_1['e_logg'+endings[0]]<0.06
    # data_1=data_1[mask]
    # mask=data_2['e_logg'+endings[1]]<0.10
    

    if converged:
        mask=[len(x)>3 for x in data_1['radial_velocities_prior']]
        data_1=data_1[mask]
        data_2=data_2[mask]
        mask=data_1['vsini'+endings[0]]<20
        data_1=data_1[mask]
        mask=data_2['vsini'+endings[0]]<20
        data_2=data_2[mask]
        for elem in ['fe_h']:
            print(elem)
            mask=data_1[elem+'_prior']<np.nanpercentile(np.array(data_1[elem+'_prior']),99)
            mask*=data_1[elem+'_prior']>np.nanpercentile(np.array(data_1[elem+'_prior']),1)
            #if data is masked then 
            masked=[np.ma.is_masked(x) for x in data_1[elem+'_prior']]
            mask=np.logical_or(mask,masked)
            mask*=data_2[elem+'_no_prior']<np.nanpercentile(np.array(data_2[elem+'_no_prior']),99)
            mask*=data_2[elem+'_no_prior']>np.nanpercentile(np.array(data_2[elem+'_no_prior']),1)
            
            masked=[np.ma.is_masked(x) for x in data_2[elem+'_no_prior']]
            mask=np.logical_or(mask,masked)

            data_1=data_1[mask]
            data_2=data_2[mask]
            print(len(data_1))
        # mask=data_1['flag_sp']==0
        # data_1=data_1[mask]
        # mask=data_2['flag_sp']==0
        # data_2=data_2[mask]
    #create bins here for both endings
    for col_dat,original_elem in zip(element_to_do,original_element_to_do):
        if '_Galah_iDR4' in endings:
            if endings[0]=='_Galah_iDR4':
                bins=np.linspace(min(min(data_1['logg']),min(data_2['logg'+endings[1]]))-0.1,max(max(data_1['logg']),max(data_2['logg'+endings[1]]))+0.1,10)
            else:
                bins=np.linspace(min(min(data_1['logg_prior']),min(data_2['logg']))-0.1,max(max(data_1['logg'+endings[0]]),max(data_2['logg']))+0.1,10)
        else:
            mask=data_1[col_dat+endings[0]]>-10
            data_temp_1=data_1[mask] 
            mask=data_2[col_dat+endings[1]]>-10
            data_temp_2=data_2[mask]
            if len(data_temp_1)==0:
                bins=np.linspace(min(data_temp_2))
            else:    
                bins=np.linspace(min(min(data_temp_1['logg'+endings[0]]),min(data_temp_2['logg'+endings[1]]))-0.3,max(max(data_temp_1['logg'+endings[0]]),max(data_temp_2['logg'+endings[1]]))+0.3,10)
        
        fig = plt.figure()
        #plt.title(col_dat)
        average_std=[]

       

        for data,colour,ending,cluster_name in zip([data_1,data_2],colours,endings,[cluster_1,cluster_2]):

            # NGC_2682_summary_dr61
                
            if len(data)!=0:
                size_font=10
                iso_type='gaia'


                if ending=="_Galah_iDR4":
                    plt.scatter(data['logg_spectroscopic'],data[col_dat.lower()],s=0.3,c=colour,alpha=0.5)
                else:
                    plt.scatter(data['logg_prior'],data[col_dat+ending],s=circle_size,marker=',',c=colour,alpha=0.5)
                average=np.array([])
                std=np.array([])
                numbers=np.array([])
                for value in range(len(bins)-1):
                    if ending=="_Galah_iDR4":
                        average=np.append(average,np.mean(data[col_dat.lower()][np.where((data['logg_spectroscopic']>bins[value]) & (data['logg_spectroscopic']<bins[value+1]))]))
                        std=np.append(std,np.std(data[col_dat.lower()][np.where((data['logg_spectroscopic']>bins[value]) & (data['logg_spectroscopic']<bins[value+1]))]))
                    else:
                        numbers=np.append(numbers,len(data[col_dat+ending][np.where((data['logg' +ending]>bins[value]) & (data['logg' +ending]<bins[value+1]))]))
                        average=np.append(average,np.mean(data[col_dat+ending][np.where((data['logg' +ending]>bins[value]) & (data['logg' +ending]<bins[value+1]))]))
                        std=np.append(std,np.std(data[col_dat+ending][np.where((data['logg' +ending]>bins[value]) & (data['logg' +ending]<bins[value+1]))]))
                plt.plot((bins[:-1]+bins[1:])/2,average,linewidth=0.5,c=colour)
                plt.fill_between((bins[:-1]+bins[1:])/2,average-std,average+std,alpha=0.5,color=colour)
            numbers-=1
            average_std=np.append(average_std,sum(std*numbers)/len(numbers))
            red_patch = mpatches.Patch(color='red', label=(cluster_1+endings[0]).replace("_", " "))
            blue_patch = mpatches.Patch(color='blue', label=(cluster_2+endings[1]).replace("_", " "))

            #fig.legend(handles=[red_patch,blue_patch],loc=(0.71,0.905))
            plt.tight_layout(w_pad=0.0,h_pad=-1.2)
            fig.set_size_inches(3.32088003321,3.32088003321/1.61)
        #x axis logg label
        plt.xlabel(r'log ($g$)')
        #y axis label
        y_label=r'['+original_elem+'/Fe]'
        plt.ylabel(y_label)
        #flip x axis
        plt.gca().invert_xaxis()
        plt.savefig('/home/kevin/Documents/Paper/new/'+cluster_1+endings[0]+'_'+cluster_2+endings[1]+'_masks_'+str(converged)+'element_'+col_dat+'_same_log_g.pdf',bbox_inches='tight')



# #load NGC_2682 data
# votable = parse('NGC_2682_no_normalization_.xml')
# data=votable.get_first_table().to_table(use_names_over_ids=True)

# fig=plt.figure()
# plt.scatter(data_1['teff_prior'],data_1['logg_prior'],c=data_1['O_Fe_prior'],s=0.8,vmin=-0.3,vmax=0.3)
# plt.gca().invert_xaxis()
# plt.gca().invert_yaxis()
# plt.colorbar()

# plt.ylim(4.5,3.5)
# plt.xlim(6250,5500)
# fig.set_size_inches(3.32088003321,3.32088003321/1.61)
# #x label and y label in latex 
# plt.xlabel(r'$T_{eff}$')
# plt.ylabel(r'$\log g$')
# plt.tight_layout()
# plt.savefig('/home/kevin/Documents/Paper/new/NGC_2682_O_Fe_prior.pdf',bbox_inches='tight')

# fig=plt.figure()
# plt.scatter(data_1['teff_prior'],data_1['logg_prior'],c=data_1['O_Fe_no_prior'],s=0.8,vmin=-0.3,vmax=0.3)
# plt.gca().invert_xaxis()
# plt.gca().invert_yaxis()
# plt.colorbar()

# plt.ylim(4.5,3.5)
# plt.xlim(6250,5500)
# fig.set_size_inches(3.32088003321,3.32088003321/1.61)
# #x label and y label in latex 
# plt.xlabel(r'$T_{eff}$')
# plt.ylabel(r'$\log g$')
# plt.tight_layout()
# plt.savefig('/home/kevin/Documents/Paper/new/NGC_2682_O_Fe_no_prior.pdf',bbox_inches='tight')


# fig=plt.figure()
# plt.scatter(data_1['teff_prior'],data_1['logg_prior'],c=data_1['Ca_Fe_prior'],s=0.8,vmin=-0.3,vmax=0.1)
# plt.gca().invert_xaxis()
# plt.gca().invert_yaxis()
# plt.colorbar()

# plt.ylim(4.5,3.5)
# plt.xlim(6250,5500)
# fig.set_size_inches(3.32088003321,3.32088003321/1.61)
# #x label and y label in latex 
# plt.xlabel(r'$T_{eff}$')
# plt.ylabel(r'$\log g$')
# plt.tight_layout()
# # plt.savefig('/home/kevin/Documents/Paper/new/NGC_2682_O_Fe_prior.pdf',bbox_inches='tight')

# fig=plt.figure()
# plt.scatter(data_1['teff_no_prior'],data_1['logg_no_prior'],c=data_1['Ca_Fe_no_prior'],s=0.8,vmin=-0.3,vmax=0.1)
# plt.gca().invert_xaxis()
# plt.gca().invert_yaxis()
# plt.colorbar()

# plt.ylim(4.5,3.5)
# plt.xlim(6250,5500)
# fig.set_size_inches(3.32088003321,3.32088003321/1.61)
# #x label and y label in latex 
# plt.xlabel(r'$T_{eff}$')
# plt.ylabel(r'$\log g$')
# plt.tight_layout()
# plt.savefig('/home/kevin/Documents/Paper/new/NGC_2682_O_Fe_no_teff_diff_prior.pdf',bbox_inches='tight')


# fig=plt.figure()
# plt.scatter(data_1['teff_prior'],data_1['logg_prior'],c=data_1['O_Fe_prior']-data_1['O_Fe_no_prior'],s=0.8,vmin=-0.1,vmax=0.1)
# plt.gca().invert_xaxis()
# plt.gca().invert_yaxis()
# plt.colorbar()

# plt.ylim(4.5,3.5)
# plt.xlim(6250,5500)
# fig.set_size_inches(3.32088003321,3.32088003321/1.61)
# #x label and y label in latex 


# plt.xlabel(r'$T_{eff}$')
# plt.ylabel(r'$\log g$')
# plt.tight_layout()