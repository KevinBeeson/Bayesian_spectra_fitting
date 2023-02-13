#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 10:15:04 2022

@author: kevin
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
def leaky_relu(z,grad=False):
    '''
    This is the activation function used by default in all our neural networks.
    '''
    limits=(z > 0)
    leaky=z*limits + 0.01*z*np.invert(limits)
    if grad:
        return leaky,limits
    return leaky
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
        
limits={'Blue':[4705,4908],'Green':[5643,5879],'Red':[6470,6743],'IR':[7577.0,7894.0]}
colour='Blue'        
tmp = np.load("NN_normalized_spectra_Blue_smaller.npz")
w_array_0 = tmp["w_array_0"]
w_array_1 = tmp["w_array_1"]
w_array_2 = tmp["w_array_2"]
b_array_0 = tmp["b_array_0"]
b_array_1 = tmp["b_array_1"]
b_array_2 = tmp["b_array_2"]
x_min = tmp["x_min"]
x_max = tmp["x_max"]
tmp.close()
NN_coeffs= (w_array_0, w_array_1, w_array_2, b_array_0, b_array_1, b_array_2, x_min, x_max)
length_synthetic=np.linspace(limits[colour][0], limits[colour][1],num=len(b_array_2))
data_training=np.load('Blue_band_5000_6000_data_training.npy',allow_pickle=True)
spectra=data_training
labels=np.load('Blue_band_5000_6000_training_labels.npy',allow_pickle=True)
spectra_synthetic=[payne_sythesize(x,x_min,x_max,NN_coeffs) for x in labels]
diff=[abs(x-y) for (x,y) in zip(spectra_synthetic,spectra)]
diff=np.vstack(diff)
diff_from_individual_mean=[abs(x-np.mean(x)) for x in spectra]
mean_individual=[np.mean(diff_from_individual_mean,axis=0)]
max_diff=[np.max(diff,axis=0)]
mean_diff=[np.mean(diff,axis=0)]
ratio=mean_diff[0]/mean_individual[0]
line, wavelength = np.loadtxt('important_lines',usecols=(0,1),unpack=True,dtype=str, comments=';')
np.append(wavelength,'4861.3230')
np.append(wavelength,'6562.7970')
np.append(line,r'H$_\beta$')
np.append(line,r'H$_\alpha$')
wavelength=wavelength.astype(float)
important_lines=np.vstack([[elem_temp,wave_temp] for elem_temp,wave_temp in zip(line,wavelength) if wave_temp>length_synthetic[0] and wave_temp<length_synthetic[-1]])
plt.figure()
# plt.plot(length_synthetic,np.transpose(max_diff))
split=4
for y in range(split):
    split_amount=len(length_synthetic)//split
    temp_synthetic_length=length_synthetic[y*split_amount:(y+1)*split_amount]
    temp_mean=mean_diff[0][y*split_amount:(y+1)*split_amount]
    
    plt.figure(figsize=(6,2))   
    plt.plot(temp_synthetic_length,np.transpose(temp_mean),label=r'Mean difference')
    plt.ylabel(r'Normalized flux')
    for individual_line in important_lines:
        if float(individual_line[1])>temp_synthetic_length[0] and float(individual_line[1])<temp_synthetic_length[-1]:
                plt.axvline(x=float(individual_line[1]),color='red',label=r'Important lines',zorder=0)
                # plt.text(float(individual_line[1]),0.02,individual_line[0][:2],fontsize=20,ha='center',color='red')
    plt.xlim((temp_synthetic_length[0],temp_synthetic_length[-1]))
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    if y==split-1:
        plt.xlabel(r'Wavelength / $\mathrm{\AA}$')
    plt.legend(by_label.values(), by_label.keys(),loc='upper right')
    plt.tight_layout()
    plt.savefig('mean difference '+ colour+str(y)+'.pdf')
for y in range(split):
    split_amount=len(length_synthetic)//split
    temp_synthetic_length=length_synthetic[y*split_amount:(y+1)*split_amount]
    temp_max=max_diff[0][y*split_amount:(y+1)*split_amount]
    plt.figure(figsize=(6,2))   
    plt.plot(temp_synthetic_length,np.transpose(temp_max),label=r'Max difference')
    plt.ylabel(r'Normalized flux')
    for individual_line in important_lines:
        if float(individual_line[1])>temp_synthetic_length[0] and float(individual_line[1])<temp_synthetic_length[-1]:
                plt.axvline(x=float(individual_line[1]),color='red',label=r'Important lines',zorder=0)
                # plt.text(float(individual_line[1]),0.02,individual_line[0][:2],fontsize=20,ha='center',color='red')
    plt.xlim((temp_synthetic_length[0],temp_synthetic_length[-1]))
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    if y==split-1:
        plt.xlabel(r'Wavelength / $\mathrm{\AA}$')
    plt.legend(by_label.values(), by_label.keys(),loc='upper right')
    plt.tight_layout()
    plt.savefig('max difference '+ colour+str(y)+'.pdf')
    
for y in range(split):
    split_amount=len(length_synthetic)//split
    temp_synthetic_length=length_synthetic[y*split_amount:(y+1)*split_amount]
    temp_mean=ratio[y*split_amount:(y+1)*split_amount]
    
    plt.figure(figsize=(6,2))   
    plt.plot(temp_synthetic_length,np.transpose(temp_mean),label=r'Mean difference')
    plt.ylabel(r'Normalized flux')
    for individual_line in important_lines:
        if float(individual_line[1])>temp_synthetic_length[0] and float(individual_line[1])<temp_synthetic_length[-1]:
                plt.axvline(x=float(individual_line[1]),color='red',label=r'Important lines',zorder=0)
                # plt.text(float(individual_line[1]),0.02,individual_line[0][:2],fontsize=20,ha='center',color='red')
    plt.xlim((temp_synthetic_length[0],temp_synthetic_length[-1]))
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    if y==split-1:
        plt.xlabel(r'Wavelength / $\mathrm{\AA}$')
    plt.legend(by_label.values(), by_label.keys(),loc='upper right')
    plt.tight_layout()
# for value,lab in enumerate(labels):
#     print(value)
#     payne_sythesize(lab,x_min,x_max,NN_coeffs)