# This fille will import all the data from the reduced open clusters and compare the variance between prior and no prior data 
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

name=names[7]


parameters=['Fe/H','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
parameters_index=['fe_h','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
parameters_new=[]
for x in parameters_index:
    if x in elements:
        parameters_new.append(x+'_Fe')
    else:
        parameters_new.append(x)
parameters_index=parameters_new

endings=['_prior','_Galah_iDR4']
numbers_all=[]
all_ratios=[]
potassium_check=[]
for name in names:
    votable = parse(name+"_no_normalization_with_new_sven.xml")
    data=votable.get_first_table().to_table(use_names_over_ids=True)
    mask=data['vsini_prior']<20
    data=data[mask]
    mask=data['flag_sp_sven']==0
    data=data[mask]
    if len(data)>1:
        if '_Galah_iDR4' in endings:
            if endings[0]=='_Galah_iDR4':
                bins=np.linspace(min(min(data['logg']),min(data['logg'+endings[1]]))-0.1,max(max(data['logg']),max(data['logg'+endings[1]]))+0.1,10)
            else:
                bins=np.linspace(min(min(data['logg'+endings[0]]),min(data['logg']))-0.1,max(max(data['logg'+endings[0]]),max(data['logg']))+0.1,10)
        else:
            bins=np.linspace(min(min(data['logg'+endings[0]]),min(data['logg'+endings[1]]))-0.1,max(max(data['logg'+endings[0]]),max(data['logg'+endings[1]]))+0.1,10)
    
        std_ending_1=[]
        std_ending_2=[]
        mean_ending_1=[]
        mean_ending_2=[]
        #checks how many potaissum is masked 
        potassium_check.append(len(data['K_Fe_prior'][[np.ma.is_masked(x) for x in data['K_Fe_no_prior']]]))
        for ending,std_all,mean_all in zip(endings,[std_ending_1,std_ending_2],[mean_ending_1,mean_ending_2]):
            temp_std=[]
            temp_mean=[]
            for value,(col_dat) in enumerate(parameters_index):
                average=np.array([])
                std=np.array([])
                numbers=np.array([])
                for value in range(len(bins)-1):
                    if ending=="_Galah_iDR4":
                        average=np.append(average,np.mean(data[col_dat.lower()+'_sven'][np.where((data['logg']>bins[value]) & (data['logg']<bins[value+1]))]))
                        std=np.append(std,np.std(data[col_dat.lower()+'_sven'][np.where((data['logg']>bins[value]) & (data['logg']<bins[value+1]))]))
                        #number os elements in each bin
                        numbers=np.append(numbers,len(data[col_dat.lower()+'_sven'][np.where((data['logg']>bins[value]) & (data['logg']<bins[value+1]))])-1)
                    else:
                        average=np.append(average,np.mean(data[col_dat+ending][np.where((data['logg' +ending]>bins[value]) & (data['logg' +ending]<bins[value+1]))]))
                        std=np.append(std,np.std(data[col_dat+ending][np.where((data['logg' +ending]>bins[value]) & (data['logg' +ending]<bins[value+1]))]))
                        #number os elements in each bin
                        numbers=np.append(numbers,len(data[col_dat+ending][np.where((data['logg' +ending]>bins[value]) & (data['logg' +ending]<bins[value+1]))])-1)
                mean_std=sum(std*numbers)/sum(numbers)
                mean_average=sum(average*numbers)/sum(numbers)
                temp_std.append(mean_std)
                temp_mean.append(mean_average)
            std_all.append(temp_std)
            mean_all.append(temp_mean)
        ratio=np.array(std_ending_1[0])/np.array(std_ending_2[0])
        all_ratios.append(ratio)
        numbers_all.append(sum(numbers))
        fig=plt.figure()
        plt.scatter(parameters,ratio)
        #plt.title('Ratio of spread  of prior and no prior data in '+name)
        #plt.xlabel('Parameters')
        plt.ylabel('Ratio of spread')
        #create vertical dashed lines for each tick
        for xtick in plt.xticks()[0]:
            plt.axvline(x=xtick, color='k', linestyle='--')
        #create a horizontal line at 1
        plt.axhline(y=1, color='k', linestyle='--')
        plt.tight_layout()
        plt.xticks(rotation=90)
        fig.set_size_inches(3.32088003321,3.32088003321/1.61)
        plt.xticks(fontsize=7)
        plt.tight_layout()
        plt.savefig('Ratio_of_spread_mask_True_'+name+'_logg.pdf')      
overall_ratio=np.zeros(len(parameters))
numbers_temp=np.zeros(len(parameters))
for rat,num in zip(np.array(all_ratios),numbers_all):
    
    for value,x in enumerate(rat):
        if not np.isnan(x) and not np.isinf(x):
            overall_ratio[value]+=x*num
            numbers_temp[value]+=num
    
overall_ratio=overall_ratio/numbers_temp

fig=plt.figure()
plt.scatter(parameters,overall_ratio)
print(np.mean(overall_ratio))
#plt.title('overall Ratio of standard deviation of prior and no prior data')
#plt.xlabel('Parameters')
plt.xticks(rotation=90) 
plt.ylabel('Ratio of standard deviation')
for xtick in plt.xticks()[0]:
    plt.axvline(x=xtick, color='k', linestyle='--')
#create a horizontal line at 1
plt.axhline(y=1, color='k', linestyle='--')    
fig.set_size_inches(3.32088003321,3.32088003321/1.61)
plt.tight_layout()
plt.xticks(fontsize=7)
plt.savefig('overall_Ratio_spread_mask_True.pdf')

#load all_clusters_no_normalization_.xml
votable = parse('all_clusters_no_normalization_.xml')
data_all=votable.get_first_table().to_table(use_names_over_ids=True)
#create a mask like above 
mask=data_all['vsini_prior']<20
data_all=data_all[mask]
mask=data_all['flag_sp']==0
data_all=data_all[mask]

#get the ratio of errors for each element
elements=['Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
other=['teff','logg','vmic','vsini']
#do a plot of elements np.ma.median(prior / no prior) 
ratio_elements=[np.ma.median(data_all['e_'+x+'_Fe_prior']/data_all['e_'+x+'_Fe_no_prior']) for x in elements]
ratio_fe=[np.ma.median(data_all['e_fe_h_prior']/data_all['e_fe_h_no_prior'])]
ratio_all=ratio_fe+ratio_elements

ratio_other=[np.ma.median(data_all['e_'+x+'_prior']/data_all['e_'+x+'_no_prior']) for x in other]
ratio_all=ratio_other+ratio_all
#plot ratio_all to overall_ratio
fig=plt.figure()
plt.scatter(ratio_all,overall_ratio)
plt.xlabel('Ratio of standard deviation of prior and no prior data')
plt.ylabel('Ratio of spread')

fig=plt.figure()
other[0]=r'$T_{\rm{eff}}$'
other[1]=r'log\tiny{(}$g$\tiny{)}'
#create bar chart of ratio of errors for each element and have it have a black box around it
plt.bar(other+parameters,ratio_all,edgecolor='black',alpha=0.5)
plt.bar(parameters,overall_ratio,edgecolor='black',alpha=0.5)

plt.xticks(rotation=90) 
plt.ylabel('Ratio of errors or variance')
#create a horizontal line every 0.25
for i in np.arange(0,2,0.25):
    plt.axhline(y=i, color='k', linestyle='--')
#shift the x axis so it starts at the first element
plt.xlim(-0.5,len(elements)-0.5)
plt.ylim(0.5,2)
fig.set_size_inches(3.32088003321,3.32088003321/1.61)
plt.xticks(fontsize=8.4)
plt.tight_layout(pad=0.2)
plt.savefig("Ratio of errors and variance mask True.pdf")



all_ratios=[]
numbers_all=[]
all_ratio_priors=[]
all_ratio_no_priors=[]
for name in names:
    votable = parse(name+"_no_normalization_.xml")
    data=votable.get_first_table().to_table(use_names_over_ids=True)
    mask=data['vsini_prior']<20
    data=data[mask]
    mask=data['flag_sp']==0
    data=data[mask]
    if len(data)>1:

        std_ending_1=[]
        std_ending_2=[]
        mean_ending_1=[]
        mean_ending_2=[]

        for ending,std_all,mean_all in zip(endings,[std_ending_1,std_ending_2],[mean_ending_1,mean_ending_2]):
            temp_std=[]
            temp_mean=[]
            numbers=[]
            for value,(col_dat) in enumerate(parameters_index):
                average=np.array([])
                std=np.array([])

                if ending=="_Galah_iDR4":
                    average=np.append(average,np.nan.mean(data[col_dat.lower()+'_sven']))
                    std=np.append(std,np.nan.std(data[col_dat.lower()+'_sven']))
                    #append numbers with the number of elements which arent masked
                    numbers=np.append(numbers,len(data[col_dat.lower()+'_sven'])-np.ma.count_masked(data[col_dat.lower()+'_sven']))

                else:
                    if np.ma.count_masked(data[col_dat+ending])!= len(data[col_dat+ending]):
                        #append average with the median of the whole data set
                        average=np.append(average,np.mean(data[col_dat+ending]))
                        #append std with the standard deviation of the whole data set
                        std=np.append(std,np.nanstd(data[col_dat+ending]))
                        numbers=np.append(numbers,len(data[col_dat+ending])-np.ma.count_masked(data[col_dat+ending]))
                    else:
                        #append average with the median of the whole data set
                        average=np.append(average,np.nan)
                        #append std with the standard deviation of the whole data set
                        std=np.append(std,np.nan)
                        numbers=np.append(numbers,np.nan)

                temp_std.append(std)
                temp_mean.append(average)
            std_all.append(temp_std)
            mean_all.append(temp_mean)
            # numbers_all.append(numbers)
        ratio=np.array(std_ending_1[0])/np.array(std_ending_2[0])
        all_ratios.append(ratio)
        numbers_all.append(len(data))
        fig=plt.figure()
        plt.scatter(parameters,ratio)
        #plt.title('Ratio of spread  of prior and no prior data in '+name)
        #plt.xlabel('Parameters')
        plt.ylabel('Ratio of spread')
        #create vertical dashed lines for each tick
        for xtick in plt.xticks()[0]:
            plt.axvline(x=xtick, color='k', linestyle='--')
        #create a horizontal line at 1
        plt.axhline(y=1, color='k', linestyle='--')
        plt.tight_layout()
        plt.xticks(rotation=90)
        fig.set_size_inches(3.32088003321,3.32088003321/1.61)
        plt.xticks(fontsize=7)
        plt.tight_layout()
        plt.title(name)
        all_ratio_priors.append(std_ending_1[0])
        all_ratio_no_priors.append(std_ending_2[0])

overall_ratio=np.zeros(len(parameters))
numbers_temp=np.zeros(len(parameters))
for rat,num in zip(np.array(all_ratios),numbers_all):
    
    for value,x in enumerate(rat):
        if not np.isnan(x) and not np.isinf(x):
            overall_ratio[value]+=x*num
            numbers_temp[value]+=num
    
overall_ratio=overall_ratio/numbers_temp

overall_variance_prior=np.zeros(len(parameters))
overall_variance_no_prior=np.zeros(len(parameters))
numbers_temp_prior=np.zeros(len(parameters))
numbers_temp_no_prior=np.zeros(len(parameters))
for rat_prior,rat_no_prior,num in zip(np.array(all_ratio_priors),np.array(all_ratio_no_priors),numbers_all):
    
    for value,(x_prior,x_no_prior) in enumerate(zip(rat_prior,rat_no_prior)):
        if not np.isnan(x_prior) and not np.isinf(x_prior):
            overall_variance_prior[value]+=x_prior*num
            numbers_temp_prior[value]+=num
        if not np.isnan(x_no_prior) and not np.isinf(x_no_prior):
            overall_variance_no_prior[value]+=x_no_prior*num
            numbers_temp_no_prior[value]+=num    
overall_variance_prior/=numbers_temp_prior
overall_variance_no_prior/=numbers_temp_no_prior
fig=plt.figure()
other[0]=r'$T_{\rm{eff}}$'
other[1]=r'log\tiny{(}$g$\tiny{)}'
#create bar chart of ratio of errors for each element and have it have a black box around it
plt.bar(other+parameters,ratio_all,edgecolor='black',alpha=0.5)
plt.bar(parameters,overall_ratio,edgecolor='black',alpha=0.5)

plt.xticks(rotation=90) 
plt.ylabel('Ratio of errors or variance')
#create a horizontal line every 0.25
for i in np.arange(0,2,0.25):
    plt.axhline(y=i, color='k', linestyle='--')
#shift the x axis so it starts at the first element
plt.xlim(-0.5,len(elements)-0.5)
plt.ylim(0.5,1.25)
fig.set_size_inches(3.32088003321,3.32088003321/1.61)
plt.xticks(fontsize=8.4)
plt.tight_layout(pad=0.2)
plt.savefig("Ratio of errors and variance mask True.pdf")

plt.figure()
plt.scatter(ratio_all[4:],overall_ratio)
#fit a 1d line to the data
m, b = np.polyfit(ratio_all[4:], overall_ratio, 1)
plt.plot(ratio_all[4:], m*np.array(ratio_all)[4:] + b)
plt.xlabel('Ratio of errors')
plt.ylabel('Ratio of variance')
plt.tight_layout()
plt.savefig("Ratio of errors and variance mask True scatter.pdf")

other=['teff','logg','vmic','vsini']
# have a list of all the median errors for each element and other in prior and no prior
median_error_prior=[np.ma.median(data_all['e_'+x+'_prior']) for x in other]
median_error_no_prior=[np.ma.median(data_all['e_'+x+'_no_prior']) for x in other]

#add e_fe to the median error lists
median_error_prior.append(np.ma.median(data_all['e_fe_h_prior']))
median_error_no_prior.append(np.ma.median(data_all['e_fe_h_no_prior']))
other.append('Fe')

for x in elements:
    median_error_prior.append(np.ma.median(data_all['e_'+x+'_Fe_prior']))
    median_error_no_prior.append(np.ma.median(data_all['e_'+x+'_Fe_no_prior']))

#save a csv text file with median error, median error no prior, overall variance prior and overall variane no prior
#first add four np.nans to the variance lists to account for the other parameters
overall_variance_prior=np.insert(overall_variance_prior,0,[np.nan,np.nan,np.nan,np.nan])
overall_variance_no_prior=np.insert(overall_variance_no_prior,0,[np.nan,np.nan,np.nan,np.nan])

#the first column is the element name 
np.savetxt('errors_overall.txt',np.column_stack((other+elements,median_error_prior,median_error_no_prior,overall_variance_prior,overall_variance_no_prior)),fmt='%s')



#Replicate the same code as above but for the mean of the errors rather than variance
# numbers=[]
# all_ratios=[]
# for name in names:
#     votable = parse(name+"_no_normalization_.xml")
#     data=votable.get_first_table().to_table(use_names_over_ids=True)
#     mask=data['vsini_prior']<20
#     data=data[mask]
#     mask=data
#     mask=data['flag_sp']==0
#     data=data[mask]
#     if len(data)>1 and not name in ['ASCC_16','Melotte_22']:
#         ratio=[]
        
#         numbers.append(len(data))
#         #cycle through parameters index and take the ratio of the standard deviation of the data with the endings
#         for x in parameters_index:
#             if x in elements:
#                 ratio.append(np.nanmedian(data['e_'+x+'_Fe'+endings[0]]/data['e_'+x+'_Fe'+endings[1]]))
#             else:
#                 temp=data['e_'+x+endings[0]]/data['e_'+x+endings[1]]
#                 mask=[not np.ma.is_masked(y) for y in temp]

#                 temp=temp[mask]
#                 if len(temp)==0:
#                     ratio.append(np.nan)
#                 else:
#                     ratio.append(np.nanmedian(temp))
#         all_ratios.append(ratio)
#         #plot the ratio having parameters on the x axis and the ratio on the y axis
#         fig=plt.figure()
#         plt.scatter(parameters,ratio)
#         #plt.title('Ratio of errors  of prior and no prior data in '+name)
#         #plt.xlabel('Parameters')
#         plt.ylabel('Ratio of Errors')
#         #create vertical dashed lines for each tick
#         for xtick in plt.xticks()[0]:
#             plt.axvline(x=xtick, color='k', linestyle='--')
#         #create a horizontal line at 1
#         plt.axhline(y=1, color='k', linestyle='--')
#         plt.tight_layout()
#         plt.xticks(rotation=90) 
#         fig.set_size_inches(3.32088003321,3.32088003321/1.61)
#         plt.xticks(fontsize=7)
#         plt.tight_layout()
#         plt.savefig('Ratio_of_errors_mask_True_'+name+'.pdf')
# overall_ratio=np.zeros(len(parameters))
# for rat,num in zip(np.array(all_ratios),numbers):
#     overall_ratio+=rat*num
# overall_ratio=overall_ratio/sum(numbers)

# fig=plt.figure()
# plt.scatter(parameters,overall_ratio)
# #plt.title('overall Ratio of standard deviation of prior and no prior data')
# #plt.xlabel('Parameters')
# plt.xticks(rotation=90) 
# plt.ylabel('Ratio of standard deviation')
# for xtick in plt.xticks()[0]:
#     plt.axvline(x=xtick, color='k', linestyle='--')
# #create a horizontal line at 1
# plt.axhline(y=1, color='k', linestyle='--')    
# fig.set_size_inches(3.32088003321,3.32088003321/1.61)
# plt.xticks(fontsize=7)
# plt.tight_layout()
# print(np.mean(overall_ratio))
# plt.savefig('overall_Ratio_errors_mask_True_.pdf')