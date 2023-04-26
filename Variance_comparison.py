# This fille will import all the data from the reduced open clusters and compare the variance between prior and no prior data 
import numpy as np
import matplotlib.pyplot as plt
import csv
from astropy.io.votable import parse

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

endings=['_prior','_no_prior']
numbers=[]
all_ratios=[]
for name in names:
    votable = parse(name+"_no_normalization_.xml")
    data=votable.get_first_table().to_table(use_names_over_ids=True)
    mask=data['vsini_prior']<20
    data=data[mask]
    mask=data['flag_sp']==0
    data=data[mask]
    if len(data)>1:
        ratio=[]
        numbers.append(len(data))
        #cycle through parameters index and take the ratio of the standard deviation of the data with the endings
        for x in parameters_index:
            if x in elements:
                ratio.append(np.std(data[x+'_Fe'+endings[0]])/np.std(data[x+'_Fe'+endings[1]]))
            else:
                ratio.append(np.std(data[x+endings[0]])/np.std(data[x+endings[1]]))
        all_ratios.append(ratio)
        #plot the ratio having parameters on the x axis and the ratio on the y axis
        plt.figure()
        plt.scatter(parameters,ratio)
        plt.title('Ratio of standard deviation of prior and no prior data in '+name)
        plt.xlabel('Parameters')
        plt.ylabel('Ratio of standard deviation')
        #create vertical dashed lines for each tick
        for xtick in plt.xticks()[0]:
            plt.axvline(x=xtick, color='k', linestyle='--')
        plt.tight_layout()
        plt.xticks(rotation=90) 
        plt.savefig('Ratio_of_variance_'+name+'.pdf')
overall_ratio=np.zeros(len(parameters))
for rat,num in zip(np.array(all_ratios),numbers):
    overall_ratio+=rat*num
overall_ratio=overall_ratio/sum(numbers)

plt.figure()
plt.scatter(parameters,overall_ratio)
plt.title('overall Ratio of standard deviation of prior and no prior data')
plt.xlabel('Parameters')
plt.xticks(rotation=90) 
plt.ylabel('Ratio of standard deviation')
for xtick in plt.xticks()[0]:
    plt.axvline(x=xtick, color='k', linestyle='--')
plt.savefig('overall_Ratio.pdf')

# #Replicate the same code as above but for the mean of the errors rather than variance
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
#     if len(data)>1:
#         ratio=[]
#         numbers.append(len(data))
#         #cycle through parameters index and take the ratio of the standard deviation of the data with the endings
#         for x in parameters_index:
#             if x in elements:
#                 ratio.append(np.nanmedian(data['e_'+x+'_Fe'+endings[0]]/data['e_'+x+'_Fe'+endings[1]]))
#             else:
#                 ratio.append(np.nanmedian(data['e_'+x+endings[0]]/data['e_'+x+endings[1]]))
#         all_ratios.append(ratio)
#         #plot the ratio having parameters on the x axis and the ratio on the y axis
#         plt.figure()
#         plt.scatter(parameters,ratio)
#         plt.title('Ratio of errors  of prior and no prior data in '+name)
#         plt.xlabel('Parameters')
#         plt.ylabel('Ratio of Errors')
#         #create vertical dashed lines for each tick
#         for xtick in plt.xticks()[0]:
#             plt.axvline(x=xtick, color='k', linestyle='--')
#         plt.tight_layout()
#         plt.xticks(rotation=90) 
#         plt.savefig('Ratio_of_errors_'+name+'.pdf')
# overall_ratio=np.zeros(len(parameters))
# for rat,num in zip(np.array(all_ratios),numbers):
#     overall_ratio+=rat*num
# overall_ratio=overall_ratio/sum(numbers)

# plt.figure()
# plt.scatter(parameters,overall_ratio)
# plt.title('overall Ratio of standard deviation of prior and no prior data')
# plt.xlabel('Parameters')
# plt.xticks(rotation=90) 
# plt.ylabel('Ratio of standard deviation')
# for xtick in plt.xticks()[0]:
#     plt.axvline(x=xtick, color='k', linestyle='--')
# plt.savefig('overall_Ratio_errors.pdf')