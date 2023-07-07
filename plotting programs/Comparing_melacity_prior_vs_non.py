#this program will cycle through the open clusters and tell me their fe/h prior and non prior content
#it will first load the names of the clusters from the dr3_clusters.txt file
#it will then load the data from the xml files
#then in a for loop print the name of the cluster and the fe/h prior and non prior content

import numpy as np
import matplotlib.pyplot as plt
import csv
from astropy.io.votable import parse,from_table,writeto
from astropy.table import Table,vstack,join

#load the names of the clusters from the dr3_clusters.txt file
cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)

metalicty_prior=[]
metalicty_no_prior=[]
#load the data from the xml files
for name in cluster_details_all[:,0]:
    votable = parse(name+"_no_normalization_.xml")
    data=votable.get_first_table().to_table(use_names_over_ids=True)
    data['cluster_name']=[name for x in range(len(data))]
    print(name)
    print(np.mean(data['fe_h_prior']))
    print(np.mean(data['fe_h_no_prior']))
    metalicty_prior.append(np.mean(data['fe_h_prior']))
    metalicty_no_prior.append(np.mean(data['fe_h_no_prior']))
meltalicty_isochrone=cluster_details_all[:,5].astype(np.float64)
metalicty_no_prior=np.array(metalicty_no_prior)
metalicty_prior=np.array(metalicty_prior)
plt.figure()
plt.scatter(metalicty_no_prior,metalicty_no_prior-metalicty_prior)
#draw y=x line use the scatter values for the min and max
x=np.linspace(min(metalicty_prior)-0.1,max(metalicty_prior)+0.1,100)
plt.plot(x,x)
plt.xlabel('fe/h no prior')
plt.ylabel('fe/h prior')
#have x lim and

plt.figure()
plt.scatter(metalicty_prior,meltalicty_isochrone)
plt.scatter(metalicty_no_prior,meltalicty_isochrone)
#plot y=x line
x=np.linspace(min(metalicty_prior)-0.1,max(metalicty_prior)+0.5,100)
plt.plot(x,x)
#plot the 1d line of best fit for both lines
fit=np.polyfit(metalicty_prior,meltalicty_isochrone,1)
fit_fn=np.poly1d(fit)
plt.plot(x,fit_fn(x),c='blue')
fit=np.polyfit(metalicty_no_prior,meltalicty_isochrone,1)
fit_fn=np.poly1d(fit)
plt.plot(x,fit_fn(x),c='orange')

plt.legend(['prior','no prior',])