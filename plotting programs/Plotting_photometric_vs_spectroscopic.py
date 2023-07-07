import numpy as np
import matplotlib.pyplot as plt
import csv
from astropy.io.votable import parse
# import names of clusters from dr3_clusters.txt first row  

cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)
print(cluster_details_all[:,0])
names=cluster_details_all[:,0]


# for name in names:
name='Melotte_22'
print(name)
votable = parse(name+"_bp_rp_run.xml")
data=votable.get_first_table().to_table(use_names_over_ids=True)
plt.figure()
plt.scatter(data['phot_g_mean_mag'],data['bp_rp'],c=data['teff'])
plt.title(name)

#checks if probability grid is inside the colnames of the data  table
if 'probability_grid' in data.colnames:
    print('probability grid in data',name)