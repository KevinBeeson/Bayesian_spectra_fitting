# will read clusters.dat a interpret the first 20 bytes of each line as the name and  412-424  btyes as radius, and  593-607 as distance. discard rest of the values
import numpy as np

file = "clusters.dat"
names = []
radius = []
distance = []
other_names = []
dist_low = []
dist_high = []
rad_low = []
rad_high = []

#open the file and for each line read the first 20 bytes as the name, 412-424 as radius, and 593-607 as distance
#in the 27-272 bytes there is a comma delimited list of other names of the cluster add it
#add 577-591 as dist_low and 609-624 as dist_high
#add 412-424 as rad_low and 440-452 as rad_high
with open(file, 'r') as f:
    for line in f:
        #get rid of the spaces 
        names.append(line[0:20].strip())
        radius.append(line[426:438].strip())
        distance.append(line[593:607].strip())
        other_names.append(line[27:272].strip())
        dist_low.append(line[577:591].strip())
        dist_high.append(line[609:624].strip())
        rad_low.append(line[412:424].strip())
        rad_high.append(line[440:452].strip())




names = np.array(names)
radius = np.array(radius)
distance = np.array(distance)
dist_low = np.array(dist_low)
dist_high = np.array(dist_high)

#seperate the commas in the other names
other_names = [x.split(',') for x in other_names]
other_names = np.array(other_names)
#make into a table
from astropy.table import Table
table = Table([names, radius, distance, other_names, dist_low, dist_high, rad_low, rad_high], names=('name', 'radius', 'distance', 'other_names', 'dist_low', 'dist_high', 'rad_low', 'rad_high'))

#open dr3_clusters.txt
cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)

#cross match the table and the first column of dr3_clusters.txt use a for loop to find the index of the matching name and then use that index to find the radius and distance
cross=[]
for x in cluster_details_all[:,0]:
    mask=table['name']==x
    #if there is no match then try to use the other names
    if np.sum(mask)==0:
        for y in table['other_names']:
            if x in y:
                mask=table['other_names']==y
                break
    cross.append(table[mask]) 
#import vstack to stack the tables
from astropy.table import vstack
cross=vstack(cross)
#delete the other names column
del cross['other_names']
#save as csv   
np.savetxt('dr3_cross_match_distances_radius.csv',cross,delimiter=',',fmt='%s')