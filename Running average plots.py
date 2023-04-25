from numba import jit

import numpy as np
import matplotlib.pyplot as plt
import csv
from astropy.io.votable import parse,from_table,writeto
from scipy.stats import kde
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kde
from pathlib import Path
import requests
import re
from os.path import exists
import csv
from scipy.interpolate import interp1d
def get_isochrone(low_age,high_age,low_metalicty,high_metalicty,a_v,age_spacing=0.01,metalicty_spacing=0.01,iso_type='gaia'):
    """
    asks the padava server for an isochrone or a list of isochrones

    Parameters
    ----------
    low_age : float
        DESCRIPTION.
    high_age : float
        DESCRIPTION.
    low_metalicty : float
        DESCRIPTION.
    high_metalicty : float
        DESCRIPTION.
    a_v : float
        extinction.
    age_spacing : float, optional
        DESCRIPTION. The default is 0.01.
    metalicty_spacing : float, optional
        DESCRIPTION. The default is 0.01.

    an array of the isochrone
    -------
    None.

    """
    #if metallicity is given as a log ([M/H]) convert it here:


#    mass=[]
#    label=[]
#    mags=[]
#    imf=[]
#    teff=[]
#    logg=[]

    #parameters other than default
    d={
    'track_parsec': 'parsec_CAF09_v1.2S',
    'track_colibri':'parsec_CAF09_v1.2S_S35',
    'track_postagb':'no',
    'n_inTPC': 10,
    'eta_reimers': 0.2,
     # 'photsys_file': 'tab_mag_odfnew/tab_mag_gaiaDR2weiler.dat',
     'photsys_file':'tab_mag_odfnew/tab_mag_gaiaEDR3.dat',
   #   'photsys_file':'tab_mag_odfnew/tab_mag_gaia.dat',
  # 'photsys_file':'tab_mag_odfnew/tab_mag_2mass_spitzer_wise.dat',
    # 'photsys_file':'tab_mag_odfnew/tab_mag_panstarrs1.dat',
     #'photsys_file':'tab_mag_odfnew/tab_mag_gaiaDR2.dat',
    'photsys_version': 'OBC',
    'dust_sourceM': 'nodustM',
    'dust_sourceC': 'nodustC',
    'extinction_av': a_v,
    'extinction_coeff':'constant',
    'extinction_curve':'cardelli',   
    'imf_file': 'tab_imf/imf_kroupa_orig.dat',
    'isoc_isagelog':'1',
    'isoc_lagelow':low_age,
    'isoc_lageupp':high_age,
    'isoc_dlage':age_spacing, #steps ages
    'isoc_ismetlog':'1',
    'isoc_metlow':low_metalicty,
    'isoc_metupp':high_metalicty,
    'isoc_dmet':metalicty_spacing, #steps M/H
    'output_kind': 0,
    'submit_form': 'Submit'}
    if iso_type=='gaia':
        d['photsys_file']='tab_mag_odfnew/tab_mag_gaiaEDR3.dat'
    elif iso_type=='panstarrs':
        d['photsys_file']='tab_mag_odfnew/tab_mag_panstarrs1.dat'
    elif iso_type=='allwise':
        d['photsys_file']='tab_mag_odfnew/tab_mag_2mass_spitzer_wise.dat'
    #Check if we already downloaded this isochrone.
    #Isochrones are saved as txt files and the filename is the hash of the dictionary values.
    webserver = 'http://stev.oapd.inaf.it'
    c = requests.get(webserver + '/cgi-bin/cmd_3.7', params=d).text
#    print(c)
    aa = re.compile('output\d+')
    fname = aa.findall(c)
    if len(fname) > 0:
        url = '{0}/tmp/{1}.dat'.format(webserver, fname[0])
        #print url
        r = requests.get(url).text
        print('successfully gotten the ischrone')
    return r
def get_and_save_ishochrones(name,low_age,high_age,low_metalicty,high_metalicty,a_v,age_spacing=0.01,metalicty_spacing=0.01,special=False,iso_type='gaia'):
    """
    You can ask the function for a isochrone it will check if we already have it saved and if not it will ask from the padava servers

    Parameters
    ----------
    name : str
        name of the cluster.
    low_age : TYPE
        DESCRIPTION.
    high_age : TYPE
        DESCRIPTION.
    low_metalicty : TYPE
        DESCRIPTION.
    high_metalicty : TYPE
        DESCRIPTION.
    a_v : TYPE
        DESCRIPTION.
    age_spacing : TYPE, optional
        DESCRIPTION. The default is 0.01.
    metalicty_spacing : TYPE, optional
        DESCRIPTION. The default is 0.01.
    special : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    Path('ALL_ISO').mkdir(parents=True,exist_ok=True)
    #gets the path name of the isochrone and checks if its already been downloaded

    if (low_age==high_age or high_age==None) and (low_metalicty==high_metalicty or high_metalicty==None):
        if special:
            path_name='ALL_ISO/'+name+'_'+special+'_age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)
        else:
            path_name='ALL_ISO/age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)
    elif low_metalicty==high_metalicty or high_metalicty==None:
        path_name='ALL_ISO/age_'+str(low_age)+'_to_'+str(high_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)

    elif low_age==high_age or high_age==None:
        path_name='ALL_ISO/age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_to_'+str(high_metalicty)+'_extinction_'+str(a_v)


    else:
        
        if special:
            path_name='ALL_ISO/'+name+'_'+special+'_mixed'
        else:
            path_name='ALL_ISO/mixed'   
        if iso_type!='gaia':
            path_name+='_'+iso_type
        path_name+='.txt'
        #asks for the isochorne
        data=get_isochrone(low_age,high_age,low_metalicty,high_metalicty,a_v,age_spacing,metalicty_spacing,iso_type=iso_type)
        data=data.split('\n')
        header=data[11]
        data=data[0:len(data)-1]
        ammendedData=''
        for x in data:
            if x[0]!='#':
                numbers=np.array(x.split(' '))
                numbers=numbers[numbers!='']
                ammendedData+=' '.join(numbers)+'\n'
            else:
                numbers=np.array(x.split(' '))
                numbers=numbers[numbers!='']
                ammendedData+=' '.join(numbers)+'\n'
                
        ammendedData=ammendedData.split('\n')
        ammendedData=[i.split(' ') for i in ammendedData]
        ammendedData=ammendedData[0:len(ammendedData)-1]
        # ammendedData=np.array(ammendedData)
    
        toWrite=''
        for x in ammendedData:
            line=''
            for y in x:
                line+=str(y)+' '
            line=line[:-1]
            toWrite+=line+'\n'    
        file=open(path_name,'w')
        file.write(toWrite)
    if iso_type!='gaia':
        path_name+='_'+iso_type
    path_name+='.txt'
    if not exists(path_name):
        
        data=get_isochrone(low_age,high_age,low_metalicty,high_metalicty,a_v,iso_type=iso_type)
        data=data.split('\n')
        header=data[11]
        data=data[0:len(data)-1]
        ammendedData=''
        for x in data:
            if x[0]!='#':
                numbers=np.array(x.split(' '))
                numbers=numbers[numbers!='']
                ammendedData+=' '.join(numbers)+'\n'
            else:
                numbers=np.array(x.split(' '))
                numbers=numbers[numbers!='']
                ammendedData+=' '.join(numbers)+'\n'
                
        ammendedData=ammendedData.split('\n')
        ammendedData=[i.split(' ') for i in ammendedData]
        ammendedData=ammendedData[0:len(ammendedData)-1]
        # ammendedData=np.array(ammendedData)
    
        toWrite=''
        for x in ammendedData:
            line=''
            for y in x:
                line+=str(y)+' '
            line=line[:-1]
            toWrite+=line+'\n'    
        file=open(path_name,'w')
        file.write(toWrite)
def interpolate(oldIsochrone,scale=5):   
    mass=oldIsochrone[:,0]
    #changes log temperature to temp.
    oldIsochrone[:,2]=10**(oldIsochrone[:,2])

    
    newIso=[[None for y in range(len(oldIsochrone[0]))] for x in range((len(oldIsochrone)-1)*scale)]
    newIso=np.array(newIso)
    newIso[:,0]=np.hstack(np.transpose(np.linspace(mass[:-1],mass[1:],num=scale,endpoint=True)))
    newIso=np.array(newIso,dtype=float)
    for x in range(len(oldIsochrone[0])-1):
        f = interp1d(mass,oldIsochrone[:,x+1])
        newIso[:,x+1]=f(newIso[:,0]) 
    newIso[:,2]=np.log10(newIso[:,2])
def iso_reader(name,low_age,low_metalicty,a_v,special=False,high_age=None, high_metalicty=None,iso_type='gaia'):
    if (low_age==high_age or high_age==None) and (low_metalicty==high_metalicty or high_metalicty==None):
        if special:
            path_name='ALL_ISO/'+name+'_'+special+'_age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)
        else:
            path_name='ALL_ISO/age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)
        if iso_type!='gaia':
            path_name+='_'+iso_type
        path_name+='.txt'
        iso = list(csv.reader(open(path_name, 'rt'), delimiter=' '))
        
        iso=iso[13:-1]
        iso=np.array(iso,dtype=float)
        #crops the isochrone down to initial ini mass, int_imf, log_T, G, G_Bp_bright, G_BP_faint, G_RP ,log_g , cur mass
        if iso_type=='gaia':
            iso=np.column_stack((iso[:,3],iso[:,4],iso[:,7],iso[:,25],iso[:,26],iso[:,27],iso[:,8],iso[:,5]))
            
                
            iso=[x for x in iso if x[3]<15]
        return [low_age,low_metalicty],np.vstack(iso)
    elif low_metalicty==high_metalicty or high_metalicty==None:
        path_name='ALL_ISO/age_'+str(low_age)+'_to_'+str(high_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)+'.txt'

    elif low_age==high_age or high_age==None:
        path_name='ALL_ISO/age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_to_'+str(high_metalicty)+'_extinction_'+str(a_v)+'.txt'

    else:
        if special:
            path_name='ALL_ISO/'+name+'_'+special+'_mixed.txt'
        else:
            path_name='ALL_ISO/mixed.txt'
    iso = list(csv.reader(open(path_name, 'rt'), delimiter=' '))
    iso_full=[]
    iso_parameters=[]
    iso_temp=[]
    
    for x in iso:
        if x[0]!='#':
            iso_temp.append(x)
        else:
            if len(iso_temp):
                iso_temp=np.array(iso_temp,dtype=float)
                age_temp=iso_temp[0][2]
                metalicity_temp=iso_temp[0][1]

                iso_temp=np.column_stack((iso_temp[:,3],iso_temp[:,4],iso_temp[:,7],iso_temp[:,25],iso_temp[:,26],iso_temp[:,27],iso_temp[:,8],iso_temp[:,5]))
                iso_temp=[x for x in iso_temp if x[3]<15]

                iso_full.append(np.vstack(iso_temp))
                iso_parameters.append([age_temp,metalicity_temp])
                iso_temp=[]
    return iso_parameters,iso_full
class isochrone:
    def __init__(self,name,age,metalicity,extinction,special=False,interpolate_scale=10,high_age=None,high_metalicity=None,limits=None,iso_type='gaia'):
        self.name=name
        self.special=special
        self.extinction=extinction
        self.iso_type=iso_type
        get_and_save_ishochrones(name,age,high_age,metalicity,high_metalicity,extinction,special=special,iso_type=iso_type)
        self.parameters,self.isochrone=iso_reader(name, age, metalicity, extinction,special=special,high_age=high_age,high_metalicty=high_metalicity,iso_type=iso_type)
        if iso_type=='allwise':
            iso_temp=isochrone(name,age,metalicity,extinction,special,interpolate_scale,high_age,high_metalicity,limits,iso_type='gaia')
            self.iso_gaia=iso_temp.isochrone
        if limits:
            if len(limits)!=2:
                raise ValueError("wrong shape for limits")
            self.limit_g=limits
        if limits:
            if isinstance(self.isochrone, np.ndarray):
                self.isochrone=np.vstack([x for x in self.isochrone if x[3]-limits[0]>0 and x[3]-limits[1]<0])  

            else:
                iso_temp=[]
                for x in self.isochrone:
                    iso_temp.append(np.vstack([x for x in isochrone if x[3]-limits[0]>0 and x[3]-limits[1]<0]))
                self.isochrone=iso_temp
        if isinstance(self.isochrone, np.ndarray) and iso_type=='gaia' and interpolate_scale!=1:
            self.isochrone=interpolate(self.isochrone,interpolate_scale)
        elif iso_type=='gaia'and interpolate_scale!=1:
            iso_temp=[]
            for x in self.isochrone:
                iso_temp.append(interpolate(x,interpolate_scale))
            self.isochrone=iso_temp
    def plot(self):
        plt.figure()
        if isinstance(self.isochrone,np.ndarray):
            plt.plot(self.isochrone[:,5]-self.isochrone[:,4],self.isochrone[:,3])
        else:
            for (x,label) in zip(self.isochrone,self.parameters):
                temp_label='age '+ str(label[0]) +' metalicty '+str(label[1])
            
                plt.plot(x[:,4]-x[:,5],x[:,3],label=temp_label)
                plt.legend(loc='best')
        plt.xlabel('bp_rp')
        plt.ylabel('G')
        plt.gca().invert_yaxis()
def converged_data(tags,tag_2=None):
    masks=data['flag_'+tags]<=convergence
    if not tag_2 is None:
        masks_2=data['flag_'+tag_2]<=convergence
        masks*masks_2
    data_temp=data[masks]
    if not tag_2 is None:
        return data_temp[tags],data_temp[tag_2]
    return data_temp[tags]





plt.rcParams['font.family']='Times New Roman'
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['text.usetex'] = True
cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)
print(cluster_details_all[:,0])
# name=input('what isochrone do you want?')

convergence=1
parameters=['teff','logg','Fe/H','vmic','vbroad','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
parameters_index=['teff','logg','fe_h','vmic','vsini','Li','C','N','O','Na','Mg','Al','Si','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Rb','Sr','Y','Zr','Mo','Ru','Ba','La','Ce','Nd','Sm','Eu']
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
circle_size=0.1
colours=['red','blue']
endings=['_prior','_Galah_iDR4']
converged=True
for name in cluster_details_all[:,0]:
    shape=(8,5)

    fig = plt.figure(figsize=shape)
    # fig.set_figwidth(shape[1])
    # fig.set_figheight(shape[0])
    ax=[]
    ax1=plt.subplot2grid(shape=shape, loc=(0, 0), colspan=3,rowspan=2)
    ax=[]
    for row in range(2):
        for col in range(3,shape[1]):
            ax.append(plt.subplot2grid(shape=shape,loc=(row,col)))
    for row in range(2,shape[0]):
        for col in range(shape[1]):
            ax.append(plt.subplot2grid(shape=shape,loc=(row,col)))
    cluster_1=cluster_2=name
    file_name_1=cluster_1+'_no_normalization_.xml'
    file_name_2=cluster_2+'_no_normalization_.xml'
    for file_name,colour,ending,cluster_name in zip([file_name_1,file_name_2],colours,endings,[cluster_1,cluster_2]):
        cluster_details=[x for x in cluster_details_all if x[0]==cluster_name][0]

        low_age=float(cluster_details[1])
        low_metalicty=float(cluster_details[2])
        low_extinction=float(cluster_details[3])
        best_age=float(cluster_details[4])
        best_metalicty=float(cluster_details[5])
        best_extinction=float(cluster_details[6])
        high_age=float(cluster_details[7])
        high_metalicty=float(cluster_details[8])
        high_extinction=float(cluster_details[9])

        interpolate_scale=1
        iso=isochrone(file_name,best_age,best_metalicty,best_extinction,special='Best_fit',high_age=best_age,high_metalicity=best_metalicty,iso_type='gaia',interpolate_scale=1)
        iso=iso.isochrone
        # NGC_2682_summary_dr61
        votable = parse(file_name)
        data=votable.get_first_table().to_table(use_names_over_ids=True)
        mask=data['vsini_prior']<20
        data=data[mask]
        if converged:
            mask=data['flag_sp']==0
            data=data[mask]
        if len(data)!=0:
            size_font=10
            iso_type='gaia'




            
            if ending=="_Galah_iDR4":
                ax1.scatter(data['teff_spectroscopic'],data['logg_spectroscopic'],s=0.3,c=colour)
            else:
                teff,logg=converged_data('teff'+ending,'logg'+ending)
                ax1.scatter(teff,logg,s=0.3,c=colour)
            ax1.set_xlabel(r'$T_{\rm{eff}}$')
            ax1.xaxis.set_label_position('top') 
            ax1.tick_params('x', top=True, labeltop=True,bottom=False,labelbottom=False)
            ax1.set_xlim(min(min(data['teff_prior']),min(data['teff_no_prior'])),max(max(data['teff_prior']),max(data['teff_no_prior'])))
            ax1.set_ylim(min(min(data['logg_prior']),min(data['logg_no_prior'])),max(max(data['logg_prior']),max(data['logg_no_prior'])))

            # x=np.array(data['teff'],dtype=float)
            # y=np.array(data['logg'],dtype=float)
            # nbins=300
            # k = kde.gaussian_kde([x,y])
            # xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
            # zi = k(np.vstack([xi.flatten(), yi.flatten()]))
            # ax1.contour(xi, yi, zi.reshape(xi.shape),levels=5)
            ax1.plot(10**iso[:,2],iso[:,-2],label=r'Best fit',c=colour,linewidth=0.8)

            ax1.text(-0.125 ,0.5,r'log(\small $g$)' ,horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes,rotation=90)

            ax1.invert_yaxis()
            ax1.invert_xaxis()
            for value,(axis,col_dat) in enumerate(zip(ax[:3],parameters_index[2:5]),2):
                if ending=="_Galah_iDR4":
                    axis.scatter(data['teff_spectroscopic'],data[col_dat.lower()],s=0.3,c=colour)
                    bins=np.linspace(min(data['teff_spectroscopic'])+200,max(data['teff_spectroscopic'])-200,10)
                else:
                    axis.scatter(data['teff'+ending],data[col_dat+ending],s=circle_size,marker=',',c=colour)
                    bins=np.linspace(min(data['teff'+ending])-200,    max(data['teff'+ending])+200,10)
                average=np.array([])
                std=np.array([])
                for value in range(len(bins)-1):
                    if ending=="_Galah_iDR4":
                        average=np.append(average,np.mean(data[col_dat.lower()][np.where((data['teff_spectroscopic']>bins[value]) & (data['teff_spectroscopic']<bins[value+1]))]))
                        std=np.append(std,np.std(data[col_dat.lower()][np.where((data['teff_spectroscopic']>bins[value]) & (data['teff_spectroscopic']<bins[value+1]))]))
                    else:
                        average=np.append(average,np.mean(data[col_dat+ending][np.where((data['teff_prior']>bins[value]) & (data['teff_prior']<bins[value+1]))]))
                        std=np.append(std,np.std(data[col_dat+ending][np.where((data['teff_prior']>bins[value]) & (data['teff_prior']<bins[value+1]))]))
                axis.plot((bins[:-1]+bins[1:])/2,average,linewidth=0.5,c=colour)
                axis.fill_between((bins[:-1]+bins[1:])/2,average-std,average+std,alpha=0.5,color=colour)
            
            for value,(axis,col_dat) in enumerate(zip(ax[3:],parameters_index[5:]),5):
                if ending=="_Galah_iDR4":
                    axis.scatter(data['teff_spectroscopic'],data[col_dat.lower()],s=0.3,c=colour)
                    bins=np.linspace(min(data['teff_spectroscopic'])-200,max(data['teff_spectroscopic'])+200,10)
                else:
                    axis.scatter(data['teff'+ending],data[col_dat+ending],s=circle_size,marker=',',c=colour)
                    bins=np.linspace(min(data['teff'+ending])-200,max(data['teff'+ending])+200,10)
                average=np.array([])
                std=np.array([])
                for value in range(len(bins)-1):
                    if ending=="_Galah_iDR4":
                        average=np.append(average,np.mean(data[col_dat.lower()][np.where((data['teff_spectroscopic']>bins[value]) & (data['teff_spectroscopic']<bins[value+1]))]))
                        std=np.append(std,np.std(data[col_dat.lower()][np.where((data['teff_spectroscopic']>bins[value]) & (data['teff_spectroscopic']<bins[value+1]))]))
                    else:
                        average=np.append(average,np.mean(data[col_dat+ending][np.where((data['teff_prior']>bins[value]) & (data['teff_prior']<bins[value+1]))]))
                        std=np.append(std,np.std(data[col_dat+ending][np.where((data['teff_prior']>bins[value]) & (data['teff_prior']<bins[value+1]))]))
                axis.plot((bins[:-1]+bins[1:])/2,average,linewidth=0.5,c=colour)
                axis.fill_between((bins[:-1]+bins[1:])/2,average-std,average+std,alpha=0.5,color=colour)

        #put names in the plots 
        for value,(axis,col_dat) in enumerate(zip(ax[3:],parameters_index[5:]),5):
            axis.text(0.15,0.7,parameters[value],horizontalalignment='center',verticalalignment='center', transform=axis.transAxes,c='black',size=size_font,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5))
        ax[0].text(0.2,0.7,parameters[2],horizontalalignment='center',verticalalignment='center', transform=ax[0].transAxes,c='black',size=size_font,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5))
        ax[1].text(0.2,0.7,parameters[3],horizontalalignment='center',verticalalignment='center', transform=ax[1].transAxes,c='black',size=size_font,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5))
        ax[2].text(0.2,0.7,parameters[4],horizontalalignment='center',verticalalignment='center', transform=ax[2].transAxes,c='black',size=size_font,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.5))

        for axis in ax[-5:]:
            axis.set_xlabel(r'$T_{\rm{eff}}$')
        for axis in ax[3:4]:
            axis.set_xlabel(r'$T_{\rm{eff}}$')
        for value,axis in enumerate(ax[4:]):
            if value%shape[1]==0:
                axis.text(-0.45,0.5,r'[x/Fe]' ,horizontalalignment='center',verticalalignment='center', transform=axis.transAxes,rotation=90)
        red_patch = mpatches.Patch(color='red', label=(cluster_1+endings[0]).replace("_", " "))
        blue_patch = mpatches.Patch(color='blue', label=(cluster_2+endings[1]).replace("_", " "))
        fig.legend(handles=[red_patch,blue_patch],loc=(0.71,0.905))
        plt.tight_layout(w_pad=0.0,h_pad=-1.2)
        fig.set_size_inches(6.97384806974,6.97384806974/7.6*8/1.6)
        plt.savefig('/home/kevin/Documents/Paper/new/'+cluster_1+endings[0]+'_'+cluster_2+endings[1]+'masks_'+str(converged)+'_running_average.pdf',bbox_inches='tight')