# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain all important figures that were created with PCA.

@author: joren
"""
import xarray as xr
import numpy as np
from functions.loading_and_processing_data import loadMITgcmData, detrend_and_average_MITgcmData, smoothDataset
import functions.pca as fc


#%%CREATE ANIMATIONS OF DEPTH INTEGRATED THETA
# input_filename='THETA_depth_integrated'
# var='THETA'
# members='all'
# window=1
# deseasonalize=True
# method='ens_mean'
# pcnumbers=10
# start_year=[1936, 1966, 1991]
# output_filename='integrated_THETA' 
# frames=48

#%%CREATE ANIMATIONS OF Curl TAU
input_filename='curlTAU'
var='curlTAU'
members=[0,1]
window=24
deseasonalize=True
method='ens_mean'
pcnumbers=5
start_year=[1936, 1966, 1991]
output_filename='curlTAU' 
frames=48

# #%%CREATE ANIMATIONS OF Curl TAU
# input_filename='oceTAUX'
# var='oceTAUX'
# members=[0,1]
# window=24
# deseasonalize=True
# method='ens_mean'
# pcnumbers=5
# start_year=[1936, 1966, 1991]
# output_filename='oceTAUX' 
# frames=48

# #%%CREATE ANIMATIONS OF Curl TAU
# input_filename='oceTAUY'
# var='oceTAUY'
# members=[0,1]
# window=24
# deseasonalize=True
# method='ens_mean'
# pcnumbers=5
# start_year=[1936, 1966, 1991]
# output_filename='oceTAUY' 
# frames=48

#%%CREATE ANIMATIONS OF Curl TAU
# input_filename='curlWIND'
# var='__xarray_dataarray_variable__'
# members=[0,1]
# window=24
# deseasonalize=True
# method='ens_mean'
# pcnumbers=5
# start_year=[1936, 1966, 1991]
# output_filename='curlWIND' 
# frames=48




def workflowPCA(input_filename, var, members, window, deseasonalize, method, pcnumbers, start_year=[1938, 1968, 1993], output_filename='', frames=36, maxy=-70, animate=True):

    if method=='ens_mean':
        try:
            print('...Directly loading '+input_filename+' ens_mean' )
            data=xr.open_dataset('../02_data/maps/'+input_filename+'_ensmean.nc')
            if 'ens' not in data.dims:
                data=data.expand_dims('ens')
        except:
            print('...or not :p' )
            data=loadMITgcmData(filename=input_filename, members=members)
    else:
        print('....Normally loading data')
        data=loadMITgcmData(filename=input_filename, members=members)
    
    #data=data.diff(dim='XC')
    
    print('...Detrending data' )
    data_dt=detrend_and_average_MITgcmData(data, var, window=window)
    if deseasonalize==True:
         data_dt = data_dt.groupby("time.month")-data_dt.groupby("time.month").mean(dim="time")
    
    if 'YC' in data_dt.dims:
        data_dt=data_dt.where(data_dt.YC<maxy, drop=True)
    else:
        data_dt=data_dt.where(data_dt.YG<maxy, drop=True)
    
    print('...Perform PCA')
    EOF_recons, PCs, tc=fc.performPCA(data=data_dt, var=var, method=method, interest=pcnumbers, filepath='../03_output/pca/', filename=output_filename)
    if animate==True:
        for year in start_year:
            print('...Create Animation for '+str(year))
            fc.makePCAAnimation(data_dt=data_dt, EOF_recons=EOF_recons, PCs=PCs, tc=tc, var=var, pcnumber=pcnumbers, frames=frames, start_year=year, filename=output_filename+'_'+str(year), filepath='../03_output/pca/')
    print('Finished!')
    return

workflowPCA(input_filename, var, members, window, deseasonalize, method, pcnumbers, start_year, output_filename, frames, animate=False)



import os
import matplotlib.pyplot as plt
import pandas as pd
from functions.loading_and_processing_data import read_all_data


#00-INITIALISATION------------------------------------------------------------

#First, read time series so that we can define events.
full='/data/hpcdata/users/grejan/mitgcm/data'
total={}
units={}
longnames={}
total, units, longnames=read_all_data(full, total, units, longnames, kind='isotherm')
total, units, longnames=read_all_data(full, total, units, longnames, kind='final')
total, units, longnames=read_all_data(full, total, units, longnames, kind='ohc')

cindfile='IPO_NOAA_unfiltered.txt' 
path='../data/'
cind=pd.read_csv(os.path.join(path, cindfile), names=range(1,13), skiprows=1, skipfooter=11, delim_whitespace=True)
cind=cind.replace(-99.000, np.nan)
cind=cind.replace(-99.99, np.nan)
cind=cind.stack()
cind.index=[pd.to_datetime(str(j)+'-'+str(i)) for i, j in cind.index]

total['IPO']=cind.to_frame()

#PCs=np.load('../02_data/principal_components/20220524_EXFvwind_PCA_componentens_mean_EXFvwind.npy')

PCs=np.load('../02_data/principal_components/20220524_curlTAU_PCA_componentens_mean_curlTAU.npy')
plt.figure(figsize=(20,10))
plt.suptitle('Comparison of different Principal Components of the Surface Stress Curl (all detrended)')
plt.subplot(2,1,1)
plt.plot(cind['1920':'2013'].index[12:-11], -1*PCs[:,0], c='m', label='curlTAU component 1')
plt.ylabel('Wind Stress Curl [N/m3]')
plt.grid()
plt.legend(loc='upper left')
plt.twinx()
plt.plot(cind.rolling(24, center=True).mean()['1920':'2013'], c='g', label='IPO')
plt.legend(loc='lower right')
plt.ylabel('IPO [deg C]')
plt.subplot(2,1,2)
plt.plot(cind['1920':'2013'].index[12:-11], PCs[:,2], c='m', label='curlTAU component 3')
plt.ylabel('Wind Stress Curl [N/m3]')
plt.grid()
plt.legend(loc='upper left')
plt.twinx()
plt.plot((total['dotson_to_cosgrove_massloss'].mean(axis=1).rolling(24, center=True).mean()-total['dotson_to_cosgrove_massloss'].mean(axis=1).rolling(25*12, min_periods=5*12, center=True).mean())['1920':'2013'], c='b',
        label='Dotson_to_Cosgrove_massloss')
plt.legend(loc='lower right')
plt.ylabel('Massloss [Gt/yr]')
from datetime import date
today = date.today()
today=today.strftime("%Y%m%d")
#plt.savefig('../03_output/pca/'+today+'_components_with_IPO.png')



import sys
sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid

gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
grid = Grid(gp)


EOFs=np.load('../02_data/principal_components/20220525_curlTAU_EOF_componentens_mean_curlTAU.npy', allow_pickle=True)

plt.figure(figsize=(20,10))

plt.subplot(3,1,1)
plt.pcolor(grid.lon_corners_1d[1:]+360, grid.lat_corners_1d[grid.lat_corners_1d<-70][1:], -1*(EOFs[0,:,:]-np.nanmean(EOFs[0,:,:]))/np.nanstd(EOFs[0,:,:]), vmax=1, vmin=-1, cmap='seismic')
plt.colorbar(label='Normalized Value')
bath=plt.contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
plt.clabel(bath, [-1000])
plt.ylim([-75.5, -70])
plt.title('curlTAU pattern corresponding to a positive IPO (1st PC)')

plt.subplot(3,1,2)
plt.pcolor(grid.lon_corners_1d[1:]+360, grid.lat_corners_1d[grid.lat_corners_1d<-70][1:], (EOFs[2,:,:]-np.nanmean(EOFs[2,:,:]))/np.nanstd(EOFs[2,:,:]), vmax=1, vmin=-1, cmap='seismic')
plt.colorbar(label='Normalized Value')
bath=plt.contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
plt.clabel(bath, [-1000])
plt.title('curlTAU pattern corresponding to a positive Melt (3rd PC)')
plt.ylim([-75.5, -70])


#combi=(EOFs[2,:,:]-np.nanmean(EOFs[2,:,:]))/np.nanstd(EOFs[2,:,:])+1*(EOFs[0,:,:]-np.nanmean(EOFs[0,:,:]))/np.nanstd(EOFs[0,:,:])
combi=abs((EOFs[2,:,:]-np.nanmean(EOFs[2,:,:]))/np.nanstd(EOFs[2,:,:]))-abs(-1*(EOFs[0,:,:]-np.nanmean(EOFs[0,:,:]))/np.nanstd(EOFs[0,:,:]))

plt.subplot(3,1,3)
plt.pcolor(grid.lon_corners_1d[1:]+360, grid.lat_corners_1d[grid.lat_corners_1d<-70][1:], combi, vmax=1, vmin=-1, cmap='seismic')
plt.colorbar(label='Normalized Value')
bath=plt.contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
plt.clabel(bath, [-1000])
plt.title('Difference between absolute components (3rd PC- 1st PC)')
plt.ylim([-75.5, -70])













# def plotPCACorrs(PCs, total, real_var, data_dt, 
#                  var_list=['IPO', 'amundsen_shelf_break_uwind_avg', 'amundsen_shelf_ohc_below_0m', 'dotson_to_cosgrove_massloss'], interest=10, 
#                  filename='', filepath='./figures/', method='all'):
#     #from scipy import signal
#     #print(cors)
#     #plt.figure(figsize=(15,10))
#     fig, ax=plt.subplots(len(var_list), 2, figsize=(15,10))
    
#     for i, var in enumerate(var_list):
#         total[var] = total[var].reindex(sorted(total[var].columns), axis=1)
#         ts=total[var]['1920':'2013']
#         ts=ts.rolling(24,center=True).mean()-ts.rolling(25*12, min_periods=5*12, center=True).mean()
#         ts=ts.iloc[12:-11]
#         #mean=ts.to_numpy().T.flatten()
#         print(np.shape(ts))
        
#         if method=='ens_mean':
#             ts=ts.mean(axis=1)
#             cors=np.zeros(interest)
#             optlag=np.zeros(interest)
#         else:
#             #print(len(mean))
#             cors=np.zeros((len(list(ts)), interest))
#             optlag=np.zeros((len(list(ts)), interest))
#         laglim=np.arange(-25, 26)
#         for j in range(0,interest):
#             PCnor=PCs[:,j]/max(abs(PCs[:,j]))
            
#             if method=='ens_mean':
#                 lagcoef=np.zeros(len(laglim))
#                 for k, l in enumerate(laglim):
#                     l=int(l)
#                     if l<0: #ts is leading
#                         lagcoef[k]=np.corrcoef(ts.iloc[:-abs(l)], PCnor[abs(l):])[0,1]
#                     elif l>0: #PCnor is leading
#                         lagcoef[k]=np.corrcoef(ts.iloc[abs(l):], PCnor[:-abs(l)])[0,1]
#                     elif l==0: #ts is leading
#                         lagcoef[k]=np.corrcoef(ts.iloc[:], PCnor[:])[0,1]
#                 optlag[j]=laglim[np.argmax(abs(lagcoef))]
#                 l=int(optlag[j])
#                 if l<0: #ts is leading
#                     cors[j]=np.corrcoef(ts.iloc[:-abs(l)], PCnor[abs(l):])[0,1]
#                 elif l>0: #PCnor is leading
#                     cors[j]=np.corrcoef(ts.iloc[abs(l):], PCnor[:-abs(l)])[0,1]
#                 elif l==0: #ts is leading
#                     cors[j]=np.corrcoef(ts.iloc[:], PCnor[:])[0,1]
#             elif method=='all':
#                 PCnor=PCnor.reshape(len(list(ts)), int(len(PCnor)/len(list(ts))))
#                 for ens in range(len(list(ts))):
#                     lagcoef=np.zeros(len(laglim))
#                     for k, l in enumerate(laglim):
#                         l=int(l)
#                         if l<0: #ts is leading
#                             lagcoef[k]=np.corrcoef(ts.iloc[:-abs(l),ens], PCnor[ens,abs(l):])[0,1]
#                         elif l>0: #PCnor is leading
#                             lagcoef[k]=np.corrcoef(ts.iloc[abs(l):,ens], PCnor[ens,:-abs(l)])[0,1]
#                         elif l==0: #ts is leading
#                             lagcoef[k]=np.corrcoef(ts.iloc[:,ens], PCnor[ens,:])[0,1]
#                     optlag[ens, j]=laglim[np.argmax(abs(lagcoef))]
#                     l=int(optlag[ens, j])
#                     if l<0: #ts is leading
#                         cors[ens, j]=np.corrcoef(ts.iloc[:-abs(l),ens], PCnor[ens,abs(l):])[0,1]
#                     elif l>0: #PCnor is leading
#                         cors[ens, j]=np.corrcoef(ts.iloc[abs(l):,ens], PCnor[ens,:-abs(l)])[0,1]
#                     elif l==0: #ts is leading
#                         cors[ens, j]=np.corrcoef(ts.iloc[:,ens], PCnor[ens,:])[0,1]
                
#                 #cors[ens, j]=np.corrcoef(ts.iloc[:,ens].shift(int(optlag[ens, j])), PCnor[ens,:])[0,1]
#                 #print(cors[ens, j])
            
# #             correlation = signal.correlate(mean, PCnor, mode="full")
            
# #             for l in laglim:
# #                 cors[j]=np.corrcoef(mean, PCnor)[0,1]
            
# #             optlag[j] = lags[np.argmax(correlation)]
# #             cors[j]=np.corrcoef(mean, PCnor)[0,1]
            
            
#         if method=='all':
#         #plt.subplot(len(var_list),1,i+1)
#             ax[i,0].errorbar(range(1,interest+1), np.mean(cors, axis=0), yerr=np.std(cors, axis=0), 
#                                  fmt='o-',  label='Ensemble Mean Pearson r')
#             ax[i,0].plot(range(1,interest+1), abs(np.mean(cors, axis=0)), 'o-', label='Absolute values')
#         else:
#             ax[i,0].plot(range(1,interest+1), cors, 
#                              'o-',  label='Ensemble Mean Pearson r')
#             ax[i,0].plot(range(1,interest+1), abs(cors), 'o-', label='Absolute values')
#         #ax[i,0].plot(range(1,interest+1), )
        
        
#         #plt.yticks([-1,0,1])
#         ax[i,0].grid(True)
#         if i==0:
#             fig.legend()
#         ax[i,0].set_ylabel('Pearson R')
#         ax[i,0].set_xlabel('Principal Component #')
#         ax[i,0].set_title('Correlation between: \n '+var+' and PCs from '+real_var)
#         if method=='all':
#             ax[i,1].errorbar(range(1,interest+1), np.mean(optlag, axis=0), yerr=np.std(optlag, axis=0), 
#                                  fmt='o-',  label='Optimal lag')
#         else:
#             ax[i,1].plot(range(1,interest+1), optlag, 
#                                  'o-',  label='Optimal lag')
#         ax[i,1].grid(True)
#         ax[i,1].legend()
#         ax[i,1].set_ylabel('Lag [months] \n (positive=PC leading)')
#         ax[i,1].set_xlabel('Principal Component #')
#         ax[i,1].set_title('Optimal lag between: \n '+var+' and PCs from '+real_var)
        
        
        
#     from datetime import date
#     today = date.today()
#     today=today.strftime("%Y%m%d")
    
#     plt.subplots_adjust(hspace=1)
#     plt.savefig(filepath+today+"_PCAcorrs_"+real_var+'_'+filename+".png", bbox_inches='tight')
    
#     return

# PCs=np.load('../02_data/principal_components/20220511_curlTAU_PCA_componentens_mean_curlTAU.npy')
# plotPCACorrs(PCs=PCs, total=total, real_var=var, filename='full_ensemble', data_dt='does this even matter?', method='ens_mean')
