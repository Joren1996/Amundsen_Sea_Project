# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should eventually contain the code to create all indices that are used in the final report.
Uncomment for the calculation of indices

@author: joren
"""
#00-INITIALISATION------------------------------------------------------------
from functions.loading_and_processing_data import loadMITgcmData, detrend_and_average_MITgcmData,  smoothDataset, createMapIndex, loadPACEData  #read_all_data,
from functions.composite_maps import createCompositeMap, createVelocityCompositeMap, defineEvents, findEvents, createCompositeMap_2, workflowCompositeMap
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

#%%-------01-PACE Indices---------------------------------------------------------------------------------------
'''
Some indices are precomputed here to safe time. Others are computed in xx_farFieldForcing.ipynb
'''
# members='all'
# PSL=loadPACEData(var='PSL', members=[0,1])
# SST=loadPACEData(var='SST', members=members) 
# PSL['time']=PSL.indexes['time'].to_datetimeindex()
# SST['time']=SST.indexes['time'].to_datetimeindex()

## Interpolate SST because otherwise an area average is difficult.
# def interpolateSST(SST, xgrid, ygrid, time, ens_list):
#     from scipy.interpolate import griddata
    
#     # data coordinates and values
#     x = SST.sel(ens=0).ULONG.to_numpy().flatten()
#     y = SST.sel(ens=0).ULAT.to_numpy().flatten()
    
#     # target grid to interpolate to
#     xi,yi = np.meshgrid(xgrid,ygrid)
    
#     zfinal=np.zeros((len(time), np.shape(xi)[0], np.shape(xi)[1], len(ens_list)))
#     #Perform this for every time step...
#     for j in ens_list:
#         print('Start with Member: '+str(j))
#         for i in range(len(time)):
#             print('Time step: '+str(i)+'; Member: '+str(j))
#             z = SST.isel(ens=j, time=i).to_numpy().flatten()
#             # interpolate
#             zfinal[i,:,:,j] = griddata((x,y),z,(xi,yi),method='linear')
        
#     # plot
#     fig = plt.figure()
#     ax = fig.add_subplot(211)
#     plt.pcolor(xi,yi,zfinal[0,:,:,0], cmap='seismic')
#     plt.xlabel('xi',fontsize=16)
#     plt.ylabel('yi',fontsize=16)

#     ax = fig.add_subplot(212)
#     img=ax.scatter(SST.ULONG.to_numpy().flatten(), SST.ULAT.to_numpy().flatten(), s=0.3,
#                           c=SST.isel(ens=0, time=0, z_t=0).to_numpy().flatten(), 
#                           cmap='seismic')
#     plt.xlabel('xi',fontsize=16)
#     plt.ylabel('yi',fontsize=16)
    
#     fig.savefig('../03_output/other/SST_interpolated.png')

#     SST_new=xr.DataArray(data=zfinal,
#                         dims=['time', 'lat', 'lon', 'ens'],
#                         coords=dict(
#                         lon=xgrid,
#                         lat=ygrid,
#                         time=time,
#                         ens=ens_list)                        
#                         )
#     return SST_new

# if members=='all':
#     ens_list=np.arange(20)
# else:
#     ens_list=np.arange(members[0], members[1])
# SST_new=interpolateSST(SST['SST'], xgrid=PSL.lon, ygrid=PSL.lat, time=PSL.time, ens_list=ens_list)

# #WPSST
# WPSST=SST_new.sel(lon=slice(120,160), lat=slice(0,15))
# WPSST=WPSST.groupby('time.month')-WPSST.sel(time=slice('1971', '2000')).groupby('time.month').mean()
# weights = np.cos(np.deg2rad(WPSST.lat))
# WPSST=WPSST.weighted(weights)
# WPSST=WPSST.mean(dim=['lat', 'lon'])

# WPSST.to_netcdf('../02_data/indices/WPSST_PACE.nc')

# #IPO
# region1=SST_new.sel(lon=slice(140, 360-140), lat=slice(25, 45))
# region1=region1.groupby('time.month')-region1.sel(time=slice('1971', '2000')).groupby('time.month').mean()
# weights = np.cos(np.deg2rad(region1.lat))
# region1=region1.weighted(weights)
# region1=region1.mean(dim=['lat', 'lon'])

# region2=SST_new.sel(lon=slice(170, 360-90), lat=slice(-10, 10))
# region2=region2.groupby('time.month')-region2.sel(time=slice('1971', '2000')).groupby('time.month').mean()
# weights = np.cos(np.deg2rad(region2.lat))
# region2=region2.weighted(weights)
# region2=region2.mean(dim=['lat', 'lon'])

# region3=SST_new.sel(lon=slice(150, 360-160), lat=slice(-50, -15))
# region3=region3.groupby('time.month')-region3.sel(time=slice('1971', '2000')).groupby('time.month').mean()
# weights = np.cos(np.deg2rad(region3.lat))
# region3=region3.weighted(weights)
# region3=region3.mean(dim=['lat', 'lon'])

# IPO=region2-(region1+region3)/2

# IPO.to_netcdf('../02_data/indices/IPO_PACE.nc')



'''
Below here is the computation of indices for the shelf break undercurrent and the PITE heat advection.
'''

#02-Shelf Break Undercurrent W115------------------------------------------------------------
# members='all'
# kind='slices'
# data=loadMITgcmData(filename='UVEL_W115', members=members, kind=kind)
# data=data.where(data!=0, drop=True)

# lon=115
# vmax=np.max([abs(data['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time').min()),abs(data['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time').max())])
# plt.figure(figsize=(20,10))
# ax=plt.subplot(1,1,1)
# ax.set_title('Longitudinal cross section of ensemble mean zonal velocity across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
# img=ax.pcolor(data.YC, data.Z, data['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
# plt.colorbar(img)
# ax.set_xlim([-71.8, -71.1])
# ax.set_ylim([-2500, 0])
# ax.set_ylabel('Depth [m]')
# ax.set_xlabel('Latitude [degree N]')
# ax.set_facecolor('grey')
# #plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/UVEL_W'+str(lon)+'_ensemble_mean_1979-2013.png')

# ybox=[-71.45, -71.40]
# hbox=[-550, -300]

# ax.plot([ybox[0], ybox[0], ybox[-1], ybox[-1], ybox[0]], [hbox[0], hbox[1], hbox[1], hbox[0], hbox[0]], '-ok')
# plt.show()

# index=data.sel(YC=slice(ybox[0], ybox[-1])).sel(Z=slice(hbox[-1], hbox[0]))
# index=index.mean(dim='YC').mean(dim='Z') #All Z differences are constant, so it is allowed? :@
# index.to_netcdf('../02_data/indices/undercurrent_W115.nc')

#DEPTH AVERAGED
#dataU=loadMITgcmData('depth_averaged_UVEL', members='all', kind='old')
#ybox=[-71.45, -71.40]
#index=dataU.sel(XG=360-115, YC=slice(ybox[0], ybox[-1])).mean(dim=['XG', 'YC'])
#index.to_netcdf('../02_data/indices/undercurrent_depth_averaged_W115.nc')

#03-PITE and Coastal Path: Advection and Velocity------------------------------------------------------------
# members='all'
# kind='maps'
# data=loadMITgcmData(filename='ADVy_TH_corrected_depth_averaged', members=members, kind=kind)
# #data=loadMITgcmData(filename='ADVy_TH_corrected_bottom100m_averaged', members=members, kind=kind)

# plt.figure()
# var='ADVy_TH'
# xcor=[m for m in data[var].coords if 'X' in m][0]
# ycor=[m for m in data[var].coords if 'Y' in m][0]

# #data=data.where(data!=0, drop=True)
# ax=plt.subplot(1,1,1)
# img=ax.pcolor(data[var][xcor], data[var][ycor], data[var].sel(ens=0).sel(time='1940').mean(dim='time'), cmap='seismic')
# plt.colorbar(img, label='[m/s]', ax=ax)
# ax.set_facecolor('grey')
# ax.set_ylim([-75.5, -70])
# #Cross sections for influx:
# #PITE
# x=[360-107,360-105]
# y=[-73, -73]
# plt.plot(x,y,'o--r')

# def find_nearest(array, value):
#     array = np.asarray(array)
#     idx = (np.abs(array - value)).argmin()
#     return array[idx]


# ynew=find_nearest(data[var].YG, y[0])
# #     flow[trough]=data_dt['VVEL'].sel(XC=slice(troughs[trough][0], troughs[trough][1]), YG=ynew).mean(dim='XC')#, skipna=True)
# index=data.sel(YG=ynew, XC=slice(x[0], x[-1]))
# index=index.mean(dim='XC')
# index.to_netcdf('../02_data/indices/'+var+'_PITE.nc')

# #Coastal
# x=[360-105,360-103]
# y=[-73, -73]
# plt.plot(x,y,'o--y')

# ynew=find_nearest(data[var].YG, y[0])
# index=data.sel(YG=ynew, XC=slice(x[0], x[-1]))
# index=index.mean(dim='XC')
# index.to_netcdf('../02_data/indices/'+var+'_coastal.nc')














