# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should eventually contain the code to create all indices that are used in the final report.

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


#--BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS------
#--BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS------
#--BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS------
#--BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS------
#--BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS------

























#---PITE FLOW------------------------------------------------------------------------------------------------------------

# #Create Index for PITE from Cross Section
# members='all'
# dataV=loadMITgcmData(filename='VVEL_S72-5', members=members, kind='slices')

# dataVmean=dataV['VVEL'].mean(dim=['ens', 'time'])
# dataVmean=dataVmean.where(dataVmean!=0)

# plt.figure(figsize=(20,10))
# ax=plt.subplot(1,1,1)
# img=plt.pcolor(dataV.XC, dataV.Z, dataVmean, cmap='seismic', vmin=-0.1, vmax=0.1)
# plt.colorbar(img, label='Velocity [m/s]')
# plt.xlim([252, 254])
# plt.ylim([-900, 0])
# plt.ylabel('Depth [m]')
# plt.xlabel('Longitude [degree E]')
# ax.set_facecolor('grey')

# ybox=[253, 253.5]
# hbox=[-750, 0]
# ax.plot([ybox[0], ybox[0], ybox[-1], ybox[-1], ybox[0]], [hbox[0], hbox[1], hbox[1], hbox[0], hbox[0]], '--og')
# hbox=[-750, -300]
# ax.plot([ybox[0], ybox[0], ybox[-1], ybox[-1], ybox[0]], [hbox[0], hbox[1], hbox[1], hbox[0], hbox[0]], '-ok')
# plt.show()

# import sys
# sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
# from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
# gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
# grid = Grid(gp)


# #Average based on overdepth (Kaitlin Utils)

# #Bottom Flow
# area=dataV['VVEL'].sel(Z=slice(-300, -800), XC=slice(253, 253.5))
# selz=area.Z
# area=area.to_numpy()

# idy=np.where(np.round(grid.lat_corners_1d, 4).astype(np.float32)==dataV.YG.values)[0][0]
# idx1=np.where((grid.lon_1d>253.0-360) & (grid.lon_1d<253.5-360))[0][0]
# idx2=np.where((grid.lon_1d>253.0-360) & (grid.lon_1d<253.5-360))[0][-1]
# idz=[i in selz for i in grid.z]

# area=np.sum(area*grid.dz[idz].reshape(1,1,sum(idz),1)*np.moveaxis(grid.hfac_s[idz, idy:idy+1, idx1:idx2+1],[0], [1]), axis=(2,3))
# area=area/np.sum(grid.dz[idz].reshape(1,1,sum(idz),1)*np.moveaxis(grid.hfac_s[idz, idy:idy+1, idx1:idx2+1],[0], [1]))
# PITE=xr.DataArray(data=area, dims=["ens", "time"], coords=dict(ens=dataV.ens, time=dataV.time))

# PITE.to_netcdf('../02_data/indices/VVEL_PITE_bottom_final.nc')

# #Barotropic Flow
# area=dataV['VVEL'].sel(Z=slice(0, -800), XC=slice(253, 253.5))
# selz=area.Z
# area=area.to_numpy()

# idy=np.where(np.round(grid.lat_corners_1d, 4).astype(np.float32)==dataV.YG.values)[0][0]
# idx1=np.where((grid.lon_1d>253.0-360) & (grid.lon_1d<253.5-360))[0][0]
# idx2=np.where((grid.lon_1d>253.0-360) & (grid.lon_1d<253.5-360))[0][-1]
# idz=[i in selz for i in grid.z]

# area=np.sum(area*grid.dz[idz].reshape(1,1,sum(idz),1)*np.moveaxis(grid.hfac_s[idz, idy:idy+1, idx1:idx2+1],[0], [1]), axis=(2,3))
# area=area/np.sum(grid.dz[idz].reshape(1,1,sum(idz),1)*np.moveaxis(grid.hfac_s[idz, idy:idy+1, idx1:idx2+1],[0], [1]))
# PITE_trop=xr.DataArray(data=area, dims=["ens", "time"], coords=dict(ens=dataV.ens, time=dataV.time))

# PITE_trop.to_netcdf('../02_data/indices/VVEL_PITE_barotropic_final.nc')

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

# ybox=[-71.46, -71.39]
# hbox=[-550, -300]

# ax.plot([ybox[0], ybox[0], ybox[-1], ybox[-1], ybox[0]], [hbox[0], hbox[1], hbox[1], hbox[0], hbox[0]], '-ok')

# ybox=[-71.46, -71.39]
# hbox=[-1500, 0]
# ax.plot([ybox[0], ybox[0], ybox[-1], ybox[-1], ybox[0]], [hbox[0], hbox[1], hbox[1], hbox[0], hbox[0]], '--og')
# plt.show()


# #Average based on overdepth (Kaitlin Utils)

# #Bottom Flow
# ybox=[-71.46, -71.39]
# hbox=[-550, -300]
# area=data['UVEL'].sel(Z=slice(hbox[-1], hbox[0]), YC=slice(ybox[0],ybox[-1]))
# selz=area.Z
# area=area.to_numpy()

# idy=np.where(360+np.round(grid.lon_corners_1d, 4).astype(np.float32)==data.XG.values)[0][0]
# idx1=np.where((grid.lat_1d>ybox[0]) & (grid.lat_1d<ybox[-1]))[0][0]
# idx2=np.where((grid.lat_1d>ybox[0]) & (grid.lat_1d<ybox[-1]))[0][-1]
# idz=[i in selz for i in grid.z]

# area=np.sum(area*grid.dz[idz].reshape(1,1,sum(idz),1)*np.moveaxis(grid.hfac_w[idz, idy:idy+1, idx1:idx2+1],[0], [1]), axis=(2,3))
# area=area/np.sum(grid.dz[idz].reshape(1,1,sum(idz),1)*np.moveaxis(grid.hfac_w[idz, idy:idy+1, idx1:idx2+1],[0], [1]))
# UW115=xr.DataArray(data=area, dims=["ens", "time"], coords=dict(ens=data.ens, time=data.time))

# UW115.to_netcdf('../02_data/indices/UW115_bottom_final.nc')


# #Barotropic Flow
# ybox=[-71.46, -71.39]
# hbox=[-1500, 0]
# area=data['UVEL'].sel(Z=slice(hbox[-1], hbox[0]), YC=slice(ybox[0],ybox[-1]))
# selz=area.Z
# area=area.to_numpy()

# idy=np.where(360+np.round(grid.lon_corners_1d, 4).astype(np.float32)==data.XG.values)[0][0]
# idx1=np.where((grid.lat_1d>ybox[0]) & (grid.lat_1d<ybox[-1]))[0][0]
# idx2=np.where((grid.lat_1d>ybox[0]) & (grid.lat_1d<ybox[-1]))[0][-1]
# idz=[i in selz for i in grid.z]

# area=np.sum(area*grid.dz[idz].reshape(1,1,sum(idz),1)*np.moveaxis(grid.hfac_w[idz, idy:idy+1, idx1:idx2+1],[0], [1]), axis=(2,3))
# area=area/np.sum(grid.dz[idz].reshape(1,1,sum(idz),1)*np.moveaxis(grid.hfac_w[idz, idy:idy+1, idx1:idx2+1],[0], [1]))
# UW115_trop=xr.DataArray(data=area, dims=["ens", "time"], coords=dict(ens=data.ens, time=data.time))

# UW115_trop.to_netcdf('../02_data/indices/UW115_barotropic_final.nc')

#---------------------------------------------------------------------------------------------------------------


#--BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS------
#--BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS------
#--BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS------
#--BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS------
#--BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS--------BONUS------


#----------------------------------------------------------------------------------------------------------------
#XX-ENSMEAN BAROCLINIC AND BAROTROPIC
# dataU=loadMITgcmData('depth_averaged_UVEL', members='all', kind='old')
# # ensmean=dataU.mean(dim='ens')
# # ensmean.to_netcdf('./02_data/maps/depth_averaged_UVEL_ensmean.nc')

# dataV=loadMITgcmData('depth_averaged_VVEL', members='all', kind='old')
# ensmean=dataV.mean(dim='ens')
# ensmean.to_netcdf('./02_data/maps/depth_averaged_VVEL_ensmean.nc')

# databottomU=loadMITgcmData('bottom100m_averaged_UVEL', members='all', kind='old')
# databottomV=loadMITgcmData('bottom100m_averaged_VVEL', members='all', kind='old')

# baroclinicU=databottomU-dataU
# baroclinicU_ensmean=baroclinicU.mean(dim='ens')
# baroclinicU.to_netcdf('./02_data/maps/baroclinic_UVEL_ensmean.nc')

# baroclinicV=databottomV-dataV
# baroclinicV_ensmean=baroclinicV.mean(dim='ens')
# baroclinicV.to_netcdf('./02_data/maps/baroclinic_VVEL_ensmean.nc')

# def computeFullField(dataU, dataV, datavar='UVEL'):
#     if 'XG' in dataU.dims:
#         dataU=dataU.rename({'XG':'XC'})
#     if 'YG' in dataU.dims:
#         dataU=dataU.rename({'YG':'YC'})
#     if 'XG' in dataV.dims:
#         dataV=dataV.rename({'XG':'XC'})
#     if 'YG' in dataV.dims:
#         dataV=dataV.rename({'YG':'YC'})

#     #Read something with the correct locations.
#     data=loadMITgcmData(filename='EXFuwind', members=[0,1], kind='maps')
#     dataU=dataU.interp(XC=data.XC, YC=data.YC)
#     dataV=dataV.interp(XC=data.XC, YC=data.YC)

#     print('...Compute the absolute combined value')
#     datavar_v=datavar.replace('U','V')
#     datavar_v=datavar_v.replace('u','v')
#     datavar_v=datavar_v.replace('x','y')
#     data=np.sqrt(dataU[datavar]**2+dataV[datavar_v]**2)
#     data=data.to_dataset(name=datavar)
#     return data


# baroclinic=computeFullField(dataU=baroclinicU, dataV=baroclinicV, datavar='UVEL')
# baroclinic=baroclinic.mean(dim='ens')
# baroclinic.to_netcdf('../02_data/maps/baroclinic_VEL_ensmean.nc')

# barotropic=computeFullField(dataU=dataU, dataV=dataV, datavar='UVEL')
# barotropic=barotropic.mean(dim='ens')
# barotropic.to_netcdf('../02_data/maps/depth_averaged_VEL_ensmean.nc')


# #01-Shelf Break Undercurrent W115------------------------------------------------------------
# # members='all'
# # kind='slices'
# # data=loadMITgcmData(filename='UVEL_W115', members=members, kind=kind)
# # data=data.where(data!=0, drop=True)

# # lon=115
# # vmax=np.max([abs(data['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time').min()),abs(data['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time').max())])
# # plt.figure(figsize=(20,10))
# # ax=plt.subplot(1,1,1)
# # ax.set_title('Longitudinal cross section of ensemble mean zonal velocity across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
# # img=ax.pcolor(data.YC, data.Z, data['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
# # plt.colorbar(img)
# # ax.set_xlim([-71.8, -71.1])
# # ax.set_ylim([-2500, 0])
# # ax.set_ylabel('Depth [m]')
# # ax.set_xlabel('Latitude [degree N]')
# # ax.set_facecolor('grey')
# # #plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/UVEL_W'+str(lon)+'_ensemble_mean_1979-2013.png')

# # ybox=[-71.45, -71.40]
# # hbox=[-550, -300]

# # ax.plot([ybox[0], ybox[0], ybox[-1], ybox[-1], ybox[0]], [hbox[0], hbox[1], hbox[1], hbox[0], hbox[0]], '-ok')
# # plt.show()

# # index=data.sel(YC=slice(ybox[0], ybox[-1])).sel(Z=slice(hbox[-1], hbox[0]))
# # index=index.mean(dim='YC').mean(dim='Z') #All Z differences are constant, so it is allowed? :@
# # index.to_netcdf('../02_data/indices/undercurrent_W115.nc')

#DEPTH AVERAGED
#dataU=loadMITgcmData('depth_averaged_UVEL', members='all', kind='old')
#ybox=[-71.45, -71.40]
#index=dataU.sel(XG=360-115, YC=slice(ybox[0], ybox[-1])).mean(dim=['XG', 'YC'])
#index.to_netcdf('../02_data/indices/undercurrent_depth_averaged_W115.nc')

# #02-PITE and Coastal Path: Advection and Velocity------------------------------------------------------------
# # members='all'
# # kind='maps'
# # data=loadMITgcmData(filename='ADVy_TH_corrected_depth_averaged', members=members, kind=kind)
# # #data=loadMITgcmData(filename='ADVy_TH_corrected_bottom100m_averaged', members=members, kind=kind)

# # plt.figure()
# # var='ADVy_TH'
# # xcor=[m for m in data[var].coords if 'X' in m][0]
# # ycor=[m for m in data[var].coords if 'Y' in m][0]

# # #data=data.where(data!=0, drop=True)
# # ax=plt.subplot(1,1,1)
# # img=ax.pcolor(data[var][xcor], data[var][ycor], data[var].sel(ens=0).sel(time='1940').mean(dim='time'), cmap='seismic')
# # plt.colorbar(img, label='[m/s]', ax=ax)
# # ax.set_facecolor('grey')
# # ax.set_ylim([-75.5, -70])
# # #Cross sections for influx:
# # #PITE
# # x=[360-107,360-105]
# # y=[-73, -73]
# # plt.plot(x,y,'o--r')

# # def find_nearest(array, value):
# #     array = np.asarray(array)
# #     idx = (np.abs(array - value)).argmin()
# #     return array[idx]


# # ynew=find_nearest(data[var].YG, y[0])
# # #     flow[trough]=data_dt['VVEL'].sel(XC=slice(troughs[trough][0], troughs[trough][1]), YG=ynew).mean(dim='XC')#, skipna=True)
# # index=data.sel(YG=ynew, XC=slice(x[0], x[-1]))
# # index=index.mean(dim='XC')
# # index.to_netcdf('../02_data/indices/'+var+'_PITE.nc')

# # #Coastal
# # x=[360-105,360-103]
# # y=[-73, -73]
# # plt.plot(x,y,'o--y')

# # ynew=find_nearest(data[var].YG, y[0])
# # index=data.sel(YG=ynew, XC=slice(x[0], x[-1]))
# # index=index.mean(dim='XC')
# # index.to_netcdf('../02_data/indices/'+var+'_coastal.nc')

# #03 SAME BUT BOTTOM 100m----------------------------------------------------
# # members='all'
# # kind='maps'
# # #data=loadMITgcmData(filename='ADVy_TH_corrected_depth_averaged', members=members, kind=kind)
# # data=loadMITgcmData(filename='ADVy_TH_corrected_bottom100m_averaged', members=members, kind=kind)

# # plt.figure()
# # var='ADVy_TH'
# # xcor=[m for m in data[var].coords if 'X' in m][0]
# # ycor=[m for m in data[var].coords if 'Y' in m][0]

# # #data=data.where(data!=0, drop=True)
# # ax=plt.subplot(1,1,1)
# # img=ax.pcolor(data[var][xcor], data[var][ycor], data[var].sel(ens=0).sel(time='1940').mean(dim='time'), cmap='seismic')
# # plt.colorbar(img, label='[m/s]', ax=ax)
# # ax.set_facecolor('grey')
# # ax.set_ylim([-75.5, -70])
# # #Cross sections for influx:
# # #PITE
# # x=[360-107,360-105]
# # y=[-73, -73]
# # plt.plot(x,y,'o--r')

# # def find_nearest(array, value):
# #     array = np.asarray(array)
# #     idx = (np.abs(array - value)).argmin()
# #     return array[idx]


# # ynew=find_nearest(data[var].YG, y[0])
# # #     flow[trough]=data_dt['VVEL'].sel(XC=slice(troughs[trough][0], troughs[trough][1]), YG=ynew).mean(dim='XC')#, skipna=True)
# # index=data.sel(YG=ynew, XC=slice(x[0], x[-1]))
# # index=index.mean(dim='XC')
# # index.to_netcdf('../02_data/indices/'+var+'_PITE_bottom100m.nc')

# # #Coastal
# # x=[360-105,360-103]
# # y=[-73, -73]
# # plt.plot(x,y,'o--y')

# # ynew=find_nearest(data[var].YG, y[0])
# # index=data.sel(YG=ynew, XC=slice(x[0], x[-1]))
# # index=index.mean(dim='XC')
# # index.to_netcdf('../02_data/indices/'+var+'_coastal_bottom100m.nc')



# #04-Wind Stress Curl
# # members='all'
# # kind='maps'
# # data=loadMITgcmData(filename='curlTAU', members=members, kind=kind)
# # data=smoothDataset(data, var='curlTAU', which_area='rAz', x=11, y=11)
# # dataind=createMapIndex(data, 'curlTAU', xlim=[250, 252], ylim=[-71.4, -71.2], which_area='rAz')
# # dataind.to_netcdf('../02_data/indices/curlTAU_smoothed.nc')




# # members='all'
# # kind='maps'
# # data=loadMITgcmData(filename='curlTAU', members=members, kind=kind)
# # ind=data['curlTAU'].sel(YG=slice(-71.4, -71.25), XG=slice(250, 252)).mean(dim='YG').mean(dim='XG')
# # ind.to_netcdf('../02_data/indices/curlTAU_for_undercurrent_new.nc')

# # members='all'
# # kind='maps'
# # data=loadMITgcmData(filename='curlWIND', members=members, kind=kind)
# # ind=data['curlWIND'].sel(YC=slice(-71.4, -71.25), XC=slice(250, 252)).mean(dim='YC').mean(dim='XC')
# # ind.to_netcdf('../02_data/indices/curlWIND_for_undercurrent.nc')



# #04-Barotropic Index W115------------------------------------------------------------
# ''' 
# This has to be done properly!!!!!!!!
# '''
# # members='all'
# # kind='slices'
# # data=loadMITgcmData(filename='THETA_W115', members=members, kind=kind)
# # data=data.where(data!=0, drop=True)

# # lon=115
# # vmax=np.max([abs(data['THETA'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time').min()),abs(data['THETA'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time').max())])
# # plt.figure(figsize=(20,10))
# # ax=plt.subplot(1,1,1)
# # ax.set_title('Longitudinal cross section of ensemble mean THETA across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
# # img=ax.pcolor(data.YC, data.Z, data['THETA'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
# # plt.colorbar(img)
# # ax.set_xlim([-71.8, -71.1])
# # ax.set_ylim([-2500, 0])
# # ax.set_ylabel('Depth [m]')
# # ax.set_xlabel('Latitude [degree N]')
# # ax.set_facecolor('grey')
# # #plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/UVEL_W'+str(lon)+'_ensemble_mean_1979-2013.png')

# # # ybox=[-71.45, -71.40]
# # # hbox=[-550, -300]

# # # ax.plot([ybox[0], ybox[0], ybox[-1], ybox[-1], ybox[0]], [hbox[0], hbox[1], hbox[1], hbox[0], hbox[0]], '-ok')
# # plt.show()


# # index=data.sel(YC=slice(-71.4, -71.3), Z=slice(0, -500)).mean(dim='Z').mean(dim='YC')-data.sel(YC=slice(-71.6, -71.5), Z=slice(0, -500)).mean(dim='Z').mean(dim='YC')
# # #All Z differences are constant, so it is allowed? :@
# # index.to_netcdf('../02_data/indices/tempgradient_W115_temp.nc')

# #Proper Barotropic Index--------------------------------------------------------------------------------------------------------------------
# #%%IMPORTING
# # import matplotlib.colors as colors

# # divnorm = colors.TwoSlopeNorm(vcenter=0)

# # plt.rcParams.update({'font.size': 14})


# # #01 MAKE CROSS SECTIONS A LA ALESSANDRO
# # lon=115

# # dataS=loadMITgcmData(filename='SALT_W'+str(lon), members='all', kind='slices')
# # dataT=loadMITgcmData(filename='THETA_W'+str(lon), members='all', kind='slices')

# # #Convert everything to density.
# # import sys
# # sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
# # from mitgcm_python_master.diagnostics import density
# # ref_depth=0
# # eosType='MDJWF'
# # dataS=dataS.where(dataS!=0, drop=True)
# # dataT=dataT.where(dataT!=0, drop=True)

# # salt=dataS['SALT'].to_numpy()
# # temp=dataT['THETA'].to_numpy()
# # q=density (eosType, salt, temp, press=ref_depth)
# # q=xr.DataArray(q, dims=('ens', 'time', "Z", "YC"), coords={'ens': dataS['SALT'].ens, 'time': dataS['SALT'].time, "YC": dataS['SALT'].YC, "Z": dataS['SALT'].Z})


# # drhody=q.differentiate(coord='YC')/1851.8 #[m per latitude degree]
# # drhodz=q.differentiate(coord='Z')
# # dzdy=drhody/drhodz

# # #-----------------------------------------------------------------------VISUALISATION
# # #ENSEMBLE MEAN

# # ##01 Ensemble Mean Density Gradient
# # plt.figure(figsize=(20,10))
# # ax=plt.subplot(1,1,1)
# # ax.set_title('Longitudinal cross section of ensemble mean density gradient \n across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
# # img=ax.pcolor(dzdy.YC, dzdy.Z, dzdy.mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time'), cmap='jet', vmin=-2, vmax=2)
# # plt.colorbar(img)
# # #ax.set_xlim([-71.8, -71.1])
# # ax.set_xlim([-71.8, -71.1])
# # ax.set_ylim([-2500, 0])
# # ax.set_ylabel('Depth [m]')
# # ax.set_xlabel('Latitude [degree N]')
# # ax.set_facecolor('grey')
# # plt.plot([-71.45, -71.4], [-215, -215], 'o-k')
# # plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/baroclinic_rhoslope_'+str(lon)+'_ensemble_mean_1979-2013.png')
# # plt.show()


# # plt.figure(figsize=(20,10))
# # ax=plt.subplot(1,1,1)
# # ax.set_title('Longitudinal Cross Section of Sine of Slope of Ensemble Mean Density Gradient \n across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
# # img=ax.pcolor(dzdy.YC, dzdy.Z, np.sin(np.arctan(dzdy.mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time'))), cmap='jet', vmin=-1, vmax=1)
# # plt.colorbar(img)
# # #ax.set_xlim([-71.8, -71.1])
# # ax.set_xlim([-71.8, -71.1])
# # ax.set_ylim([-2500, 0])
# # ax.set_ylabel('Depth [m]')
# # ax.set_xlabel('Latitude [degree N]')
# # ax.set_facecolor('grey')
# # plt.plot([-71.45, -71.4], [-215, -215], 'o-k')
# # plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/baroclinic_sineofrhoslope_'+str(lon)+'_ensemble_mean_1979-2013.png')
# # plt.show()

# # #CREATE THE INDEX!
# # index=np.sin(np.arctan(dzdy.sel(YC=slice(-71.45, -71.4), Z=slice(-210, -220)).mean(dim='YC').mean(dim='Z'))) #Take Index at dzdy to confuse confusion with average of angle....
# # index.to_netcdf('../02_data/indices/baroclinicity_W115.nc')
# # #CREATE THE INDEX!--> JUST THE SLOPE
# # index=dzdy.sel(YC=slice(-71.45, -71.4), Z=slice(-210, -220)).mean(dim='YC').mean(dim='Z') #Take Index at dzdy to confuse confusion with average of angle....
# # index.to_netcdf('../02_data/indices/baroclinicity_slope_W115.nc')

# # #Create Lower Baroclinicity Index
# # index=dzdy.sel(YC=slice(-71.41, -71.37), Z=slice(-550, -560)).mean(dim='YC').mean(dim='Z') #Take Index at dzdy to confuse confusion with average of angle....
# # index.to_netcdf('../02_data/indices/baroclinicity_slope_lower_W115.nc')


# # p=total['amundsen_shelf_break_uwind_avg']['1920':'2013']
# # q=adv['UW115'].to_dataframe()['UVEL'].unstack(0)['1920':'2013']
# # q.columns=['ens'+str(i+1).zfill(2) for i in range(20)]

# # p=(p-p.mean())/p.std()
# # q=(q-q.mean())/q.std()

# # diff=p.rolling(24, center=True).mean()-q.rolling(24, center=True).mean()   


# # xrdf=xr.DataArray(data=diff.values, dims=["time","ens"],

# #     coords=dict(

# #         time=diff.index, ens=diff.columns))
# # xrdf.to_netcdf('./02_data/indices/normalized_diff_wind_undercurrent.nc')

# # index=xr.open_dataarray('../02_data/indices/normalized_diff_wind_undercurrent.nc')
# # dzdy['time']=dzdy.indexes['time'].to_datetimeindex()
# # res=xr.corr(dzdy.mean(dim='ens'), index.mean(dim='ens'), dim='time')
# # plt.figure(figsize=(20,10))
# # ax=plt.subplot(1,1,1)
# # ax.set_title('Longitudinal cross section of ensemble mean density gradient correlation map \n across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
# # img=ax.pcolor(res.YC, res.Z, res, cmap='jet')
# # plt.colorbar(img)
# # cs=ax.contour(res.YC, res.Z, res, np.arange(-1, 1, 0.1))
# # ax.clabel(cs, cs.levels, inline=True, fontsize=10)
# # ax.set_xlim([-74.5, -70])
# # ax.set_ylim([-2500, 0])
# # ax.set_ylabel('Depth [m]')
# # ax.set_xlabel('Latitude [degree N]')
# # ax.set_facecolor('grey')
# # plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/baroclinic_rhoslope_correlationmap_'+str(lon)+'_ensemble_mean_1979-2013.png')
# # plt.plot([-71, -70.8], [-555, -555], 'o-k')
# # plt.show()


# #BAROTROPIC BAROCLINIC:
# print('Everything is fine! Let us go!')
# members='all'
# kind='old'
# # dataTrop=loadMITgcmData(filename='depth_averaged_UVEL', members=members, kind=kind)
# dataBot=loadMITgcmData(filename='bottom100m_averaged_UVEL', members=members, kind=kind)

# # tropleft=dataTrop.sel(XG=slice(360-115.1, 360-114.9), YC=slice(-71.6, -71.3)).mean(dim=['YC', 'XG'])
# # tropleft.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/indices/tropleft.nc')

# # botleft=dataBot.sel(XG=slice(360-115.1, 360-114.9), YC=slice(-71.6, -71.3)).mean(dim=['YC', 'XG'])
# # botleft.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/indices/botleft.nc')

# # tropleft['time']=tropleft.indexes['time'].to_datetimeindex()
# # botleft['time']=botleft.indexes['time'].to_datetimeindex()
# # clinleft=botleft-tropleft
# # clinleft.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/indices/clinleft.nc')

# # tropright=dataTrop.sel(XG=slice(360-108.1, 360-107.9), YC=slice(-71.2, -70.9)).mean(dim=['YC', 'XG'])
# # tropright.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/indices/tropright.nc')

# botright=dataBot.sel(XG=slice(360-108.1, 360-107.9), YC=slice(-71.2, -70.9)).mean(dim=['YC', 'XG'])
# botright.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/indices/botright.nc')

# # tropright['time']=tropright.indexes['time'].to_datetimeindex()
# # botright['time']=botright.indexes['time'].to_datetimeindex()
# # clinright=botright-tropright
# # clinright.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/indices/clinright.nc')

# print('Finished!')








