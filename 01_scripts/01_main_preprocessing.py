# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain the exact pathway with which the data was saved that is used in the final analysis.


LACKING: 
- HOW DID I ADD THE TIMESERIES THAT KAITLIN HAS.
- HOW DID I READ THE COORDINATES ALESSANDRO.


@author: joren
"""
import numpy as np
from functions.preprocessing import readMITgcmData, readMITgcmData_depth, readMITgcmData_Transect, readMITgcmData_advection, readMITgcmData_Slice, makeEverythingSmaller
from functions.loading_and_processing_data import loadMITgcmData

#01-MAPS--------------------------------------------------------------------------------
'''
Reading MAPS:
- Zonal and Meridional Wind
- Zonal and Meridional Wind Stress
- Wind Stress Curl
- Depth Integrated THETA
'''
#---WIND-------------------------
# members='all'
# var='EXFuwind'
# #READ ZONAL WIND
# dataU=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', save=True, filename=None)
# umean=dataU.mean(dim='ens')
# umean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/EXFuwind_ensmean.nc')
# var='EXFvwind'
# #READ MERIDIONAL WIND
# dataV=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', save=True, filename=None)
# vmean=dataV.mean(dim='ens')
# vmean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/EXFvwind_ensmean.nc')
#----------------------------------------------------------------


#---Surface Wind Stress-------------------------
# members='all'
# var='oceTAUX'
# #READ ZONAL SURFACE STRESS
# dataU=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=True, filename=None)
# umean=dataU.mean(dim='ens')
# umean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/oceTAUX_ensmean.nc')
# var='oceTAUY'
# #READ MERIDIONAL SURFACE STRESS
# dataV=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=True, filename=None)
# vmean=dataV.mean(dim='ens')
# vmean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/oceTAUY_ensmean.nc')
#----------------------------------------------------------------


#---Surface Wind Stress Curl-------------------------
'''
KAN WAARSCHIJNLIJK VEEL BETER WANT DATASET HEEFT AL dyC/dyG and dxC/dxG. Dat is al de afstand in meter vgm...
'''
# members='all'
# varx='oceTAUX'
# vary='oceTAUY'
# dataU=loadMITgcmData(filename=varx, members=members, kind='maps')
# dataV=loadMITgcmData(filename=vary, members=members, kind='maps')

# R=6371000 #[m]; Radius Earth

# dudy=dataU[varx].diff('YC')/(2*np.pi*R/360*dataU['YC'].diff('YC'))
# dvdx=dataV[vary].diff('XC')/(2*np.pi*R*np.cos(np.deg2rad(dataV['YG']))*dataV['XC'].diff('XC')/360)

# dvdx['XC']=dataU.XG[1:].values
# dvdx=dvdx.rename({'XC':'XG'})

# dudy['YC']=dataV.YG[1:].values
# dudy=dudy.rename({'YC':'YG'})

# curlTAU=dvdx-dudy

# curlTAU=curlTAU.to_dataset()
# for i in range(20):
#     curlTAU.sel(ens=i).to_netcdf('../02_data/maps/curlTAU_ens'+str(i+1)+'.nc')
#curlTAU.mean(dim='ens').to_netcdf('../02_data/maps/curlTAU_ensmean.nc')
#----------------------------------------------------------------


#---Wind Curl---------------------------------------------------
# members='all'
# varx='EXFuwind'
# vary='EXFvwind'
# dataU=loadMITgcmData(filename=varx, members=members, kind='maps')
# dataV=loadMITgcmData(filename=vary, members=members, kind='maps')

# R=6371000 #[m]; Radius Earth

# dudy=dataU[varx].diff('YC')/(2*np.pi*R/360*dataU['YC'].diff('YC'))
# dvdx=dataV[vary].diff('XC')/(2*np.pi*R*np.cos(np.deg2rad(dataV['YC']))*dataV['XC'].diff('XC')/360)

# dvdx['XC']=dataU.XC[1:].values
# #dvdx=dvdx.rename({'XC':'XG'})

# dudy['YC']=dataV.YC[1:].values
# #dudy=dudy.rename({'YC':'YG'})

# curlWIND=dvdx-dudy

# curlWIND=curlWIND.to_dataset(name='curlWIND')
# for i in range(20):
#     curlWIND.sel(ens=i).to_netcdf('../02_data/maps/curlWIND_ens'+str(i+1)+'.nc')
# curlWIND.mean(dim='ens').to_netcdf('../02_data/maps/curlWIND_ensmean.nc')
#------------------------------------------------------------------------------------------------------------


#---Surface Heat Flux--------------------------------------------
# members='all'
# var='oceQnet'
# #READ SURFACE HEAT FLUX
# data=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=True, filename=None)
# mean=data.mean(dim='ens')
# mean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/oceQnet_ensmean.nc')
#----------------------------------------------------------------


#---SIheff--------------------------------------------
# members='all'
# var='SIheff'
# #READ SURFACE Snow Height
# data=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=False, filename=None)
# mean=data.mean(dim='ens')
# mean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/SIheff_ensmean.nc')
#----------------------------------------------------------------


#---SIhsnow--------------------------------------------
# members='all'
# var='SIhsnow'
# #READ SURFACE Snow Height
# data=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=False, filename=None)
# mean=data.mean(dim='ens')
# mean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/SIhsnow_ensmean.nc')
#----------------------------------------------------------------


#---EXFpress--------------------------------------------
# members='all'
# var='EXFpress'
# #READ SURFACE Snow Height
# data=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=False, filename=None)
# mean=data.mean(dim='ens')
# mean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/EXFpress_ensmean.nc')
#----------------------------------------------------------------



#--THETA--------------------------------------------------------------
# members='all'
# var='THETA'
# #READ DEPTH INTEGRATED THETA
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='integrate', members=members, z0=None, save=True, bottom=False, zbottom=100, output=False)
#----------------------------------------------------------------


#--ADV--------------------------------------------------------------
# members='all'
# var='ADVx_TH' #The next code reads both the x and y component of the advection, and it also reads both depth averaged and bottom 100 m.
# # #01-READ DEPTH INTEGRATED AND BOTTOM ADVECTION
# readMITgcmData_advection(var_name=var, members=members, bd='/data/oceans_output/shelf/kaight/mitgcm',  gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/')
#----------------------------------------------------------------


#-WVEL at 200m------------------------------------------------------------------------------------------------
# members='all'
# var='WVEL'
# #READ WVEL at 200m
# readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', save=True, filename='Z195', at_z=-195)
# mean=data.mean(dim='ens')
# mean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/Z195_WVEL_ensmean.nc')
#----------------------------------------------------------------

#-SALT at -5m------------------------------------------------------------------------------------------------
members=[2,20]
var='SALT'
#READ SALT at -5m
readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=True, filename='Z5', at_z=-5)
#mean=data.mean(dim='ens')
#mean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/Z195_WVEL_ensmean.nc')
#----------------------------------------------------------------




#---Surface Fresh Water Flux--------------------------------------------
# members='all'
# var='oceFWflx'
# #READ SURFACE HEAT FLUX
# data=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=True, filename=None)
# mean=data.mean(dim='ens')
# mean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/oceFWflx_ensmean.nc')
#----------------------------------------------------------------


# #---Sea Surface Height--------------------------------------------
# members='all'
# var='ETAN'
# #READ SURFACE HEAT FLUX
# data=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=True, filename=None)
# mean=data.mean(dim='ens')
# mean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/ETAN_ensmean.nc')
# #----------------------------------------------------------------


#02-TRANSECTS--------------------------------------------------------------------------------
'''
Reading TRANSECTS:
- PITE THETA
'''

#--THETA-PITE-------------------------------------------------------------
# members='all'
# var='THETA'
# readMITgcmData_Transect(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, 
#                            x=[-98, -105.5, -106, -106.6, -106, -102], 
#                            y=[-75.5, -73.8, -72.8, -72.5, -71.7, -71], transname='PITE_full', save=True)
#----------------------------------------------------------------



#03-SLICES-----------------------------------------------------------------------------------
'''
Reading Slices:
- U, THETA, and SALT for 115W
- U, THETA, and SALT for 123W
'''

#----U, THETA, and SALT for 115W------------------------------------------------------------
# members='all'
# var='THETA'
# readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, method='lon', x=115, y=-72.5)
# var='SALT'
# readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, method='lon', x=115, y=-72.5)
# var='UVEL's
# readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, method='lon', x=115, y=-72.5)
#----------------------------------------------------------------

#-----U, THETA, and SALT for 123W-----------------------------------------------------------
# members='all'
# var='THETA'
# readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, method='lon', x=123, y=-72.5)
# var='SALT'
# readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, method='lon', x=123, y=-72.5)
# var='UVEL'
# readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, method='lon', x=123, y=-72.5)
#----------------------------------------------------------------





# #04-Make Everything Smaller (to save space)----------------------------------------------------
makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='SALT')
#-------------------------------------------------------------------------------------------