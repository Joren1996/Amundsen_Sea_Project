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
from functions.preprocessing import readMITgcmData, readMITgcmData_depth, readMITgcmData_Transect, readMITgcmData_advection, readMITgcmData_Slice, makeEverythingSmaller, readGEOData, readGEOData_depth, readGEOData_Slice
from functions.loading_and_processing_data import loadMITgcmData

#%%01-MAPS--------------------------------------------------------------------------------
'''
Reading MAPS:
- Zonal and Meridional Wind
- Surface Freshwater Flux
- Salinity (0-200m)
- Bottom 100m Averaged Velocity
- Depth Averaged Velocity
- Bottom 100m Averaged and Depth Averaged Heat Advection

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

#---Surface Fresh Water Flux--------------------------------------------
# members='all'
# var='oceFWflx'
# #READ SURFACE HEAT FLUX
# data=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=True, filename=None)
#----------------------------------------------------------------

#---Salinity-----------------------------------------------------
# var='SALT'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', 
#                      z0=[-200, 0], bottom=False, zbottom=100, save=True, output=False)
#----------------------------------------------------------------

#---Bottom 100m Current-------------------------------------------
# var='UVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=None, bottom=True, zbottom=100, save=True, output=False, gtype='u')

# var='VVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=None, bottom=True, zbottom=100, save=True, output=False, gtype='v')

#-----------------------------------------------------------------

#---Depth Averaged Current-------------------------------------------
# var='UVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=None, bottom=False, zbottom=100, save=True, output=False, gtype='u')

# var='VVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=None, bottom=False, zbottom=100, save=True, output=False, gtype='v')

#-----------------------------------------------------------------

#--ADV------------------------------------------------------------
# members='all'
# var='ADVx_TH' #The next code reads both the x and y component of the advection, and it also reads both depth averaged and bottom 100 m.
# # #01-READ DEPTH INTEGRATED AND BOTTOM ADVECTION
# readMITgcmData_advection(var_name=var, members=members, bd='/data/oceans_output/shelf/kaight/mitgcm',  gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/')
#----------------------------------------------------------------

# makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='all')


'''
Maps for New Geometry
- Salinity (200-700m)
- Temperature (200-700m)
- Bottom 100m Velocity
'''


#---SHIfwFlx--------------------------------------------

#for mem in [4,6,10,11,12]:
#     members=[mem-1,mem]
#     var='SHIfwFlx'
#     #READ Shelf FW flx
#     data=readGEOData(var, bd='/data/oceans_output/shelf/grejan/archer2_mitgcm_GEO', members=members, save=True, filename=None)
#     #READ Shelf FW flx
#     data=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=True, filename=None)

#---Temperature (200-700m)----------------------------------------
# for mem in [4,6,10,11,12]:
#     members=[mem-1,mem]
#     var='THETA'
#     data=readGEOData_depth(var, method='average', members=members, z0=[-700, -200], bottom=False, zbottom=100, save=True, output=False)
#     data=readMITgcmData_depth(var, method='average', members=members, z0=[-700, -200], bottom=False, zbottom=100, save=True, output=False)

#---Salinity (200-700m)----------------------------------------
# for mem in [4,6,10,11,12]:
#     members=[mem-1,mem]
#     var='SALT'
#     data=readGEOData_depth(var, method='average', members=members, z0=[-700, -200], bottom=False, zbottom=100, save=True, output=False)
#     data=readMITgcmData_depth(var, method='average', members=members, z0=[-700, -200], bottom=False, zbottom=100, save=True, output=False)
#----------------------------------------------------------

#---BOTTOM-100M-VELOCITY----------------------------------------    
# for mem in [4,6,10,11,12]:
#     members=[mem-1,mem]
#     var='UVEL'
#     readGEOData_depth(var, bd='/data/oceans_output/shelf/grejan/archer2_mitgcm_GEO', method='average', members=members, z0=None, bottom=True, zbottom=100, save=True, output=False, gtype='u')
    
#     var='VVEL'
#     readGEOData_depth(var, bd='/data/oceans_output/shelf/grejan/archer2_mitgcm_GEO', method='average', members=members, z0=None, bottom=True, zbottom=100, save=True, output=False, gtype='v')





#02-SLICES-----------------------------------------------------------------------------------
'''
Reading Slices:
- U, THETA, and SALT for 115W
- U, THETA, and SALT for 72.5S
'''

#-----U, THETA, and SALT for 115W-----------------------------------------------------------
# members='all'
# var='THETA'
# readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, method='lon', x=115, y=-72.5)
# var='SALT'
# readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, method='lon', x=115, y=-72.5)
# var='UVEL'
# readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, method='lon', x=115, y=-72.5)
#----------------------------------------------------------------

#----U, THETA, and SALT for 72.5S------------------------------------------------------------

# members='all'
# var='THETA'
# readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, method='lat', x=115, y=-72.5)
# var='SALT'
# readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, method='lat', x=115, y=-72.5)
# var='VVEL'
# readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, method='lat', x=115, y=-72.5)

# for mem in [4,6,10,11,12]:
#     members=[mem-1,mem]
#     var='THETA'
#     readGEOData_Slice(var, members=members, method='lat', x=115, y=-72.5)
#     var='SALT'
#     readMITgcmData_Slice(var, members=members, method='lat', x=115, y=-72.5)
#     var='VVEL'
#     readMITgcmData_Slice(var, members=members, method='lat', x=115, y=-72.5)
    
#----------------------------------------------------------------

# makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/slices/', newtype=np.float32, kind='all')



#03-TIME-SERIES
'''
Time series used for analysis: most are just taken from Naughten et al. Here, 

'''

#----DERIVED FROM NAUGHTEN ET AL (2022)--------------------------
# import sys
# import os
# sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
# from mitgcm_python_master.postprocess import precompute_timeseries_coupled

# segment_dir=[i for i in os.listdir('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output') if (i[0]=='1') | (i[0]=='2')]
# segment_dir=[i for i in segment_dir if (int(i[:4])>=1920)]
##segment_dir=[i for i in segment_dir if (int(i[:4])<=1995)]
# print(segment_dir)
# for i in range(1,21):
#     print('Start with member: '+str(i))
#     precompute_timeseries_coupled (output_dir='/data/oceans_output/shelf/grejan/archer2_mitgcm_GEO/PAS_{}_GEO/output'.format(str(i).zfill(2)), timeseries_file='/data/hpcdata/users/grejan/mitgcm/02_data/timeseries/timeseries_full_GEO{}.nc'.format(str(i).zfill(2)), 
#                                    hovmoller_file='hovmoller.nc', file_name='output.nc', segment_dir=segment_dir, timeseries_types=['amundsen_shelf_break_uwind_avg', 'all_massloss', 'dotson_to_cosgrove_massloss', 'dotson_massloss', 'pig_massloss', 'thwaites_massloss', 'cosgrove_massloss', 'crosson_massloss', 'amundsen_shelf_ohc_below_0m', 'amundsen_shelf_temp_btw_200_700m', 'amundsen_shelf_salt_below_0m', 'all_melting', 'dotson_to_cosgrove_melting', 'dotson_melting', 'pig_melting', 'thwaites_melting', 'cosgrove_melting', 'crosson_melting', 'amundsen_shelf_seaice_melt', 'amundsen_shelf_seaice_freeze'], 
#                                    hovmoller_loc=[], key='PAS', time_average=False)

#     print('Finished with member: '+str(i))
# print('Finished with *EVERYTHING*!')

#----NEW GEOMETRY-------------------------------------------------
# import sys
# import os
# sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
# from mitgcm_python_master.postprocess import precompute_timeseries_coupled

# segment_dir=[i for i in os.listdir('/data/oceans_output/shelf/grejan/archer2_mitgcm_GEO/PAS_06_GEO/output') if (i[0]=='1') | (i[0]=='2')]
# segment_dir=[i for i in segment_dir if (int(i[:4])>=1920)]
# segment_dir=[i for i in segment_dir if (int(i[:4])<=1995)]
# print(segment_dir)
# for i in [4,6,10,11,12]:
#     print('Start with member: '+str(i))
#     precompute_timeseries_coupled (output_dir='/data/oceans_output/shelf/grejan/archer2_mitgcm_GEO/PAS_{}_GEO/output'.format(str(i).zfill(2)), timeseries_file='/data/hpcdata/users/grejan/mitgcm/02_data/timeseries/timeseries_full_GEO{}.nc'.format(str(i).zfill(2)), 
#                                    hovmoller_file='hovmoller.nc', file_name='output.nc', segment_dir=segment_dir, timeseries_types=['amundsen_shelf_break_uwind_avg', 'all_massloss', 'dotson_to_cosgrove_massloss', 'dotson_massloss', 'pig_massloss', 'thwaites_massloss', 'cosgrove_massloss', 'crosson_massloss', 'amundsen_shelf_ohc_below_0m', 'amundsen_shelf_temp_btw_200_700m', 'amundsen_shelf_salt_below_0m', 'all_melting', 'dotson_to_cosgrove_melting', 'dotson_melting', 'pig_melting', 'thwaites_melting', 'cosgrove_melting', 'crosson_melting', 'amundsen_shelf_seaice_melt', 'amundsen_shelf_seaice_freeze'], 
#                                    hovmoller_loc=[], key='PAS', time_average=False)

#     print('Finished with member: '+str(i))
# print('Finished with *EVERYTHING*!') 


#--------BONUS----------------------------#--------BONUS----------------------------#--------BONUS----------------------------
#--------BONUS----------------------------#--------BONUS----------------------------#--------BONUS----------------------------
#--------BONUS----------------------------#--------BONUS----------------------------#--------BONUS----------------------------
#--------BONUS----------------------------#--------BONUS----------------------------#--------BONUS----------------------------
#--------BONUS----------------------------#--------BONUS----------------------------#--------BONUS----------------------------








# makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='UVEL')

# for mem in [4,6,10,11,12]:
#     members=[mem-1,mem]

#makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='VVEL')

# var='VVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=None, bottom=True, zbottom=100, save=True, output=False, gtype='v')

# makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='VVEL_bottom')

# var='VVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=None, bottom=False, zbottom=100, save=True, output=False, gtype='v')

# makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='VVEL_depth')

# var='UVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=None, bottom=False, zbottom=100, save=True, output=False, gtype='u')


# members='all'
# var='ADVx_TH' #The next code reads both the x and y component of the advection, and it also reads both depth averaged and bottom 100 m.
# # #01-READ DEPTH INTEGRATED AND BOTTOM ADVECTION
# readMITgcmData_advection(var_name=var, members=members, bd='/data/oceans_output/shelf/kaight/mitgcm',  gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/')



#02-Maps for New Geometry Runs--------------------------
'''
- Shelfice Freshwater Flux
- Salinity (0-200m; 200-700m)
- Temperature (0-200m; 200-700m)

'''




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


#---Bottom Current-------------------------------------------
# var='UVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=None, bottom=True, zbottom=100, save=True, output=False, gtype='u')

# makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='UVEL')

# for mem in [4,6,10,11,12]:
#     members=[mem-1,mem]

#makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='VVEL')

# var='VVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=None, bottom=True, zbottom=100, save=True, output=False, gtype='v')

# makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='VVEL_bottom')

# var='VVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=None, bottom=False, zbottom=100, save=True, output=False, gtype='v')

# makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='VVEL_depth')

# var='UVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=None, bottom=False, zbottom=100, save=True, output=False, gtype='u')


#---Surface Current-------------------------------------------
# var='UVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=[-50, 0], bottom=False, zbottom=100, save=True, output=False)

# var='VVEL'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=[-50, 0], bottom=False, zbottom=100, save=True, output=False)

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


# #---SIfwfrz--------------------------------------------
# members='all'
# var='SIfwmelt'
# #READ SURFACE Snow Height
# data=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=True, filename=None)
# #----------------------------------------------------------------

# # #04-Make Everything Smaller (to save space)----------------------------------------------------
#makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='THETA')
# #-------------------------------------------------------------------------------------------

# #---SIfwfrz--------------------------------------------
# members='all'
# var='SIfwfrz'
# #READ SURFACE Snow Height
# data=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=True, filename=None)
# #----------------------------------------------------------------

# # #04-Make Everything Smaller (to save space)----------------------------------------------------
#makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='VEL_depth_averaged')
#makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/maps/', newtype=np.float32, kind='VEL_averaged_-5500')
# #-------------------------------------------------------------------------------------------


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
# members=[2,20]
# var='SALT'
# #READ SALT at -5m
# readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members=members, save=True, filename='Z5', at_z=-5)
#mean=data.mean(dim='ens')
#mean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/Z195_WVEL_ensmean.nc')

# var='SALT'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=[-700, -200], bottom=False, zbottom=100, save=True, output=False)

# var='SALT'
# readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', 
#                      z0=[-200, 0], bottom=False, zbottom=100, save=True, output=False)

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





#----------------------------------------------------------------------------------------
#%%READ TIME SERIES PAS
# import sys
# import os
# sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
# from mitgcm_python_master.postprocess import precompute_timeseries_coupled

# segment_dir=[i for i in os.listdir('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output') if (i[0]=='1') | (i[0]=='2')]
# segment_dir=[i for i in segment_dir if (int(i[:4])>=1920)]

# print(segment_dir)

# for i in range(1,21):
#     print('Start with member: '+str(i))
#     precompute_timeseries_coupled (output_dir='/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE{}/output'.format(str(i).zfill(2)), timeseries_file='/data/hpcdata/users/grejan/mitgcm/02_data/timeseries/timeseries_vwind_PACE{}.nc'.format(str(i).zfill(2)), 
#                                    hovmoller_file='hovmoller.nc', file_name='output.nc', segment_dir=segment_dir, timeseries_types=['inner_amundsen_shelf_vwind_avg', 'pine_island_bay_vwind_avg'], 
#                                    hovmoller_loc=[], key='PAS', time_average=False)

#     print('Finished with member: '+str(i))
# print('Finished with *EVERYTHING*!')
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
#%% READ TIME SERIES NEW GEOMETRY

#makeEverythingSmaller(path='/data/hpcdata/users/grejan/mitgcm/02_data/timeseries/', newtype=np.float32, kind='all')

# import sys
# import os
# sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
# from mitgcm_python_master.postprocess import precompute_timeseries_coupled

# segment_dir=[i for i in os.listdir('/data/oceans_output/shelf/grejan/archer2_mitgcm_GEO/PAS_06_GEO/output') if (i[0]=='1') | (i[0]=='2')]
# segment_dir=[i for i in segment_dir if (int(i[:4])>=1920)]
# segment_dir=[i for i in segment_dir if (int(i[:4])<=1995)]


# print(segment_dir)

# for i in [4,6,10,11,12]:
#     print('Start with member: '+str(i))
#     precompute_timeseries_coupled (output_dir='/data/oceans_output/shelf/grejan/archer2_mitgcm_GEO/PAS_{}_GEO/output'.format(str(i).zfill(2)), timeseries_file='/data/hpcdata/users/grejan/mitgcm/02_data/timeseries/timeseries_full_GEO{}.nc'.format(str(i).zfill(2)), 
#                                    hovmoller_file='hovmoller.nc', file_name='output.nc', segment_dir=segment_dir, timeseries_types=['amundsen_shelf_break_uwind_avg', 'all_massloss', 'dotson_to_cosgrove_massloss', 'dotson_massloss', 'pig_massloss', 'thwaites_massloss', 'cosgrove_massloss', 'crosson_massloss', 'amundsen_shelf_ohc_below_0m', 'amundsen_shelf_temp_btw_200_700m', 'amundsen_shelf_salt_below_0m', 'all_melting', 'dotson_to_cosgrove_melting', 'dotson_melting', 'pig_melting', 'thwaites_melting', 'cosgrove_melting', 'crosson_melting', 'amundsen_shelf_seaice_melt', 'amundsen_shelf_seaice_freeze'], 
#                                    hovmoller_loc=[], key='PAS', time_average=False)

#     print('Finished with member: '+str(i))
# print('Finished with *EVERYTHING*!') 
#----------------------------------------------------------------------------------------
# import sys
# import os
# sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
# from mitgcm_python_master.postprocess import precompute_timeseries_coupled

# segment_dir=[i for i in os.listdir('/data/oceans_output/shelf/grejan/archer2_mitgcm_GEO/PAS_06_GEO/output') if (i[0]=='1') | (i[0]=='2')]
# segment_dir=[i for i in segment_dir if (int(i[:4])>=1920)]
# #segment_dir=[i for i in segment_dir if (int(i[:4])<=1952)]

# print(segment_dir)

# for i in range(6,7):
#     print('Start with member: '+str(i))
#     precompute_timeseries_coupled (output_dir='/data/oceans_output/shelf/grejan/archer2_mitgcm_GEO/PAS_{}_GEO/output'.format(str(i).zfill(2)), timeseries_file='/data/hpcdata/users/grejan/mitgcm/02_data/timeseries/timeseries_melting_GEO{}.nc'.format(str(i).zfill(2)), 
#                                    hovmoller_file='hovmoller.nc', file_name='output.nc', segment_dir=segment_dir, timeseries_types=[ 'all_melting', 'dotson_to_cosgrove_melting', 'dotson_melting', 'pig_melting', 'thwaites_melting', 'cosgrove_melting'], 
#                                    hovmoller_loc=[], key='PAS', time_average=False)

#     print('Finished with member: '+str(i))
# print('Finished with *EVERYTHING*!')

#-----------------------------------------------------------------------------------------------------
# #%%READ TIME SERIES PAS
# import sys
# import os
# sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
# from mitgcm_python_master.postprocess import precompute_timeseries_coupled

# segment_dir=[i for i in os.listdir('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output') if (i[0]=='1') | (i[0]=='2')]
# segment_dir=[i for i in segment_dir if (int(i[:4])>=1920)]

# print(segment_dir)

# for i in [4,10,11,12]:
#     print('Start with member: '+str(i))
#     precompute_timeseries_coupled (output_dir='/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE{}/output'.format(str(i).zfill(2)), timeseries_file='/data/hpcdata/users/grejan/mitgcm/02_data/timeseries/timeseries_melting_PACE{}.nc'.format(str(i).zfill(2)), 
#                                    hovmoller_file='hovmoller.nc', file_name='output.nc', segment_dir=segment_dir, timeseries_types=['all_melting', 'dotson_to_cosgrove_melting', 'dotson_melting', 'pig_melting', 'thwaites_melting', 'cosgrove_melting'], 
#                                    hovmoller_loc=[], key='PAS', time_average=False)

#     print('Finished with member: '+str(i))
# print('Finished with *EVERYTHING*!')