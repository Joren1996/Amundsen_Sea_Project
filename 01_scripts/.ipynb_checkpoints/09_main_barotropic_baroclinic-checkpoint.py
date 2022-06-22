# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain test!.

@author: joren
"""

#%%IMPORTING
import xarray as xr
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import pandas as pd
import datetime
import matplotlib.colors as colors

import sys
sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.

# def selectCoordinates(var_name, lat_list, lon_list, gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/', gtype='u', filename=''):
#     from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
#     from mitgcm_python_master.file_io import read_netcdf, netcdf_time
#     from mitgcm_python_master.interpolation import interp_bilinear
    
#     grid = Grid(gp)
    
#     #%% MAIN CODE
#     #%% INITIALISATION
#     segment_dir=[i for i in os.listdir('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output') if (i[0]=='1') | (i[0]=='2')]
#     segment_dir.sort()
#     for i in range(1,21):
#         print('Start with member: '+str(i))

#         output_dir='/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE{}/output'.format(str(i).zfill(2))
#         for k, seg in enumerate(segment_dir):
#             print('Start with '+seg)
#             file_path=os.path.join(output_dir, seg, 'MITgcm/output.nc')
#             data=read_netcdf(file_path, var_name)
#             time=netcdf_time (file_path)
#             z=grid.z
            
#             lats=np.unique(lat_list[0])
#             lons=np.unique(lon_list[0])
            
            
#             new=xr.DataArray(np.full((len(time),len(lats), len(lons), len(z)), np.nan), dims=('time', 'lat', 'lon',"Z"), coords={ 'time': time, 'lat': lats, 'lon': lons, 'Z':z})
#             hfacnew=xr.DataArray(np.full((len(lats), len(lons), len(z)), np.nan), dims=('lat', 'lon',"Z"), coords={'lat': lats, 'lon': lons, 'Z':z})

            
#             interp=np.zeros((len(lat_list[0]), 12, 50))
#             hfac=np.zeros((len(lat_list[0]), 50))
#             for j in range(len(lat_list[0])):
#                 if (i==1) & (k==0):
#                     interp[j,:,:], hfac[j,:]=interp_bilinear (data=data, lon0=lon_list[0][j], lat0=lat_list[0][j], grid=grid, gtype=gtype, return_hfac=True)
#                     hfacnew.loc[lat_list[0][j], lon_list[0][j],:]=hfac[j,:]
#                 else:
#                     interp[j,:,:]=interp_bilinear (data=data, lon0=lon_list[0][j], lat0=lat_list[0][j], grid=grid, gtype=gtype, return_hfac=False)
#                 new.loc[:, lat_list[0][j], lon_list[0][j],:]=interp[j,:,:]
                
#             #print(np.shape(interp))
#             #new_array=xr.DataArray(interp, dims=('space', 'time', "Z"), coords={'space':np.arange(len(lat_list[0])), 'lat':(('space'), lat_list[0]), 'lon':(('space'), lon_list[0]), 'time': time, 'Z':grid.z}) 
                
            
#             print(k)
#             if k==0:
#                 final_array=new
#             else:
#                 final_array=xr.concat([final_array, new], dim='time')
        
#             #                                                              "Z": (trans[3][:,0]+trans[4][:,0])/2, 
#             #                                                               "Dist": trans[2][0,:]})
#         final_array.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/coordinates/'+var_name+'_'+filename+'_ens'+str(i).zfill(2)+'.nc')
#     hfacnew.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/coordinates/hFac_'+var_name+'_'+filename+'.nc')
#     return

# import scipy.io
# coords=scipy.io.loadmat('/data/hpcdata/users/grejan/mitgcm/coord.mat')
# data=selectCoordinates(var_name='UVEL', lat_list=coords['lat_undercurrent'], lon_list=coords['lon_undercurrent'], gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/', gtype='u', filename='undercurrentAll')
# data=selectCoordinates(var_name='VVEL', lat_list=coords['lat_undercurrent'], lon_list=coords['lon_undercurrent'], gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/', gtype='v', filename='undercurrentAll')


#--------------------------------------------------------------------------------------------------------------------
# #Compute the ensemble mean.
# from functions.loading_and_processing_data import loadMITgcmData
# dataU=loadMITgcmData(filename='UVEL_undercurrentAll', members='all', kind='coordinates')
# dataV=loadMITgcmData(filename='VVEL_undercurrentAll', members='all', kind='coordinates')

# umean=dataU.mean(dim='ens')
# vmean=dataV.mean(dim='ens')
# umean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/coordinates/UVEL_undercurrentAll_ensmean.nc')
# vmean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/coordinates/VVEL_undercurrentAll_ensmean.nc')                


#-------------------------------------------------------------------------------------------
# import pandas as pd
# import xarray as xr
# from functions.loading_and_processing_data import loadMITgcmData
# #dataU=loadMITgcmData(filename='UVEL_undercurrentA', members='all', kind='coordinates')
# #dataV=loadMITgcmData(filename='VVEL_undercurrentA', members='all', kind='coordinates')


# for ens in ['ens01', 'ens19', 'ensmean']:
#     b=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/coordinates/UVEL_undercurrentA_'+ens+'.nc')
#     c=b.copy()
#     c=c.drop_vars('lat')
#     c=c.drop_vars('lon')
#     c.coords['space'] = pd.MultiIndex.from_tuples([(b.lat.values[i], b.lon.values[i]) for i in b.space], names=['lat', 'lon'])
#     c.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/coordinates/UVEL_undercurrentA_multiindex_'+ens+'.nc')
    
#     b=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/coordinates/VVEL_undercurrentA_'+ens+'.nc')
#     c=b.copy()
#     c=c.drop_vars('lat')
#     c=c.drop_vars('lon')
#     c.coords['space'] = pd.MultiIndex.from_tuples([(b.lat.values[i], b.lon.values[i]) for i in b.space], names=['lat', 'lon'])
#     c.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/coordinates/VVEL_undercurrentA_multiindex_'+ens+'.nc')



# time=umean.time
# z=umean.Z
# lat_list=coords['lat_undercurrent']
# lon_list=coords['lon_undercurrent']
# lats=np.unique(lat_list[0])
# lons=np.unique(lon_list[0])

# new=xr.DataArray(np.full((len(time),len(lats), len(lons), len(z)), np.nan), dims=('time', 'lat', 'lon',"Z"), coords={ 'time': time, 'lat': lats, 'lon': lons, 'Z':z})

# for i in range(len(lat_list[0]))[:2]:
#     new.loc[:, lat_list[0][i], lon_list[0][i],:]=umean.sel(space=i).to_array().to_numpy()[0,:,:]



#-------------------------------------------------------------------------------------------------------------
#Loading Winds
from functions.loading_and_processing_data import loadMITgcmData
from functions.preprocessing import readMITgcmData

var='EXFuwind'
dataU=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', save=True, filename=None)
umean=dataU.mean(dim='ens')
umean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/EXFuwind_ensmean.nc')


var='EXFvwind'
dataV=readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', save=True, filename=None)
vmean=dataV.mean(dim='ens')
vmean.to_netcdf('/data/hpcdata/users/grejan/mitgcm/02_data/maps/EXFvwind_ensmean.nc')
