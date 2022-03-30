# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 13:30:03 2022

@author: joren
"""


import unixFunctionsMITgcm_final as fc

# members='all'

# var='UVEL'
# data=fc.readMITgcmData_depth_averaged(var,save=True, members=members, bottom=True)

# var='VVEL'
# data=fc.readMITgcmData_depth_averaged(var,save=True, members=members, bottom=True)



# #COMPUTING TIME SERIES OF SEA ICE
# from mitgcm_python_master.postprocess import precompute_timeseries_coupled
# import os

# segment_dir=[i for i in os.listdir('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output') if (i[0]=='1') | (i[0]=='2')]
# segment_dir=[i for i in segment_dir if (int(i[:4])>=1920)]

# print(segment_dir)

# for i in range(1,21):
#     print('Start with member: '+str(i))
#     precompute_timeseries_coupled (output_dir='/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE{}/output'.format(str(i).zfill(2)), timeseries_file='./data/timeseries_seaice_PACE{}.nc'.format(str(i).zfill(2)), 
#                                    hovmoller_file='hovmoller.nc', file_name='output.nc', segment_dir=segment_dir, timeseries_types=['amundsen_shelf_seaice_freeze', 'pine_island_bay_seaice_freeze'], 
#                                    hovmoller_loc=[], key='PAS', time_average=False)

#     print('Finished with member: '+str(i))
# print('Finished with *EVERYTHING*!')


#COMPUTING OCEAN HEAT CONTENT
from mitgcm_python_master.timeseries import read_data_xy
import os
import numpy as np
import pandas as pd
import xarray as xr
from mitgcm_python_master.utils import mask_2d_to_3d, apply_mask, select_bottom
from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
from mitgcm_python_master.calculus import vertical_average





#%% MAIN CODE
#%% INITIALISATION
segment_dir=[i for i in os.listdir('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output') if (i[0]=='1') | (i[0]=='2')]
#segment_dir=[i for i in segment_dir if (int(i[:4])>=1920)]
var_name='ADVx_TH'
var_y = var_name.replace('x', 'y')

print(segment_dir)

#% PREPARE ALL
gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
grid = Grid(gp)
#Bottom Flux
print('Start with bottom mask!')
data_xr=xr.open_dataset('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output/194001/MITgcm/output.nc')
depth=data_xr[var_name].Z.to_numpy()
depth = np.stack([depth]*np.shape(grid.bathy)[0], axis=1)
depth = np.stack([depth]*np.shape(grid.bathy)[1], axis=2)
numdat=data_xr[var_name].to_numpy()
mask=(numdat[-1,:,:,:]==select_bottom(numdat[-1,:,:,:], grid=grid, masked=False)) | (depth<=(grid.bathy+100)) & (numdat[-1,:,:]!=0)

print('Created a bottom mask!')


for i in range(1,2):
    print('Start with member: '+str(i))
    
    output_dir='/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE{}/output'.format(str(i).zfill(2))
    for k, seg in enumerate(segment_dir):
        print('Start with '+seg)
        file_path=os.path.join(output_dir, seg, 'MITgcm/output.nc')
        data_x, data_y=read_data_xy (file_path, var_name, time_index=None, t_start=None, t_end=None, time_average=False)
        
        print('Reading complete, taking averages now')
        #Bottom Average
        numdatx = apply_mask(data_x, np.invert(mask), time_dependent=True)
        numdatx=vertical_average (numdatx, grid, gtype='t', time_dependent=True)
        ax =xr.DataArray(numdatx, dims=('time', "YC", "XG"), coords={'time': data_xr[var_name].time, "YC": data_xr[var_name].YC, "XG": data_xr[var_name].XG})
        
        numdaty = apply_mask(data_y, np.invert(mask), time_dependent=True)
        numdaty=vertical_average (numdaty, grid, gtype='t', time_dependent=True)
        ay =xr.DataArray(numdaty, dims=('time', "YG", "XC"), coords={'time': data_xr[var_y].time, "YG": data_xr[var_y].YG, "XC": data_xr[var_y].XC})
        
        #Depth Average
        data_x=vertical_average (data_x, grid, gtype='t', time_dependent=True)
        bx =xr.DataArray(data_x, dims=('time', "YC", "XG"), coords={'time': data_xr[var_name].time, "YC": data_xr[var_name].YC, "XG": data_xr[var_name].XG})
    
        data_y=vertical_average (data_y, grid, gtype='t', time_dependent=True)
        by =xr.DataArray(data_y, dims=('time', "YG", "XC"), coords={'time': data_xr[var_y].time, "YG": data_xr[var_y].YG, "XC": data_xr[var_y].XC})
        
        bx=bx.to_dataset(name=var_name)
        by=by.to_dataset(name=var_y)
        ax=ax.to_dataset(name=var_name)
        ay=ay.to_dataset(name=var_y)
        
        if k==0:
            fullx=bx
            fully=by
            bottomx=ax
            bottomy=ay

        else:
            fullx=xr.concat([fullx,bx], dim='time')
            fully=xr.concat([fully,by], dim='time')
            bottomx=xr.concat([bottomx,ax], dim='time')
            bottomy=xr.concat([bottomy,ay], dim='time')
            
    fullx.to_netcdf('./data/new_depth_averaged_'+var_name+'_ens'+str(i+1)+'.nc')
    fully.to_netcdf('./data/new_depth_averaged_'+var_y+'_ens'+str(i+1)+'.nc')    
        
    bottomx.to_netcdf('./data/new_bottom100m_averaged_'+var_name+'_ens'+str(i+1)+'.nc')
    bottomy.to_netcdf('./data/new_bottom100m_averaged_'+var_y+'_ens'+str(i+1)+'.nc')   
        
        
        
    #precompute_timeseries_coupled (output_dir='/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE{}/output'.format(str(i).zfill(2)), timeseries_file='./data/timeseries_seaice_PACE{}.nc'.format(str(i).zfill(2)), 
    #                               hovmoller_file='hovmoller.nc', file_name='output.nc', segment_dir=segment_dir, timeseries_types=['amundsen_shelf_seaice_freeze', 'pine_island_bay_seaice_freeze'], 
    #                               hovmoller_loc=[], key='PAS', time_average=False)

    print('Finished with member: '+str(i))
print('Finished with *EVERYTHING*!')


















# import xarray as xr
# import os
# import shutil


# real=True

# goal='hovmoller_adv.nc'
# #goal='timeseries_final.nc'

# if real==True:
#     mainpath=r'/data/oceans_output/shelf/kaight/mitgcm'
# else:
#     mainpath=r'C:\Users\joren\OneDrive\Documenten\Master Thesis\Project\MITgcm\data_MITgcm\timeseries'

# for folder in os.listdir(mainpath):
#     if 'PAS_PACE' in folder:
#         if 'hb' in folder:
#             print('Skip!')
#         else:
#             if real==True:
#                 file=os.path.join(mainpath, folder, 'output', goal)
#             else:
#                 file=os.path.join(mainpath, folder, goal)
#             print('Start with {}!'.format(file))
#             shutil.copyfile(file, goal[:-3]+folder[-7:]+'.nc')
#             print('Finished with {}!'.format(file))
