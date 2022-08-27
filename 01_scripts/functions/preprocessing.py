# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain all functions to read data from the output of MITgcm to smaller files that are useful for the analysis.

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

#%% READING AND LOADING DATA

def readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', save=True, filename=None, at_z=False):
    '''
    Read dataset with the variable of interest directly from MITgcm output nc-files.
    
    INPUT:
    var (string): name of variable as it is named in MITgcm.
    bd (string): base directory, with PAS_PACE folders.
    members (string or list): members to read ('all' or list with boundaries (e.g. [0, 1] for member 1)).
    save (bool):  indicate if the processed data should be saved.
    filename (None or string): if data is saved, you can give the nc-file some extra name.
    
    OUPUT:
    data (xr.Dataset): dataset containing the variable of interest and the ensemble members along dimension 'ens'.
    
    '''
    
    #Creating list with files to read.
    if members=='all':
        pdir=[os.path.join(bd, f) for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)]
        number=20
    else:
        pdir=[os.path.join(bd, f) for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)][members[0]:members[-1]]
        number=members[-1]-members[0]
    
    if 'SI' in var:
        files=[os.path.join('output', f, 'MITgcm/output.nc') for f in os.listdir('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output') if ((f[0]==str(1)) & (f[1]!=str(8)) & (f[2]!=str(0)) & (f[2]!=str(1))) | (f[0]==str(2))]
    else:
        files=[os.path.join('output', f, 'MITgcm/output.nc') for f in os.listdir('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output') if (f[0]==str(1)) | (f[0]==str(2))]
    
    #Reading the actual data.
    data=None

    for i, p in enumerate(pdir):

        fd=[os.path.join(p, f) for f in files]
        fd.sort()
        b=xr.open_dataset(fd[0]) #Read the first year
        #print(b)
        b=b[[var]] #Select the variable of interest
        
        if at_z!=False:
            if 'Z' in b.dims:
                b=b.sel(Z=at_z)
            else:
                #print(b.dims)
                b=b.sel(Zl=at_z)
        
#        if var=='THETA':
#            b=b.sel(Z=slice(-200, -700)).mean(dim='Z')
#        if var=='ADVx_TH':
#            b=b.sel(Z=-405)
#        if var=='ADVy_TH':
#            b=b.sel(Z=-405)
            
        for f in fd[1:]: #Read all other years and add them to the dataset
            a=xr.open_dataset(f)
            a=a[[var]]
            
            if at_z!=False:
                if 'Z' in a.dims:
                    a=a.sel(Z=at_z)
                else:
                    a=a.sel(Zl=at_z)
            
#             if var=='THETA':
#                 a=a.sel(Z=slice(-200, -700)).mean(dim='Z')
#             if var=='ADVx_TH':
#                 a=a.sel(Z=-405)
#             if var=='ADVy_TH':
#                 a=a.sel(Z=-405)

            b=xr.concat([b,a], dim='time')

        #Save the ensemble member: do it now so the members are not mixed up.
        if save==True:
            if number==1:
                i=members[0]
            if type(filename)==type(None):
                b.to_netcdf('../02_data/'+var+'_ens'+str(i+1)+'.nc')
            else:
                b.to_netcdf('../02_data/'+filename+'_'+var+'_ens'+str(i+1)+'.nc') #in case we want to give it an extra name.
            if number==1:
                i=0
            
        if i==0:
            data=b
        else:
            data=xr.concat([data,b], dim='ens')
        
    #Making sure that the dataset has the dimension 'ens'
    if number==1:
        data=data.expand_dims("ens")
            
    return data


def readMITgcmData_depth(var, bd='/data/oceans_output/shelf/kaight/mitgcm', method='average', members='all', z0=None, bottom=False, zbottom=100, save=False, output=True):
    '''
    Read dataset with the variable of interest from MITgcm output nc-files and take a vertical average or sum.
    Makes use of the Kaitlin's code.
    
    INPUT
    var: name of variable (string).
    bd: base directory, with PAS_PACE folders (string).
    method: indicate if we want to integrate or average (string: 'average' or 'integrate')
    members: members to read ('all' or list with boundaries (e.g. [0, 1] for member 1)).
    z0: depth range of interest (None or list with boundaries (e.g. [-700, -200] for range between 700m and 200m depth with respect to surface)). If not given, we take a full depth average.
    save (bool):  indicate if the processed data should be saved.
    output (bool):  indicate if the processed data should be returned.
    bottom (bool):  indicate if we are interested in a specific depth interval with respect to the bottom. NOTE: OVERRIDES DEPTH RANGE!
    zbottom (int): in case we are looking at a specific area above the bottom, define the upper limit (this case 100m above seafloor). Is only used when bottom==True. NOTE: IS 100m ABOVE THE GRID BATHYMETRY: HOWEVER, SOME LOCATIONS HAVE VALUES BELOW THIS BATHYMETRY!
    
    OUTPUT:
    data (xr.Dataset): dataset containing the variable of interest and the ensemble members along dimension 'ens'.
    '''

    #Import Kaitlin's code
    from mitgcm_python_master.utils import mask_2d_to_3d, apply_mask, select_bottom
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from mitgcm_python_master.calculus import vertical_average, vertical_integral
    
    
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    if bottom==False:
        mask = np.ones([grid.ny, grid.nx]).astype(bool)
        if type(z0)!=type(None):
            mask = mask_2d_to_3d(mask, grid, zmin=z0[0], zmax=z0[1])
    
    #Creating list with files to read.
    if members=='all':
        pdir=[os.path.join(bd, f) for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)]
        number=20
    else:
        pdir=[os.path.join(bd, f) for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)][members[0]:members[-1]]
        number=members[-1]-members[0]
        
    files=[os.path.join('output', f, 'MITgcm/output.nc') for f in os.listdir('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output') if (f[0]==str(1)) | (f[0]==str(2))]
    
    #Reading the actual data.
    data=None

    for i, p in enumerate(pdir):
        print(i)
        print(p)
        fd=[os.path.join(p, f) for f in files]
        fd.sort()
        b=xr.open_dataset(fd[0])
        b=b[[var]]
        
        
        numdat=b[var].to_numpy()
        
        if bottom==True:
            print('Apply Bottom Mask')
            depth=b[var].Z.to_numpy()
            depth = np.stack([depth]*np.shape(grid.bathy)[0], axis=1)
            depth = np.stack([depth]*np.shape(grid.bathy)[1], axis=2)
            mask=(numdat[-1,:,:,:]==select_bottom(numdat[-1,:,:,:], grid=grid, masked=False)) | (depth<=(grid.bathy+zbottom)) & (numdat[-1,:,:]!=0)
            numdat = apply_mask(numdat, np.invert(mask), time_dependent=True)
        
        elif type(z0)!=type(None):
            print('Apply Depth Range Mask')
            numdat = apply_mask(numdat, np.invert(mask), time_dependent=True)
        
        if method=='integrate':
            print('Method: integrate')
            numdat=vertical_integral (numdat, grid, gtype='t', time_dependent=True)
        else:
            print('Method: average')
            numdat=vertical_average (numdat, grid, gtype='t', time_dependent=True)
        
        #CONVERT BACK TO DATASET
        if 'XC' in b.variables:
            if 'YC' in b.variables:
                b =xr.DataArray(numdat, dims=('time', "YC", "XC"), coords={'time': b.time, "YC": b.YC, "XC": b.XC})
            else:
                b =xr.DataArray(numdat, dims=('time', "YG", "XC"), coords={'time': b.time, "YG": b.YG, "XC": b.XC})
        else:
            b =xr.DataArray(numdat, dims=('time', "YC", "XG"), coords={'time': b.time, "YC": b.YC, "XG": b.XG})
        b=b.to_dataset(name=var)
        
#         if var=='THETA':
#             numdat=b[var].to_numpy()
#             numdat = apply_mask(numdat, np.invert(mask), time_dependent=True)
#             numdat=vertical_average (numdat, grid, gtype='t', time_dependent=True)
#             b =xr.DataArray(numdat, dims=('time', "YC", "XC"), coords={'time': b.time, "YC": b.YC, "XC": b.XC})
#             b=b.to_dataset(name='THETA')
        
#         if ('VEL' in var) |  ('ADV' in var):
#             numdat=b[var].to_numpy()
            
#             if bottom==True:
#                 print('Apply Bottom Mask')
#                 depth=b[var].Z.to_numpy()
#                 depth = np.stack([depth]*np.shape(grid.bathy)[0], axis=1)
#                 depth = np.stack([depth]*np.shape(grid.bathy)[1], axis=2)
#                 mask=(numdat[-1,:,:,:]==select_bottom(numdat[-1,:,:,:], grid=grid, masked=False)) | (depth<=(grid.bathy+100)) & (numdat[-1,:,:]!=0)
#                 numdat = apply_mask(numdat, np.invert(mask), time_dependent=True)
#             numdat=vertical_average (numdat, grid, gtype='t', time_dependent=True)
#             if 'U' in var:
#                 b =xr.DataArray(numdat, dims=('time', "YC", "XG"), coords={'time': b.time, "YC": b.YC, "XG": b.XG})
#             else:
#                 b =xr.DataArray(numdat, dims=('time', "YG", "XC"), coords={'time': b.time, "YG": b.YG, "XC": b.XC})
#             b=b.to_dataset(name=var)
            #b=b.sel(Z=slice(-200, -700)).mean(dim='Z')

        for f in fd[1:]:
            a=xr.open_dataset(f)
            a=a[[var]]
            
            numdat=a[var].to_numpy()
        
            if bottom==True:
                print('Apply Bottom Mask')
                depth=a[var].Z.to_numpy()
                depth = np.stack([depth]*np.shape(grid.bathy)[0], axis=1)
                depth = np.stack([depth]*np.shape(grid.bathy)[1], axis=2)
                mask=(numdat[-1,:,:,:]==select_bottom(numdat[-1,:,:,:], grid=grid, masked=False)) | (depth<=(grid.bathy+zbottom)) & (numdat[-1,:,:]!=0)
                numdat = apply_mask(numdat, np.invert(mask), time_dependent=True)

            elif type(z0)!=type(None):
                print('Apply Depth Range Mask')
                numdat = apply_mask(numdat, np.invert(mask), time_dependent=True)
            
            if method=='integrate':
                print('Method: integrate')
                numdat=vertical_integral(numdat, grid, gtype='t', time_dependent=True)
            else:
                print('Method: average')
                numdat=vertical_average (numdat, grid, gtype='t', time_dependent=True)

            #CONVERT BACK TO DATASET
            if 'XC' in a.variables:
                if 'YC' in a.variables:
                    a =xr.DataArray(numdat, dims=('time', "YC", "XC"), coords={'time': a.time, "YC": a.YC, "XC": a.XC})
                else:
                    a =xr.DataArray(numdat, dims=('time', "YG", "XC"), coords={'time': a.time, "YG": a.YG, "XC": a.XC})
            else:
                a =xr.DataArray(numdat, dims=('time', "YC", "XG"), coords={'time': a.time, "YC": a.YC, "XG": a.XG})
            a=a.to_dataset(name=var)
            
#             if var=='THETA':
#                 numdat=a[var].to_numpy()
#                 numdat = apply_mask(numdat, np.invert(mask), time_dependent=True)
#                 numdat=vertical_average (numdat, grid, gtype='t', time_dependent=True)
#                 a = xr.DataArray(numdat, dims=('time', "YC", "XC"), coords={'time': a.time, "YC": a.YC, "XC": a.XC})
#                 a=a.to_dataset(name='THETA')
#             if ('VEL' in var) |  ('ADV' in var):
#                 numdat=a[var].to_numpy()
#                 if bottom==True:
#                     numdat = apply_mask(numdat, np.invert(mask), time_dependent=True)
#                 numdat=vertical_average (numdat, grid, gtype='t', time_dependent=True)
#                 if 'U' in var:
#                     a = xr.DataArray(numdat, dims=('time', "YC", "XG"), coords={'time': a.time, "YC": a.YC, "XG": a.XG})
#                 else:
#                     a = xr.DataArray(numdat, dims=('time', "YG", "XC"), coords={'time': a.time, "YG": a.YG, "XC": a.XC})
#                 a=a.to_dataset(name=var)
                #a=a.sel(Z=slice(-200, -700)).mean(dim='Z')
            b=xr.concat([b,a], dim='time')
        if save==True:
            if method=='average':
                if bottom==True:
                    b.to_netcdf('../02_data/maps/'+var+'_bottom100m_averaged_ens'+str(i+1)+'.nc')
                else:
                    if type(z0)==type(None):
                        b.to_netcdf('../02_data/maps/'+var+'_depth_averaged_ens'+str(i+1)+'.nc')
                    else:
                        b.to_netcdf('../02_data/maps/'+var+'_averaged_'+str(z0[0])+'to'+str(z0[-1])+'_ens'+str(i+1)+'.nc')
            elif method=='integrate':
                if bottom==True:
                    b.to_netcdf('../02_data/maps/'+var+'_bottom100m_integrated_ens'+str(i+1)+'.nc')
                else:
                    if type(z0)==type(None):
                        b.to_netcdf('../02_data/maps/'+var+'_depth_integrated_ens'+str(i+1)+'.nc')
                    else:
                        b.to_netcdf('../02_data/maps/'+var+'_integrated_'+str(z0[0])+'to'+str(z0[-1])+'_ens'+str(i+1)+'.nc')
        if output==True:
            if i==0:
                data=b
            else:
                data=xr.concat([data,b], dim='ens')
#    for i in range(len(data.ens)):
#        
    if output==True:
        #Making sure that the dataset has the dimension 'ens'
        if number==1:
            data=data.expand_dims("ens")
        return data
    else:
        return



def readMITgcmData_advection(var_name='ADVx_TH', members='all', bd='/data/oceans_output/shelf/kaight/mitgcm',  gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'):
    from mitgcm_python_master.timeseries import read_data_xy #MAKE SURE THE PART WHERE THEY CORRECT THE FLUXES IS UNCOMMENTED!
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
    
    var_y = var_name.replace('x', 'y')

    #print(segment_dir)

    #% PREPARE ALL
    
    grid = Grid(gp)
    #Bottom Flux
    print('Start with bottom mask!')
    data_xr=xr.open_dataset('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output/194001/MITgcm/output.nc')
    depth=data_xr[var_name].Z.to_numpy()
    depth = np.stack([depth]*np.shape(grid.bathy)[0], axis=1)
    depth = np.stack([depth]*np.shape(grid.bathy)[1], axis=2)
    numdat=data_xr[var_name].to_numpy()
    mask=(numdat[-1,:,:,:]==select_bottom(numdat[-1,:,:,:], grid=grid, masked=False)) | (depth<=(grid.bathy+100)) & (numdat[-1,:,:]!=0)
    
    #Creating list with files to read.
    if members=='all':
        pdir=[os.path.join(bd, f, 'output') for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)]
        number=20
    else:
        pdir=[os.path.join(bd, f, 'output') for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)][members[0]:members[-1]]
        number=members[-1]-members[0]
    
    print('Created a bottom mask!')

    pdir.sort()
    segment_dir.sort()
    for i, p in enumerate(pdir):
        print('Start with member: '+str(i+1))

        output_dir=p
        for k, seg in enumerate(segment_dir):
            print('Start with '+seg)
            file_path=os.path.join(output_dir, seg, 'MITgcm/output.nc')
            
            data_xr=xr.open_dataset(file_path) #Make sure that the time axis is correct...
            
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

        fullx.to_netcdf('../02_data/maps/'+var_name+'_corrected_depth_averaged_ens'+str(i+1)+'.nc')
        fully.to_netcdf('../02_data/maps/'+var_y+'_corrected_depth_averaged_ens'+str(i+1)+'.nc')    

        bottomx.to_netcdf('../02_data/maps/'+var_name+'_corrected_bottom100m_averaged_ens'+str(i+1)+'.nc')
        bottomy.to_netcdf('../02_data/maps/'+var_y+'_corrected_bottom100m_averaged_ens'+str(i+1)+'.nc')   
        
        print('Finished with member: '+str(i))
    print('Finished with *EVERYTHING*!')
    
    return
    
    

def readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', method='lon', x=106, y=-72.5, save=True, filename=None):
    '''
    var: name of variable (string).
    bd: base directory, with PAS_PACE folders.
    members: members to read ('all' or list with boundaries (e.g. [0, 1] for member 1)).
    '''
    #Creating list with files to read.
    if members=='all':
        pdir=[os.path.join(bd, f) for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)]
        number=20
    else:
        pdir=[os.path.join(bd, f) for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)][members[0]:members[-1]]
        number=members[-1]-members[0]
        
    files=[os.path.join('output', f, 'MITgcm/output.nc') for f in os.listdir('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output') if (f[0]==str(1)) | (f[0]==str(2))]
    
    #Reading the actual data.
    data=None

    for i, p in enumerate(pdir):
        print(p)
        fd=[os.path.join(p, f) for f in files]
        fd.sort()
        b=xr.open_dataset(fd[0])
        b=b[[var]]
        
        def find_nearest(array, value):
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return array[idx]
        
        if method=='lon':
            if 'XC' in b.variables:
                xnew=find_nearest(b.XC.to_numpy(), 360-x)
                b=b.sel(XC=xnew)
            else:
                xnew=find_nearest(b.XG.to_numpy(), 360-x)
                b=b.sel(XG=xnew)
        elif method=='lat':
            if 'YC' in b.variables:
                ynew=find_nearest(b.YC.to_numpy(), y)
                b=b.sel(YC=ynew)
            else:
                ynew=find_nearest(b.YG.to_numpy(), y)
                b=b.sel(YG=ynew)
        
            
        for f in fd[1:]:
            print(f)
            a=xr.open_dataset(f)
            a=a[[var]]
            
            if method=='lon':
                if 'XC' in a.variables:
                    xnew=find_nearest(a.XC.to_numpy(), 360-x)
                    a=a.sel(XC=xnew)
                else:
                    xnew=find_nearest(a.XG.to_numpy(), 360-x)
                    a=a.sel(XG=xnew)
            elif method=='lat':
                if 'YC' in a.variables:
                    ynew=find_nearest(a.YC.to_numpy(), y)
                    a=a.sel(YC=ynew)
                else:
                    ynew=find_nearest(a.YG.to_numpy(), y)
                    a=a.sel(YG=ynew)

            
            b=xr.concat([b,a], dim='time')
            
        #Save the ensemble member: do it now so the members are not mixed up.
        if save==True:
            if method=='lon':
                if type(filename)==type(None):
                    b.to_netcdf('../02_data/slices/'+var+'_W'+str(abs(x)).replace(".", "-")+'_ens'+str(i+1)+'.nc')
                else:
                    b.to_netcdf('../02_data/slices/'+filename+'_'+var+'_W'+str(abs(x)).replace(".", "-")+'_ens'+str(i+1)+'.nc') #in case we want to give it an extra name.
            elif method=='lat':
                if type(filename)==type(None):
                    b.to_netcdf('../02_data/slices/'+var+'_S'+str(abs(y)).replace(".", "-")+'_ens'+str(i+1)+'.nc')
                else:
                    b.to_netcdf('../02_data/slices/'+filename+'_'+var+'_S'+str(abs(y)).replace(".", "-")+'_ens'+str(i+1)+'.nc') #in case we want to give it an extra name.
                    
                    
        if i==0:
            data=b
        else:
            data=xr.concat([data,b], dim='ens')
            
    #Making sure that the dataset has the dimension 'ens'
    if number==1:
        data=data.expand_dims("ens")
    
    return data


def readMITgcmData_Transect(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', 
                            x=[-98, -105.5, -106, -106.6, -106, -102], 
                            y=[-75.5, -73.8, -72.8, -72.5, -71.7, -71], transname='test', save=True):
    '''
    var: name of variable (string).
    bd: base directory, with PAS_PACE folders.
    members: members to read ('all' or list with boundaries (e.g. [0, 1] for member 1)).
    '''
    
    from mitgcm_python_master.utils import mask_2d_to_3d, apply_mask, select_bottom
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from mitgcm_python_master.calculus import vertical_average, vertical_integral
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    from mitgcm_python_master.plot_utils.slices import get_transect
    #Creating list with files to read.
    if members=='all':
        pdir=[os.path.join(bd, f) for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)]
        number=20
    else:
        pdir=[os.path.join(bd, f) for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)][members[0]:members[-1]]
        number=members[-1]-members[0]
        
    files=[os.path.join('output', f, 'MITgcm/output.nc') for f in os.listdir('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output') if (f[0]==str(1)) | (f[0]==str(2))]
    
    #Reading the actual data.
    data=None

    for i, p in enumerate(pdir):
        print(p)
        fd=[os.path.join(p, f) for f in files]
        fd.sort()
        b=xr.open_dataset(fd[0])
        b=b[var]
        
        datab={}
        for iota in range(len(x)-1):
            trans=get_transect (b.to_numpy(), grid, 
                                point0=(x[iota], y[iota]), 
                                point1=(x[iota+1], y[iota+1]), 
                                gtype='t', 
                                return_grid_vars=True, 
                                time_dependent=True)
            
            datab[iota]=xr.DataArray(trans[0], dims=('time', "Z", "Dist"), coords={'time': b.time, 
                                                                       "Z": (trans[3][:,0]+trans[4][:,0])/2, 
                                                                       "Dist": trans[2][0,:]})
            datab[iota]=datab[iota].to_dataset(name=var)
            
        for f in fd[1:]:
            print(f)
            a=xr.open_dataset(f)
            a=a[var]
            
            dataa={}
            for iota in range(len(x)-1):
                trans=get_transect (a.to_numpy(), grid, 
                                    point0=(x[iota], y[iota]), 
                                    point1=(x[iota+1], y[iota+1]), 
                                    gtype='t', 
                                    return_grid_vars=True, 
                                    time_dependent=True)

                dataa[iota]=xr.DataArray(trans[0], dims=('time', "Z", "Dist"), coords={'time': a.time, 
                                                                           "Z": (trans[3][:,0]+trans[4][:,0])/2, 
                                                                           "Dist": trans[2][0,:]})

                
                dataa[iota]=dataa[iota].to_dataset(name=var)
                datab[iota]=xr.concat([datab[iota],dataa[iota]], dim='time')
        
        
        if save==True:
            for iota in range(len(x)-1):
                datab[iota].to_netcdf('../02_data/transects/'+var+'_transect_'+transname+'_'+str(iota)+'_ens'+str(i+1)+'.nc')
        
#         if i==0:
#             data=b
#         else:
#             data=xr.concat([data,b], dim='ens')
            
#     #Making sure that the dataset has the dimension 'ens'
#     if number==1:
#         data=data.expand_dims("ens")
            
    return



def readMITgcmTimeSeries(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', tf='timeseries_final.nc'):
    '''
    INPUT:
    var: name of variable (string).
    bd: base directory, with PAS_PACE folders.
    members: members to read ('all' or list with boundaries (e.g. [0, 1] for member 1)).
    tf: timeseries file name; has to be a netcdf file (string)
    
    OUTPUT:
    ts: timeseries file with  (xr.DataArray)
    
    '''
    
    if members=='all':
        pdir=[os.path.join(bd, f) for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)]
        number=20
    else:
        pdir=[os.path.join(bd, f) for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)][members[0]:members[-1]]
        number=members[-1]-members[0]
        
    tsfiles=[os.path.join(p, 'output', tf) for p in pdir]
    
    ts=None

    #print(p)
    fd=tsfiles
    fd.sort()
    b=xr.open_dataset(fd[0])
    b=b[[var]]

    for i, f in enumerate(fd[1:]):
        a=xr.open_dataset(f)
        a=a[[var]]
        b=xr.concat([b,a], dim='ens')
    
    ts=b
    if number==1:
        ts=ts.expand_dims("ens")
    
    return ts



def selectCoordinates(var_name, lat, lon, gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/', gtype='u'):
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from mitgcm_python_master.file_io import read_netcdf
    from mitgcm_python_master.interpolation import interp_bilinear
    
    grid = Grid(gp)
    
    #%% MAIN CODE
    #%% INITIALISATION
    segment_dir=[i for i in os.listdir('/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01/output') if (i[0]=='1') | (i[0]=='2')]
    segment_dir.sort()
    for i in range(1,2):
        print('Start with member: '+str(i))

        output_dir='/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE{}/output'.format(str(i).zfill(2))
        for k, seg in enumerate(segment_dir):
            print('Start with '+seg)
            file_path=os.path.join(output_dir, seg, 'MITgcm/output.nc')
            data=read_netcdf(file_path, var_name)
            interp=np.zeros((len(lat), 50))
            for j in range(len(lat)):
                interp[j,:]=interp_bilinear (data=data, lon0=lon[j], lat0=lat[j], grid=grid, gtype=gtype, return_hfac=False)
                print(shape)
            #new_array=xr.DataArray(interp, dims=('time', "Z", "Dist"), coords={'time': a.time, 
            #                                                              "Z": (trans[3][:,0]+trans[4][:,0])/2, 
            #                                                               "Dist": trans[2][0,:]})
            
    return data




def makeEverythingSmaller(path, newtype, kind=None):
    '''
    Convert the data in netcdf file to new datatype, in order to safe space.
    
    INPUT:
    path: (string) path indicating the folder in which we want to convert the files.
    newtype: (type) type to which we want to convert the data
    
    OUTPUT:
    -
    
    '''
    print('Working on directory:'+path)
    for file in os.listdir(path):
        temppath=os.path.join(path, 'temp.nc')
        if '.nc' in file:
            if (kind==None) | (kind in file):
                print('...Start with '+file)
                fullpath=os.path.join(path, file)
                a=xr.open_dataset(fullpath)
                b=a.astype(newtype)
                b.to_netcdf(temppath)
                os.remove(fullpath)
                os.rename(temppath, fullpath)
                print('...Finished with '+file)
    print('Finished!')
    return
        
    

