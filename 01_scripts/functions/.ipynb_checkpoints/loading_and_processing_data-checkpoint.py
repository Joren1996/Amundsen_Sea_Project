# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain all functions to load and detrend data that is already preprocessed.

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

def loadMITgcmData(filename='zonal_winds', members='all', kind='maps'):
    '''
    INPUT:
    filename (string): name of the file that has to be loaded. This should be followed by '_ensXX'.
    members (string or list): members to read ('all' or list with boundaries (e.g. [0, 1] for member 1)).
    kind: which subfolder we want to use
    
    RETURN:
    data (xr.Dataset): dataset containing the variable of interest and the ensemble members along dimension 'ens'.
    '''
    if members=='all':
        ens_list=np.arange(20)
    else:
        ens_list=np.arange(members[0], members[-1])
        
    for j,i in enumerate(ens_list):
        print('Loading number:'+str(i))
        if kind=='old':
            print('../data/'+filename+'_ens'+str(i+1)+'.nc')
            a=xr.open_dataset('../data/'+filename+'_ens'+str(i+1)+'.nc')
        else:
            try:
                print('../02_data/'+kind+'/'+filename+'_ens'+str(i+1)+'.nc')
                a=xr.open_dataset('../02_data/'+kind+'/'+filename+'_ens'+str(i+1)+'.nc')
            except:
                print('../02_data/'+kind+'/'+filename+'_ens'+str(i+1).zfill(2)+'.nc')
                a=xr.open_dataset('../02_data/'+kind+'/'+filename+'_ens'+str(i+1).zfill(2)+'.nc')
        if j==0:
            data=a
        else:
            data=xr.concat([data,a], dim='ens')   
    
    if len(ens_list)==1:
        data=data.expand_dims("ens")
    
    return data





def loadTimeSeriesData(filename='timeseries_final', members='all', kind='old'):
    '''
    INPUT:
    filename (string): name of the file that has to be loaded. This should be followed by '_ensXX'.
    members (string or list): members to read ('all' or list with boundaries (e.g. [0, 1] for member s1)).
    kind: which subfolder we want to use
    
    RETURN:
    data (xr.Dataset): dataset containing the variable of interest and the ensemble members along dimension 'ens'.
    '''
    if members=='all':
        ens_list=np.arange(20)
    else:
        ens_list=np.arange(members[0], members[-1])
        
    for j,i in enumerate(ens_list):
        print('Loading number:'+str(i))
        if kind=='old':
            print('../data/'+filename+'_PACE'+str(i+1).zfill(2)+'.nc')
            a=xr.open_dataset('../data/'+filename+'_PACE'+str(i+1).zfill(2)+'.nc')
        else:
            try:
                print('../02_data/'+kind+'/'+filename+'_PACE'+str(i+1).zfill(2)+'.nc')
                a=xr.open_dataset('../02_data/'+kind+'/'+filename+'_PACE'+str(i+1).zfill(2)+'.nc')
            except:
                print('../02_data/'+kind+'/'+filename+'_ens'+str(i+1).zfill(2)+'.nc')
                a=xr.open_dataset('../02_data/'+kind+'/'+filename+'_PACE'+str(i+1).zfill(2)+'.nc')
        if j==0:
            data=a
        else:
            data=xr.concat([data,a], dim='ens')   
    
    if len(ens_list)==1:
        data=data.expand_dims("ens")
    
    return data




def loadPACEData(var='PSL', members='all'):
    '''
    INPUT:
    filename (string): name of the variable that has to be loaded. This should be followed by '_ensXX'.
    members (string or list): members to read ('all' or list with boundaries (e.g. [0, 1] for member 1)).
    
    RETURN:
    data (xr.Dataset): dataset containing the variable of interest and the ensemble members along dimension 'ens'.
    '''
    if members=='all':
        ens_list=np.arange(20)
    else:
        ens_list=np.arange(members[0], members[-1])
        
    for j,i in enumerate(ens_list):
        print('Loading number:'+str(i))
        if var=='PSL':
            a=xr.open_dataset('/data/oceans_input/raw_input_data/CESM/PPACE/monthly/'+var+'/b.e11.B20TRLENS.f09_g16.SST.restoring.ens'+str(i+1).zfill(2)+'.cam.h0.'+var+'.192001-200512.nc')[var]
            b=xr.open_dataset('/data/oceans_input/raw_input_data/CESM/PPACE/monthly/'+var+'/b.e11.BRCP85LENS.f09_g16.SST.restoring.ens'+str(i+1).zfill(2)+'.cam.h0.'+var+'.200601-201312.nc')[var]
            
            a=xr.concat([a,b], dim='time')
            a=a.to_dataset()
        elif var=='SST':
            a=xr.open_dataset('/data/oceans_input/raw_input_data/CESM/PPACE/monthly/'+var+'/b.e11.B20TRLENS.f09_g16.SST.restoring.ens'+str(i+1).zfill(2)+'.pop.h.SST.192001-200512.nc')[var]
            b=xr.open_dataset('/data/oceans_input/raw_input_data/CESM/PPACE/monthly/'+var+'/b.e11.BRCP85LENS.f09_g16.SST.restoring.ens'+str(i+1).zfill(2)+'.pop.h.SST.200601-201312.nc')[var]
            a=xr.concat([a,b], dim='time')
            a=a.to_dataset()
            
        else:
            print('Not coded yet...')
    
        if j==0:
            data=a
        else:
            data=xr.concat([data,a], dim='ens')   
    
    if len(ens_list)==1:
        data=data.expand_dims("ens")
    
    return data


def detrend_and_average_MITgcmData(data, var, window=24, method='mean', start='1920', end='2013', longwindow=12*25, min_periods=5*12):
    '''
    data: (xarray dataset)
    var: variable name (string)
    window: rolling mean window [months] (int)
    method: mean, linear or quadratic (string)
    start: start year (string)
    end: end year (string)
    longwindow: rolling mean window for detrening, 
                only during mean detrending[months] (int)
    min_periods: minimum number of datapoints for rolling mean window for detrending
                 only during mean detrending[months] (int)
    '''


    def detrend_dim(da, dim, deg=1):
        # detrend along a single dimension
        p = da.polyfit(dim=dim, deg=deg, skipna=True)
        #return p
        x=da[dim]
        fit = xr.polyval(x, p.polyfit_coefficients)
        return da - fit

    def detrend(da, dims, deg=1):
        # detrend along multiple dimensions
        # only valid for linear detrending (deg=1)
        da_detrended = da
        for dim in dims:
            da_detrended = detrend_dim(da_detrended, dim, deg=deg)
        return da_detrended
    
    data=data.sel(time=slice(start, end)).rolling(time=window, center=True).mean()
    
    if method=='linear':
        deg=1
        data=detrend(data[var], dims=['time'], deg=deg)
    elif method=='quadratic':
        deg=2
        data=detrend(data[var], dims=['time'], deg=deg)
    elif method=='mean':
        long_mean=data.rolling(time=longwindow, center=True, min_periods=min_periods).mean()
        data=data-long_mean
        data=data[var]
    
    return data


def read_all_data(full, total, units, longnames, kind):
    file_paths=[os.path.join(full, file) for file in os.listdir(full) if kind in file]
    
    if kind!='trough':    
        names=['ens'+str(file[-5:-3]) for file in file_paths]

        #print(file_paths)
        varlist=list(xr.open_dataset(file_paths[1]).variables)

        for var in varlist:
            total[var]=pd.DataFrame()    
        for i, file in enumerate(file_paths):
            data=xr.open_dataset(file)

            timeind=netcdf_time(file, monthly=False)

            for var in varlist:
                total[var][names[i]]=pd.Series(data[var], index=timeind)
                if (i==0) & (var!='time'):
                    units[var]=data[var].units
                    longnames[var]=data[var].long_name
    
    else:
        varlist=['PITE', 'Burke', 'Dotson', 'BRE', 'BRW']
        for var in varlist:
            ts=xr.open_dataarray('./data/'+var+'_depth_averaged_vflow.nc')
            ts=ts.drop('YG')
            names=['ens'+str(i+1).zfill(2) for i in ts.ens.values]
            total[var]=pd.DataFrame(data=ts.values.T, columns=names, index=ts.indexes['time'].to_datetimeindex())
            units[var]='[m/s]'
            longnames[var]='Depth Averaged Meridional Flow through Cross Section of '+var
            
            varb=var+'B'
            ts=xr.open_dataarray('./data/'+var+'_bottom100m_vflow.nc')
            ts=ts.drop('YG')
            names=['ens'+str(i+1).zfill(2) for i in ts.ens.values]
            total[varb]=pd.DataFrame(data=ts.values.T, columns=names, index=ts.indexes['time'].to_datetimeindex())
            units[varb]='[m/s]'
            longnames[varb]='Bottom 100m Meridional Flow (depth averaged flow subtracted) through Cross Section of '+var
            
            varb=var+'_full'
            ts=xr.open_dataarray('./data/'+var+'_depth_averaged_vflow_full.nc')
            ts=ts.drop('YG')
            names=['ens'+str(i+1).zfill(2) for i in ts.ens.values]
            total[varb]=pd.DataFrame(data=ts.values.T, columns=names, index=ts.indexes['time'].to_datetimeindex())
            units[varb]='[m/s]'
            longnames[varb]='Depth Averaged Meridional Flow through Cross Section of '+var+' (Not Detrended)'
            
            varb=var+'B_full'
            ts=xr.open_dataarray('./data/'+var+'_bottom100m_vflow_full.nc')
            ts=ts.drop('YG')
            names=['ens'+str(i+1).zfill(2) for i in ts.ens.values]
            total[varb]=pd.DataFrame(data=ts.values.T, columns=names, index=ts.indexes['time'].to_datetimeindex())
            units[varb]='[m/s]'
            longnames[varb]='Bottom 100m Meridional Flow (depth averaged flow subtracted) through Cross Section of '+var
            
            
    return total, units, longnames

def netcdf_time (file_path, var_name='time', t_start=None, t_end=None, return_date=True, monthly=True, return_units=False):

    import netCDF4 as nc

    # Open the file and get the length of the record
    id = nc.Dataset(file_path, 'r')
    time_id = id.variables[var_name]
    units = time_id.units
    try:
        calendar = time_id.calendar
    except(AttributeError):
        calendar = 'standard'
    num_time = time_id.size

    # Choose range of time values to consider
    # If t_start and/or t_end are already set, use those bounds
    # Otherwise, start at the first time_index and/or end at the last time_index in the file
    if t_start is None:
        t_start = 0
    if t_end is None:
        t_end = num_time

    # Read the variable
    if return_date:
        # Return as handy Date objects
        time = nc.num2date(time_id[t_start:t_end], units=units, calendar=calendar)
    else:
        # Return just as scalar values
        time = time_id[t_start:t_end]
    id.close()

    if return_date:
        # Want to convert to a datetime object
        if monthly:
            # Back up to previous month
            for t in range(time.size):
                month = time[t].month-1
                year = time[t].year
                if month < 1:
                    month += 12
                    year -= 1
                time[t] = datetime.datetime(year, month, 1)
        else:
            for t in range(time.size):
                time[t] = datetime.datetime(time[t].year, time[t].month, time[t].day)             

    if return_units:
        return time, units, calendar
    else:
        return time
    
    
def smoothDataset(data, var, x=21, y=21, which_area='rA'):
    import xarray as xr
    data=data.where(data.sel(time=data.time[-30])!=0, drop=True)
    ra=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/coordinates/'+which_area+'.nc')
    if which_area=='rA':
        data=(data[var]*ra).rolling(YC=y, center=True).sum().rolling(XC=x, center=True).sum()/ra.rolling(YC=y, center=True).sum().rolling(XC=x, center=True).sum()
    elif which_area=='rAz':
        data=(data[var]*ra).rolling(YG=y, center=True).sum().rolling(XG=x, center=True).sum()/ra.rolling(YG=y, center=True).sum().rolling(XG=x, center=True).sum()
    return data.to_dataset(name=var)


def createMapIndex(data, var, xlim, ylim, which_area='rA'):
    ra=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/coordinates/'+which_area+'.nc')
    if which_area=='rA':
        data=(data*ra).sel(XC=slice(xlim[0], xlim[1]), YC=slice(ylim[0], ylim[1])).sum(dim=['XC', 'YC'])/ra.sel(XC=slice(xlim[0], xlim[1]), YC=slice(ylim[0], ylim[1])).sum()
    elif which_area=='rAz':
            data=(data*ra).sel(XG=slice(xlim[0], xlim[1]), YG=slice(ylim[0], ylim[1])).sum(dim=['XG', 'YG'])/ra.sel(XG=slice(xlim[0], xlim[1]), YG=slice(ylim[0], ylim[1])).sum()
    return data