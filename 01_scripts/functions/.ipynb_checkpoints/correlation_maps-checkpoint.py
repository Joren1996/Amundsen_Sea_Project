# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain all functions to create the different kinds of correlation maps used in this project.

@author: joren
"""
#%%IMPORTING
import xarray as xr
import xskillscore as xs
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import pandas as pd
import datetime
import matplotlib.colors as colors
import sys

from .loading_and_processing_data import loadMITgcmData, loadPACEData, detrend_and_average_MITgcmData, read_all_data, loadTimeSeriesData
sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid

def simpleCorrelationMap(members='all',
                         data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                         ind=None, indfn='undercurrent_W115', indvar='UVEL',
                         detrend=True, window=24, title='Correlation Map', ymax=-70, save=False, savefn='test', cumulative=False,
                         draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                         inverse=False, ygradient=False, xgradient=False,
                         give_back=False):
    print('Start preparing for correlation map!')
    plt.rcParams.update({'font.size': 14})
    if type(data)==type(None):
        if members=='all':
            print('...Try reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_ensmean.nc')
            try:
                data=xr.open_dataset('/data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_ensmean.nc')
                data=data.expand_dims('ens')
                print('...Took a shortcut! Read Ensmean file.')
            except:
                try:
                    data=xr.open_dataset('/data/hpcdata/users/grejan/mitgcm/data/'+datafn+'_ensmean.nc')
                    data=data.expand_dims('ens')
                    print('...Took a shortcut! Read Ensmean file.')
                except:
                    print('...Shortcut failed, so...')
                    data=loadMITgcmData(filename=datafn, members=members, kind=datakind)
        else:
            data=loadMITgcmData(filename=datafn, members=members, kind=datakind)
    
    if xgradient==True:
        print('...NOTE: Taking XGRADIENT!')
        if 'XC' in data.dims:
            data=data.diff('XC')
        else:
            data=data.diff('XG')
    if ygradient==True:
        print('...NOTE: Taking YGRADIENT!')
        if 'YC' in data.dims:
            data=data.diff('YC')
        else:
            data=data.diff('YG')
    
    
    if type(ind)==type(None):
        print('...Reading: '+indfn)
        ind=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+indfn+'.nc')
        if members!='all':
            print('...Selecting member '+str(members[0]+1))
            ind=ind.sel(ens=[members[0]])
        
    if detrend==True:
        print('...Detrending data')
        data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
        ind=ind.to_dataset()
        ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)
        
    if window!=1:
        print('...Applying rolling mean')
        data=data.rolling(time=window, center=True).mean()
        ind=ind.rolling(time=window, center=True).mean()
    
    print('...Computing correlation coefficients')
    if detrend==True:
        if cumulative==True:
            print('...NOTE: Taking CUMULATIVE!')
            res=xr.corr(data.mean(dim='ens').cumsum(), ind.mean(dim='ens'), dim='time')
        else:
            res=xr.corr(data.mean(dim='ens'), ind.mean(dim='ens'), dim='time')
    else:
        res=xr.corr(data[datavar].mean(dim='ens'), ind.mean(dim='ens'), dim='time')
    
    print('...Plotting')
    fig, ax=plt.subplots(1,1, figsize=(20,10))
    
    if 'XG' in res.dims:
        res=res.rename({'XG':'XC'})
    if 'YG' in res.dims:
        res=res.rename({'YG':'YC'})
    
    if 'lon' in res.dims:
        res=res.rename({'lon':'XC'})
    if 'lat' in res.dims:
        res=res.rename({'lat':'YC'})
        
    if 'nlon' in res.dims:
        res=res.rename({'nlon':'XC'})
    if 'nlat' in res.dims:
        res=res.rename({'nlat':'YC'})
    
    if type(ymax)!=type(None):
        res=res.sel(YC=slice(-75.5, ymax))
        
    
    if members=='all':
        ax.set_title(title+'\n Ensemble Mean')
    else:
        ax.set_title(title+'\n Ensemble Member '+ str(members[-1]))
    
    if inverse==True:
        print('...NOTE: Taking the INVERSE of the correlation!')
        res=-1*res
    
    img=ax.pcolor(res.XC, res.YC, res, cmap='seismic', vmin=-1, vmax=1)
    plt.colorbar(img, label='Pearson Correlation Coefficient')
    CS=ax.contour(res.XC, res.YC, res, np.arange(-1, 1, 0.25))
    ax.clabel(CS, CS.levels, inline=True,fontsize=10)
    
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    bath=ax.contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax.set_ylim([-75.5, ymax])
    ax.set_facecolor('grey')
    ax.set_xlabel('Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')
    
    if draw_box==True:
        ax.plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        if members=='all':
            fig.savefig('../03_output/correlation/'+today+'_cormap_'+savefn+'_ensmean.png')
        else:
            fig.savefig('../03_output/correlation/'+today+'_cormap_'+savefn+'_member'+str(members[-1])+'.png')

    if give_back==False:
        return
    else:
        return data, ind, res

    
    
def correlationMap_2D(members='all',
                      data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      detrend=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False):
    '''
    Input: Always give U component as input!
    
    '''
    
    print('Start preparing for correlation map!')
    plt.rcParams.update({'font.size': 14})
    
    #First, we read the data!
    if type(data)==type(None):
            print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_....')
            dataU=loadMITgcmData(filename=datafn, members=members, kind=datakind)
            datafn_v=datafn.replace('U','V')
            datafn_v=datafn_v.replace('u','v')
            datafn_v=datafn_v.replace('x','y')
            print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn_v+'_....')
            dataV=loadMITgcmData(filename=datafn_v, members=members, kind=datakind)
    
    else:
        dataU=data[0]
        dataV=data[1]
    
    print('...Making sure interpolate everything to the center of the cells')
    if 'XG' in dataU.dims:
        dataU=dataU.rename({'XG':'XC'})
    if 'YG' in dataU.dims:
        dataU=dataU.rename({'YG':'YC'})
    if 'XG' in dataV.dims:
        dataV=dataV.rename({'XG':'XC'})
    if 'YG' in dataV.dims:
        dataV=dataV.rename({'YG':'YC'})
    
    #Read something with the correct locations.
    data=loadMITgcmData(filename='EXFuwind', members=[0,1], kind='maps')
    
    dataU=dataU.interp(XC=data.XC, YC=data.YC)
    dataV=dataV.interp(XC=data.XC, YC=data.YC)
    
    print('...Compute the absolute combined value')
    datavar_v=datavar.replace('U','V')
    datavar_v=datavar_v.replace('u','v')
    datavar_v=datavar_v.replace('x','y')
    data=np.sqrt(dataU[datavar]**2+dataV[datavar_v]**2)
    data=data.to_dataset(name=datavar)
    
    
    
    if type(ind)==type(None):
        print('...Reading: '+indfn)
        if 'timeseries' in indfn:
            ind=loadTimeSeriesData(filename=indfn, members=members, kind='old')
        else:        
            ind=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+indfn+'.nc')
            if members!='all':
                print('...Selecting member '+str(members[0]+1))
                ind=ind.sel(ens=[members[0]])
            ind=ind.to_dataset()
        
    if detrend==True:
        print('...Detrending data')
        data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
        ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)
        
    if window!=1:
        print('...Applying rolling mean')
        data=data.rolling(time=window, center=True).mean()
        ind=ind.rolling(time=window, center=True).mean()
    
    print('...Computing correlation coefficients')
    
    data=data.stack(ens_time=['ens', 'time'])
    for z, lag in enumerate(lags):
        print('...Looking at lag: '+str(lag))
        ind_temp=ind.shift(time=lag) #Positive means index is leading!
        ind_temp=ind_temp.stack(ens_time=['ens', 'time'])
        
        
        if detrend==True:
            res_temp=xr.corr(data, ind_temp, dim='ens_time')
        else:
            res_temp=xr.corr(data[datavar], ind_temp[indvar], dim='ens_time')

        if z==0:
            res=res_temp
        else:
            res=xr.concat([res, res_temp], dim='lag')
    
    print('...Plotting')
    fig, ax=plt.subplots(2,1, figsize=(20,12))
    
    if type(ymax)!=type(None):
        res=res.sel(YC=slice(-75.5, ymax))
    
    if members=='all':
        fig.suptitle(title+'\n All members')
    else:
        if members[1]-members[0]==1:
            fig.suptitle(title+'\n Ensemble Member '+ str(members[-1]))
        else:
            fig.suptitle(title+'\n Ensemble Members '+ str(members[0]+1)+' to '+str(members[-1]))
    
    if inverse==True:
        print('...NOTE: Taking the INVERSE of the correlation!')
        res=-1*res
    
    maxcor=res.max(dim='lag')
    #PLOT CORRELATION COEFFICIENTS
    img=ax[0].pcolor(res.XC, res.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
    plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
    CS=ax[0].contour(res.XC, res.YC, maxcor, np.arange(-1, 1, 0.25))
    ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)
    
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax[0].set_ylim([-75.5, ymax])
    ax[0].set_facecolor('grey')
    ax[0].set_xlabel('Longitude [deg]')
    ax[0].set_ylabel('Latitude [deg]')
    ax[0].set_title('Optimal Correlation Coefficient')
    
    if draw_box==True:
        ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
    
    
    print('...Finding optimal lag')
    idx=np.argmax([i=='lag' for i in res.dims])
    masked=np.ma.masked_invalid(res.to_numpy())
    opt_lag=np.nanargmax(masked, axis=idx)
    opt_lag=opt_lag.astype(float)
    opt_lag[masked.mask.max(axis=idx)]=np.nan
    opt_lag=opt_lag+lags[0]
    
    
    #PLOT OPTIMAL LAG
    img=ax[1].pcolor(res.XC, res.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
    plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
    #CS=ax[1].contour(res.XC, res.YC, opt_lag, np.arange(-1, 1, 0.25))
    #ax[1].clabel(CS, CS.levels, inline=True,fontsize=10)
    
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax[1].set_ylim([-75.5, ymax])
    ax[1].set_facecolor('grey')
    ax[1].set_xlabel('Longitude [deg]')
    ax[1].set_ylabel('Latitude [deg]')
    
    if draw_box==True:
        ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')
    
    fig.subplots_adjust(hspace=0.4)
    
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        if members=='all':
            fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_all.png')
        else:
            fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_member'+str(members[-1])+'.png')

    if give_back==False:
        return
    else:
        return data, ind, res
    
    
def correlationMap_1D(members='all',
                      data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      detrend=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False):
    '''
    Input: Only one map
    
    '''
    
    print('Start preparing for correlation map!')
    plt.rcParams.update({'font.size': 14})
    
    #First, we read the data!
    if type(data)==type(None):
            print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_....')
            data=loadMITgcmData(filename=datafn, members=members, kind=datakind)
    
    print('...Making sure interpolate everything to the center of the cells')
    if 'XG' in data.dims:
        data=data.rename({'XG':'XC'})
    if 'YG' in data.dims:
        data=data.rename({'YG':'YC'})

    #Read something with the correct locations.
    data1=loadMITgcmData(filename='EXFuwind', members=[0,1], kind='maps')
    data=data.interp(XC=data1.XC, YC=data1.YC)
    
    
    if type(ind)==type(None):
        print('...Reading: '+indfn)
        if 'timeseries' in indfn:
            ind=loadTimeSeriesData(filename=indfn, members=members, kind='old')
        else:        
            ind=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+indfn+'.nc')
            if members!='all':
                print('...Selecting member '+str(members[0]+1))
                ind=ind.sel(ens=[members[0]])
            ind=ind.to_dataset()
        
    if detrend==True:
        print('...Detrending data')
        data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
        ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)
        
    if window!=1:
        print('...Applying rolling mean')
        data=data.rolling(time=window, center=True).mean()
        ind=ind.rolling(time=window, center=True).mean()
    
    print('...Computing correlation coefficients')
    
    data=data.stack(ens_time=['ens', 'time'])
    for z, lag in enumerate(lags):
        print('...Looking at lag: '+str(lag))
        ind_temp=ind.shift(time=lag) #Positive means index is leading!
        ind_temp=ind_temp.stack(ens_time=['ens', 'time'])
        
        
        if detrend==True:
            res_temp=xr.corr(data, ind_temp, dim='ens_time')
        else:
            res_temp=xr.corr(data[datavar], ind_temp[indvar], dim='ens_time')

        if z==0:
            res=res_temp
        else:
            res=xr.concat([res, res_temp], dim='lag')
    
    print('...Plotting')
    fig, ax=plt.subplots(2,1, figsize=(20,12))
    
    if type(ymax)!=type(None):
        res=res.sel(YC=slice(-75.5, ymax))
    
    if members=='all':
        fig.suptitle(title+'\n All members')
    else:
        if members[1]-members[0]==1:
            fig.suptitle(title+'\n Ensemble Member '+ str(members[-1]))
        else:
            fig.suptitle(title+'\n Ensemble Members '+ str(members[0]+1)+' to '+str(members[-1]))
    
    if inverse==True:
        print('...NOTE: Taking the INVERSE of the correlation!')
        res=-1*res
    
    maxcor=res.max(dim='lag')
    #PLOT CORRELATION COEFFICIENTS
    img=ax[0].pcolor(res.XC, res.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
    plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
    CS=ax[0].contour(res.XC, res.YC, maxcor, np.arange(-1, 1, 0.25))
    ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)
    
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax[0].set_ylim([-75.5, ymax])
    ax[0].set_facecolor('grey')
    ax[0].set_xlabel('Longitude [deg]')
    ax[0].set_ylabel('Latitude [deg]')
    ax[0].set_title('Optimal Correlation Coefficient')
    
    if draw_box==True:
        ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
    
    
    print('...Finding optimal lag')
    idx=np.argmax([i=='lag' for i in res.dims])
    masked=np.ma.masked_invalid(res.to_numpy())
    opt_lag=np.nanargmax(masked, axis=idx)
    opt_lag=opt_lag.astype(float)
    opt_lag[masked.mask.max(axis=idx)]=np.nan
    opt_lag=opt_lag+lags[0]
    
    
    #PLOT OPTIMAL LAG
    img=ax[1].pcolor(res.XC, res.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
    plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
    #CS=ax[1].contour(res.XC, res.YC, opt_lag, np.arange(-1, 1, 0.25))
    #ax[1].clabel(CS, CS.levels, inline=True,fontsize=10)
    
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax[1].set_ylim([-75.5, ymax])
    ax[1].set_facecolor('grey')
    ax[1].set_xlabel('Longitude [deg]')
    ax[1].set_ylabel('Latitude [deg]')
    
    if draw_box==True:
        ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')
    
    fig.subplots_adjust(hspace=0.4)
    
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        if members=='all':
            fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_all.png')
        else:
            fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_member'+str(members[-1])+'.png')

    if give_back==False:
        return
    else:
        return data, ind, res
    
    
def correlationMap_1D_v2(members='all',
                      data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      detrend=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False):
    '''
    Input: Only one map
    New: using xskillscore
    
    '''
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 2!')
    plt.rcParams.update({'font.size': 14})
    
    #First, we read the data!
    if type(data)==type(None):
            print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_....')
            data=loadMITgcmData(filename=datafn, members=members, kind=datakind)
    
    print('...Making sure interpolate everything to the center of the cells')
    if 'XG' in data.dims:
        data=data.rename({'XG':'XC'})
    if 'YG' in data.dims:
        data=data.rename({'YG':'YC'})

    #Read something with the correct locations.
    data1=loadMITgcmData(filename='EXFuwind', members=[0,1], kind='maps')
    data=data.interp(XC=data1.XC, YC=data1.YC)
    
    
    if type(ind)==type(None):
        print('...Reading: '+indfn)
        if 'timeseries' in indfn:
            ind=loadTimeSeriesData(filename=indfn, members=members, kind='old')
        else:        
            ind=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+indfn+'.nc')
            if members!='all':
                print('...Selecting member '+str(members[0]+1))
                ind=ind.sel(ens=[members[0]])
            ind=ind.to_dataset()
        
    if detrend==True:
        print('...Detrending data')
        data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
        ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)
        
    if window!=1:
        print('...Applying rolling mean')
        data=data.rolling(time=window, center=True).mean()
        ind=ind.rolling(time=window, center=True).mean()
    
    print('...Computing correlation coefficients')
    
    data=data.stack(ens_time=['ens', 'time'])
    for z, lag in enumerate(lags):
        print('...Looking at lag: '+str(lag))
        ind_temp=ind.shift(time=lag) #Positive means index is leading!
        ind_temp=ind_temp.stack(ens_time=['ens', 'time'])
        
        
#         if detrend==True:
#             res_temp=xr.corr(data, ind_temp, dim='ens_time')
#         else:
#             res_temp=xr.corr(data[datavar], ind_temp[indvar], dim='ens_time')
        if detrend==True:
            res_temp=xs.pearson_r(data, ind_temp, dim='ens_time', skipna=True)
            p_temp=xs.pearson_r_p_value(data, ind_temp, dim='ens_time', skipna=True)
        else:
            res_temp=xs.pearson_r(data[datavar], ind_temp[indvar], dim='ens_time', skipna=True)
            p_temp=xs.pearson_r_p_value(data[datavar], ind_temp[indvar], dim='ens_time', skipna=True)
            
        if z==0:
            res=res_temp
            p=p_temp
            
        else:
            res=xr.concat([res, res_temp], dim='lag')
            p=xr.concat([p, p_temp], dim='lag')
            
    print('...Plotting')
    fig, ax=plt.subplots(2,1, figsize=(20,12))
    
    if type(ymax)!=type(None):
        res=res.sel(YC=slice(-75.5, ymax))
        p=p.sel(YC=slice(-75.5, ymax))
        
    if members=='all':
        fig.suptitle(title+'\n All members; p<0.05')
    else:
        if members[1]-members[0]==1:
            fig.suptitle(title+'\n Ensemble Member '+ str(members[-1])+'; p<0.05')
        else:
            fig.suptitle(title+'\n Ensemble Members '+ str(members[0]+1)+' to '+str(members[-1])+'; p<0.05')
    
    if inverse==True:
        print('...NOTE: Taking the INVERSE of the correlation!')
        res=-1*res
    
#     print('Now: p')
#     print(p)
#     print('Now: cor')
#     print(res)
    
    maxcor=res.max(dim='lag')
    mincor=res.min(dim='lag')
    maxp=p.min(dim='lag')
    
#     print('Now: maxp')
    print(maxp)
    
    maxcor=maxcor.where((maxp<0.05) | (np.isnan(maxp)), other=0)
    mincor=mincor.where((maxp<0.05) | (np.isnan(maxp)), other=0)
    
    maxcor.to_numpy()[maxcor<abs(mincor)]=mincor.to_numpy()[maxcor<abs(mincor)]
    
    #PLOT CORRELATION COEFFICIENTS
    img=ax[0].pcolor(res.XC, res.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
    plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
    CS=ax[0].contour(res.XC, res.YC, maxcor, np.arange(-1, 1, 0.25))
    ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)
    
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax[0].set_ylim([-75.5, ymax])
    ax[0].set_facecolor('grey')
    ax[0].set_xlabel('Longitude [deg]')
    ax[0].set_ylabel('Latitude [deg]')
    ax[0].set_title('Optimal Correlation Coefficient')
    
    if draw_box==True:
        ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
    
    
#     cmap = mpl.cm.seismic
#     bounds = np.arange(lags[0], lags[-1], 3)
#     norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
    
    
    print('...Finding optimal lag')
    idx=np.argmax([i=='lag' for i in res.dims])
    masked=np.ma.masked_invalid(res.where((maxp<0.05) & (maxp>0)).to_numpy())
    opt_lag=np.nanargmax(abs(masked), axis=idx)
    opt_lag=opt_lag.astype(float)
    opt_lag[masked.mask.max(axis=idx)]=np.nan
    opt_lag=opt_lag+lags[0]
    
    #PLOT OPTIMAL LAG
    img=ax[1].pcolor(res.XC, res.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
    plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
    img=ax[1].pcolor(maxp.where(maxp>0.05).XC, maxp.where(maxp>0.05).YC, maxp.where(maxp>0.05), color=['green'])
    #CS=ax[1].contour(res.XC, res.YC, opt_lag, np.arange(-1, 1, 0.25))
    #ax[1].clabel(CS, CS.levels, inline=True,fontsize=10)
    
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax[1].set_ylim([-75.5, ymax])
    ax[1].set_facecolor('grey')
    ax[1].set_xlabel('Longitude [deg]')
    ax[1].set_ylabel('Latitude [deg]')
    
    if draw_box==True:
        ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')
    
    fig.subplots_adjust(hspace=0.4)
    
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        if members=='all':
            fig.savefig('../03_output/correlation/'+today+'_1Dcormap_'+savefn+'_all.png')
        else:
            fig.savefig('../03_output/correlation/'+today+'_1Dcormap_'+savefn+'_member'+str(members[-1])+'.png')

    if give_back==False:
        return
    else:
        return data, ind, res
    
    
def correlationMap_2D_v2(members='all',
                      data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      detrend=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False):
    '''
    Input: Always give U component as input!
    
    '''
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 2!')
    plt.rcParams.update({'font.size': 14})
    
    #First, we read the data!
    if type(data)==type(None):
            print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_....')
            dataU=loadMITgcmData(filename=datafn, members=members, kind=datakind)
            datafn_v=datafn.replace('U','V')
            datafn_v=datafn_v.replace('u','v')
            datafn_v=datafn_v.replace('x','y')
            print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn_v+'_....')
            dataV=loadMITgcmData(filename=datafn_v, members=members, kind=datakind)
    
    else:
        dataU=data[0]
        dataV=data[1]
    
    print('...Making sure interpolate everything to the center of the cells')
    if 'XG' in dataU.dims:
        dataU=dataU.rename({'XG':'XC'})
    if 'YG' in dataU.dims:
        dataU=dataU.rename({'YG':'YC'})
    if 'XG' in dataV.dims:
        dataV=dataV.rename({'XG':'XC'})
    if 'YG' in dataV.dims:
        dataV=dataV.rename({'YG':'YC'})
    
    #Read something with the correct locations.
    data=loadMITgcmData(filename='EXFuwind', members=[0,1], kind='maps')
    
    dataU=dataU.interp(XC=data.XC, YC=data.YC)
    dataV=dataV.interp(XC=data.XC, YC=data.YC)
    
    print('...Compute the absolute combined value')
    datavar_v=datavar.replace('U','V')
    datavar_v=datavar_v.replace('u','v')
    datavar_v=datavar_v.replace('x','y')
    data=np.sqrt(dataU[datavar]**2+dataV[datavar_v]**2)
    data=data.to_dataset(name=datavar)
    
    
    
    if type(ind)==type(None):
        print('...Reading: '+indfn)
        if 'timeseries' in indfn:
            ind=loadTimeSeriesData(filename=indfn, members=members, kind='old')
        else:        
            ind=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+indfn+'.nc')
            if members!='all':
                print('...Selecting member '+str(members[0]+1))
                ind=ind.sel(ens=[members[0]])
            ind=ind.to_dataset()
        
    if detrend==True:
        print('...Detrending data')
        data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
        ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)
        
    if window!=1:
        print('...Applying rolling mean')
        data=data.rolling(time=window, center=True).mean()
        ind=ind.rolling(time=window, center=True).mean()
    
    print('...Computing correlation coefficients')
    
    data=data.stack(ens_time=['ens', 'time'])
    for z, lag in enumerate(lags):
        print('...Looking at lag: '+str(lag))
        ind_temp=ind.shift(time=lag) #Positive means index is leading!
        ind_temp=ind_temp.stack(ens_time=['ens', 'time'])
        
        if detrend==True:
            res_temp=xs.pearson_r(data, ind_temp, dim='ens_time', skipna=True)
            p_temp=xs.pearson_r_p_value(data, ind_temp, dim='ens_time', skipna=True)
        else:
            res_temp=xs.pearson_r(data[datavar], ind_temp[indvar], dim='ens_time', skipna=True)
            p_temp=xs.pearson_r_p_value(data[datavar], ind_temp[indvar], dim='ens_time', skipna=True)
            
        if z==0:
            res=res_temp
            p=p_temp
            
        else:
            res=xr.concat([res, res_temp], dim='lag')
            p=xr.concat([p, p_temp], dim='lag')
            
    print('...Plotting')
    fig, ax=plt.subplots(2,1, figsize=(20,12))
    
    if type(ymax)!=type(None):
        res=res.sel(YC=slice(-75.5, ymax))
        p=p.sel(YC=slice(-75.5, ymax))
        
    if members=='all':
        fig.suptitle(title+'\n All members; p<0.05')
    else:
        if members[1]-members[0]==1:
            fig.suptitle(title+'\n Ensemble Member '+ str(members[-1])+'; p<0.05')
        else:
            fig.suptitle(title+'\n Ensemble Members '+ str(members[0]+1)+' to '+str(members[-1])+'; p<0.05')
    
    if inverse==True:
        print('...NOTE: Taking the INVERSE of the correlation!')
        res=-1*res
    
#     print('Now: p')
#     print(p)
#     print('Now: cor')
#     print(res)
    
    maxcor=res.max(dim='lag')
    mincor=res.min(dim='lag')
    maxp=p.min(dim='lag')
    
#     print('Now: maxp')
    print(maxp)
    
    maxcor=maxcor.where((maxp<0.05) | (np.isnan(maxp)), other=0)
    mincor=mincor.where((maxp<0.05) | (np.isnan(maxp)), other=0)
    
    maxcor.to_numpy()[maxcor<abs(mincor)]=mincor.to_numpy()[maxcor<abs(mincor)]
    
    #PLOT CORRELATION COEFFICIENTS
    img=ax[0].pcolor(res.XC, res.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
    plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
    CS=ax[0].contour(res.XC, res.YC, maxcor, np.arange(-1, 1, 0.25))
    ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)
    
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax[0].set_ylim([-75.5, ymax])
    ax[0].set_facecolor('grey')
    ax[0].set_xlabel('Longitude [deg]')
    ax[0].set_ylabel('Latitude [deg]')
    ax[0].set_title('Optimal Correlation Coefficient')
    
    if draw_box==True:
        ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
    
    
#     cmap = mpl.cm.seismic
#     bounds = np.arange(lags[0], lags[-1], 3)
#     norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
    
    
    print('...Finding optimal lag')
    idx=np.argmax([i=='lag' for i in res.dims])
    masked=np.ma.masked_invalid(res.where((maxp<0.05) & (maxp>0)).to_numpy())
    opt_lag=np.nanargmax(abs(masked), axis=idx)
    opt_lag=opt_lag.astype(float)
    opt_lag[masked.mask.max(axis=idx)]=np.nan
    opt_lag=opt_lag+lags[0]
    
    #PLOT OPTIMAL LAG
    img=ax[1].pcolor(res.XC, res.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
    plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
    img=ax[1].pcolor(maxp.where(maxp>0.05).XC, maxp.where(maxp>0.05).YC, maxp.where(maxp>0.05), color=['green'])
    #CS=ax[1].contour(res.XC, res.YC, opt_lag, np.arange(-1, 1, 0.25))
    #ax[1].clabel(CS, CS.levels, inline=True,fontsize=10)
    
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax[1].set_ylim([-75.5, ymax])
    ax[1].set_facecolor('grey')
    ax[1].set_xlabel('Longitude [deg]')
    ax[1].set_ylabel('Latitude [deg]')
    
    if draw_box==True:
        ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')
    
    fig.subplots_adjust(hspace=0.4)
    
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        if members=='all':
            fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_all.png')
        else:
            fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_member'+str(members[-1])+'.png')

    if give_back==False:
        return
    else:
        return data, ind, res
    
    

def correlationMap_1D_v3(members='all',
                      data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      mask='land', years=['1920', '2013'],
                      detrend=True, deseasonalize=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False, save_results=True, give_back=False):
    '''
    Input: Only one map
    New: using xskillscore
    
    '''
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 14})
    if members=='all':
        members=[0,20]
    
    data_og=data
    ind_og=ind
    for ens in range(members[0]+1,members[-1]+1):
        #First, we read the data!
        if type(data_og)==type(None):
                print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_....')
                data=loadMITgcmData(filename=datafn, members=[ens-1, ens], kind=datakind)

        print('...Making sure interpolate everything to the center of the cells')
        if 'XG' in data.dims:
            data=data.rename({'XG':'XC'})
        if 'YG' in data.dims:
            data=data.rename({'YG':'YC'})

        #Read something with the correct locations.
        data1=loadMITgcmData(filename='EXFuwind', members=[0,1], kind='maps')
        data=data.interp(XC=data1.XC, YC=data1.YC)

        #Reading Grid
        gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
        grid = Grid(gp)

        if mask=='land':
            print('...Applying land mask.')
            data=data.where(grid.land_mask==False)
        elif mask=='ice':
            print('...Applying ice mask.')
            data=data.where(grid.ice_mask==False)
        elif mask=='land_and_ice':
            print('...Applying land AND ice mask.')
            data=data.where(grid.land_mask==False)
            data=data.where(grid.ice_mask==False)    

        if type(ind_og)==type(None):
            print('...Reading: '+indfn)
            if 'timeseries' in indfn:
                ind=loadTimeSeriesData(filename=indfn, members=[ens-1, ens], kind='old')
            else:        
                ind=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+indfn+'.nc')
                if members!='all':
                    print('...Selecting member '+str(members[0]+1))
                    ind=ind.sel(ens=[members[0]])
                ind=ind.to_dataset()

        #Selecting the proper years:
        print('...Selecting times between: '+str(years[0])+' and '+str(years[-1]))
        data=data.sel(time=slice(years[0], years[-1]))
        ind=ind.sel(time=slice(years[0], years[-1]))


        if detrend==True:
            print('...Detrending data')
            data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
            ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)
            
        if deseasonalize==True:
            print('...Deseasonalizing the data')
            monthlymean=data.groupby('time.month').mean("time")
            data=data.groupby('time.month')-monthlymean
            
            monthlymean=ind.groupby('time.month').mean("time")
            ind=ind.groupby('time.month')-monthlymean

        if window!=1:
            print('...Applying rolling mean')
            data=data.rolling(time=window, center=True).mean()
            ind=ind.rolling(time=window, center=True).mean()

        print('...Computing correlation coefficients')

        data=data.stack(ens_time=['ens', 'time'])
        for z, lag in enumerate(lags):
            print('...Looking at lag: '+str(lag))
            ind_temp=ind.shift(time=lag) #Positive means index is leading!
            ind_temp=ind_temp.stack(ens_time=['ens', 'time'])


    #         if detrend==True:
    #             res_temp=xr.corr(data, ind_temp, dim='ens_time')
    #         else:
    #             res_temp=xr.corr(data[datavar], ind_temp[indvar], dim='ens_time')
            if detrend==True:
                res_temp=xs.pearson_r(data, ind_temp, dim='ens_time', skipna=True)
                p_temp=xs.pearson_r_p_value(data, ind_temp, dim='ens_time', skipna=True)
            else:
                res_temp=xs.pearson_r(data[datavar], ind_temp[indvar], dim='ens_time', skipna=True)
                p_temp=xs.pearson_r_p_value(data[datavar], ind_temp[indvar], dim='ens_time', skipna=True)

            if z==0:
                res=res_temp
                p=p_temp

            else:
                res=xr.concat([res, res_temp], dim='lag')
                p=xr.concat([p, p_temp], dim='lag')

        print('...Plotting')
        fig, ax=plt.subplots(2,1, figsize=(20,12))

        if type(ymax)!=type(None):
            res=res.sel(YC=slice(-75.5, ymax))
            p=p.sel(YC=slice(-75.5, ymax))

        #if members=='all':
        #    fig.suptitle(title+'\n All members; p<0.05')
        #else:
        #    if members[1]-members[0]==1:
        fig.suptitle(title+'\n Ensemble Member '+ str(ens)+'; p<0.05')
        #    else:
        #        fig.suptitle(title+'\n Ensemble Members '+ str(members[0]+1)+' to '+str(members[-1])+'; p<0.05')

        if inverse==True:
            print('...NOTE: Taking the INVERSE of the correlation!')
            res=-1*res

    #     print('Now: p')
    #     print(p)
    #     print('Now: cor')
    #     print(res)

        maxcor=res.max(dim='lag')
        mincor=res.min(dim='lag')
        maxp=p.min(dim='lag')

    #     print('Now: maxp')
    #    print(maxp)

        maxcor=maxcor.where((maxp<=0.05) | (np.isnan(maxp)), other=0)
        mincor=mincor.where((maxp<=0.05) | (np.isnan(maxp)), other=0)

        maxcor.to_numpy()[maxcor<abs(mincor)]=mincor.to_numpy()[maxcor<abs(mincor)]

        #PLOT CORRELATION COEFFICIENTS
        img=ax[0].pcolor(res.XC, res.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
        plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
        CS=ax[0].contour(res.XC, res.YC, maxcor, np.arange(-1, 1, 0.25))
        ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)

        bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
        plt.clabel(bath, [-1000])
        if type(ymax)!=type(None):
            ax[0].set_ylim([-75.5, ymax])
        if 'ice' in mask:
            ax[0].pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
        ax[0].set_facecolor('grey')
        ax[0].set_xlabel('Longitude [deg]')
        ax[0].set_ylabel('Latitude [deg]')
        ax[0].set_title('Optimal Correlation Coefficient')

        if draw_box==True:
            ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')



    #     cmap = mpl.cm.seismic
    #     bounds = np.arange(lags[0], lags[-1], 3)
    #     norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')


        print('...Finding optimal lag')
        idx=np.argmax([i=='lag' for i in res.dims])
        masked=np.ma.masked_invalid(res.where((maxp<=0.05) & (maxp>=0)).to_numpy())
        opt_lag=np.nanargmax(abs(masked), axis=idx)
        opt_lag=opt_lag.astype(float)
        opt_lag[masked.mask.max(axis=idx)]=np.nan
        opt_lag=opt_lag+lags[0]

        #PLOT OPTIMAL LAG
        img=ax[1].pcolor(res.XC, res.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
        plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
        img=ax[1].pcolor(maxp.where(maxp>0.05).XC, maxp.where(maxp>0.05).YC, maxp.where(maxp>0.05)*0, cmap='ocean')
        #CS=ax[1].contour(res.XC, res.YC, opt_lag, np.arange(-1, 1, 0.25))
        #ax[1].clabel(CS, CS.levels, inline=True,fontsize=10)

    #     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    #     grid = Grid(gp)
        bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
        plt.clabel(bath, [-1000])
        if type(ymax)!=type(None):
            ax[1].set_ylim([-75.5, ymax])

        if 'ice' in mask:
             ax[1].pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
        ax[1].set_facecolor('grey')
        ax[1].set_xlabel('Longitude [deg]')
        ax[1].set_ylabel('Latitude [deg]')

        if draw_box==True:
            ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
        ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')

        fig.subplots_adjust(hspace=0.4)

        if save==True:
            print('...Saving figure')
            from datetime import date
            today = date.today()
            today=today.strftime("%Y%m%d")
            #if members=='all':
            #    fig.savefig('../03_output/correlation/'+today+'_1Dcormap_'+savefn+'_all.png')
            #else:
            fig.savefig('../03_output/correlation/'+today+'_1Dcormap_'+savefn+'_member'+str(ens)+'.png')

        if save_results==True:
            np.save('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy', opt_lag)
            maxp.to_netcdf('../02_data/correlation/maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            maxcor.to_netcdf('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        
    return
  
    
    
def correlationMap_2D_v3(members='all',
                      data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      detrend=True, deseasonalize=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      mask='land', years=['1920', '2013'],
                      inverse=False, save_results=True,
                      give_back=False):
    '''
    Input: Always give U component as input!
    
    '''
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 14})
    
    if members=='all':
        members=[0,20]
    
    data_og=data
    ind_og=ind
    for ens in range(members[0]+1,members[-1]+1):
        #First, we read the data!
        if type(data_og)==type(None):
                print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_....')
                dataU=loadMITgcmData(filename=datafn, members=[ens-1, ens], kind=datakind)
                datafn_v=datafn.replace('U','V')
                datafn_v=datafn_v.replace('u','v')
                datafn_v=datafn_v.replace('x','y')
                print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn_v+'_....')
                dataV=loadMITgcmData(filename=datafn_v, members=[ens-1, ens], kind=datakind)

        else:
            dataU=data[0]
            dataV=data[1]

        print('...Making sure interpolate everything to the center of the cells')
        if 'XG' in dataU.dims:
            dataU=dataU.rename({'XG':'XC'})
        if 'YG' in dataU.dims:
            dataU=dataU.rename({'YG':'YC'})
        if 'XG' in dataV.dims:
            dataV=dataV.rename({'XG':'XC'})
        if 'YG' in dataV.dims:
            dataV=dataV.rename({'YG':'YC'})

        #Read something with the correct locations.
        data=loadMITgcmData(filename='EXFuwind', members=[0,1], kind='maps')

        dataU=dataU.interp(XC=data.XC, YC=data.YC)
        dataV=dataV.interp(XC=data.XC, YC=data.YC)

        print('...Compute the absolute combined value')
        datavar_v=datavar.replace('U','V')
        datavar_v=datavar_v.replace('u','v')
        datavar_v=datavar_v.replace('x','y')
        data=np.sqrt(dataU[datavar]**2+dataV[datavar_v]**2)
        data=data.to_dataset(name=datavar)
        
        #Reading Grid
        gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
        grid = Grid(gp)

        if mask=='land':
            print('...Applying land mask.')
            data=data.where(grid.land_mask==False)
        elif mask=='ice':
            print('...Applying ice mask.')
            data=data.where(grid.ice_mask==False)
        elif mask=='land_and_ice':
            print('...Applying land AND ice mask.')
            data=data.where(grid.land_mask==False)
            data=data.where(grid.ice_mask==False)    


        if type(ind_og)==type(None):
            print('...Reading: '+indfn)
            if 'timeseries' in indfn:
                ind=loadTimeSeriesData(filename=indfn, members=[ens-1, ens], kind='old')
            else:        
                ind=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+indfn+'.nc')
                if members!='all':
                    print('...Selecting member '+str(members[0]+1))
                    ind=ind.sel(ens=[members[0]])
                ind=ind.to_dataset()
        
        #Selecting the proper years:
        print('...Selecting times between: '+str(years[0])+' and '+str(years[-1]))
        data=data.sel(time=slice(years[0], years[-1]))
        ind=ind.sel(time=slice(years[0], years[-1]))

        if detrend==True:
            print('...Detrending data')
            data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
            ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)

        if deseasonalize==True:
            print('...Deseasonalizing the data')
            monthlymean=data.groupby('time.month').mean("time")
            data=data.groupby('time.month')-monthlymean
            
            monthlymean=ind.groupby('time.month').mean("time")
            ind=ind.groupby('time.month')-monthlymean
            
        if window!=1:
            print('...Applying rolling mean')
            data=data.rolling(time=window, center=True).mean()
            ind=ind.rolling(time=window, center=True).mean()

        print('...Computing correlation coefficients')

        data=data.stack(ens_time=['ens', 'time'])
        for z, lag in enumerate(lags):
            print('...Looking at lag: '+str(lag))
            ind_temp=ind.shift(time=lag) #Positive means index is leading!
            ind_temp=ind_temp.stack(ens_time=['ens', 'time'])

            if detrend==True:
                res_temp=xs.pearson_r(data, ind_temp, dim='ens_time', skipna=True)
                p_temp=xs.pearson_r_p_value(data, ind_temp, dim='ens_time', skipna=True)
            else:
                res_temp=xs.pearson_r(data[datavar], ind_temp[indvar], dim='ens_time', skipna=True)
                p_temp=xs.pearson_r_p_value(data[datavar], ind_temp[indvar], dim='ens_time', skipna=True)

            if z==0:
                res=res_temp
                p=p_temp

            else:
                res=xr.concat([res, res_temp], dim='lag')
                p=xr.concat([p, p_temp], dim='lag')

        print('...Plotting')
        fig, ax=plt.subplots(2,1, figsize=(20,12))

        if type(ymax)!=type(None):
            res=res.sel(YC=slice(-75.5, ymax))
            p=p.sel(YC=slice(-75.5, ymax))

        fig.suptitle(title+'\n Ensemble Member '+ str(ens)+'; p<0.05')
        #    else:
        #        fig.suptitle(title+'\n Ensemble Members '+ str(members[0]+1)+' to '+str(members[-1])+'; p<0.05')

        if inverse==True:
            print('...NOTE: Taking the INVERSE of the correlation!')
            res=-1*res

    #     print('Now: p')
    #     print(p)
    #     print('Now: cor')
    #     print(res)

        maxcor=res.max(dim='lag')
        mincor=res.min(dim='lag')
        maxp=p.min(dim='lag')

    #     print('Now: maxp')
    #    print(maxp)

        maxcor=maxcor.where((maxp<=0.05) | (np.isnan(maxp)), other=0)
        mincor=mincor.where((maxp<=0.05) | (np.isnan(maxp)), other=0)

        maxcor.to_numpy()[maxcor<abs(mincor)]=mincor.to_numpy()[maxcor<abs(mincor)]

        #PLOT CORRELATION COEFFICIENTS
        img=ax[0].pcolor(res.XC, res.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
        plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
        CS=ax[0].contour(res.XC, res.YC, maxcor, np.arange(-1, 1, 0.25))
        ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)

        bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
        plt.clabel(bath, [-1000])
        if type(ymax)!=type(None):
            ax[0].set_ylim([-75.5, ymax])
        if 'ice' in mask:
            ax[0].pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
        ax[0].set_facecolor('grey')
        ax[0].set_xlabel('Longitude [deg]')
        ax[0].set_ylabel('Latitude [deg]')
        ax[0].set_title('Optimal Correlation Coefficient')

        if draw_box==True:
            ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')



    #     cmap = mpl.cm.seismic
    #     bounds = np.arange(lags[0], lags[-1], 3)
    #     norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')


        print('...Finding optimal lag')
        idx=np.argmax([i=='lag' for i in res.dims])
        masked=np.ma.masked_invalid(res.where((maxp<=0.05) & (maxp>=0)).to_numpy())
        opt_lag=np.nanargmax(abs(masked), axis=idx)
        opt_lag=opt_lag.astype(float)
        opt_lag[masked.mask.max(axis=idx)]=np.nan
        opt_lag=opt_lag+lags[0]

        #PLOT OPTIMAL LAG
        img=ax[1].pcolor(res.XC, res.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
        plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
        img=ax[1].pcolor(maxp.where(maxp>0.05).XC, maxp.where(maxp>0.05).YC, maxp.where(maxp>0.05)*0, cmap='ocean')
        #CS=ax[1].contour(res.XC, res.YC, opt_lag, np.arange(-1, 1, 0.25))
        #ax[1].clabel(CS, CS.levels, inline=True,fontsize=10)

    #     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    #     grid = Grid(gp)
        bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
        plt.clabel(bath, [-1000])
        if type(ymax)!=type(None):
            ax[1].set_ylim([-75.5, ymax])

        if 'ice' in mask:
             ax[1].pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
        ax[1].set_facecolor('grey')
        ax[1].set_xlabel('Longitude [deg]')
        ax[1].set_ylabel('Latitude [deg]')

        if draw_box==True:
            ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
        ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')

        fig.subplots_adjust(hspace=0.4)

        if save==True:
            print('...Saving figure')
            from datetime import date
            today = date.today()
            today=today.strftime("%Y%m%d")
            #if members=='all':
            #    fig.savefig('../03_output/correlation/'+today+'_1Dcormap_'+savefn+'_all.png')
            #else:
            fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_member'+str(ens)+'.png')

        if save_results==True:
            np.save('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy', opt_lag)
            maxp.to_netcdf('../02_data/correlation/2D_maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            maxcor.to_netcdf('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')

    if give_back==False:
        return
    else:
        return data, ind, res
    
    
#%%---------------------------------------------------------------------------------

def showCombinedCorrelationMap_2D_v3(members='all',
                      data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      detrend=True, deseasonalize=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      mask='land', years=['1920', '2013'],
                      inverse=False, save_results=True,
                      give_back=False):
    '''
    Input: Always give U component as input!
    
    '''
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 22})
    
    if members=='all':
        members=[0,20]
    
    for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
        if deseasonalize==True:
            try:
                temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            except:
                temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxp=xr.open_dataarray('../02_data/correlation/maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        else:
            try:
                temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            except:
                temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxp=xr.open_dataarray('../02_data/correlation/maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        if i==0:
            opt_lag=temp_opt_lag
            maxp=temp_maxp
            maxcor=temp_maxcor
        else:
            opt_lag=np.dstack((opt_lag, temp_opt_lag))
            maxp=xr.concat([maxp, temp_maxp], dim='ens')
            maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
    if i==0:
        maxp=maxp.expand_dims('ens')
        maxcor=maxcor.expand_dims('ens')
    if members[-1]!=21:
        std_lag=np.nanstd(opt_lag, axis=2)
        opt_lag=np.nanmean(opt_lag, axis=2)
    else:
        std_lag=np.zeros(np.shape(opt_lag))
    
    stdp=maxp.mean(dim='ens')
    maxp=maxp.mean(dim='ens')
    
    stdcor=maxcor.std(dim='ens')
    maxcor=maxcor.mean(dim='ens')
    
    #Reading Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    print('...Plotting')
    fig, ax=plt.subplots(2,1, figsize=(20,14))
    
    if members[-1]==21:
        fig.suptitle(title+'\n Ensemble Mean; '+str(window)+' month rolling mean; p<0.05')
    else:
        fig.suptitle(title+'\n Mean over members '+ str(members[0]+1)+' to '+ str(members[-1])+'; '+str(window)+' month rolling mean; p<0.05')


    #PLOT CORRELATION COEFFICIENTS
    img=ax[0].pcolor(maxcor.XC, maxcor.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
    plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
    CS=ax[0].contour(maxcor.XC, maxcor.YC, maxcor, [-0.5, 0.5], colors='yellow')
    ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)

    bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax[0].set_ylim([-75.5, ymax])
    if 'ice' in mask:
        ax[0].pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
    ax[0].set_facecolor('grey')
    ax[0].set_xlabel('Longitude [deg]')
    ax[0].set_ylabel('Latitude [deg]')
    ax[0].set_title('Optimal Correlation Coefficient')

    if draw_box==True:
        ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')

    #PLOT OPTIMAL LAG
    img=ax[1].pcolor(maxcor.XC, maxcor.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
    plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
    img=ax[1].pcolor(maxp.where(maxp>0.05).XC, maxp.where(maxp>0.05).YC, maxp.where(maxp>0.05)*0, cmap='ocean')

    bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax[1].set_ylim([-75.5, ymax])

    if 'ice' in mask:
        ax[1].pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
    ax[1].set_facecolor('grey')
    ax[1].set_xlabel('Longitude [deg]')
    ax[1].set_ylabel('Latitude [deg]')

    if draw_box==True:
        ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')

    fig.subplots_adjust(hspace=0.5, top=0.85)

    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
    fig=plt.figure(figsize=(20,10))
    plt.suptitle('Extra statistics over ensemble members for : \n'+title)
    plt.subplot(2,2,1)
    img=plt.pcolor(maxp.XC, maxp.YC, maxp, vmax=0.05, cmap='seismic')
    plt.colorbar(img, label='Pvalue')
    #CS=plt.contour(maxp.XC, maxp.YC, maxp, [1e-4], colors='yellow')
    #plt.clabel(CS, CS.levels, inline=True,fontsize=10)
    plt.title('Mean P value with highest lag')# \n Yellow line indicates 1e-4')
    
    plt.subplot(2,2,2)
    img=plt.pcolor(stdp.XC, stdp.YC, stdp, vmax=0.05, cmap='seismic')
    plt.colorbar(img, label='Std of Pvalue')
    #CS=plt.contour(stdp.XC, stdp.YC, stdp, [1e-4], colors='yellow')
    #plt.clabel(CS, CS.levels, inline=True,fontsize=10)
    plt.title('Standard deviation of P value')# \n Yellow line indicates 1e-4')
    
    plt.subplot(2,2,3)
    img=plt.pcolor(stdcor.XC, stdcor.YC, stdcor, cmap='seismic', vmax=0.2)
    plt.colorbar(img, label='Std of cor')
    #CS=plt.contour(stdcor.XC, stdcor.YC, stdcor, [0.1], colors='yellow')
    #plt.clabel(CS, CS.levels, inline=True,fontsize=10)
    plt.title('Standard deviation of correlation coefficient')#  \n Yellow line indicates 0.1')
    
    plt.subplot(2,2,4)
    img=plt.pcolor(stdcor.XC, stdcor.YC, std_lag, cmap='seismic', vmax=10)
    plt.colorbar(img, label='Std of lag')
    #CS=plt.contour(stdcor.XC, stdcor.YC, std_lag, [4], colors='yellow')
    #plt.clabel(CS, CS.levels, inline=True,fontsize=10)
    plt.title('Standard deviation of lag')# \n Yellow line indicates 4 months')
    
    fig.subplots_adjust(hspace=0.5, wspace=0.3, top=0.8)
    
    if save==True:
        print('...Saving second figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_BONUS_2Dcormap_'+savefn+'_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
    return


