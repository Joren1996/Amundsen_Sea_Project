# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain all functions to related to checking for convection. This is shown in the appendix. Note that functions are virtually similar to the normal correlation maps.....

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
from functions.loading_and_processing_data import loadMITgcmData,loadTimeSeriesData, detrend_and_average_MITgcmData, read_all_data
from functions.composite_maps import createCompositeMap, createVelocityCompositeMap, defineEvents, findEvents
import sys
import xskillscore as xs
from functions.loading_and_processing_data import loadMITgcmData, loadPACEData, detrend_and_average_MITgcmData, read_all_data, loadTimeSeriesData
sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid

divnorm = colors.TwoSlopeNorm(vcenter=0)
plt.rcParams.update({'font.size': 18})

def find_convection_times(isotherm, method='full', recov=60, fname='', save=True):
    '''
    Find times when convection takes place.
    
    INPUT:
    isotherm (xr.DataArray()): contains isotherms
    method (string): How do we want to compute convection 
                    (full: all definitions, 
                     pine: only using isotherms in Pine Island Bay
                     dotson: only using isotherms in Dotson Bay
                     amundsen: only using isotherms in amundsen shelf)
    recov (int): recovery period
    fname (str): file name of figure to be saved when save==True
    save (boolean): Do we want to save the figure?
    
    OUTPUT:
    convection (xr.DataArray): boolean indicating when convection takes place according to proposed definition.
    '''
    if method=='full':
        convection=(isotherm['amundsen_shelf_isotherm_0.5C_below_100m']\
                    .rolling(time=24, center=True).mean()<-440) &\
                    (isotherm['pine_island_bay_isotherm_0C_below_100m']\
                     .rolling(time=24, center=True).mean()<-470) & \
                    (isotherm['dotson_bay_isotherm_-1C_below_100m']\
                     .rolling(time=24, center=True).mean()<-430)
    elif method=='pine':
        convection=isotherm['pine_island_bay_isotherm_0C_below_100m']\
                     .rolling(time=24, center=True).mean()<-470
    elif method=='dotson':
        convection=isotherm['dotson_bay_isotherm_-1C_below_100m']\
                     .rolling(time=24, center=True).mean()<-430
        
    elif method=='amundsen':
        convection=isotherm['amundsen_shelf_isotherm_0.5C_below_100m']\
                    .rolling(time=24, center=True).mean()<-440
    
    plt.figure(figsize=(10,7.5))
    ax=plt.subplot(1,1,1)
    a=convection
    a=0.5*a+a.ens
    a.plot.line(x='time')
    plt.yticks(range(0,20), range(1,21))
    plt.title('Deep Convection in the MITgcm model')
    plt.ylabel('Member')
    ax.get_legend().remove()
    
    if recov!=0:
        convection=convection.to_numpy()
        recov_convection=np.zeros(np.shape(convection))
        for i in range(1, recov+1):
            recov_convection[:,recov:]+=convection[:,recov-i:-i]
        recov_convection=recov_convection.astype(bool)
        
        convection=xr.DataArray(data=recov_convection, dims=["ens", "time"],
                                coords=dict(ens=isotherm.ens, time=isotherm.time))
        
        a=convection
        a=0.5*a+a.ens
        a.plot.line(x='time', ls='--')
        ax.get_legend().remove()
        
        
    plt.yticks(range(0,20), range(1,21))
    plt.title('Selected Deep Convection in the MITgcm model')
    plt.ylabel('Member')
    plt.xlim([convection.sel(time='1920-01').time.values[0], \
              convection.sel(time='2013-12').time.values[0]])
        
    convection=convection.sel(time=slice('1920', '2013'))
    if save==True:
        plt.savefig('../03_output/convection/conv_times_'+method+'_rec'+str(recov)+'_'+fname+'.png')
    return convection


def plot_isotherm_convection(isotherm, convection, members=[0,20], recov=60, fname='', save=True):
    '''
    Plot isotherms without selected convection periods.
    
    INPUT:
    isotherm (xr.DataArray()): contains isothermsa
    convection (xr.DataArray): boolean indicating when convection takes place according to proposed definition.
    members (list): members to choose [0,1] is member 1; [0,20] is all 20 members.
    recov (int): recovery period [months]
    fname (str): file name of figure to be saved when save==True
    save (boolean): Do we want to save the figure?
    
    '''
    
    plt.figure(figsize=(15,15))
    isotherm=isotherm.sel(time=slice('1920', '2013'))
    iso_conv=isotherm.where(convection==False)
    
    for i, var in enumerate(['amundsen_shelf_isotherm_0.5C_below_100m',
                             'pine_island_bay_isotherm_0C_below_100m',
                             'dotson_bay_isotherm_-1C_below_100m']):
        plt.subplot(3,1,i+1)
        for ens in range(members[0], members[-1]):
            isotherm[var].sel(ens=ens).rolling(time=24,
            center=True).mean().plot(lw=0.2, c='grey')
        for ens in range(members[0], members[-1]):
            iso_conv[var].sel(ens=ens).rolling(time=24,
            center=True).mean().plot(lw=1, c='blue')
        plt.title(isotherm[var].long_name)
    plt.subplots_adjust(hspace=0.8)
    plt.suptitle('Isotherms without Selected Convection Periods')
    if save==True:
        plt.savefig('../03_output/convection/conv_isoth_rec'+str(recov)+'_'+fname+'.png')
    return


def convectionCorMap(members='all',
    data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
    ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
    mask='land', years=['1920', '2013'],
    detrend=True, deseasonalize=True, window=24, title='Convection Correlation Map', ymax=-70, save=False, savefn='test',
    draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
    inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False):
    
    '''
    Create a correlation map with removal of convective periods 
    with one spatial variable with one component (e.g. Temperature, not velocity).
    
    INPUT:
    members (string or list): members to use ('all' or [0,1] for member 1)
    data (None or xr.DataSet): Spatial 2D data field. If 'None', then it will be read within the function.
    datafn (string): filename of the spatial 2D data field to be read.
    datavar (string): variable to be read.
    datakind (string): subfolder of 02_data to be read.
    ind (None or xr.DataSet): 1D time series to correlate with.
    indfn (string): filename of 1D time series to be read.
    indvar (string): variable in index to be used.
    lags (list): lags to be used for the cross correlation
    mask (string): whether we want to mask out the 'land', 'ice', or 'land and ice' in the figure.
    years (list): years that are at the edges of the time interval we are interested in.
    detrend (boolean): do we want to detrend both the index and the spatial map
    deseasonalize (boolean): do we want to deseasonalize the data.
    window (int): what kind of running mean do we want to apply? [months]
    title (string): titel of the correlation map to use
    ymax (float): northern limit of the map that we are interested in. [deg N]
    save (boolean): doe we want to save the figure?
    savefn (string): extra text for the filename of the figure.
    draw_box (boolean): do we want to draw a box in the domian.
    box_x (list): x limits of the box [deg E]
    box_y (list): y limits of the box [deg N]
    inverse (boolean): do we want to invert the correlation?
    save_results (boolean): do we want to save the correlation coefficients and the lag
    give_back (boolean): whether we want to return the data.
    cumsum_map (boolean): make the spatial variable cumulative
    cumsum_ind (boolean): make the 1D time series cumulative

    OUTPUT:
    '''
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 14})
    if members=='all':
        members=[0,20]
        
    print('Start Preparing for Convection....')
    print('...Read Isotherms')
    isotherm=loadTimeSeriesData(filename='timeseries_isotherm', members=members, kind='old')
    time=isotherm.time   
    isotherm['time']=isotherm.indexes['time'].to_datetimeindex()
     
    print('...Select Convective Events')
    conv=find_convection_times(isotherm, method='full', recov=60, save=save,\
                              fname=datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window))
    plot_isotherm_convection(isotherm, convection=conv, members=members, save=save,\
                             fname=datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window))
    print('Convection done!')
    conv['time']=time.sel(time=slice(years[0], years[-1]))
    conv=conv.rename({'ens':'mem'})
    
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
                if ('final' in indfn) or ('forcing' in indfn):
                    ind=loadTimeSeriesData(filename=indfn, members=[ens-1, ens], kind='old')
                else:
                    ind=loadTimeSeriesData(filename=indfn, members=[ens-1, ens], kind='timeseries')
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

        fig=plt.figure(figsize=(15,7.5))
        plt.subplot(2,1,1)
        plt.plot(data[datavar].sel(ens=0).rolling(time=24, center=True).mean().values[:,190,300], c='grey', lw=0.5)
        plt.subplot(2,1,2)
        plt.plot(ind[indvar].sel(ens=0).rolling(time=24, center=True).mean().values, c='grey', lw=0.5)
        
        #Remove Convective Periods
        print('...Remove Convective Periods!')
        data=data.where(conv.sel(mem=ens-1)==False)
        ind=ind.where(conv.sel(mem=ens-1)==False)
        
        plt.subplot(2,1,1)
        plt.plot(data[datavar].sel(ens=0).rolling(time=24, center=True).mean().values[:,190,300], 'b', lw=2)
        plt.title(datavar)
        plt.grid(True)
        plt.subplot(2,1,2)
        plt.plot(ind[indvar].sel(ens=0).rolling(time=24, center=True).mean().values, 'r', lw=2)
        plt.title(indvar)
        plt.suptitle('Removed Convective Periods for: '+datafn+' vs '+indvar+'\n Member '+str(ens))
        plt.grid(True)
        plt.subplots_adjust(hspace=0.5, wspace=0.5)
        if save==True:
            plt.savefig('../03_output/convection/convmember_'+str(ens)+'_'\
                +datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'.png')
        
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
    
        if cumsum_map==True:
            data=data.cumsum(dim='time')
        if cumsum_ind==True:
            ind=ind.cumsum(dim='time')
        
        if window!=1:
            print('...Applying rolling mean')
            data=data.rolling(time=window, center=True, min_periods=window).mean()
            ind=ind.rolling(time=window, center=True, min_periods=window).mean()

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

        fig.suptitle(title+'(Ex. Convection) \n Ensemble Member '+ str(ens)+'; p<0.05')

        if inverse==True:
            print('...NOTE: Taking the INVERSE of the correlation!')
            res=-1*res

        maxcor=res.max(dim='lag')
        mincor=res.min(dim='lag')
        maxp=p.min(dim='lag')

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
            fig.savefig('../03_output/convection/'+today+'_1Dcormap_conv_'+savefn+'_member'+str(ens)+'.png')

        if save_results==True:
            if cumsum_map==True:
                print('Saving cumsum')
                np.save('../02_data/convection/conv_opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy', opt_lag)
                maxp.to_netcdf('../02_data/convection/conv_maxp_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                maxcor.to_netcdf('../02_data/convection/conv_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')

            else:
                np.save('../02_data/convection/conv_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy', opt_lag)
                maxp.to_netcdf('../02_data/convection/conv_maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                maxcor.to_netcdf('../02_data/convection/conv_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')

    return


def showConvectionCorMap(members='all',
    data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
    ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
    mask='land', years=['1920', '2013'],
    detrend=True, deseasonalize=True, window=24, title='Convection Correlation Map', ymax=-70, save=False, savefn='test',
    draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
    inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False):
'''
    Create a correlation map without convection for multiple ensemble members, after the correlations have been performed!
    
    INPUT:
    members (string or list): members to use ('all' or [0,1] for member 1); [20 21] for ensemble mean.
    data (None or list): List with 2 spatial 2D data fields: index 0 must be x component, index 1 must by y component. If 'None', then it will be read within the function.
    datafn (string): filename of the spatial 2D data field to be read.
    datavar (string): variable to be read.
    datakind (string): subfolder of 02_data to be read.
    ind (None or xr.DataSet): 1D time series to correlate with.
    indfn (string): filename of 1D time series to be read.
    indvar (string): variable in index to be used.
    lags (list): lags to be used for the cross correlation
    mask (string): whether we want to mask out the 'land', 'ice', or 'land and ice' in the figure.
    years (list): years that are at the edges of the time interval we are interested in.
    detrend (boolean): do we want to detrend both the index and the spatial map
    deseasonalize (boolean): do we want to deseasonalize the data.
    window (int): what kind of running mean do we want to apply? [months]
    title (string): titel of the correlation map to use
    ymax (float): northern limit of the map that we are interested in. [deg N]
    save (boolean): doe we want to save the figure?
    savefn (string): extra text for the filename of the figure.
    draw_box (boolean): do we want to draw a box in the domian.
    box_x (list): x limits of the box [deg E]
    box_y (list): y limits of the box [deg N]
    inverse (boolean): do we want to invert the correlation?
    save_results (boolean): do we want to save the correlation coefficients and the lag
    give_back (boolean): whether we want to return the data.
    cumsum_map (boolean): make the spatial variable cumulative
    cumsum_ind (boolean): make the 1D time series cumulative
    
    OUTPUT:
    '''
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from scipy import stats
    
    print('Start preparing for Convection Correlation Map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 22})
    
    if members=='all':
        members=[0,20]
    
    for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
        if deseasonalize==True:
            if cumsum_map==True:
                try:
                    temp_opt_lag=np.load('../02_data/convection/conv_2D_opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_2D_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/convection/conv_opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            else:
                try:
                    temp_opt_lag=np.load('../02_data/convection/conv_2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/convection/conv_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        else:
            try:
                temp_opt_lag=np.load('../02_data/convection/conv_2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxcor=xr.open_dataarray('../02_data/convection/conv_2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            except:
                try:
                    temp_opt_lag=np.load('../02_data/convection/conv_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/convection/conv_2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    
        if i==0:
            opt_lag=temp_opt_lag
            maxcor=temp_maxcor
        else:
            opt_lag=np.dstack((opt_lag, temp_opt_lag))
            maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
    if i==0:
        maxcor=maxcor.expand_dims('ens')

    if members[-1]!=21:
        maxp=stats.ttest_1samp(maxcor, popmean=0)[1]
        opt_lag_p_pos=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='greater')[1]
        opt_lag_p_neg=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='less')[1]
    else:
        maxp=np.zeros(np.shape(opt_lag))
        opt_lag_p_neg=np.zeros(np.shape(opt_lag))
        opt_lag_p_pos=np.zeros(np.shape(opt_lag))
        
    stdcor=maxcor.std(dim='ens')
    maxcor=maxcor.mean(dim='ens')
    if members[-1]!=21:
        std_lag=np.nanstd(opt_lag, axis=2)
        opt_lag=np.nanmean(opt_lag, axis=2)
    else:
        std_lag=np.zeros(np.shape(opt_lag))
    
    maxcor=maxcor.where((maxp<0.05) | maxcor.isnull(), other=0)
    opt_lag[std_lag>8]=0
    
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
        fig.savefig('../03_output/convection/conv_'+today+'_2Dcormap_'+savefn+'_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
    fig=plt.figure(figsize=(20,10))
    plt.suptitle('Extra statistics over ensemble members for : \n'+title)
    plt.subplot(2,3,1)
    img=plt.pcolor(maxcor.XC, maxcor.YC, maxp, vmax=0.05, cmap='seismic')
    plt.colorbar(img, label='P value')
    plt.title('P value, corr!=0'))
    
    plt.subplot(2,3,2)
    img=plt.pcolor(maxcor.XC, maxcor.YC, opt_lag_p_pos, vmin=0, vmax=0.05, cmap='seismic')
    plt.colorbar(img, label='P value')
    plt.title('P value, lag>0')
    
    plt.subplot(2,3,3)
    img=plt.pcolor(maxcor.XC, maxcor.YC, opt_lag_p_neg, vmin=0, vmax=0.05, cmap='seismic')
    plt.colorbar(img, label='P value')
    plt.title('P value, lag<0')
    
    plt.subplot(2,2,3)
    img=plt.pcolor(stdcor.XC, stdcor.YC, stdcor, cmap='seismic', vmax=0.3)
    plt.colorbar(img, label='Std of cor')
    plt.title('Standard deviation of correlation coefficient')#  \n Yellow line indicates 0.1')
    
    plt.subplot(2,2,4)
    img=plt.pcolor(stdcor.XC, stdcor.YC, std_lag, cmap='seismic', vmax=10)
    plt.colorbar(img, label='Std of lag')
    plt.title('Standard deviation of lag')
    
    fig.subplots_adjust(hspace=0.5, wspace=0.3, top=0.8)
    
    if save==True:
        print('...Saving second figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/convection/conv_'+today+'_BONUS_2Dcormap_'+savefn+'_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
    return



def showConvectionCorMapAxis(ax, members='all',
    data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
    ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
    mask='land', years=['1920', '2013'],
    detrend=True, deseasonalize=True, window=24, title='Convection Correlation Map', ymax=-70, save=False, savefn='test',
    draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
    inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False):
'''
    Same as showConvectionCorMap, but only the correlation map and only for one figure axis.
    
    Create a correlation map for multiple ensemble members, after the correlations have been performed!
    
    INPUT:
    ax (Matplotlib figure axis): axis on which the correlation map should be plotted. 
    members (string or list): members to use ('all' or [0,1] for member 1); [20 21] for ensemble mean.
    data (None or list): List with 2 spatial 2D data fields: index 0 must be x component, index 1 must by y component. If 'None', then it will be read within the function.
    datafn (string): filename of the spatial 2D data field to be read.
    datavar (string): variable to be read.
    datakind (string): subfolder of 02_data to be read.
    ind (None or xr.DataSet): 1D time series to correlate with.
    indfn (string): filename of 1D time series to be read.
    indvar (string): variable in index to be used.
    lags (list): lags to be used for the cross correlation
    mask (string): whether we want to mask out the 'land', 'ice', or 'land and ice' in the figure.
    years (list): years that are at the edges of the time interval we are interested in.
    detrend (boolean): do we want to detrend both the index and the spatial map
    deseasonalize (boolean): do we want to deseasonalize the data.
    window (int): what kind of running mean do we want to apply? [months]
    title (string): titel of the correlation map to use
    ymax (float): northern limit of the map that we are interested in. [deg N]
    save (boolean): doe we want to save the figure?
    savefn (string): extra text for the filename of the figure.
    draw_box (boolean): do we want to draw a box in the domian.
    box_x (list): x limits of the box [deg E]
    box_y (list): y limits of the box [deg N]
    inverse (boolean): do we want to invert the correlation?
    save_results (boolean): do we want to save the correlation coefficients and the lag
    give_back (boolean): whether we want to return the data.
    cumsum_map (boolean): make the spatial variable cumulative
    cumsum_ind (boolean): make the 1D time series cumulative
    
    OUTPUT:
    ax (Matplotlib figure axis): figure axis with plot
    (maxcor): correlation coefficients
    '''
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from scipy import stats
    
    print('Start preparing for Convection Correlation Map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 22})
    
    if members=='all':
        members=[0,20]
    
    for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
        if deseasonalize==True:
            if cumsum_map==True:
                try:
                    temp_opt_lag=np.load('../02_data/convection/conv_2D_opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_2D_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/convection/conv_opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            else:
                try:
                    temp_opt_lag=np.load('../02_data/convection/conv_2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/convection/conv_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        else:
            try:
                temp_opt_lag=np.load('../02_data/convection/conv_2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxcor=xr.open_dataarray('../02_data/convection/conv_2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            except:
                try:
                    temp_opt_lag=np.load('../02_data/convection/conv_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/convection/conv_2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    
        if i==0:
            opt_lag=temp_opt_lag
            maxcor=temp_maxcor
        else:
            opt_lag=np.dstack((opt_lag, temp_opt_lag))
            maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
    if i==0:
        maxcor=maxcor.expand_dims('ens')

    if members[-1]!=21:
        maxp=stats.ttest_1samp(maxcor, popmean=0)[1]
        opt_lag_p_pos=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='greater')[1]
        opt_lag_p_neg=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='less')[1]
    else:
        maxp=np.zeros(np.shape(opt_lag))
        opt_lag_p_neg=np.zeros(np.shape(opt_lag))
        opt_lag_p_pos=np.zeros(np.shape(opt_lag))
        
    stdcor=maxcor.std(dim='ens')
    maxcor=maxcor.mean(dim='ens')
    if members[-1]!=21:
        std_lag=np.nanstd(opt_lag, axis=2)
        opt_lag=np.nanmean(opt_lag, axis=2)
    else:
        std_lag=np.zeros(np.shape(opt_lag))

    maxcor=maxcor.where((maxp<0.05) | maxcor.isnull(), other=0)
    opt_lag[std_lag>8]=0
    
    #Reading Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    print('...Plotting')
    #PLOT CORRELATION COEFFICIENTS
    img=ax.pcolor(maxcor.XC, maxcor.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
    plt.colorbar(img, ax=ax, label='Correlation Coefficient')
    CS=ax.contour(maxcor.XC, maxcor.YC, maxcor, [-0.5, 0.5], colors='yellow')
    ax.clabel(CS, CS.levels, inline=True,fontsize=10)

    bath=ax.contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax.set_ylim([-75.5, ymax])
    if 'ice' in mask:
        ax.pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
    ax.set_facecolor('grey')
    ax.set_xlabel('Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')
    ax.set_title(title)

    if draw_box==True:
        ax.plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')

    return ax
        
def showConvectionLagMapAxis(ax, members='all',
    data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
    ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
    mask='land', years=['1920', '2013'],
    detrend=True, deseasonalize=True, window=24, title='Convection Correlation Map', ymax=-70, save=False, savefn='test',
    draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
    inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False):
'''
    Same as showConvectionCorMap, but only the lag map and only for one figure axis.
    
    Create a lag map for multiple ensemble members, after the correlations have been performed!
    
    INPUT:
    ax (Matplotlib figure axis): axis on which the correlation map should be plotted. 
    members (string or list): members to use ('all' or [0,1] for member 1); [20 21] for ensemble mean.
    data (None or list): List with 2 spatial 2D data fields: index 0 must be x component, index 1 must by y component. If 'None', then it will be read within the function.
    datafn (string): filename of the spatial 2D data field to be read.
    datavar (string): variable to be read.
    datakind (string): subfolder of 02_data to be read.
    ind (None or xr.DataSet): 1D time series to correlate with.
    indfn (string): filename of 1D time series to be read.
    indvar (string): variable in index to be used.
    lags (list): lags to be used for the cross correlation
    mask (string): whether we want to mask out the 'land', 'ice', or 'land and ice' in the figure.
    years (list): years that are at the edges of the time interval we are interested in.
    detrend (boolean): do we want to detrend both the index and the spatial map
    deseasonalize (boolean): do we want to deseasonalize the data.
    window (int): what kind of running mean do we want to apply? [months]
    title (string): titel of the correlation map to use
    ymax (float): northern limit of the map that we are interested in. [deg N]
    save (boolean): doe we want to save the figure?
    savefn (string): extra text for the filename of the figure.
    draw_box (boolean): do we want to draw a box in the domian.
    box_x (list): x limits of the box [deg E]
    box_y (list): y limits of the box [deg N]
    inverse (boolean): do we want to invert the correlation?
    save_results (boolean): do we want to save the correlation coefficients and the lag
    give_back (boolean): whether we want to return the data.
    cumsum_map (boolean): make the spatial variable cumulative
    cumsum_ind (boolean): make the 1D time series cumulative
    
    OUTPUT:
    ax (Matplotlib figure axis): figure axis with plot
    (maxcor): correlation coefficients
    '''
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from scipy import stats
    
    print('Start preparing for Convection Correlation Map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 22})
    
    if members=='all':
        members=[0,20]
    
    for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
        if deseasonalize==True:
            if cumsum_map==True:
                try:
                    temp_opt_lag=np.load('../02_data/convection/conv_2D_opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_2D_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/convection/conv_opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            else:
                try:
                    temp_opt_lag=np.load('../02_data/convection/conv_2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/convection/conv_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        else:
            try:
                temp_opt_lag=np.load('../02_data/convection/conv_2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxcor=xr.open_dataarray('../02_data/convection/conv_2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            except:
                try:
                    temp_opt_lag=np.load('../02_data/convection/conv_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/convection/conv_2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/convection/conv_2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    
        if i==0:
            opt_lag=temp_opt_lag
            maxcor=temp_maxcor
        else:
            opt_lag=np.dstack((opt_lag, temp_opt_lag))
            maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
    if i==0:
        maxcor=maxcor.expand_dims('ens')

    
    if members[-1]!=21:
        maxp=stats.ttest_1samp(maxcor, popmean=0)[1]
        opt_lag_p_pos=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='greater')[1]
        opt_lag_p_neg=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='less')[1]
    else:
        maxp=np.zeros(np.shape(opt_lag))
        opt_lag_p_neg=np.zeros(np.shape(opt_lag))
        opt_lag_p_pos=np.zeros(np.shape(opt_lag))
        
    stdcor=maxcor.std(dim='ens')
    maxcor=maxcor.mean(dim='ens')
    if members[-1]!=21:
        std_lag=np.nanstd(opt_lag, axis=2)
        opt_lag=np.nanmean(opt_lag, axis=2)
    else:
        std_lag=np.zeros(np.shape(opt_lag))
    
    maxcor=maxcor.where((maxp<0.05) | maxcor.isnull(), other=0)
    opt_lag[std_lag>8]=0
    
    #Reading Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    print('...Plotting')

    if lags[-1]<0:
        divnorm=MidpointNormalize(vmin=lags[0], vmax=0.1, midpoint=0)  
    else:
        divnorm=MidpointNormalize(vmin=lags[0], vmax=lags[-1], midpoint=0)
    
    #PLOT OPTIMAL LAG
    img=ax.pcolor(maxcor.XC, maxcor.YC, opt_lag, cmap='seismic', norm=divnorm)
    plt.colorbar(img, ax=ax, label='Lag [months]', pad=0.10, ticks=np.arange(lags[0], lags[-1]+1, 6))
    
    bath=ax.contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax.set_ylim([-75.5, ymax])

    if 'ice' in mask:
        ax.pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
    ax.set_facecolor('grey')
    ax.set_xlabel('Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')

    if draw_box==True:
        ax.plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    ax.set_title(title)
    return ax


def createConvectionMaps(save=True):
    '''
    Show showing correlation maps for different variables in the mechanism.
    '''
    
    window=60
    fig=plt.figure(figsize=(20,15))
    ax=plt.subplot(2,2,1)
    ax=showConvectionCorMapAxis(ax=ax, members='all',
                  data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Zonal Wind', ymax=None, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=True, cumsum_map=False)

    ax=plt.subplot(2,2,2)
    showConvectionCorMapAxis(ax=ax, members='all',
                      data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                      detrend=True, window=window, title='Meridional Wind', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False, deseasonalize=True, cumsum_map=False)
    
    ax=plt.subplot(2,2,3)
    ax=showConvectionCorMapAxis(ax=ax, members='all',
                  data=None, datafn='oceFWflx', datavar='oceFWflx', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  detrend=True, window=window, title='Surface Freshwater Flux', ymax=None, save=True, savefn='oceFWflx_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=True, cumsum_map=False)
    
    ax=plt.subplot(2,2,4)
    ax=showConvectionCorMapAxis(ax=ax, members='all',
                  data=None, datafn='SALT_averaged_-200to0', datavar='SALT', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', 
                  lags=range(-24, 24),
                  detrend=True, window=window, 
                  title='Average Salinity (-200 to 0m)', ymax=-70, save=True, 
                  savefn='SALT-200to0_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=True, cumsum_map=False)


    fig.suptitle('Correlations without Convection with Dotson to Cosgrove Mass Loss')
    #plt.rcParams.update({'font.size': 18})
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/convection/'+today+'_convection_corrcoeff.png')
        
    #----LAGS-----
    fig=plt.figure(figsize=(20,15))
    
    ax=plt.subplot(2,2,1)
    ax=showConvectionLagMapAxis(ax=ax, members='all',
                  data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Zonal Wind', ymax=None, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=True, cumsum_map=False)

    ax=plt.subplot(2,2,2)
    showConvectionLagMapAxis(ax=ax, members='all',
                      data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                      detrend=True, window=window, title='Meridional Wind', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False, deseasonalize=True, cumsum_map=False)
    
    ax=plt.subplot(2,2,3)
    ax=showConvectionLagMapAxis(ax=ax, members='all',
                  data=None, datafn='oceFWflx', datavar='oceFWflx', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  detrend=True, window=window, title='Surface Freshwater Flux', ymax=None, save=True, savefn='oceFWflx_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=True, cumsum_map=False)
    
    ax=plt.subplot(2,2,4)
    ax=showConvectionLagMapAxis(ax=ax, members='all',
                  data=None, datafn='SALT_averaged_-200to0', datavar='SALT', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', 
                  lags=range(-24, 24),
                  detrend=True, window=window, 
                  title='Average Salinity (-200 to 0m)', ymax=-70, save=True, 
                  savefn='SALT-200to0_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=True, cumsum_map=False)

    fig.suptitle('Lags without Convection with Dotson to Cosgrove Mass Loss \n positive = mass loss is leading')
    plt.rcParams.update({'font.size': 18})
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/convection/'+today+'_convection_lags.png')
    return

import matplotlib as mpl
class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return np.ma.masked_array(np.interp(value, x, y))