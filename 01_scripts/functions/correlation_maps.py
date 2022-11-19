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
import string

from .loading_and_processing_data import loadMITgcmData, loadPACEData, detrend_and_average_MITgcmData, loadTimeSeriesData
sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
from mitgcm_python_master.plot_utils.labels import latlon_axes

def correlationMap_1D_v3(members='all',
                      data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      mask='land', years=['1920', '2013'],
                      detrend=True, deseasonalize=True, window=24, title='2D Correlation Map', ymax=-70, save=False,
                      savefn='test',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False):
    '''
    Create a correlation map with one spatial variable with one component (e.g. Temperature, not velocity).
    
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
    
    #Save original data
    if type(ind)!=type(None):
        ind_og=ind.copy()
    else:
        ind_og=None
    if type(data)!=type(None):
        data_og=data.copy()
    else:
        data_og=None
    #data_og=data_og.rename({'ens':'ens1'})
    #ind_og=ind_og.rename({'ens':'ens1'})
    
    #Loop over ensemble members!
    for ens in range(members[0]+1,members[-1]+1):
        #First, we read the data!
        if type(data_og)==type(None):
            print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_....')
            data=loadMITgcmData(filename=datafn, members=[ens], kind=datakind)
        else:
            data=data_og.sel(ens=[ens-1])
        
        #Interpolate everything to the center of the cell.
        print('...Making sure interpolate everything to the center of the cells')
        
        #Rename the axes
        if 'XG' in data.dims:
            data=data.rename({'XG':'XC'})
        if 'YG' in data.dims:
            data=data.rename({'YG':'YC'})

        #Read something with the correct locations.
        data1=loadMITgcmData(filename='EXFuwind', members=[1], kind='maps')
        
        #Interpolate
        data=data.interp(XC=data1.XC, YC=data1.YC)

        #Reading Grid
        gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
        grid = Grid(gp)

        #Mask out land and/or ice
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

        #Read the 1D time series if necessary
        if type(ind_og)==type(None):
            print('...Reading: '+indfn)
            if 'timeseries' in indfn:
                if 'GEO' in datafn:
                    ind=loadTimeSeriesData(filename=indfn, members=[ens], kind='timeseries', PACE=False)
                    _, index = np.unique(ind['time'], return_index=True)
                    ind=ind.isel(time=index)
                elif ('final' in indfn) or ('forcing' in indfn):
                    ind=loadTimeSeriesData(filename=indfn, members=[ens], kind='old')
                else:
                    ind=loadTimeSeriesData(filename=indfn, members=[ens], kind='timeseries')
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

        #Detrend data
        if detrend==True:
            print('...Detrending data')
            data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
            ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)
            
        #Deseasonalize data
        if deseasonalize==True:
            print('...Deseasonalizing the data')
            monthlymean=data.groupby('time.month').mean("time")
            data=data.groupby('time.month')-monthlymean
            
            monthlymean=ind.groupby('time.month').mean("time")
            ind=ind.groupby('time.month')-monthlymean
    
        #Create cumulative data
        if cumsum_map==True:
            data=data.cumsum(dim='time')
        if cumsum_ind==True:
            ind=ind.cumsum(dim='time')
        
        #Apply rolling mean
        if window!=1:
            print('...Applying rolling mean')
            data=data.rolling(time=window, center=True).mean()
            ind=ind.rolling(time=window, center=True).mean()

        print('...Computing correlation coefficients')

        #Create a new axis: ens_time!
        data=data.stack(ens_time=['ens', 'time'])
        
        #Find out optimal lag
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
      
        #Take inverse of correlation
        if inverse==True:
            print('...NOTE: Taking the INVERSE of the correlation!')
            res=-1*res
    
        #Find maximum absolute correlation
        maxcor=res.max(dim='lag')
        mincor=res.min(dim='lag')
        
        #Find minimal p value.
        maxp=p.min(dim='lag')
    
        #Remove not significant cells
        maxcor=maxcor.where((maxp<=0.05) | (np.isnan(maxp)), other=0)
        mincor=mincor.where((maxp<=0.05) | (np.isnan(maxp)), other=0)

        #Cells get the maximum absolute value for the correlation coefficinents
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

        
        #Compute optimal lag
        print('...Finding optimal lag')
        
        #Select which axis is lag
        idx=np.argmax([i=='lag' for i in res.dims])
        
        #Mask out values where not significant.
        masked=np.ma.masked_invalid(res.where((maxp<=0.05) & (maxp>=0)).to_numpy())
        
        #Select the lag where the absolute correlation coefficient is maximum
        opt_lag=np.nanargmax(abs(masked), axis=idx)
        opt_lag=opt_lag.astype(float)
        #Remove cells where any of the values are invalid
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
            fig.savefig('../03_output/correlation/'+today+'_1Dcormap_'+savefn+'_member'+str(ens)+'.png')

        if save_results==True:
            if cumsum_map==True:
                print('Saving cumsum')
                np.save('../02_data/correlation/opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy', opt_lag)
                maxp.to_netcdf('../02_data/correlation/maxp_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                maxcor.to_netcdf('../02_data/correlation/maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')

            else:
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
                      give_back=False, cumsum_map=False, cumsum_ind=False):
    '''
    Create a correlation map with 2 spatial variables (x and y component) (e.g. Velocity, not Temperature).
    Give the correlation of absolute value, not direction.
    NOTE: ALWAYS GIVE X COMPONENT!
    
    
    INPUT:
    members (string or list): members to use ('all' or [0,1] for member 1)
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
    
    OUTPUT:
    '''
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 14})
    
    if members=='all':
        members=[0,20]
    
    #Original data
    data_og=data
    ind_og=ind
    
    for ens in range(members[0]+1,members[-1]+1):
        #First, we read the data! Both x and y compoenet
        if type(data_og)==type(None):
                print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_....')
                dataU=loadMITgcmData(filename=datafn, members=[ens], kind=datakind)
                datafn_v=datafn.replace('U','V')
                datafn_v=datafn_v.replace('u','v')
                datafn_v=datafn_v.replace('x','y')
                print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn_v+'_....')
                dataV=loadMITgcmData(filename=datafn_v, members=[ens], kind=datakind)

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
        data=loadMITgcmData(filename='EXFuwind', members=[1], kind='maps')

        #Interpolate it to the correct location
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

        #Read 1D time series
        if type(ind_og)==type(None):
            print('...Reading: '+indfn)
            if 'timeseries' in indfn:
                if 'GEO' in datafn:
                    ind=loadTimeSeriesData(filename=indfn, members=[ens], kind='timeseries', PACE=False)
                    _, index = np.unique(ind['time'], return_index=True)
                    ind=ind.isel(time=index)
                elif ('final' in indfn) or ('forcing' in indfn):
                    ind=loadTimeSeriesData(filename=indfn, members=[ens], kind='old')
                else:
                    ind=loadTimeSeriesData(filename=indfn, members=[ens], kind='timeseries')
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
        
        #Create cumulative data
        if cumsum_map==True:
            data=data.cumsum(dim='time')
        if cumsum_ind==True:
            ind=ind.cumsum(dim='time')
        
        if window!=1:
            print('...Applying rolling mean')
            data=data.rolling(time=window, center=True).mean()
            ind=ind.rolling(time=window, center=True).mean()

        print('...Computing correlation coefficients')
        
        #Create new axis: ens_time.
        data=data.stack(ens_time=['ens', 'time'])
        
        #Look at different lags.
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
        
        #Give northern edge of plot
        if type(ymax)!=type(None):
            res=res.sel(YC=slice(-75.5, ymax))
            p=p.sel(YC=slice(-75.5, ymax))

        fig.suptitle(title+'\n Ensemble Member '+ str(ens)+'; p<0.05')
        
        #Take inverse of correlation
        if inverse==True:
            print('...NOTE: Taking the INVERSE of the correlation!')
            res=-1*res

        #Show maximum absolute correlation
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
        
        #Find optimal lag corresponding to maximum correlation.
        idx=np.argmax([i=='lag' for i in res.dims])
        masked=np.ma.masked_invalid(res.where((maxp<=0.05) & (maxp>=0)).to_numpy())
        opt_lag=np.nanargmax(abs(masked), axis=idx)
        opt_lag=opt_lag.astype(float)

        #Remove insignificant data points
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
                      give_back=False, cumsum_map=False, cumsum_ind=False):
    '''
    Create a correlation map for multiple ensemble members, after the correlations have been performed!
    NOTE: ALWAYS GIVE X COMPONENT!
    
    
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
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 22})
    
    if members=='all':
        members=[0,20]
    
    for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
        #Read Data: complicated because things were not saved in a proper format....
        if deseasonalize==True:
            if cumsum_map==True:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            
            else:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        
        else:
            try:
                temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            
            except:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
           

        if i==0:
            opt_lag=temp_opt_lag
            maxcor=temp_maxcor
        else:
            opt_lag=np.dstack((opt_lag, temp_opt_lag))
            maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
    if i==0:
        maxcor=maxcor.expand_dims('ens')

    #Perform t-tests if not is ensemble mean.
    if members[-1]!=21: 
        maxp=stats.ttest_1samp(maxcor, popmean=0)[1]
        opt_lag_p_pos=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='greater')[1]
        opt_lag_p_neg=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='less')[1]
    else:
        maxp=np.zeros(np.shape(opt_lag))
        opt_lag_p_neg=np.zeros(np.shape(opt_lag))
        opt_lag_p_pos=np.zeros(np.shape(opt_lag))
        
    #Compute standard deviations and mean values of lag and corr
    stdcor=maxcor.std(dim='ens')
    maxcor=maxcor.mean(dim='ens')
    if members[-1]!=21:
        std_lag=np.nanstd(opt_lag, axis=2)
        opt_lag=np.nanmean(opt_lag, axis=2)
    else:
        std_lag=np.zeros(np.shape(opt_lag))
    
    #Select only data points where p value of correlation is smaller than 0.05
    maxcor=maxcor.where((maxp<0.05) | maxcor.isnull(), other=0)
    
    #Select only lags where lag is smaller than 8.
    opt_lag[std_lag>8]=0
    
    #Reading Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    print('...Plotting')
    fig, ax=plt.subplots(2,1, figsize=(20,14))
    fig.suptitle(title)


    #PLOT CORRELATION COEFFICIENTS
    img=ax[0].pcolor(maxcor.XC-360, maxcor.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
    plt.colorbar(img, ax=ax[0], label='Correlation Coefficient')
    CS=ax[0].contour(maxcor.XC-360, maxcor.YC, maxcor, [-0.8, -0.5, 0.5, 0.8], linewidths=[2], colors=['white', 'yellow', 'yellow', 'white'])
    ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)
    
    bathy=grid.bathy
    bathy[grid.ice_mask | ((grid.lon_2d>-119) & (grid.lat_2d<-72.5))]=np.nan
    bath=ax[0].contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], linewidths=[2], colors=['black'])
    plt.clabel(bath, [-1000])
    
    
    if type(ymax)!=type(None):
        ax[0].set_ylim([-75.5, ymax])
        
    if 'ice' in mask:
        ax[0].pcolor(grid.lon_2d, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool', zorder=10)
    ax[0].set_facecolor('grey')
    ax[0].set_xlabel('Longitude')
    ax[0].set_ylabel('Latitude')
    ax[0].set_title('Correlation Coefficient')
    
    if type(ymax)!=type(None):
        latlon_axes(ax=ax[0], x=maxcor.XC-360, y=maxcor.YC[maxcor.YC<ymax])
    else:
        latlon_axes(ax=ax[0], x=maxcor.XC-360, y=maxcor.YC)
    
    if draw_box==True:
        ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
    
    
    #PLOT OPTIMAL LAG
    if lags[-1]==0:
        divnorm=MidpointNormalize(vmin=lags[0], vmax=lags[-1], midpoint=0.001)
    else:
        divnorm=MidpointNormalize(vmin=lags[0], vmax=lags[-1], midpoint=0)
    
    img=ax[1].pcolor(maxcor.XC-360, maxcor.YC, opt_lag, cmap='seismic', norm=divnorm)
    plt.colorbar(img, ax=ax[1], label='Lag [months]', ticks=np.arange(lags[0], lags[-1]+1, 6))
    bath=ax[1].contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], linewidths=[2], colors=['black'])
    plt.clabel(bath, [-1000])
    
    
    if type(ymax)!=type(None):
        ax[1].set_ylim([-75.5, ymax])

    if 'ice' in mask:
        ax[1].pcolor(grid.lon_2d, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool', zorder=10)
    ax[1].set_facecolor('grey')
    ax[1].set_xlabel('Longitude')
    ax[1].set_ylabel('Latitude')
    
    if type(ymax)!=type(None):
        latlon_axes(ax=ax[1], x=maxcor.XC-360, y=maxcor.YC[maxcor.YC<ymax])
    else:
        latlon_axes(ax=ax[1], x=maxcor.XC-360, y=maxcor.YC)
    
    if draw_box==True:
        ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
    
    
    if indvar=='dotson_to_cosgrove_massloss':
        timeseries_name='mass loss'
    elif indvar=='amundsen_shelf_break_uwind_avg':
        timeseries_name='wind'
    ax[1].set_title('Lag \n (positive = '+timeseries_name+' leading)')
    
    import string
    for n, ax in enumerate(fig.axes[:2]):
        #Add letters to subplot
        ax.text(-0.1, 0.9, string.ascii_uppercase[n]+'.', transform=ax.transAxes, 
                size=25, weight='bold')
    
    fig.subplots_adjust(hspace=0.5, top=0.85)

    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
    fig=plt.figure(figsize=(20,10))
    plt.suptitle('Extra statistics over ensemble members for : \n'+title)
    ax1=plt.subplot(2,3,1)
    img=plt.pcolor(maxcor.XC-360, maxcor.YC, maxp, vmin=0, vmax=0.1, cmap='jet')
    plt.colorbar(img, label='P value', extend='max')
    plt.title('P value, corr!=0')
    plt.xticks(rotation = 25)
    
    ax2=plt.subplot(2,3,2)
    img=plt.pcolor(maxcor.XC-360, maxcor.YC, opt_lag_p_pos, vmin=0, vmax=0.1, cmap='jet')
    plt.colorbar(img, label='P value', extend='max')
    plt.title('P value, lag>0')
    plt.xticks(rotation = 25)
    
    ax3=plt.subplot(2,3,3)
    img=plt.pcolor(maxcor.XC-360, maxcor.YC, opt_lag_p_neg, vmin=0, vmax=0.1, cmap='jet')
    plt.colorbar(img, label='P value', extend='max')
    plt.title('P value, lag<0')
    plt.xticks(rotation = 25)
    
    ax4=plt.subplot(2,2,3)
    img=plt.pcolor(stdcor.XC-360, stdcor.YC, stdcor, cmap='jet', vmin=0, vmax=0.5)
    plt.colorbar(img, label='Std of cor', extend='max')
    plt.title('Standard deviation of correlation coefficient')
    
    ax5=plt.subplot(2,2,4)
    img=plt.pcolor(stdcor.XC-360, stdcor.YC, std_lag, cmap='jet', vmin=0, vmax=12)
    plt.colorbar(img, label='Std of lag', extend='max')
    plt.title('Standard deviation of lag')
    
    
    if 'ice' in mask:
        for ax in fig.axes[::2]:
            ax.pcolor(grid.lon_2d, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)),\
                                                            np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')

    import string
    for n, ax in enumerate(fig.axes[::2]):
        ax.set_facecolor('grey')
        ax.text(-0.25, 0.9, string.ascii_uppercase[n]+'.', transform=ax.transAxes, 
                size=25, weight='bold')
    
    fig.subplots_adjust(hspace=0.5, wspace=0.5, top=0.8)
    
    for ax in fig.axes[::2]:
        latlon_axes(ax=ax, x=maxcor.XC-360, y=maxcor.YC)
    
    if save==True:
        print('...Saving second figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_BONUS_2Dcormap_'+savefn+'_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
    return

def showCorMapAxis(ax, members='all',
                  data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  detrend=True, deseasonalize=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  mask='land', years=['1920', '2013'],
                  inverse=False, save_results=True,
                  give_back=False, cumsum_map=False):
    
    '''
    Same as showCombinedCorrelationMap_2D_v3, but only the correlation map and only for one figure axis.
    
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
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 22})
    
    if members=='all':
        members=[0,20]
    
    for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
        #Read the data
        if deseasonalize==True:
            if cumsum_map==True:
                try:
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                     temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            else:
                try:
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        else:
            try:
                temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            except:
                try:
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
           
        
        if i==0:
            maxcor=temp_maxcor
        else:
            maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
    if i==0:
        maxcor=maxcor.expand_dims('ens')

    #Perform statistical test
    if members[1]-members[0]==1:
        maxp=np.zeros(np.shape(maxcor.mean(dim='ens')))
    else:
        maxp=stats.ttest_1samp(maxcor, popmean=0)[1]
    
    stdcor=maxcor.std(dim='ens')
    maxcor=maxcor.mean(dim='ens')
    
    maxcor=maxcor.where((maxp<0.05) | maxcor.isnull(), other=0)
    
    #Reading Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    print('...Plotting')

    #PLOT CORRELATION COEFFICIENTS
    img=ax.pcolor(maxcor.XC-360, maxcor.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
    plt.colorbar(img, ax=ax, label='Correlation Coefficient')
    CS=ax.contour(maxcor.XC-360, maxcor.YC, maxcor, [-0.8, -0.5, 0.5, 0.8], colors=['white', 'yellow','yellow','white'])
    ax.clabel(CS, CS.levels, inline=True,fontsize=10)

    bathy=grid.bathy
    bathy[grid.ice_mask | ((grid.lon_2d>-119) & (grid.lat_2d<-72.5))]=np.nan
    bath=ax.contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], linewidths=[2], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax.set_ylim([-75.5, ymax])
    if 'ice' in mask:
        ax.pcolor(grid.lon_2d, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
    ax.set_facecolor('grey')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(title)
    #latlon_axes(ax=ax, x=grid.lon_2d, y=grid.lat_2d)
    
    if draw_box==True:
        ax.plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
    if give_back==True:
        return ax, maxcor
    else:
        return ax



def showLagMapAxis(ax, members='all',
                  data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  detrend=True, deseasonalize=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  mask='land', years=['1920', '2013'],
                  inverse=False, save_results=True,
                  give_back=False, cumsum_map=False):
    
    '''
    Same as showCombinedCorrelationMap_2D_v3, but only the lag map and only for one figure axis.
    
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
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 22})
    
    if members=='all':
        members=[0,20]
    
    for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
        if deseasonalize==True:
            if cumsum_map==True:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            else:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        else:
            try:
                temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            except:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    
        if i==0:
            opt_lag=temp_opt_lag
            maxcor=temp_maxcor
        else:
            opt_lag=np.dstack((opt_lag, temp_opt_lag))
            maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
    if i==0:
        maxcor=maxcor.expand_dims('ens')

    #Perform statistical tests
    maxp=stats.ttest_1samp(maxcor, popmean=0)[1]
    #opt_lag_p_pos=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='greater')[1]
    #opt_lag_p_neg=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='less')[1]
    
    stdcor=maxcor.std(dim='ens')
    maxcor=maxcor.mean(dim='ens')
    
    #Find values if not ensemble mean.
    if (members[-1]-members[0])==1:
        std_lag=np.zeros(np.shape(opt_lag))
    else:
        std_lag=np.nanstd(opt_lag, axis=2)
        opt_lag=np.nanmean(opt_lag, axis=2)
        
    
    maxcor=maxcor.where((maxp<0.05) | maxcor.isnull(), other=0)
    opt_lag[std_lag>8]=0
    
    #Reading Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    print('...Plotting')

    #PLOT LAGS
    if lags[-1]<0:
        divnorm=MidpointNormalize(vmin=lags[0], vmax=0.1, midpoint=0)  
    else:
        divnorm=MidpointNormalize(vmin=lags[0], vmax=lags[-1], midpoint=0)    
    img=ax.pcolor(maxcor.XC-360, maxcor.YC, opt_lag, cmap='seismic', norm=divnorm)
    plt.colorbar(img, ax=ax, label='Lag [months]', ticks=np.arange(lags[0], lags[-1]+1, 6))
 
    bathy=grid.bathy
    bathy[grid.ice_mask | ((grid.lon_2d>-119) & (grid.lat_2d<-72.5))]=np.nan
    bath=ax.contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], linewidths=[2], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax.set_ylim([-75.5, ymax])

    if 'ice' in mask:
        ax.pcolor(grid.lon_2d, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
    ax.set_facecolor('grey')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(title)

    if draw_box==True:
        ax.plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    #latlon_axes(ax=ax, x=grid.lon_2d, y=grid.lat_2d)
    return ax






def createWindCorrelationMaps(save=True, window=60):
    '''
    Create correlation maps for wind components; also create lag maps.
    INPUT:
    save (boolean): save figures
    window (int): use correlation coefficients with rolling mean [months]
    '''
    #Import Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)

    fig=plt.figure(figsize=(20,10))

    ax=plt.subplot(2,2,1)
    showCorMapAxis(ax,members='all',
                  data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Zonal Wind', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=False )
    ax=plt.subplot(2,2,2)
    showCorMapAxis(ax,members='all',
                  data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Meridional Wind', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=False )

    ax=plt.subplot(2,2,3)
    showCorMapAxis(ax,members='all',
                  data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Cumulative Zonal Wind', ymax=None, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=True)

    ax=plt.subplot(2,2,4)
    showCorMapAxis(ax,members='all',
                  data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Cumulative Meridional Wind', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=True)
    fig.suptitle('Correlations between the Wind and Dotson to Cosgrove Mass Loss')
    plt.rcParams.update({'font.size': 18})
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
    for n, ax in enumerate(fig.axes[::2]):
        ax.set_xticks([-130, -115, -100, -85])
        ax.set_yticks([-75, -70, -65])
        latlon_axes(ax=ax, x=grid.lon_2d, y=grid.lat_2d)
        ax.text(-0.25, 0.9, string.ascii_uppercase[n]+'.', transform=ax.transAxes, 
            size=25, weight='bold')
        
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_all_wind_corrcoeff.png')
        
    #----LAGS-----
    fig=plt.figure(figsize=(20,10))

    ax=plt.subplot(2,2,1)
    showLagMapAxis(ax,members='all',
                  data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Zonal Wind', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=False )
    ax=plt.subplot(2,2,2)
    showLagMapAxis(ax,members='all',
                  data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Meridional Wind', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=False )

    ax=plt.subplot(2,2,3)
    showLagMapAxis(ax,members='all',
                  data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Cumulative Zonal Wind', ymax=None, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=True)

    ax=plt.subplot(2,2,4)
    showLagMapAxis(ax,members='all',
                  data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Cumulative Meridional Wind', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=True)
    fig.suptitle('Lags between the Wind and Dotson to Cosgrove Mass Loss \n negative = wind is leading')
    plt.rcParams.update({'font.size': 18})
    fig.subplots_adjust(hspace=0.5, wspace=0.5, top=0.8)
    
    for n, ax in enumerate(fig.axes[::2]):
        ax.set_xticks([-130, -115, -100, -85])
        ax.set_yticks([-75, -70, -65])
        latlon_axes(ax=ax, x=grid.lon_2d, y=grid.lat_2d)
        ax.text(-0.25, 0.9, string.ascii_uppercase[n]+'.', transform=ax.transAxes, 
            size=25, weight='bold')
        
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_all_wind_lags.png')
        
    return


def createFlowCorrelationMaps(save=True, window=60):
    '''
    Create correlation maps for baroclinic and barotropic flow; also create lag maps.
    INPUT:
    save (boolean): save figures
    window (int): use correlation coefficients with rolling mean [months]
    '''
    #Reading Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    fig=plt.figure(figsize=(20,15))
    ax=plt.subplot(3,1,1)
    ax, corClin=showCorMapAxis(ax=ax, members='all',
                   data=None, datafn='baroclinic_UVEL', datavar='UVEL', datakind='old',
                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                   detrend=True, window=window, title='Baroclinic Flow (lowest 100m-depth averaged)', ymax=-70, save=True, 
                   savefn='baroclinic_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                   inverse=False,
                   give_back=True, deseasonalize=False)
    ax=plt.subplot(3,1,2)
    ax, corTrop=showCorMapAxis(ax,members='all',
                  data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  mask='land', years=['1920', '2013'],
                  detrend=True, deseasonalize=False, window=window, title='Barotropic Flow', ymax=-70, 
                save=True, savefn='barotropic_VEL_dotson_to_cosgrove_massloss_ensmean_'+str(window)+'month_dt',
                draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                inverse=False, save_results=True, give_back=True, cumsum_map=False)
    
    #Show difference
    ax=plt.subplot(3,1,3)
    difCor=abs(corClin)-abs(corTrop)
    img=ax.pcolor(corTrop.XC-360, corTrop.YC, difCor, cmap='seismic', vmin=-1, vmax=1)
    plt.colorbar(img, ax=ax, label='Correlation Coefficient')
    bathy=grid.bathy
    bathy[grid.ice_mask | ((grid.lon_2d>-119) & (grid.lat_2d<-72.5))]=np.nan
    bath=ax.contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], linewidths=[2], colors=['black'])
    plt.clabel(bath, [-1000])
    ymax=-70
    if type(ymax)!=type(None):
        ax.set_ylim([-75.5, ymax])
    ax.set_facecolor('grey')
    ax.set_xlabel('Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')
    ax.set_title('Difference Absolute Correlation Coefficients (Baroclinic-Barotropic)')

    fig.suptitle('Correlations between the Flow Components and Dotson to Cosgrove Mass Loss')
    #plt.rcParams.update({'font.size': 18})
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
    for n, ax in enumerate(fig.axes[::2]):
        ax.set_xticks([-130, -115, -100, -85])
        ax.set_yticks([-75, -73, -71])
        latlon_axes(ax=ax, x=grid.lon_1d, y=grid.lat_1d[grid.lat_1d<-70])
        ax.text(-0.1, 0.9, string.ascii_uppercase[n]+'.', transform=ax.transAxes, 
            size=25, weight='bold')
    
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_flow_corrcoeff.png')
        
    #----LAGS-----
    fig=plt.figure(figsize=(20,15))
    
    ax=plt.subplot(2,1,1)
    showLagMapAxis(ax=ax, members='all',
                   data=None, datafn='baroclinic_UVEL', datavar='UVEL', datakind='old',
                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                   detrend=True, window=window, title='Baroclinic Flow (lowest 100m-depth averaged)', ymax=-70, save=True, 
                   savefn='baroclinic_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                   inverse=False,
                   give_back=False, deseasonalize=False)
    ax=plt.subplot(2,1,2)
    showLagMapAxis(ax,members='all',
                  data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  mask='land', years=['1920', '2013'],
                  detrend=True, deseasonalize=False, window=window, title='Barotropic Flow', ymax=-70, 
                save=True, savefn='barotropic_VEL_dotson_to_cosgrove_massloss_ensmean_'+str(window)+'month_dt',
                draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                inverse=False, save_results=True, give_back=False, cumsum_map=False)

    fig.suptitle('Lags between the Flow Components and Dotson to Cosgrove Mass Loss \n negative = wind is leading')
    plt.rcParams.update({'font.size': 18})
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
    for n, ax in enumerate(fig.axes[::2]):
        ax.set_xticks([-130, -115, -100, -85])
        ax.set_yticks([-75, -73, -71])
        latlon_axes(ax=ax, x=grid.lon_1d, y=grid.lat_1d[grid.lat_1d<-70])
        ax.text(-0.1, 0.9, string.ascii_uppercase[n]+'.', transform=ax.transAxes, 
            size=25, weight='bold')
        
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_flow_lags.png')
        
    return




def correlationMapPACE(members='all',
                      data=None, datavar='PSL',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                      years=['1920', '2013'],
                      detrend=True, deseasonalize=True, window=24, title='2D Correlation Map', ymax=-70, 
                      save=False, savefn='test',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False):
    '''
    Create a correlation map with one spatial variable with one component (e.g. PSL, SST) from PACE Ensemble.
    
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
    
    data_og=data
    ind_og=ind
    for ens in range(members[0]+1,members[-1]+1):
        #First, we read the data!
        if type(data_og)==type(None):
                print('...Reading: '+datavar+'...')
                data=loadPACEData(datavar, members=[ens]) 
        if datavar=='SST':
            data=data.isel(z_t=0)
        
        if type(ind_og)==type(None):
            print('...Reading: '+indfn)
            if 'timeseries' in indfn:
                ind=loadTimeSeriesData(filename=indfn, members=[ens], kind='old')
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
        
        #This step was necessary because for some reason 1920-01 was not contained in the PACE output (?)
        if len(ind.time)!=len(data.time):
            ind=ind[:,1:]
        
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
        fig=plt.figure(figsize=(20,12))


        fig.suptitle(title+'\n Ensemble Member '+ str(ens)+'; ('+str(window)+' month running mean); p<0.05')

        if inverse==True:
            print('...NOTE: Taking the INVERSE of the correlation!')
            res=-1*res

        maxcor=res.max(dim='lag')
        mincor=res.min(dim='lag')
        maxp=p.min(dim='lag')

        maxcor=maxcor.where((maxp<=0.05) | (np.isnan(maxp)), other=0)
        mincor=mincor.where((maxp<=0.05) | (np.isnan(maxp)), other=0)        
        
        maxcor.to_numpy()[maxcor<abs(mincor)]=mincor.to_numpy()[maxcor<abs(mincor)]

        import cartopy.crs as ccrs
        ax = plt.subplot(211, projection=ccrs.PlateCarree(central_longitude=180.0))
        if 'PSL' in datavar:
            img=ax.pcolormesh(res.lon, res.lat, 
                          maxcor, 
                          cmap='seismic', transform=ccrs.PlateCarree(), vmin=-1, vmax=1)
            ax.contour(res.lon, res.lat, 
                          maxcor, levels=[-0.5, 0.5], c=['yellow', 'yellow'],
                          transform=ccrs.PlateCarree())
        elif ('SST' in datavar): 
            img=ax.scatter(res.TLONG.to_numpy().flatten(), res.TLAT.to_numpy().flatten(), s=0.3,
                          c=maxcor.to_numpy().flatten(), 
                          cmap='seismic', transform=ccrs.PlateCarree(), vmin=-1, vmax=1)
        elif ('aice' in datavar):
                img=ax.scatter(res.TLON.to_numpy().flatten(), res.TLAT.to_numpy().flatten(), s=0.3,
                          c=maxcor.to_numpy().flatten(), 
                          cmap='seismic', transform=ccrs.PlateCarree(), vmin=-1, vmax=1)

        ax.set_facecolor('grey')
        plt.colorbar(img, label='Correlation Coefficient', pad=0.10, ax=ax)
        ax.coastlines()
        ax.set_title('Optimal correlation coefficient')
        ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree())
        

        print('...Finding optimal lag')
        idx=np.argmax([i=='lag' for i in res.dims])
        masked=np.ma.masked_invalid(res.where((maxp<=0.05) & (maxp>=0)).to_numpy())
        opt_lag=np.nanargmax(abs(masked), axis=idx)
        opt_lag=opt_lag.astype(float)
        opt_lag[masked.mask.max(axis=idx)]=np.nan
        opt_lag=opt_lag+lags[0]

        ax = plt.subplot(212, projection=ccrs.PlateCarree(central_longitude=180.0))
        if 'PSL' in datavar:
            img=ax.pcolormesh(res.lon, res.lat, 
                          opt_lag, 
                          cmap='seismic', transform=ccrs.PlateCarree(), vmin=lags[0], vmax=lags[-1])
        elif ('SST' in datavar):
            img=ax.scatter(res.TLONG.to_numpy().flatten(), res.TLAT.to_numpy().flatten(), s=0.3,
                          c=opt_lag.flatten(), 
                          cmap='seismic', transform=ccrs.PlateCarree(), vmin=lags[0], vmax=lags[-1])
        elif ('aice' in datavar):
            img=ax.scatter(res.TLON.to_numpy().flatten(), res.TLAT.to_numpy().flatten(), s=0.3,
                          c=opt_lag.flatten(), 
                          cmap='seismic', transform=ccrs.PlateCarree(), vmin=lags[0], vmax=lags[-1])
        ax.set_facecolor('grey')
        plt.colorbar(img, label='Optimal Lag [months]', pad=0.10, ax=ax)
        ax.coastlines()
        ax.set_title('Optimal lag \n (positive=time series leading)')
        ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree())
    
        fig.subplots_adjust(hspace=0.4)

        if save==True:
            print('...Saving figure')
            from datetime import date
            today = date.today()
            today=today.strftime("%Y%m%d")
            fig.savefig('../03_output/correlation/'+today+'_PACEcormap_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_member'+str(ens)+'.png')

        if save_results==True:
            if cumsum_map==True:
                print('Saving cumsum')
                np.save('../02_data/correlation/opt_lag_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy', opt_lag)
                maxp.to_netcdf('../02_data/correlation/maxp_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                maxcor.to_netcdf('../02_data/correlation/maxcor_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')

            else:
                np.save('../02_data/correlation/opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy', opt_lag)
                maxp.to_netcdf('../02_data/correlation/maxp_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                maxcor.to_netcdf('../02_data/correlation/maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')

    return
  

def showMapPACE(members='all',
              data=None, datavar='PSL',
              ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
              years=['1920', '2013'],
              detrend=True, deseasonalize=True, window=24, title='2D Correlation Map', ymax=-70, 
              save=False, savefn='test',
              draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
              inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False):
    '''
    Show a correlation map with one spatial variable with one component (e.g. PSL, SST) from PACE Ensemble.
    The correlation already has to be performed.
    
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
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from scipy import stats
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 22})
    
    if members=='all':
        members=[0,20]
    
    #Read data
    for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
        if deseasonalize==True:
            if cumsum_map==True:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            else:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        else:
            try:
                temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            except:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    
        if i==0:
            opt_lag=temp_opt_lag
            maxcor=temp_maxcor
        else:
            opt_lag=np.dstack((opt_lag, temp_opt_lag))
            #maxp=xr.concat([maxp, temp_maxp], dim='ens')
            maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
    if i==0:
        maxcor=maxcor.expand_dims('ens')

    #Perform statistical tests
    maxp=stats.ttest_1samp(maxcor, popmean=0)[1]
    opt_lag_p_pos=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='greater')[1]
    opt_lag_p_neg=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='less')[1]
    
    #Compute the mean and standard deviations
    stdcor=maxcor.std(dim='ens')
    maxcor=maxcor.mean(dim='ens')
    if members[-1]!=21:
        std_lag=np.nanstd(opt_lag, axis=2)
        opt_lag=np.nanmean(opt_lag, axis=2)
    else:
        std_lag=np.zeros(np.shape(opt_lag))
    
    #Select only signficant data
    maxcor=maxcor.where((maxp<0.05) | maxcor.isnull(), other=0)
    opt_lag[(opt_lag_p_pos>0.05)&(opt_lag_p_neg>0.05)]=0
        
    #Reading Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    #Plot Figure
    print('...Plotting')
    fig=plt.figure(figsize=(20,14))
    
    if members[-1]==21:
        fig.suptitle(title+'\n Ensemble Mean; '+str(window)+' month rolling mean; p<0.05')
    else:
        fig.suptitle(title+'\n Mean over members '+ str(members[0]+1)+' to '+ str(members[-1])+'; '+str(window)+' month rolling mean; p<0.05')


    #PLOT CORRELATION COEFFICIENTS
    import cartopy.crs as ccrs
    ax = plt.subplot(211, projection=ccrs.PlateCarree(central_longitude=180.0))
    
    
    #Interpolate the correlation coefficients when the variable is SST
    if ('SST' in datavar) | ('aice' in datavar):
        def regridSST(z, x, y):
            from scipy.interpolate import griddata
            PSL=loadPACEData(var='PSL', members=[0,1])
            xgrid=PSL.lon
            ygrid=PSL.lat
            xi,yi = np.meshgrid(xgrid,ygrid)
            
            zfinal=griddata((x,y),z,(xi,yi),method='linear')
            var_new=xr.DataArray(data=zfinal,
                                dims=['lat', 'lon'],
                                coords=dict(
                                lon=xgrid,
                                lat=ygrid)                        
                                    )
            return var_new
        if 'SST' in datavar:
            opt_lag=regridSST(z=opt_lag.flatten(),
                             x=maxcor.ULONG.to_numpy().flatten(),
                             y=maxcor.ULAT.to_numpy().flatten())

            maxcor=regridSST(z=maxcor.to_numpy().flatten(),
                             x=maxcor.ULONG.to_numpy().flatten(),
                             y=maxcor.ULAT.to_numpy().flatten())
        else:
            opt_lag=regridSST(z=opt_lag.flatten(),
                             x=maxcor.ULON.to_numpy().flatten(),
                             y=maxcor.ULAT.to_numpy().flatten())

            maxcor=regridSST(z=maxcor.to_numpy().flatten(),
                             x=maxcor.ULON.to_numpy().flatten(),
                             y=maxcor.ULAT.to_numpy().flatten())
        
    #Show the correlation coefficients.
    img=ax.pcolormesh(maxcor.lon, maxcor.lat, 
                      maxcor, 
                      cmap='seismic', transform=ccrs.PlateCarree(), vmin=-1, vmax=1)
    CS=ax.contour(maxcor.lon, maxcor.lat, 
                      maxcor, levels=[-0.5, 0.5], colors='yellow',
                      transform=ccrs.PlateCarree())
    ax.clabel(CS, CS.levels, inline=True,fontsize=10)
        

    ax.set_facecolor('grey')
    plt.colorbar(img, label='Correlation Coefficient', pad=0.10, ax=ax)
    ax.coastlines()
    ax.set_title('Optimal correlation coefficient')
    ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree())
  

    #PLOT OPTIMAL LAG
    ax = plt.subplot(212, projection=ccrs.PlateCarree(central_longitude=180.0))
    #if 'PSL' in datavar:
    img=ax.pcolormesh(maxcor.lon, maxcor.lat, 
                      opt_lag, 
                      cmap='seismic', transform=ccrs.PlateCarree(), vmin=lags[0], vmax=lags[-1])
    ax.set_facecolor('grey')
    plt.colorbar(img, label='Optimal Lag [months]', pad=0.10, ax=ax)
    ax.coastlines()
    ax.set_title('Optimal lag \n (positive=time series leading)')
    ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree())
    
    fig.subplots_adjust(hspace=0.5, top=0.85)
    
    #Save the data
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+datavar+'_'+indvar+'_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
    #Show statistics
    if 'PSL' in datavar:
        fig=plt.figure(figsize=(20,10))
        plt.suptitle('Extra statistics over ensemble members for : \n'+title)
        plt.subplot(2,3,1)
        img=plt.pcolor(maxcor.lon, maxcor.lat, maxp, vmax=0.05, cmap='seismic')
        plt.colorbar(img, label='P value')
        plt.title('P value, corr!=0')

        plt.subplot(2,3,2)
        img=plt.pcolor(maxcor.lon, maxcor.lat, opt_lag_p_pos, vmin=0, vmax=0.05, cmap='seismic')
        plt.colorbar(img, label='P value')
        plt.title('P value, lag>0')

        plt.subplot(2,3,3)
        img=plt.pcolor(maxcor.lon, maxcor.lat, opt_lag_p_neg, vmin=0, vmax=0.05, cmap='seismic')
        plt.colorbar(img, label='P value')
        plt.title('P value, lag<0')

        plt.subplot(2,2,3)
        img=plt.pcolor(stdcor.lon, stdcor.lat, stdcor, cmap='seismic', vmax=0.3)
        plt.colorbar(img, label='Std of cor')
        plt.title('Standard deviation of correlation coefficient')#  \n Yellow line indicates 0.1')

        plt.subplot(2,2,4)
        img=plt.pcolor(stdcor.lon, stdcor.lat, std_lag, cmap='seismic', vmax=10)
        plt.colorbar(img, label='Std of lag')
        plt.title('Standard deviation of lag')

        fig.subplots_adjust(hspace=0.5, wspace=0.3, top=0.8)

    if save==True:
        print('...Saving second figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_BONUS_2Dcormap_'+datavar+'_'+indvar+'_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
    return



def showCorMapPACEAxis(ax, members='all',
              data=None, datavar='PSL',
              ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
              years=['1920', '2013'],
              detrend=True, deseasonalize=True, window=24, title='2D Correlation Map', ymax=-70, 
              save=False, savefn='test',
              draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
              inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False):
    '''
    Plot correlation map in specific axis for PACE. Use spatial variable with one component (e.g. PSL, SST) from PACE Ensemble.
    The correlation already has to be performed.
    
    INPUT:
    ax (mpl Figure axis): axis on which the correlation map has to be plotted.
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
    ax (mpl Figure axis): axis with correlation map
    '''
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from scipy import stats
    import cartopy.crs as ccrs
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 22})
    
    if members=='all':
        members=[0,20]
    
    #Read the precomputed correlation coefficients
    for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
        if deseasonalize==True:
            if cumsum_map==True:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            else:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        else:
            try:
                temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            except:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    
        if i==0:
            opt_lag=temp_opt_lag
            maxcor=temp_maxcor
        else:
            opt_lag=np.dstack((opt_lag, temp_opt_lag))
            maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
    if i==0:
        maxcor=maxcor.expand_dims('ens')

    
    #Compute statistics
    maxp=stats.ttest_1samp(maxcor, popmean=0)[1]
    opt_lag_p_pos=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='greater')[1]
    opt_lag_p_neg=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='less')[1]
    
    #Compute standard deviation and mean
    stdcor=maxcor.std(dim='ens')
    maxcor=maxcor.mean(dim='ens')
    #If only ensemble mean, than only show 0s.
    if members[-1]!=21:
        std_lag=np.nanstd(opt_lag, axis=2)
        opt_lag=np.nanmean(opt_lag, axis=2)
    else:
        std_lag=np.zeros(np.shape(opt_lag))
    
    #Only show significant data
    maxcor=maxcor.where((maxp<0.05) | maxcor.isnull(), other=0)
    opt_lag[(opt_lag_p_pos>0.05)&(opt_lag_p_neg>0.05)]=0
        
    #Reading Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    print('...Plotting')

    #PLOT CORRELATION COEFFICIENTS
    
    #Regrid the SST correlation coefficients.
    if ('SST' in datavar) | ('aice' in datavar):
        def regridSST(z, x, y):
            from scipy.interpolate import griddata
            PSL=xr.open_dataarray('../02_data/correlation/maxcor_PSL_dotson_to_cosgrove_massloss_dtTrue_dSTrue_w24_ens1.nc')#loadPACEData(var='PSL', members=[0,1])
            xgrid=PSL.lon
            ygrid=PSL.lat
            xi,yi = np.meshgrid(xgrid,ygrid)
            
            zfinal=griddata((x,y),z,(xi,yi),method='linear')
            var_new=xr.DataArray(data=zfinal,
                                dims=['lat', 'lon'],
                                coords=dict(
                                lon=xgrid,
                                lat=ygrid)                        
                                    )
            return var_new
        if 'SST' in datavar:
            maxcor=regridSST(z=maxcor.to_numpy().flatten(),
                             x=maxcor.ULONG.to_numpy().flatten(),
                             y=maxcor.ULAT.to_numpy().flatten())
        else:
            maxcor=regridSST(z=maxcor.to_numpy().flatten(),
                             x=maxcor.ULON.to_numpy().flatten(),
                             y=maxcor.ULAT.to_numpy().flatten())
    
    #Plot axis
    img=ax.pcolormesh(maxcor.lon, maxcor.lat, 
                      maxcor, 
                      cmap='seismic', transform=ccrs.PlateCarree(), vmin=-1, vmax=1)
    CS=ax.contour(maxcor.lon, maxcor.lat, 
                      maxcor, levels=[-0.5, 0.5], colors='yellow',
                      transform=ccrs.PlateCarree())
    ax.clabel(CS, CS.levels, inline=True,fontsize=10)

    ax.set_facecolor('grey')
    plt.colorbar(img, label='Correlation Coefficient', pad=0.10, ax=ax)
    ax.coastlines()
    ax.set_title(title)
    ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree())
  
    return ax


def showLagMapPACEAxis(ax, members='all',
              data=None, datavar='PSL',
              ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
              years=['1920', '2013'],
              detrend=True, deseasonalize=True, window=24, title='2D Correlation Map', ymax=-70, 
              save=False, savefn='test',
              draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
              inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False):
    '''
    Plot lag map in specific axis for PACE. Use spatial variable with one component (e.g. PSL, SST) from PACE Ensemble.
    The correlation already has to be performed.
    
    INPUT:
    ax (mpl Figure axis): axis on which the lag map has to be plotted.
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
    ax (mpl Figure axis): axis with lag map
    '''

    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from scipy import stats
    import cartopy.crs as ccrs
    
    print('Start preparing for correlation map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 22})
    
    if members=='all':
        members=[0,20]
    
    #Read precomputed correlation maps
    for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
        if deseasonalize==True:
            if cumsum_map==True:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_cum'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            else:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        else:
            try:
                temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            except:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datavar+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    
        if i==0:
            opt_lag=temp_opt_lag
            maxcor=temp_maxcor
        else:
            opt_lag=np.dstack((opt_lag, temp_opt_lag))
            maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
    if i==0:
        maxcor=maxcor.expand_dims('ens')

    
    #Compute statistics
    maxp=stats.ttest_1samp(maxcor, popmean=0)[1]
    opt_lag_p_pos=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='greater')[1]
    opt_lag_p_neg=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='less')[1]
    
    #Compute standard deviation and mean
    stdcor=maxcor.std(dim='ens')
    maxcor=maxcor.mean(dim='ens')
    if members[-1]!=21:
        std_lag=np.nanstd(opt_lag, axis=2)
        opt_lag=np.nanmean(opt_lag, axis=2)
    else:
        std_lag=np.zeros(np.shape(opt_lag))
    
    #Only show significant data.
    maxcor=maxcor.where((maxp<0.05) | maxcor.isnull(), other=0)
    opt_lag[std_lag>8]=0
    
    #Reading Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    print('...Plotting')


    #PLOT CORRELATION COEFFICIENTS
    
    #Regrid SST lags.
    if ('SST' in datavar) | ('aice' in datavar):
        def regridSST(z, x, y):
            from scipy.interpolate import griddata
            PSL=xr.open_dataarray('../02_data/correlation/maxcor_PSL_dotson_to_cosgrove_massloss_dtTrue_dSTrue_w24_ens1.nc')
            xgrid=PSL.lon
            ygrid=PSL.lat
            xi,yi = np.meshgrid(xgrid,ygrid)
            
            zfinal=griddata((x,y),z,(xi,yi),method='linear')
            var_new=xr.DataArray(data=zfinal,
                                dims=['lat', 'lon'],
                                coords=dict(
                                lon=xgrid,
                                lat=ygrid)                        
                                    )
            return var_new
        if 'SST' in datavar:
            opt_lag=regridSST(z=opt_lag.flatten(),
                             x=maxcor.ULONG.to_numpy().flatten(),
                             y=maxcor.ULAT.to_numpy().flatten())

            maxcor=regridSST(z=maxcor.to_numpy().flatten(),
                             x=maxcor.ULONG.to_numpy().flatten(),
                             y=maxcor.ULAT.to_numpy().flatten())
        else:
            opt_lag=regridSST(z=opt_lag.flatten(),
                             x=maxcor.ULON.to_numpy().flatten(),
                             y=maxcor.ULAT.to_numpy().flatten())

            maxcor=regridSST(z=maxcor.to_numpy().flatten(),
                             x=maxcor.ULON.to_numpy().flatten(),
                             y=maxcor.ULAT.to_numpy().flatten())

    #Use proper colormap
    if lags[-1]<0:
        divnorm=MidpointNormalize(vmin=lags[0], vmax=0.1, midpoint=0)  
    else:
        divnorm=MidpointNormalize(vmin=lags[0], vmax=lags[-1], midpoint=0)
    
    img=ax.pcolormesh(maxcor.lon, maxcor.lat, 
                      opt_lag, 
                      cmap='seismic', transform=ccrs.PlateCarree(), norm=divnorm)

    plt.colorbar(img, ax=ax, label='Lag [months]', pad=0.10, ticks=np.arange(lags[0], lags[-1]+1, 6))
    
    ax.set_facecolor('grey')
    ax.coastlines()
    ax.set_title(title)
    ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree())
    
    return ax


def createPACECorrelationMaps(save=True):
    '''
    Create correlation map with SST and PSL for main report.
    
    INPUT:
    save (boolean): whether to save the figure
    '''
    window=60
    import cartopy.crs as ccrs
    fig=plt.figure(figsize=(20,10))

    ax = plt.subplot(223, projection=ccrs.PlateCarree(central_longitude=180.0))
    showCorMapPACEAxis(ax,members='all',
                   data=None, datavar='PSL',
                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 0),
                   years=['1920', '2013'],
                   detrend=True, deseasonalize=True, window=window, title='SLP vs Amundsen \n Shelf Break Zonal Wind', ymax=-70, 
                   save=True, savefn='PACE_PSL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                   inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False)
    ax.text(-0.1, 1.1, string.ascii_uppercase[2], transform=ax.transAxes, 
            size=20, weight='bold')
    
    ax = plt.subplot(224, projection=ccrs.PlateCarree(central_longitude=180.0))
    showCorMapPACEAxis(ax,members='all',
                   data=None, datavar='PSL',
                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                   years=['1920', '2013'],
                   detrend=True, deseasonalize=True, window=window, title='SLP vs Dotson \n to Cosgrove Mass Loss', ymax=-70, 
                   save=True, savefn='PACE_PSL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                   inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False)
    ax.text(-0.1, 1.1, string.ascii_uppercase[3], transform=ax.transAxes, 
            size=20, weight='bold')
    
    ax = plt.subplot(221, projection=ccrs.PlateCarree(central_longitude=180.0))
    showCorMapPACEAxis(ax,members='all',
                   data=None, datavar='SST',
                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 0),
                   years=['1920', '2013'],
                   detrend=True, deseasonalize=True, window=window, title='SST vs Amundsen \n Shelf Break Zonal Wind', ymax=-70, 
                   save=True, savefn='PACE_SST_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                   inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False)
    ax.text(-0.1, 1.1, string.ascii_uppercase[0], transform=ax.transAxes, 
            size=20, weight='bold')
    
    ax = plt.subplot(222, projection=ccrs.PlateCarree(central_longitude=180.0))
    showCorMapPACEAxis(ax,members='all',
                   data=None, datavar='SST',
                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                   years=['1920', '2013'],
                   detrend=True, deseasonalize=True, window=window, title='SST vs Ice Shelf Basal Mass \n Loss (Dotson to Cosgrove)', ymax=-70, 
                   save=True, savefn='PACE_PSL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                   inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False)
    ax.text(-0.1, 1.1, string.ascii_uppercase[1], transform=ax.transAxes, 
            size=20, weight='bold')
    
    fig.suptitle('Correlation of Amundsen Sea Zonal Wind and Mass Loss with Global SST and SLP')
                 #'\n Mean over Members 1-20; 60 month running mean; p<0.05')
    plt.rcParams.update({'font.size': 18})
    fig.subplots_adjust(hspace=0.5, wspace=0.5, top=0.85)
    
#     for n, ax in enumerate(fig.axes[::2]):
#         ax.text(-0.25, 0.9, string.ascii_uppercase[n]+'.', transform=ax.transAxes, 
#             size=25, weight='bold')
        
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_all_PACE_corrcoeff.png')
        
    #----LAGS-----

    window=60
    fig=plt.figure(figsize=(20,10))
    
    ax = plt.subplot(223, projection=ccrs.PlateCarree(central_longitude=180.0))
    showLagMapPACEAxis(ax,members='all',
                   data=None, datavar='PSL',
                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 0),
                   years=['1920', '2013'],
                   detrend=True, deseasonalize=True, window=window, title='SLP vs Amundsen \n Shelf Break Zonal Wind', ymax=-70, 
                   save=True, savefn='PACE_PSL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                   inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False)
    
    ax = plt.subplot(224, projection=ccrs.PlateCarree(central_longitude=180.0))
    showLagMapPACEAxis(ax,members='all',
                   data=None, datavar='PSL',
                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                   years=['1920', '2013'],
                   detrend=True, deseasonalize=True, window=window, title='SLP vs Ice Shelf Basal \n Mass Loss (Dotson to Cosgrove)', ymax=-70, 
                   save=True, savefn='PACE_PSL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                   inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False)

    ax = plt.subplot(221, projection=ccrs.PlateCarree(central_longitude=180.0))
    showLagMapPACEAxis(ax,members='all',
                   data=None, datavar='SST',
                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 0),
                   years=['1920', '2013'],
                   detrend=True, deseasonalize=True, window=window, title='SST vs Amundsen \n Shelf Break Zonal Wind', ymax=-70, 
                   save=True, savefn='PACE_SST_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                   inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False)
    
    ax = plt.subplot(222, projection=ccrs.PlateCarree(central_longitude=180.0))
    showLagMapPACEAxis(ax,members='all',
                   data=None, datavar='SST',
                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                   years=['1920', '2013'],
                   detrend=True, deseasonalize=True, window=window, title='SST vs Dotson \n to Cosgrove Mass Loss', ymax=-70, 
                   save=True, savefn='PACE_PSL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                   inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False)
    fig.suptitle('Lags of Amundsen Sea Zonal Wind and Mass Loss with Global SST and SLP \n negative = global forcing is leading')  #\n Mean over Members 1-20; 60 month running mean; p<0.05')
    plt.rcParams.update({'font.size': 18})
    fig.subplots_adjust(hspace=0.5, wspace=0.5, top=0.8)
    
    for n, ax in enumerate(fig.axes[::2]):
        ax.text(-0.25, 0.9, string.ascii_uppercase[n]+'.', transform=ax.transAxes, 
            size=25, weight='bold')
    
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_all_PACE_lags.png')
        
    return


#Derived from: https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
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
    

    

def createMapsOfProcess(save=True):
    window=12
    
    #Import Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    fig=plt.figure(figsize=(20,15))
    ax=plt.subplot(3,2,1)
    ax=showCorMapAxis(ax=ax, members='all',
                  data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Zonal Wind', ymax=None, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=True, cumsum_map=False)

    ax=plt.subplot(3,2,2)
    showCorMapAxis(ax=ax, members='all',
                      data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                      detrend=True, window=window, title='Meridional Wind', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False, deseasonalize=True, cumsum_map=False)
    
    ax=plt.subplot(3,2,3)
    ax=showCorMapAxis(ax=ax, members=[0,14],
                  data=None, datafn='oceFWflx', datavar='oceFWflx', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  detrend=True, window=window, title='Surface Freshwater Flux', ymax=-70, save=True, savefn='oceFWflx_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=True, cumsum_map=False)
    
    ax=plt.subplot(3,2,4)
    ax=showCorMapAxis(ax=ax, members='all',
                  data=None, datafn='SALT_averaged_-200to0', datavar='SALT', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', 
                  lags=range(-24, 24),
                  detrend=True, window=window, 
                  title='Average Salinity (-200 to 0m)', ymax=-70, save=True, 
                  savefn='SALT-200to0_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=True, cumsum_map=False)
    
    
    ax=plt.subplot(3,2,5)
    ax=showCorMapAxis(ax=ax, members='all',
                  data=None, datafn='ADVx_TH_corrected_bottom100m_averaged', datavar='ADVx_TH',\
                  datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss',\
                  lags=range(-24, 24),
                  detrend=True, window=window, title='Heat Advection in Bottom 100m', ymax=-70,\
                  save=True,
                  savefn='bottom100m_ADV_TH_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=False)

    fig.suptitle('Correlations with Dotson to Cosgrove Mass Loss \n {} months running mean'.format(str(window)))
    #plt.rcParams.update({'font.size': 18})
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
    for n, ax in enumerate(fig.axes[::2]):
        ax.text(-0.25, 0.9, string.ascii_uppercase[n]+'.', transform=ax.transAxes, 
            size=25, weight='bold')
    for n, ax in enumerate(fig.axes[::2][:2]):
        ax.set_xticks([-130, -115, -100, -85])
        ax.set_yticks([-75, -70, -65])
        latlon_axes(ax=ax, x=grid.lon_1d, y=grid.lat_1d)
    for n, ax in enumerate(fig.axes[::2][2:]):
        ax.set_xticks([-130, -115, -100, -85])
        ax.set_yticks([-75, -73, -71])
        latlon_axes(ax=ax, x=grid.lon_1d, y=grid.lat_1d[grid.lat_1d<-70])
    
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_process_{}months_corrcoeff.png'.format(str(window)))
        
    #----LAGS-----
    fig=plt.figure(figsize=(20,15))
    
    ax=plt.subplot(3,2,1)
    ax=showLagMapAxis(ax=ax, members='all',
                  data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Zonal Wind', ymax=None, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=True, cumsum_map=False)

    ax=plt.subplot(3,2,2)
    showLagMapAxis(ax=ax, members='all',
                      data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                      detrend=True, window=window, title='Meridional Wind', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False, deseasonalize=True, cumsum_map=False)
    
    ax=plt.subplot(3,2,3)
    ax=showLagMapAxis(ax=ax, members=[0,14],
                  data=None, datafn='oceFWflx', datavar='oceFWflx', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  detrend=True, window=window, title='Surface Freshwater Flux', ymax=-70, save=True, savefn='oceFWflx_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=True, cumsum_map=False)
    
    ax=plt.subplot(3,2,4)
    ax=showLagMapAxis(ax=ax, members='all',
                  data=None, datafn='SALT_averaged_-200to0', datavar='SALT', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', 
                  lags=range(-24, 24),
                  detrend=True, window=window, 
                  title='Average Salinity (-200 to 0m)', ymax=-70, save=True, 
                  savefn='SALT-200to0_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=True, cumsum_map=False)
    
    ax=plt.subplot(3,2,5)
    ax=showLagMapAxis(ax=ax, members='all',
                  data=None, datafn='ADVx_TH_corrected_bottom100m_averaged', datavar='ADVx_TH',\
                  datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss',\
                  lags=range(-24, 24),
                  detrend=True, window=window, title='Heat Advection in Bottom 100m', ymax=-70,\
                  save=True,
                  savefn='bottom100m_ADV_TH_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, deseasonalize=False)

    
    
    fig.suptitle('Lags with Dotson to Cosgrove Mass Loss \n {} months running mean'.format(str(window)))
    plt.rcParams.update({'font.size': 18})
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    
    for n, ax in enumerate(fig.axes[::2]):
        ax.text(-0.25, 0.9, string.ascii_uppercase[n]+'.', transform=ax.transAxes, 
            size=25, weight='bold')
    for n, ax in enumerate(fig.axes[::2][:2]):
        ax.set_xticks([-130, -115, -100, -85])
        ax.set_yticks([-75, -70, -65])
        latlon_axes(ax=ax, x=grid.lon_1d, y=grid.lat_1d)
    for n, ax in enumerate(fig.axes[::2][2:]):
        ax.set_xticks([-130, -115, -100, -85])
        ax.set_yticks([-75, -73, -71])
        latlon_axes(ax=ax, x=grid.lon_1d, y=grid.lat_1d[grid.lat_1d<-70])
        
    
    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_process_{}months_lags.png'.format(str(window)))
    return    
    
    
#def simpleCorrelationMap(members='all',
#                          data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
#                          ind=None, indfn='undercurrent_W115', indvar='UVEL',
#                          detrend=True, window=24, title='Correlation Map', ymax=-70, save=False, savefn='test', cumulative=False,
#                          draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                          inverse=False, ygradient=False, xgradient=False,
#                          give_back=False):
    
    
    
#     print('Start preparing for correlation map!')
#     plt.rcParams.update({'font.size': 14})
#     if type(data)==type(None):
#         if members=='all':
#             print('...Try reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_ensmean.nc')
#             try:
#                 data=xr.open_dataset('/data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_ensmean.nc')
#                 data=data.expand_dims('ens')
#                 print('...Took a shortcut! Read Ensmean file.')
#             except:
#                 try:
#                     data=xr.open_dataset('/data/hpcdata/users/grejan/mitgcm/data/'+datafn+'_ensmean.nc')
#                     data=data.expand_dims('ens')
#                     print('...Took a shortcut! Read Ensmean file.')
#                 except:
#                     print('...Shortcut failed, so...')
#                     data=loadMITgcmData(filename=datafn, members=members, kind=datakind)
#         else:
#             data=loadMITgcmData(filename=datafn, members=members, kind=datakind)
    
#     if xgradient==True:
#         print('...NOTE: Taking XGRADIENT!')
#         if 'XC' in data.dims:
#             data=data.diff('XC')
#         else:
#             data=data.diff('XG')
#     if ygradient==True:
#         print('...NOTE: Taking YGRADIENT!')
#         if 'YC' in data.dims:
#             data=data.diff('YC')
#         else:
#             data=data.diff('YG')
    
    
#     if type(ind)==type(None):
#         print('...Reading: '+indfn)
#         ind=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+indfn+'.nc')
#         if members!='all':
#             print('...Selecting member '+str(members[0]+1))
#             ind=ind.sel(ens=[members[0]])
        
#     if detrend==True:
#         print('...Detrending data')
#         data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
#         ind=ind.to_dataset()
#         ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)
        
#     if window!=1:
#         print('...Applying rolling mean')
#         data=data.rolling(time=window, center=True).mean()
#         ind=ind.rolling(time=window, center=True).mean()
    
#     print('...Computing correlation coefficients')
#     if detrend==True:
#         if cumulative==True:
#             print('...NOTE: Taking CUMULATIVE!')
#             res=xr.corr(data.mean(dim='ens').cumsum(), ind.mean(dim='ens'), dim='time')
#         else:
#             res=xr.corr(data.mean(dim='ens'), ind.mean(dim='ens'), dim='time')
#     else:
#         res=xr.corr(data[datavar].mean(dim='ens'), ind.mean(dim='ens'), dim='time')
    
#     print('...Plotting')
#     fig, ax=plt.subplots(1,1, figsize=(20,10))
    
#     if 'XG' in res.dims:
#         res=res.rename({'XG':'XC'})
#     if 'YG' in res.dims:
#         res=res.rename({'YG':'YC'})
    
#     if 'lon' in res.dims:
#         res=res.rename({'lon':'XC'})
#     if 'lat' in res.dims:
#         res=res.rename({'lat':'YC'})
        
#     if 'nlon' in res.dims:
#         res=res.rename({'nlon':'XC'})
#     if 'nlat' in res.dims:
#         res=res.rename({'nlat':'YC'})
    
#     if type(ymax)!=type(None):
#         res=res.sel(YC=slice(-75.5, ymax))
        
    
#     if members=='all':
#         ax.set_title(title+'\n Ensemble Mean')
#     else:
#         ax.set_title(title+'\n Ensemble Member '+ str(members[-1]))
    
#     if inverse==True:
#         print('...NOTE: Taking the INVERSE of the correlation!')
#         res=-1*res
    
#     img=ax.pcolor(res.XC, res.YC, res, cmap='seismic', vmin=-1, vmax=1)
#     plt.colorbar(img, label='Pearson Correlation Coefficient')
#     CS=ax.contour(res.XC, res.YC, res, np.arange(-1, 1, 0.25))
#     ax.clabel(CS, CS.levels, inline=True,fontsize=10)
    
#     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
#     grid = Grid(gp)
#     bath=ax.contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#     plt.clabel(bath, [-1000])
#     if type(ymax)!=type(None):
#         ax.set_ylim([-75.5, ymax])
#     ax.set_facecolor('grey')
#     ax.set_xlabel('Longitude [deg]')
#     ax.set_ylabel('Latitude [deg]')
    
#     if draw_box==True:
#         ax.plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
#     if save==True:
#         print('...Saving figure')
#         from datetime import date
#         today = date.today()
#         today=today.strftime("%Y%m%d")
#         if members=='all':
#             fig.savefig('../03_output/correlation/'+today+'_cormap_'+savefn+'_ensmean.png')
#         else:
#             fig.savefig('../03_output/correlation/'+today+'_cormap_'+savefn+'_member'+str(members[-1])+'.png')

#     if give_back==False:
#         return
#     else:
#         return data, ind, res

    
    
# def correlationMap_2D(members='all',
#                       data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
#                       ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                       detrend=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False):
#     '''
#     Input: Always give U component as input!
    
#     '''
    
#     print('Start preparing for correlation map!')
#     plt.rcParams.update({'font.size': 14})
    
#     #First, we read the data!
#     if type(data)==type(None):
#             print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_....')
#             dataU=loadMITgcmData(filename=datafn, members=members, kind=datakind)
#             datafn_v=datafn.replace('U','V')
#             datafn_v=datafn_v.replace('u','v')
#             datafn_v=datafn_v.replace('x','y')
#             print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn_v+'_....')
#             dataV=loadMITgcmData(filename=datafn_v, members=members, kind=datakind)
    
#     else:
#         dataU=data[0]
#         dataV=data[1]
    
#     print('...Making sure interpolate everything to the center of the cells')
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
    
    
    
#     if type(ind)==type(None):
#         print('...Reading: '+indfn)
#         if 'timeseries' in indfn:
#             ind=loadTimeSeriesData(filename=indfn, members=members, kind='old')
#         else:        
#             ind=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+indfn+'.nc')
#             if members!='all':
#                 print('...Selecting member '+str(members[0]+1))
#                 ind=ind.sel(ens=[members[0]])
#             ind=ind.to_dataset()
        
#     if detrend==True:
#         print('...Detrending data')
#         data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
#         ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)
        
#     if window!=1:
#         print('...Applying rolling mean')
#         data=data.rolling(time=window, center=True).mean()
#         ind=ind.rolling(time=window, center=True).mean()
    
#     print('...Computing correlation coefficients')
    
#     data=data.stack(ens_time=['ens', 'time'])
#     for z, lag in enumerate(lags):
#         print('...Looking at lag: '+str(lag))
#         ind_temp=ind.shift(time=lag) #Positive means index is leading!
#         ind_temp=ind_temp.stack(ens_time=['ens', 'time'])
        
        
#         if detrend==True:
#             res_temp=xr.corr(data, ind_temp, dim='ens_time')
#         else:
#             res_temp=xr.corr(data[datavar], ind_temp[indvar], dim='ens_time')

#         if z==0:
#             res=res_temp
#         else:
#             res=xr.concat([res, res_temp], dim='lag')
    
#     print('...Plotting')
#     fig, ax=plt.subplots(2,1, figsize=(20,12))
    
#     if type(ymax)!=type(None):
#         res=res.sel(YC=slice(-75.5, ymax))
    
#     if members=='all':
#         fig.suptitle(title+'\n All members')
#     else:
#         if members[1]-members[0]==1:
#             fig.suptitle(title+'\n Ensemble Member '+ str(members[-1]))
#         else:
#             fig.suptitle(title+'\n Ensemble Members '+ str(members[0]+1)+' to '+str(members[-1]))
    
#     if inverse==True:
#         print('...NOTE: Taking the INVERSE of the correlation!')
#         res=-1*res
    
#     maxcor=res.max(dim='lag')
#     #PLOT CORRELATION COEFFICIENTS
#     img=ax[0].pcolor(res.XC, res.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
#     plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
#     CS=ax[0].contour(res.XC, res.YC, maxcor, np.arange(-1, 1, 0.25))
#     ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)
    
#     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
#     grid = Grid(gp)
#     bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#     plt.clabel(bath, [-1000])
#     if type(ymax)!=type(None):
#         ax[0].set_ylim([-75.5, ymax])
#     ax[0].set_facecolor('grey')
#     ax[0].set_xlabel('Longitude [deg]')
#     ax[0].set_ylabel('Latitude [deg]')
#     ax[0].set_title('Optimal Correlation Coefficient')
    
#     if draw_box==True:
#         ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
    
    
#     print('...Finding optimal lag')
#     idx=np.argmax([i=='lag' for i in res.dims])
#     masked=np.ma.masked_invalid(res.to_numpy())
#     opt_lag=np.nanargmax(masked, axis=idx)
#     opt_lag=opt_lag.astype(float)
#     opt_lag[masked.mask.max(axis=idx)]=np.nan
#     opt_lag=opt_lag+lags[0]
    
    
#     #PLOT OPTIMAL LAG
#     img=ax[1].pcolor(res.XC, res.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
#     plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
#     #CS=ax[1].contour(res.XC, res.YC, opt_lag, np.arange(-1, 1, 0.25))
#     #ax[1].clabel(CS, CS.levels, inline=True,fontsize=10)
    
#     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
#     grid = Grid(gp)
#     bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#     plt.clabel(bath, [-1000])
#     if type(ymax)!=type(None):
#         ax[1].set_ylim([-75.5, ymax])
#     ax[1].set_facecolor('grey')
#     ax[1].set_xlabel('Longitude [deg]')
#     ax[1].set_ylabel('Latitude [deg]')
    
#     if draw_box==True:
#         ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
#     ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')
    
#     fig.subplots_adjust(hspace=0.4)
    
#     if save==True:
#         print('...Saving figure')
#         from datetime import date
#         today = date.today()
#         today=today.strftime("%Y%m%d")
#         if members=='all':
#             fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_all.png')
#         else:
#             fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_member'+str(members[-1])+'.png')

#     if give_back==False:
#         return
#     else:
#         return data, ind, res
    
    
# def correlationMap_1D(members='all',
#                       data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
#                       ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                       detrend=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False):
#     '''
#     Input: Only one map
    
#     '''
    
#     print('Start preparing for correlation map!')
#     plt.rcParams.update({'font.size': 14})
    
#     #First, we read the data!
#     if type(data)==type(None):
#             print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_....')
#             data=loadMITgcmData(filename=datafn, members=members, kind=datakind)
    
#     print('...Making sure interpolate everything to the center of the cells')
#     if 'XG' in data.dims:
#         data=data.rename({'XG':'XC'})
#     if 'YG' in data.dims:
#         data=data.rename({'YG':'YC'})

#     #Read something with the correct locations.
#     data1=loadMITgcmData(filename='EXFuwind', members=[0,1], kind='maps')
#     data=data.interp(XC=data1.XC, YC=data1.YC)
    
    
#     if type(ind)==type(None):
#         print('...Reading: '+indfn)
#         if 'timeseries' in indfn:
#             ind=loadTimeSeriesData(filename=indfn, members=members, kind='old')
#         else:        
#             ind=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+indfn+'.nc')
#             if members!='all':
#                 print('...Selecting member '+str(members[0]+1))
#                 ind=ind.sel(ens=[members[0]])
#             ind=ind.to_dataset()
        
#     if detrend==True:
#         print('...Detrending data')
#         data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
#         ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)
        
#     if window!=1:
#         print('...Applying rolling mean')
#         data=data.rolling(time=window, center=True).mean()
#         ind=ind.rolling(time=window, center=True).mean()
    
#     print('...Computing correlation coefficients')
    
#     data=data.stack(ens_time=['ens', 'time'])
#     for z, lag in enumerate(lags):
#         print('...Looking at lag: '+str(lag))
#         ind_temp=ind.shift(time=lag) #Positive means index is leading!
#         ind_temp=ind_temp.stack(ens_time=['ens', 'time'])
        
        
#         if detrend==True:
#             res_temp=xr.corr(data, ind_temp, dim='ens_time')
#         else:
#             res_temp=xr.corr(data[datavar], ind_temp[indvar], dim='ens_time')

#         if z==0:
#             res=res_temp
#         else:
#             res=xr.concat([res, res_temp], dim='lag')
    
#     print('...Plotting')
#     fig, ax=plt.subplots(2,1, figsize=(20,12))
    
#     if type(ymax)!=type(None):
#         res=res.sel(YC=slice(-75.5, ymax))
    
#     if members=='all':
#         fig.suptitle(title+'\n All members')
#     else:
#         if members[1]-members[0]==1:
#             fig.suptitle(title+'\n Ensemble Member '+ str(members[-1]))
#         else:
#             fig.suptitle(title+'\n Ensemble Members '+ str(members[0]+1)+' to '+str(members[-1]))
    
#     if inverse==True:
#         print('...NOTE: Taking the INVERSE of the correlation!')
#         res=-1*res
    
#     maxcor=res.max(dim='lag')
#     #PLOT CORRELATION COEFFICIENTS
#     img=ax[0].pcolor(res.XC, res.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
#     plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
#     CS=ax[0].contour(res.XC, res.YC, maxcor, np.arange(-1, 1, 0.25))
#     ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)
    
#     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
#     grid = Grid(gp)
#     bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#     plt.clabel(bath, [-1000])
#     if type(ymax)!=type(None):
#         ax[0].set_ylim([-75.5, ymax])
#     ax[0].set_facecolor('grey')
#     ax[0].set_xlabel('Longitude [deg]')
#     ax[0].set_ylabel('Latitude [deg]')
#     ax[0].set_title('Optimal Correlation Coefficient')
    
#     if draw_box==True:
#         ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
    
    
#     print('...Finding optimal lag')
#     idx=np.argmax([i=='lag' for i in res.dims])
#     masked=np.ma.masked_invalid(res.to_numpy())
#     opt_lag=np.nanargmax(masked, axis=idx)
#     opt_lag=opt_lag.astype(float)
#     opt_lag[masked.mask.max(axis=idx)]=np.nan
#     opt_lag=opt_lag+lags[0]
    
    
#     #PLOT OPTIMAL LAG
#     img=ax[1].pcolor(res.XC, res.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
#     plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
#     #CS=ax[1].contour(res.XC, res.YC, opt_lag, np.arange(-1, 1, 0.25))
#     #ax[1].clabel(CS, CS.levels, inline=True,fontsize=10)
    
#     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
#     grid = Grid(gp)
#     bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#     plt.clabel(bath, [-1000])
#     if type(ymax)!=type(None):
#         ax[1].set_ylim([-75.5, ymax])
#     ax[1].set_facecolor('grey')
#     ax[1].set_xlabel('Longitude [deg]')
#     ax[1].set_ylabel('Latitude [deg]')
    
#     if draw_box==True:
#         ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
#     ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')
    
#     fig.subplots_adjust(hspace=0.4)
    
#     if save==True:
#         print('...Saving figure')
#         from datetime import date
#         today = date.today()
#         today=today.strftime("%Y%m%d")
#         if members=='all':
#             fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_all.png')
#         else:
#             fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_member'+str(members[-1])+'.png')

#     if give_back==False:
#         return
#     else:
#         return data, ind, res
    
    
# def correlationMap_1D_v2(members='all',
#                       data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
#                       ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                       detrend=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False):
#     '''
#     Input: Only one map
#     New: using xskillscore
    
#     '''
    
#     print('Start preparing for correlation map!')
#     print('NOTE: this is version 2!')
#     plt.rcParams.update({'font.size': 14})
    
#     #First, we read the data!
#     if type(data)==type(None):
#             print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_....')
#             data=loadMITgcmData(filename=datafn, members=members, kind=datakind)
    
#     print('...Making sure interpolate everything to the center of the cells')
#     if 'XG' in data.dims:
#         data=data.rename({'XG':'XC'})
#     if 'YG' in data.dims:
#         data=data.rename({'YG':'YC'})

#     #Read something with the correct locations.
#     data1=loadMITgcmData(filename='EXFuwind', members=[0,1], kind='maps')
#     data=data.interp(XC=data1.XC, YC=data1.YC)
    
    
#     if type(ind)==type(None):
#         print('...Reading: '+indfn)
#         if 'timeseries' in indfn:
#             ind=loadTimeSeriesData(filename=indfn, members=members, kind='old')
#         else:        
#             ind=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+indfn+'.nc')
#             if members!='all':
#                 print('...Selecting member '+str(members[0]+1))
#                 ind=ind.sel(ens=[members[0]])
#             ind=ind.to_dataset()
        
#     if detrend==True:
#         print('...Detrending data')
#         data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
#         ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)
        
#     if window!=1:
#         print('...Applying rolling mean')
#         data=data.rolling(time=window, center=True).mean()
#         ind=ind.rolling(time=window, center=True).mean()
    
#     print('...Computing correlation coefficients')
    
#     data=data.stack(ens_time=['ens', 'time'])
#     for z, lag in enumerate(lags):
#         print('...Looking at lag: '+str(lag))
#         ind_temp=ind.shift(time=lag) #Positive means index is leading!
#         ind_temp=ind_temp.stack(ens_time=['ens', 'time'])
        
        
# #         if detrend==True:
# #             res_temp=xr.corr(data, ind_temp, dim='ens_time')
# #         else:
# #             res_temp=xr.corr(data[datavar], ind_temp[indvar], dim='ens_time')
#         if detrend==True:
#             res_temp=xs.pearson_r(data, ind_temp, dim='ens_time', skipna=True)
#             p_temp=xs.pearson_r_p_value(data, ind_temp, dim='ens_time', skipna=True)
#         else:
#             res_temp=xs.pearson_r(data[datavar], ind_temp[indvar], dim='ens_time', skipna=True)
#             p_temp=xs.pearson_r_p_value(data[datavar], ind_temp[indvar], dim='ens_time', skipna=True)
            
#         if z==0:
#             res=res_temp
#             p=p_temp
            
#         else:
#             res=xr.concat([res, res_temp], dim='lag')
#             p=xr.concat([p, p_temp], dim='lag')
            
#     print('...Plotting')
#     fig, ax=plt.subplots(2,1, figsize=(20,12))
    
#     if type(ymax)!=type(None):
#         res=res.sel(YC=slice(-75.5, ymax))
#         p=p.sel(YC=slice(-75.5, ymax))
        
#     if members=='all':
#         fig.suptitle(title+'\n All members; p<0.05')
#     else:
#         if members[1]-members[0]==1:
#             fig.suptitle(title+'\n Ensemble Member '+ str(members[-1])+'; p<0.05')
#         else:
#             fig.suptitle(title+'\n Ensemble Members '+ str(members[0]+1)+' to '+str(members[-1])+'; p<0.05')
    
#     if inverse==True:
#         print('...NOTE: Taking the INVERSE of the correlation!')
#         res=-1*res
    
# #     print('Now: p')
# #     print(p)
# #     print('Now: cor')
# #     print(res)
    
#     maxcor=res.max(dim='lag')
#     mincor=res.min(dim='lag')
#     maxp=p.min(dim='lag')
    
# #     print('Now: maxp')
#     print(maxp)
    
#     maxcor=maxcor.where((maxp<0.05) | (np.isnan(maxp)), other=0)
#     mincor=mincor.where((maxp<0.05) | (np.isnan(maxp)), other=0)
    
#     maxcor.to_numpy()[maxcor<abs(mincor)]=mincor.to_numpy()[maxcor<abs(mincor)]
    
#     #PLOT CORRELATION COEFFICIENTS
#     img=ax[0].pcolor(res.XC, res.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
#     plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
#     CS=ax[0].contour(res.XC, res.YC, maxcor, np.arange(-1, 1, 0.25))
#     ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)
    
#     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
#     grid = Grid(gp)
#     bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#     plt.clabel(bath, [-1000])
#     if type(ymax)!=type(None):
#         ax[0].set_ylim([-75.5, ymax])
#     ax[0].set_facecolor('grey')
#     ax[0].set_xlabel('Longitude [deg]')
#     ax[0].set_ylabel('Latitude [deg]')
#     ax[0].set_title('Optimal Correlation Coefficient')
    
#     if draw_box==True:
#         ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
    
    
# #     cmap = mpl.cm.seismic
# #     bounds = np.arange(lags[0], lags[-1], 3)
# #     norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
    
    
#     print('...Finding optimal lag')
#     idx=np.argmax([i=='lag' for i in res.dims])
#     masked=np.ma.masked_invalid(res.where((maxp<0.05) & (maxp>0)).to_numpy())
#     opt_lag=np.nanargmax(abs(masked), axis=idx)
#     opt_lag=opt_lag.astype(float)
#     opt_lag[masked.mask.max(axis=idx)]=np.nan
#     opt_lag=opt_lag+lags[0]
    
#     #PLOT OPTIMAL LAG
#     img=ax[1].pcolor(res.XC, res.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
#     plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
#     img=ax[1].pcolor(maxp.where(maxp>0.05).XC, maxp.where(maxp>0.05).YC, maxp.where(maxp>0.05), color=['green'])
#     #CS=ax[1].contour(res.XC, res.YC, opt_lag, np.arange(-1, 1, 0.25))
#     #ax[1].clabel(CS, CS.levels, inline=True,fontsize=10)
    
#     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
#     grid = Grid(gp)
#     bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#     plt.clabel(bath, [-1000])
#     if type(ymax)!=type(None):
#         ax[1].set_ylim([-75.5, ymax])
#     ax[1].set_facecolor('grey')
#     ax[1].set_xlabel('Longitude [deg]')
#     ax[1].set_ylabel('Latitude [deg]')
    
#     if draw_box==True:
#         ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
#     ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')
    
#     fig.subplots_adjust(hspace=0.4)
    
#     if save==True:
#         print('...Saving figure')
#         from datetime import date
#         today = date.today()
#         today=today.strftime("%Y%m%d")
#         if members=='all':
#             fig.savefig('../03_output/correlation/'+today+'_1Dcormap_'+savefn+'_all.png')
#         else:
#             fig.savefig('../03_output/correlation/'+today+'_1Dcormap_'+savefn+'_member'+str(members[-1])+'.png')

#     if give_back==False:
#         return
#     else:
#         return data, ind, res
    
    
# def correlationMap_2D_v2(members='all',
#                       data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
#                       ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                       detrend=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False):
#     '''
#     Input: Always give U component as input!
    
#     '''
    
#     print('Start preparing for correlation map!')
#     print('NOTE: this is version 2!')
#     plt.rcParams.update({'font.size': 14})
    
#     #First, we read the data!
#     if type(data)==type(None):
#             print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn+'_....')
#             dataU=loadMITgcmData(filename=datafn, members=members, kind=datakind)
#             datafn_v=datafn.replace('U','V')
#             datafn_v=datafn_v.replace('u','v')
#             datafn_v=datafn_v.replace('x','y')
#             print('...Reading: /data/hpcdata/users/grejan/mitgcm/02_data/'+datakind+'/'+datafn_v+'_....')
#             dataV=loadMITgcmData(filename=datafn_v, members=members, kind=datakind)
    
#     else:
#         dataU=data[0]
#         dataV=data[1]
    
#     print('...Making sure interpolate everything to the center of the cells')
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
    
    
    
#     if type(ind)==type(None):
#         print('...Reading: '+indfn)
#         if 'timeseries' in indfn:
#             ind=loadTimeSeriesData(filename=indfn, members=members, kind='old')
#         else:        
#             ind=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+indfn+'.nc')
#             if members!='all':
#                 print('...Selecting member '+str(members[0]+1))
#                 ind=ind.sel(ens=[members[0]])
#             ind=ind.to_dataset()
        
#     if detrend==True:
#         print('...Detrending data')
#         data=detrend_and_average_MITgcmData(data, var=datavar, window=1)
#         ind=detrend_and_average_MITgcmData(ind, var=indvar, window=1)
        
#     if window!=1:
#         print('...Applying rolling mean')
#         data=data.rolling(time=window, center=True).mean()
#         ind=ind.rolling(time=window, center=True).mean()
    
#     print('...Computing correlation coefficients')
    
#     data=data.stack(ens_time=['ens', 'time'])
#     for z, lag in enumerate(lags):
#         print('...Looking at lag: '+str(lag))
#         ind_temp=ind.shift(time=lag) #Positive means index is leading!
#         ind_temp=ind_temp.stack(ens_time=['ens', 'time'])
        
#         if detrend==True:
#             res_temp=xs.pearson_r(data, ind_temp, dim='ens_time', skipna=True)
#             p_temp=xs.pearson_r_p_value(data, ind_temp, dim='ens_time', skipna=True)
#         else:
#             res_temp=xs.pearson_r(data[datavar], ind_temp[indvar], dim='ens_time', skipna=True)
#             p_temp=xs.pearson_r_p_value(data[datavar], ind_temp[indvar], dim='ens_time', skipna=True)
            
#         if z==0:
#             res=res_temp
#             p=p_temp
            
#         else:
#             res=xr.concat([res, res_temp], dim='lag')
#             p=xr.concat([p, p_temp], dim='lag')
            
#     print('...Plotting')
#     fig, ax=plt.subplots(2,1, figsize=(20,12))
    
#     if type(ymax)!=type(None):
#         res=res.sel(YC=slice(-75.5, ymax))
#         p=p.sel(YC=slice(-75.5, ymax))
        
#     if members=='all':
#         fig.suptitle(title+'\n All members; p<0.05')
#     else:
#         if members[1]-members[0]==1:
#             fig.suptitle(title+'\n Ensemble Member '+ str(members[-1])+'; p<0.05')
#         else:
#             fig.suptitle(title+'\n Ensemble Members '+ str(members[0]+1)+' to '+str(members[-1])+'; p<0.05')
    
#     if inverse==True:
#         print('...NOTE: Taking the INVERSE of the correlation!')
#         res=-1*res
    
# #     print('Now: p')
# #     print(p)
# #     print('Now: cor')
# #     print(res)
    
#     maxcor=res.max(dim='lag')
#     mincor=res.min(dim='lag')
#     maxp=p.min(dim='lag')
    
# #     print('Now: maxp')
#     print(maxp)
    
#     maxcor=maxcor.where((maxp<0.05) | (np.isnan(maxp)), other=0)
#     mincor=mincor.where((maxp<0.05) | (np.isnan(maxp)), other=0)
    
#     maxcor.to_numpy()[maxcor<abs(mincor)]=mincor.to_numpy()[maxcor<abs(mincor)]
    
#     #PLOT CORRELATION COEFFICIENTS
#     img=ax[0].pcolor(res.XC, res.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
#     plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
#     CS=ax[0].contour(res.XC, res.YC, maxcor, np.arange(-1, 1, 0.25))
#     ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)
    
#     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
#     grid = Grid(gp)
#     bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#     plt.clabel(bath, [-1000])
#     if type(ymax)!=type(None):
#         ax[0].set_ylim([-75.5, ymax])
#     ax[0].set_facecolor('grey')
#     ax[0].set_xlabel('Longitude [deg]')
#     ax[0].set_ylabel('Latitude [deg]')
#     ax[0].set_title('Optimal Correlation Coefficient')
    
#     if draw_box==True:
#         ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
    
    
# #     cmap = mpl.cm.seismic
# #     bounds = np.arange(lags[0], lags[-1], 3)
# #     norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
    
    
#     print('...Finding optimal lag')
#     idx=np.argmax([i=='lag' for i in res.dims])
#     masked=np.ma.masked_invalid(res.where((maxp<0.05) & (maxp>0)).to_numpy())
#     opt_lag=np.nanargmax(abs(masked), axis=idx)
#     opt_lag=opt_lag.astype(float)
#     opt_lag[masked.mask.max(axis=idx)]=np.nan
#     opt_lag=opt_lag+lags[0]
    
#     #PLOT OPTIMAL LAG
#     img=ax[1].pcolor(res.XC, res.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
#     plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
#     img=ax[1].pcolor(maxp.where(maxp>0.05).XC, maxp.where(maxp>0.05).YC, maxp.where(maxp>0.05), color=['green'])
#     #CS=ax[1].contour(res.XC, res.YC, opt_lag, np.arange(-1, 1, 0.25))
#     #ax[1].clabel(CS, CS.levels, inline=True,fontsize=10)
    
#     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
#     grid = Grid(gp)
#     bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#     plt.clabel(bath, [-1000])
#     if type(ymax)!=type(None):
#         ax[1].set_ylim([-75.5, ymax])
#     ax[1].set_facecolor('grey')
#     ax[1].set_xlabel('Longitude [deg]')
#     ax[1].set_ylabel('Latitude [deg]')
    
#     if draw_box==True:
#         ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
#     ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')
    
#     fig.subplots_adjust(hspace=0.4)
    
#     if save==True:
#         print('...Saving figure')
#         from datetime import date
#         today = date.today()
#         today=today.strftime("%Y%m%d")
#         if members=='all':
#             fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_all.png')
#         else:
#             fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_member'+str(members[-1])+'.png')

#     if give_back==False:
#         return
#     else:
#         return data, ind, res
    
    
#    def showCombinedCorrelationMap_2D_v3_old(members='all',
#                       data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
#                       ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                       detrend=True, deseasonalize=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       mask='land', years=['1920', '2013'],
#                       inverse=False, save_results=True,
#                       give_back=False, cumsum_map=False):
#     '''
#     Input: Always give U component as input!
    
#     '''
#     from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
#     print('Start preparing for correlation map!')
#     print('NOTE: this is version 3!')
#     plt.rcParams.update({'font.size': 22})
    
#     if members=='all':
#         members=[0,20]
    
#     for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
#         if deseasonalize==True:
#             if cumsum_map==True:
#                 try:
#                     temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                     temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                     temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                 except:
#                     temp_opt_lag=np.load('../02_data/correlation/opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                     temp_maxp=xr.open_dataarray('../02_data/correlation/maxp_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                     temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#             else:
#                 try:
#                     temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                     temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                     temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                 except:
#                     temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                     temp_maxp=xr.open_dataarray('../02_data/correlation/maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                     temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#         else:
#             try:
#                 temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                 temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                 temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#             except:
#                 try:
#                     temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                     temp_maxp=xr.open_dataarray('../02_data/correlation/maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                     temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                 except:
#                     temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                     temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                     temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    
#         if i==0:
#             opt_lag=temp_opt_lag
#             maxp=temp_maxp
#             maxcor=temp_maxcor
#         else:
#             opt_lag=np.dstack((opt_lag, temp_opt_lag))
#             maxp=xr.concat([maxp, temp_maxp], dim='ens')
#             maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
#     if i==0:
#         maxp=maxp.expand_dims('ens')
#         maxcor=maxcor.expand_dims('ens')
#     if members[-1]!=21:
#         std_lag=np.nanstd(opt_lag, axis=2)
#         opt_lag=np.nanmean(opt_lag, axis=2)
#     else:
#         std_lag=np.zeros(np.shape(opt_lag))
    
#     stdp=maxp.mean(dim='ens')
#     maxp=maxp.mean(dim='ens')
    
#     stdcor=maxcor.std(dim='ens')
#     maxcor=maxcor.mean(dim='ens')
    
#     #Reading Grid
#     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
#     grid = Grid(gp)
    
#     print('...Plotting')
#     fig, ax=plt.subplots(2,1, figsize=(20,14))
    
#     if members[-1]==21:
#         fig.suptitle(title+'\n Ensemble Mean; '+str(window)+' month rolling mean; p<0.05')
#     else:
#         fig.suptitle(title+'\n Mean over members '+ str(members[0]+1)+' to '+ str(members[-1])+'; '+str(window)+' month rolling mean; p<0.05')


#     #PLOT CORRELATION COEFFICIENTS
#     img=ax[0].pcolor(maxcor.XC, maxcor.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
#     plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
#     CS=ax[0].contour(maxcor.XC, maxcor.YC, maxcor, [-0.5, 0.5], colors='yellow')
#     ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)

#     bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#     plt.clabel(bath, [-1000])
#     if type(ymax)!=type(None):
#         ax[0].set_ylim([-75.5, ymax])
#     if 'ice' in mask:
#         ax[0].pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
#     ax[0].set_facecolor('grey')
#     ax[0].set_xlabel('Longitude [deg]')
#     ax[0].set_ylabel('Latitude [deg]')
#     ax[0].set_title('Optimal Correlation Coefficient')

#     if draw_box==True:
#         ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')

#     #PLOT OPTIMAL LAG
#     img=ax[1].pcolor(maxcor.XC, maxcor.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
#     plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
#     img=ax[1].pcolor(maxp.where(maxp>0.05).XC, maxp.where(maxp>0.05).YC, maxp.where(maxp>0.05)*0, cmap='ocean')

#     bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#     plt.clabel(bath, [-1000])
#     if type(ymax)!=type(None):
#         ax[1].set_ylim([-75.5, ymax])

#     if 'ice' in mask:
#         ax[1].pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
#     ax[1].set_facecolor('grey')
#     ax[1].set_xlabel('Longitude [deg]')
#     ax[1].set_ylabel('Latitude [deg]')

#     if draw_box==True:
#         ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
#     ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')

#     fig.subplots_adjust(hspace=0.5, top=0.85)

#     if save==True:
#         print('...Saving figure')
#         from datetime import date
#         today = date.today()
#         today=today.strftime("%Y%m%d")
#         fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
#     fig=plt.figure(figsize=(20,10))
#     plt.suptitle('Extra statistics over ensemble members for : \n'+title)
#     plt.subplot(2,2,1)
#     img=plt.pcolor(maxp.XC, maxp.YC, maxp, vmax=0.05, cmap='seismic')
#     plt.colorbar(img, label='Pvalue')
#     #CS=plt.contour(maxp.XC, maxp.YC, maxp, [1e-4], colors='yellow')
#     #plt.clabel(CS, CS.levels, inline=True,fontsize=10)
#     plt.title('Mean P value with highest lag')# \n Yellow line indicates 1e-4')
    
#     plt.subplot(2,2,2)
#     img=plt.pcolor(stdp.XC, stdp.YC, stdp, vmax=0.05, cmap='seismic')
#     plt.colorbar(img, label='Std of Pvalue')
#     #CS=plt.contour(stdp.XC, stdp.YC, stdp, [1e-4], colors='yellow')
#     #plt.clabel(CS, CS.levels, inline=True,fontsize=10)
#     plt.title('Standard deviation of P value')# \n Yellow line indicates 1e-4')
    
#     plt.subplot(2,2,3)
#     img=plt.pcolor(stdcor.XC, stdcor.YC, stdcor, cmap='seismic', vmax=0.2)
#     plt.colorbar(img, label='Std of cor')
#     #CS=plt.contour(stdcor.XC, stdcor.YC, stdcor, [0.1], colors='yellow')
#     #plt.clabel(CS, CS.levels, inline=True,fontsize=10)
#     plt.title('Standard deviation of correlation coefficient')#  \n Yellow line indicates 0.1')
    
#     plt.subplot(2,2,4)
#     img=plt.pcolor(stdcor.XC, stdcor.YC, std_lag, cmap='seismic', vmax=10)
#     plt.colorbar(img, label='Std of lag')
#     #CS=plt.contour(stdcor.XC, stdcor.YC, std_lag, [4], colors='yellow')
#     #plt.clabel(CS, CS.levels, inline=True,fontsize=10)
#     plt.title('Standard deviation of lag')# \n Yellow line indicates 4 months')
    
#     fig.subplots_adjust(hspace=0.5, wspace=0.3, top=0.8)
    
#     if save==True:
#         print('...Saving second figure')
#         from datetime import date
#         today = date.today()
#         today=today.strftime("%Y%m%d")
#         fig.savefig('../03_output/correlation/'+today+'_BONUS_2Dcormap_'+savefn+'_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
#     return


#def showCombinedCorMap_Barotropic_Baroclinic(save=True):
#     '''
#     INPUT:
#     save (boolean): do we want to save the map.
    
#     '''
#     import sys
#     sys.path.append('/data/hpcdata/users/grejan/mitgcm/')
#     from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
#     from scipy import stats
#     print('Start preparing for correlation map!')
#     print('NOTE: this is version 3!')
#     plt.rcParams.update({'font.size': 22})
    
#     datavar='UVEL'
#     members='all'
#     data=None 
#     datakind='old',
#     ind=None
#     indfn='timeseries_final'
#     indvar='dotson_to_cosgrove_massloss'
#     lags=range(-24, 24),
#     detrend=True
#     window=60
#     title='Correlation between Baroclinic and Barotropic Flow and Dotson to Cosgrove Mass Loss'
#     ymax=-70
#     save=True
#     savefn='surface_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
#     draw_box=False
#     box_x=[360-115, 360-102]
#     box_y=[-71.8, -70.2]
#     inverse=False
#     give_back=False 
#     deseasonalize=False
#     mask='land'
#     #Reading Grid
#     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
#     grid = Grid(gp)
    
 
#     print('...Plotting')
#     fig, ax=plt.subplots(3,1, figsize=(20,14))
#     for count, datafn in enumerate(['baroclinic_UVEL', 'depth_averaged_UVEL']):
#         if members=='all':
#             members=[0,20]

#         for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
#             if deseasonalize==True:
#                 if cumsum_map==True:
#                     try:
#                         temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                         temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                         temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                     except:
#                         temp_opt_lag=np.load('../02_data/correlation/opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                         temp_maxp=xr.open_dataarray('../02_data/correlation/maxp_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                         temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                 else:
#                     try:
#                         temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                         temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                         temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                     except:
#                         temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                         temp_maxp=xr.open_dataarray('../02_data/correlation/maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                         temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#             else:
#                 try:
#                     temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                     temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                     temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                 except:
#                     temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
#                     temp_maxp=xr.open_dataarray('../02_data/correlation/maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#                     temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
#             if i==0:
#                 opt_lag=temp_opt_lag
#                 maxp=temp_maxp
#                 maxcor=temp_maxcor
#             else:
#                 opt_lag=np.dstack((opt_lag, temp_opt_lag))
#                 maxp=xr.concat([maxp, temp_maxp], dim='ens')
#                 maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
#         if i==0:
#             maxp=maxp.expand_dims('ens')
#             maxcor=maxcor.expand_dims('ens')


#         maxp=stats.ttest_1samp(maxcor, popmean=0)[1]
#         opt_lag_p_pos=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='greater')[1]
#         opt_lag_p_neg=stats.ttest_1samp(opt_lag, axis=2, popmean=0, alternative='less')[1]

#         stdcor=maxcor.std(dim='ens')
#         maxcor=maxcor.mean(dim='ens')
#         if members[-1]!=21:
#             std_lag=np.nanstd(opt_lag, axis=2)
#             opt_lag=np.nanmean(opt_lag, axis=2)
#         else:
#             std_lag=np.zeros(np.shape(opt_lag))

#         #maxcor=maxcor.where(maxcor!=0, other=np.nan)
#         maxcor=maxcor.where((maxp<0.05) | maxcor.isnull(), other=0)
#         opt_lag[(opt_lag_p_pos>0.05)&(opt_lag_p_neg>0.05)]=0
        

#         if members[-1]==21:
#             fig.suptitle(title+'\n Ensemble Mean; '+str(window)+' month rolling mean; p<0.05')
#         else:
#             fig.suptitle(title+'\n Mean over members '+ str(members[0]+1)+' to '+ str(members[-1])+'; '+str(window)+' month rolling mean; p<0.05')


#         #PLOT CORRELATION COEFFICIENTS
#         img=ax[count].pcolor(maxcor.XC, maxcor.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
#         plt.colorbar(img, ax=ax[count], label='Max Pearson \n Correlation Coefficient')
#         CS=ax[count].contour(maxcor.XC, maxcor.YC, maxcor, [-0.5, 0.5], colors='yellow')
#         ax[count].clabel(CS, CS.levels, inline=True,fontsize=10)

#         bath=ax[count].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#         plt.clabel(bath, [-1000])
#         if type(ymax)!=type(None):
#             ax[count].set_ylim([-75.5, ymax])
#         if 'ice' in mask:
#             ax[count].pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
#         ax[count].set_facecolor('grey')
#         ax[count].set_xlabel('Longitude [deg]')
#         ax[count].set_ylabel('Latitude [deg]')
#         if count==0:
#             ax[count].set_title('Baroclinic flow: optimal correlation coefficient')
#         if count==1:
#             ax[count].set_title('Barotropic flow: optimal correlation coefficient')

#         if draw_box==True:
#             ax[count].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    
#         if count==0:
#             oldcor=maxcor.copy()
    
    
#     count+=1
#     #PLOT CORRELATION COEFFICIENTS
#     difcor=abs(oldcor)-abs(maxcor)
#     img=ax[count].pcolor(maxcor.XC, maxcor.YC, difcor, cmap='seismic', vmin=-1, vmax=1)
#     plt.colorbar(img, ax=ax[count], label='Max Pearson \n Correlation Coefficient')
#     #CS=ax[count].contour(maxcor.XC, maxcor.YC, difcor, [-0.5, 0.5], colors='yellow')
#     #ax[count].clabel(CS, CS.levels, inline=True,fontsize=10)

#     bath=ax[count].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
#     plt.clabel(bath, [-1000])
#     if type(ymax)!=type(None):
#         ax[count].set_ylim([-75.5, ymax])
#     if 'ice' in mask:
#         ax[count].pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
#     ax[count].set_facecolor('grey')
#     ax[count].set_xlabel('Longitude [deg]')
#     ax[count].set_ylabel('Latitude [deg]')
#     ax[count].set_title('Baroclinic-Barotropic (difference in absolute correlation coeficient)')

#     fig.subplots_adjust(hspace=0.5, top=0.85)

#     if save==True:
#         print('...Saving figure')
#         from datetime import date
#         today = date.today()
#         today=today.strftime("%Y%m%d")
#         fig.savefig('../03_output/correlation/'+today+'_baroclinic_barotropic_diff_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
#    return