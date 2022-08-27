# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain all functions that are needed to create composite maps and animations.

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

#------------------------------------------------------------------------------------------------------------

def makeCompositeAnimation(data, months=[-25,5], show_box=True, title='Composite map or not',
                           filename='map',
                           vmin_pct=0.2, vmax_pct=0.8,
                           gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/',
                           unit='[m/s]',
                           region='amundsen_sea', wvar='dotson_to_cosgrove_massloss', pct=90, larger_than=True):
    
    #From my code to Kaitlin's code
    from mitgcm_python_master.plot_latlon import latlon_plot
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from mitgcm_python_master.utils import real_dir, daily_to_monthly, fix_lon_range, split_longitude, mask_land_ice, moving_average, index_year_start, index_year_end, index_period, mask_2d_to_3d, days_per_month, add_time_dim, z_to_xyz, select_bottom, convert_ismr, mask_except_ice, xy_to_xyz, apply_mask, var_min_max, mask_3d, average_12_months, depth_of_isoline, mask_land, axis_edges
    from mitgcm_python_master.plot_utils.latlon import shade_mask, overlay_vectors
    from mitgcm_python_master.plot_utils.windows import finished_plot, set_panels
    from matplotlib.animation import FuncAnimation
    import matplotlib.animation as animation
    
    fig, ax=plt.subplots(figsize=(14,5))
    moi=months[0]
    
    #if len(np.shape(data))==4:
    #    print('hoi')
    cors=data.sel(time=moi).mean(dim='events')
    #else:
    #    cors=data.sel(time=moi).mean(dim='time')
    try:
        cors=cors.to_array().to_numpy()[:,:]
    except:
        cors=cors.to_numpy()[:,:]
    
    grid = Grid(gp)
    region=region
    mask = grid.get_region_mask(region)
    mask=mask.astype(bool)

    vmax=data.quantile(vmax_pct).to_numpy()
    vmin=data.quantile(vmin_pct).to_numpy()

    #Prepare everything for the map!
    r_data = np.zeros([grid.ny, grid.nx])
    r_data[mask] = cors[mask]
    r_data = mask_land(r_data, grid)
    r_data[np.isnan(r_data)]=0
    
    title_new=title+' (time since time of detection: {} months)'.format(str(moi))
    #if larger_than==True:
    #    title_new=title+'\n (time since time of detection: {} months, threshold: above {} percentile of {})'.format(str(moi), str(pct), wvar)
    #else:
    #    title_new=title+'\n (time since time of detection: {} months, threshold: below {} percentile of {})'.format(str(moi), str(pct), wvar)

    #Plot the Map
    #if vmax==0:
    #    vmax=0.1
    img=latlon_plot(r_data, grid, ctype='plusminus', vmin=vmin, vmax=vmax, title=title_new, titlesize=14, figsize=(14,5), ax=ax, make_cbar=False)
    cbar = plt.colorbar(img)

    cbar.ax.set_ylabel(unit, rotation=270)

    bath=ax.contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])

    show_box=True
    if show_box==True:
        x=[-115, -102, -102, -115, -115]
        y=[-71.8, -71.8, -70.2,  -70.2, -71.8] 
        ax.plot(x, y, c='g', marker='o')


    # animation function
    def animate(i):
        print(i)
        moi=i-25


        cors=data.sel(time=moi).mean(dim='events')
        try:
            cors=cors.to_array().to_numpy()[:,:]
        except:
            cors=cors.to_numpy()[:,:]

        #Prepare everything for the map!
        r_data = np.zeros([grid.ny, grid.nx])
        r_data[mask] = cors[mask]
        r_data = mask_land(r_data, grid)
        r_data[np.isnan(r_data)]=0

        #if larger_than==True:
        title_new=title+' (time since time of detection: {} months)'.format(str(moi))
        #else:
        #    title_new=title+' (time since time of detection: {} months, threshold: below {} percentile of {})'.format(str(moi), str(pct), wvar)

        #Plot the Map
        img=latlon_plot(r_data, grid, ctype='plusminus', vmin=vmin, vmax=vmax, title=title_new, titlesize=14, figsize=(14,5), ax=ax, make_cbar=False)

        bath=ax.contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
        plt.clabel(bath, [-1000])

        show_box=True
        if show_box==True:
            x=[-115, -102, -102, -115, -115]
            y=[-71.8, -71.8, -70.2,  -70.2, -71.8] 
            ax.plot(x, y, c='g', marker='o')

        return 

    
    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")
    
    anim = animation.FuncAnimation(fig, animate, frames=months[-1]-months[0], repeat=False, interval=5000)
    writergif = animation.PillowWriter(fps=2)
    anim.save('../03_output/composite/animation_composite_'+filename+'_'+today+'.gif', writer=writergif)
    return



def createCompositeMap(data, selected, title='Composite map', 
                       filename='map', window=24, method='mean', 
                       longwindow=12*25, 
                       unit='[m/s]',
                       members='all', larger_than=True, pct=90, startyr='1920', endyr='2013', months=[-25, 10]): 
    
    #from unixFunctionsMITgcm_final import makeCompositeAnimation, makeCompositeMap
    
    def month_diff(a, b):
        return round(((a - b)/np.timedelta64(1, 'M')))
    
    #Find the start date per ensemble member.
    start_list=pd.DataFrame()
    for ens in selected.columns:
        start=np.array([False,]*len(selected[ens]))
        for i in range(len(selected[ens])):
            if (selected[ens].iloc[i]==True) & (selected[ens].iloc[i-1]==False):
                start[i]=True
        start_list[ens]=pd.Series(start, index=selected[ens].index)
    
    
    #Load the correct data and give it a relative time name.
    columns=start_list.columns
    columns=columns.sort_values()
    
    
    months=np.arange(months[0], months[-1])
    
    data_new=None
    for month in months:
        b=None
        print('Month is:'+str(month))
        for ensnum, col in enumerate(columns):
            #print(col)
            #print(ensnum)
            true_list=start_list[col][start_list[col]==True]
            for i in range(len(true_list)):
                timestring=(true_list.index[i]+pd.DateOffset(months=month))

                if (timestring<pd.to_datetime(startyr)) | (timestring>pd.to_datetime(endyr)):
                    timestring=true_list.index[i].strftime('%Y-%m-%d')
                    a=data.sel(ens=ensnum).sel(time=timestring)*np.nan
                else:
                    timestring=timestring.strftime('%Y-%m-%d')
                    a=data.sel(ens=ensnum).sel(time=timestring)
                a['time']=[month]
                if type(b)==type(None):
                    b=a
                else:
                    b=xr.concat([b,a], dim='events')
        if type(data_new)==type(None):
            data_new=b
        else:
            data_new=xr.concat([data_new, b], dim='time')
            
            
    #Completing the title
    if members=='all':
        title+=' (full ensemble)'
    else:
        title+=' ({} members)'.format(str(members[-1]-members[0]))
        
    title=title+'\n (both {} month moving average,'.format(str(window))
    if method=='linear':
        title+=' linearly detrended)'
    elif method=='quadratic':
        title+=' quadratically detrended)'
    elif method=='mean':
        title+='  detrended with {} month rolling mean)'.format(str(longwindow))
    
    
    print('Starting with plotting')
    
    makeCompositeAnimation(data=data_new, months=[months[0],months[-1]], title=title, filename=filename, unit=unit, 
                           pct=pct, larger_than=larger_than)
    
    print('Finished!')
    return






def createVelocityCompositeMap(dataU, dataV, selected, title='Composite map', 
                       filename='map', window=24, method='mean', 
                       longwindow=12*25, 
                       unit='[m/s]',
                       members='all', larger_than=True, pct=90, startyr='1920', endyr='2013', months=[-25, 10], vmax=None): 
    
    #from unixFunctionsMITgcm_final import makeCompositeAnimation, makeCompositeMap
    
    def month_diff(a, b):
        return round(((a - b)/np.timedelta64(1, 'M')))
    
    #Find the start date per ensemble member.
    start_list=pd.DataFrame()
    for ens in selected.columns:
        start=np.array([False,]*len(selected[ens]))
        for i in range(len(selected[ens])):
            if (selected[ens].iloc[i]==True) & (selected[ens].iloc[i-1]==False):
                start[i]=True
        start_list[ens]=pd.Series(start, index=selected[ens].index)
    
    
    #Load the correct data and give it a relative time name.
    columns=start_list.columns
    columns=columns.sort_values()
    
    
    months=np.arange(months[0], months[-1])
    
    data=dataU
    data_new=None
    for month in months:
        b=None
        print('Month is:'+str(month))
        for ensnum, col in enumerate(columns):
            #print(col)
            #print(ensnum)
            true_list=start_list[col][start_list[col]==True]
            for i in range(len(true_list)):
                timestring=(true_list.index[i]+pd.DateOffset(months=month))

                if (timestring<pd.to_datetime(startyr)) | (timestring>pd.to_datetime(endyr)):
                    timestring=true_list.index[i].strftime('%Y-%m-%d')
                    a=data.sel(ens=ensnum).sel(time=timestring)*np.nan
                else:
                    timestring=timestring.strftime('%Y-%m-%d')
                    a=data.sel(ens=ensnum).sel(time=timestring)
                a['time']=[month]
                if type(b)==type(None):
                    b=a
                    b=b.expand_dims("events") #added
                else:
                    b=xr.concat([b,a], dim='events')
        if type(data_new)==type(None):
            data_new=b
        else:
            data_new=xr.concat([data_new, b], dim='time')
    
    data_newU=data_new.copy()
    
    data=dataV
    data_new=None
    for month in months:
        b=None
        print('Month is:'+str(month))
        for ensnum, col in enumerate(columns):
            #print(col)
            #print(ensnum)
            true_list=start_list[col][start_list[col]==True]
            for i in range(len(true_list)):
                timestring=(true_list.index[i]+pd.DateOffset(months=month))

                if (timestring<pd.to_datetime(startyr)) | (timestring>pd.to_datetime(endyr)):
                    timestring=true_list.index[i].strftime('%Y-%m-%d')
                    a=data.sel(ens=ensnum).sel(time=timestring)*np.nan
                else:
                    timestring=timestring.strftime('%Y-%m-%d')
                    a=data.sel(ens=ensnum).sel(time=timestring)
                a['time']=[month]
                if type(b)==type(None):
                    b=a
                    b=b.expand_dims("events") #added
                else:
                    b=xr.concat([b,a], dim='events')
        if type(data_new)==type(None):
            data_new=b
        else:
            data_new=xr.concat([data_new, b], dim='time')
    
    data_newV=data_new
    
            
    #Completing the title
    if members=='all':
        title+=' (full ensemble)'
    else:
        title+=' ({} members)'.format(str(members[-1]-members[0]))
        
    title=title+'\n (both {} month moving average,'.format(str(window))
    if method=='linear':
        title+=' linearly detrended)'
    elif method=='quadratic':
        title+=' quadratically detrended)'
    elif method=='mean':
        title+='  detrended with {} month rolling mean)'.format(str(longwindow))
    
    
    print('Starting with plotting')
    
    makeVelocityCompositeAnimation(dataU=data_newU, dataV=data_newV, months=[months[0],months[-1]], title=title, filename=filename, unit=unit, 
                           pct=pct, larger_than=larger_than, vmax=vmax)
    
    print('Finished!')
    return

def makeVelocityCompositeAnimation(dataU, dataV, months=[-25,5], show_box=True, title='Composite map or not',
                           filename='map',
                           vmin_pct=0.3, vmax_pct=0.7,
                           gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/',
                           unit='[m/s]',
                           region='amundsen_sea', wvar='dotson_to_cosgrove_massloss', pct=90, larger_than=True, vmax=None):
    
    #From my code to Kaitlin's code
    from mitgcm_python_master.plot_latlon import latlon_plot
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from mitgcm_python_master.utils import real_dir, daily_to_monthly, fix_lon_range, split_longitude, mask_land_ice, moving_average, index_year_start, index_year_end, index_period, mask_2d_to_3d, days_per_month, add_time_dim, z_to_xyz, select_bottom, convert_ismr, mask_except_ice, xy_to_xyz, apply_mask, var_min_max, mask_3d, average_12_months, depth_of_isoline, mask_land, axis_edges
    from mitgcm_python_master.plot_utils.latlon import shade_mask, overlay_vectors
    from mitgcm_python_master.plot_utils.windows import finished_plot, set_panels
    from matplotlib.animation import FuncAnimation
    import matplotlib.animation as animation
    
    
    
    vfield=np.sqrt(dataU.values**2+dataV.values**2)
    
    print(np.shape(vfield))
    print(np.shape(dataU.time))
    print(np.shape(dataU.events))
    print(np.shape(dataV.XC))
    
    data = xr.DataArray(vfield, coords=[dataU.events, dataU.time, dataU.YC, dataV.XC], dims=[ "events", "time", "YC", "XC"])
    
    fig, ax=plt.subplots(figsize=(14,5))
    moi=months[0]
    
    #if len(np.shape(data))==4:
    #    print('hoi')
    cors=data.sel(time=moi).mean(dim='events')
    #else:
    #    cors=data.sel(time=moi).mean(dim='time')
    try:
        cors=cors.to_array().to_numpy()[:,:]
    except:
        cors=cors.to_numpy()[:,:]
    
    grid = Grid(gp)
    region=region
    mask = grid.get_region_mask(region)
    mask=mask.astype(bool)

    if type(vmax)==type(None):
        vmax=data.quantile(vmax_pct).to_numpy()
        vmin=data.quantile(vmin_pct).to_numpy()
    else:
        vmax=vmax
        vmin=0
    
    #Prepare everything for the map!
    r_data = np.zeros([grid.ny, grid.nx])
    r_data[mask] = cors[mask]
    r_data = mask_land_ice(r_data, grid)
    r_data[np.isnan(r_data)]=0
    
    title_new=title+' (time since time of detection: {} months)'.format(str(moi))
    #if larger_than==True:
    #    title_new=title+'\n (time since time of detection: {} months, threshold: above {} percentile of {})'.format(str(moi), str(pct), wvar)
    #else:
    #    title_new=title+'\n (time since time of detection: {} months, threshold: below {} percentile of {})'.format(str(moi), str(pct), wvar)

    #Plot the Map
    img=latlon_plot(r_data, grid, vmin=vmin, vmax=vmax, title=title_new, ymax=-70, titlesize=14, figsize=(14,5), ax=ax, make_cbar=False)
    cbar = plt.colorbar(img)

    cbar.ax.set_ylabel(unit, rotation=270)

    bath=ax.contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    
    step=6
    ax.quiver(grid.lon_2d[::step, ::step], grid.lat_2d[::step, ::step], 
              dataU.sel(time=moi).mean(dim='events').values[::step, ::step], 
              dataV.sel(time=moi).mean(dim='events').values[::step, ::step], edgecolor='white', linewidth = 0.5,
              color='black')
        
    
#     step=5    
#     overlay_vectors(ax, dataU.sel(time=moi).mean(dim='events').values[::step, ::step], dataV.sel(time=moi).mean(dim='events').values[::step, ::step], grid)#, chunk_x=40, chunk_y=40, scale=1e2, headwidth=2, headlength=3)
    
    show_box=True
    if show_box==True:
        x=[-115, -102, -102, -115, -115]
        y=[-71.8, -71.8, -70.2,  -70.2, -71.8] 
        ax.plot(x, y, c='g', marker='o')


    # animation function
    def animate(i):
        print(i)
        moi=i+months[0]


        cors=data.sel(time=moi).mean(dim='events')
        try:
            cors=cors.to_array().to_numpy()[:,:]
        except:
            cors=cors.to_numpy()[:,:]

        #Prepare everything for the map!
        r_data = np.zeros([grid.ny, grid.nx])
        r_data[mask] = cors[mask]
        r_data = mask_land_ice(r_data, grid)
        r_data[np.isnan(r_data)]=0

        #if larger_than==True:
        title_new=title+' (time since time of detection: {} months)'.format(str(moi))
        #else:
        #    title_new=title+' (time since time of detection: {} months, threshold: below {} percentile of {})'.format(str(moi), str(pct), wvar)

        #Plot the Map
        img=latlon_plot(r_data, grid, vmin=vmin, vmax=vmax, ymax=-70, title=title_new, titlesize=14, figsize=(14,5), ax=ax, make_cbar=False)

        bath=ax.contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
        plt.clabel(bath, [-1000])
        
        
        step=6
        ax.quiver(grid.lon_2d[::step, ::step], grid.lat_2d[::step, ::step], 
                   dataU.sel(time=moi).mean(dim='events').values[::step, ::step], 
                   dataV.sel(time=moi).mean(dim='events').values[::step, ::step], edgecolor='white', linewidth = 0.5,
                   color='black')
        
#         step=5    
#         overlay_vectors(ax, dataU.sel(time=moi).mean(dim='events').values[::step, ::step], dataV.sel(time=moi).mean(dim='events').values[::step, ::step], grid)
        
        show_box=True
        if show_box==True:
            x=[-115, -102, -102, -115, -115]
            y=[-71.8, -71.8, -70.2,  -70.2, -71.8] 
            ax.plot(x, y, c='g', marker='o')

        return 

    
    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")
    
    anim = animation.FuncAnimation(fig, animate, frames=months[-1]-months[0], repeat=False, interval=5000)
    writergif = animation.PillowWriter(fps=2)
    anim.save('../03_output/composite/animation_composite_'+filename+'_'+today+'.gif', writer=writergif)
    return




def defineEvents(total, units, longnames, pct=90, wvar='amundsen_shelf_ohc_below_0m', window=24,
               larger_than=True, method='mean', longwindow=12*25, min_periods=5*12,
               start='1920', end='1958', detrend=True, filename=''):
    
    #Create a long period running mean for detrending. We are interested in decadal variability.
    long=total[wvar][start:end].rolling(longwindow, min_periods=min_periods, center=True).mean()
    
    if detrend==True:
        threshold=np.nanpercentile(total[wvar][start:end].rolling(window, center=True).mean()-long, pct)
    else:
        threshold=np.nanpercentile(total[wvar][start:end].rolling(window, center=True).mean(), pct)
    #windows=[window]
    
    #Creat a dataframe with the selected time periods.
    selected=pd.DataFrame() 
    
    #Create arrays for the time axis
    a=[pd.to_datetime(str(j)) for j in np.arange(1920, 2011, 20)] 
    b=[str(j) for j in np.arange(1920, 2011, 20)]
    xlim=[1920,2020]
    
    lim=total[wvar][start:end].rolling(window, center=True).mean().std().mean()
    
    #Create figure
    plt.figure(figsize=(20,10))
    for i in range(len(total[wvar].columns)):    
        ax=plt.subplot(5,4,i+1)
        
        
        col=total[wvar].columns.sort_values()[i]
        if detrend==True:
            series_dt=total[wvar][col][start:end].rolling(window, center=True).mean()-total[wvar][col][start:end].rolling(longwindow, min_periods=min_periods, center=True).mean()
        else:
            series_dt=total[wvar][col][start:end].rolling(window, center=True).mean()
        series_dt.plot()
        
        if larger_than==True:
            sel=series_dt>threshold
        elif larger_than==False:
            sel=series_dt<threshold
        
        selected[col]=sel
        for j in range(len(sel)):
            if (sel.iloc[j]==True) & (sel.iloc[j-1]==False):
                c=0
                while sel.iloc[j+c]==True:
                    c+=1
                    if j+c==len(sel):
                        break
                plt.axvspan(sel.index[j], sel.index[j+c-1], color='red', alpha=0.3)
        #plt.hlines(0, pd.to_datetime('1910'), pd.to_datetime('2020'), colors='k', linestyles='dashed')
        plt.grid()
        plt.xlim([pd.to_datetime(str(k)) for k in xlim])
        
        plt.title(col)
        plt.ylabel(units[wvar])    
        #plt.ylim([-lim*2.5, lim*2.5])
        plt.xticks(a,b)
        ax.grid(True)


    plt.subplots_adjust(hspace=0.4, wspace=0.4)

    plt.suptitle(longnames[wvar]+';\n {} month rolling mean (centered)'.format(window)+'\n Threshold: {} in: '.format(str(threshold))+longnames[wvar])
    
    
    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")
    
    plt.savefig('../03_output/composite/defined_periods_'+filename+'_'+today+'.png')
    
    return threshold


def findEvents(total, units, longnames, threshold=50, wvar='amundsen_shelf_ohc_below_0m', window=24,
               larger_than=True, method='mean', longwindow=12*25, min_periods=5*12,
               start='1920', end='2013', detrend=True, sep=24, filename=''):
    
    #Create a long period running mean for detrending. We are interested in decadal variability.
    long=total[wvar][start:end].rolling(longwindow, min_periods=min_periods, center=True).mean()
    
    
    #threshold=np.nanpercentile(total[wvar][start:end].rolling(window, center=True).mean()-long, pct)
    #windows=[window]
    
    #Creat a dataframe with the selected time periods.
    selected=pd.DataFrame() 
    
    #Create arrays for the time axis
    a=[pd.to_datetime(str(j)) for j in np.arange(1920, 2011, 20)] 
    b=[str(j) for j in np.arange(1920, 2011, 20)]
    xlim=[1920,2020]
    
    lim=total[wvar][start:end].rolling(window, center=True).mean().std().mean()
    
    #Create figure
    plt.figure(figsize=(20,10))
    for i in range(len(total[wvar].columns)):    
        ax=plt.subplot(5,4,i+1)
        
        
        col=total[wvar].columns.sort_values()[i]
        if detrend==True:
            series_dt=total[wvar][col][start:end].rolling(window, center=True).mean()-total[wvar][col][start:end].rolling(longwindow, min_periods=min_periods, center=True).mean()
        else:
            series_dt=total[wvar][col][start:end].rolling(window, center=True).mean()
        series_dt.plot()
        
        if larger_than==True:
            sel=series_dt>threshold
        elif larger_than==False:
            sel=series_dt<threshold
    
        for j in range(len(sel)):
            if (sel.iloc[j]==True) & (sel.iloc[j-1]==False):
                c=0
                while sel.iloc[j+c]==True:
                    c+=1
                    if j+c==len(sel):
                        break
                plt.axvspan(sel.index[j], sel.index[j+c-1], color='red', alpha=0.3)
                
        for j in range(len(sel)):
            if (sel.iloc[j]==False) & (sel.iloc[j-1]==True):
                c=0
                while (c<sep) and (c+j+1)<len(sel):
                    c+=1
                    if sel.iloc[j+c]==True:
                        print(col)
                        print('Sir, we have one...')
                        plt.axvspan(sel.index[j], sel.index[j+c], color='blue', alpha=0.3)                   
                        sel.iloc[j:j+c]=True
                        break
            selected[col]=sel
        
        
            
        #plt.hlines(0, pd.to_datetime('1910'), pd.to_datetime('2020'), colors='k', linestyles='dashed')
        plt.grid()
        plt.xlim([pd.to_datetime(str(k)) for k in xlim])
        
        plt.title(col)
        plt.ylabel(units[wvar])    
        #plt.ylim([-lim*2.5, lim*2.5])
        plt.xticks(a,b)
        ax.grid(True)


    plt.subplots_adjust(hspace=0.4, wspace=0.4)

    plt.suptitle(longnames[wvar]+';\n {} month rolling mean (centered)'.format(window)+'\n Threshold: {} in: '.format(str(round(threshold)))+longnames[wvar])
    
    
    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")
    
    plt.savefig('../03_output/composite/selected_periods_'+filename+'_'+today+'.png')
    
    return selected, units




def makeCompositeAnimation_2(data, months=[-25,5], show_box=True, title='Composite map or not',
                           filename='map',
                           vmin_pct=0.2, vmax_pct=0.8,
                           gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/',
                           unit='[m/s]',
                           region='amundsen_sea', wvar='dotson_to_cosgrove_massloss', pct=90, larger_than=True):
    
    #From my code to Kaitlin's code
    from mitgcm_python_master.plot_latlon import latlon_plot
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from mitgcm_python_master.utils import real_dir, daily_to_monthly, fix_lon_range, split_longitude, mask_land_ice, moving_average, index_year_start, index_year_end, index_period, mask_2d_to_3d, days_per_month, add_time_dim, z_to_xyz, select_bottom, convert_ismr, mask_except_ice, xy_to_xyz, apply_mask, var_min_max, mask_3d, average_12_months, depth_of_isoline, mask_land, axis_edges
    from mitgcm_python_master.plot_utils.latlon import shade_mask, overlay_vectors
    from mitgcm_python_master.plot_utils.windows import finished_plot, set_panels
    from matplotlib.animation import FuncAnimation
    import matplotlib.animation as animation
    
    fig, ax=plt.subplots(figsize=(14,5))
    moi=months[0]
    
    #if len(np.shape(data))==4:
    #    print('hoi')
    cors=data.sel(time=moi).mean(dim='events')
    #else:
    #    cors=data.sel(time=moi).mean(dim='time')
    try:
        cors=cors.to_array().to_numpy()[:,:]
    except:
        cors=cors.to_numpy()[:,:]
    
    grid = Grid(gp)
    region=region
    mask = grid.get_region_mask(region)
    mask=mask.astype(bool)

    vmax=data.quantile(vmax_pct).to_numpy()
    vmin=data.quantile(vmin_pct).to_numpy()

    #Prepare everything for the map!
    r_data = np.zeros([grid.ny, grid.nx])
    r_data[mask] = cors[mask]
    r_data = mask_land(r_data, grid)
    r_data[np.isnan(r_data)]=0
    
    title_new=title+' (time since time of detection: {} months)'.format(str(moi))
    #if larger_than==True:
    #    title_new=title+'\n (time since time of detection: {} months, threshold: above {} percentile of {})'.format(str(moi), str(pct), wvar)
    #else:
    #    title_new=title+'\n (time since time of detection: {} months, threshold: below {} percentile of {})'.format(str(moi), str(pct), wvar)

    #Plot the Map
    #if vmax==0:
    #    vmax=0.1
    img=latlon_plot(r_data, grid, ctype='plusminus', vmin=vmin, vmax=vmax, title=title_new, titlesize=14, figsize=(14,5), ax=ax, make_cbar=False, ymax=-70)
    cbar = plt.colorbar(img)

    cbar.ax.set_ylabel(unit, rotation=270)

    bath=ax.contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])

    show_box=True
    if show_box==True:
        x=[-115, -102, -102, -115, -115]
        y=[-71.8, -71.8, -70.2,  -70.2, -71.8] 
        ax.plot(x, y, c='g', marker='o')


    # animation function
    def animate(i):
        print(i)
        moi=i-25


        cors=data.sel(time=moi).mean(dim='events')
        try:
            cors=cors.to_array().to_numpy()[:,:]
        except:
            cors=cors.to_numpy()[:,:]

        #Prepare everything for the map!
        r_data = np.zeros([grid.ny, grid.nx])
        r_data[mask] = cors[mask]
        r_data = mask_land(r_data, grid)
        r_data[np.isnan(r_data)]=0

        #if larger_than==True:
        title_new=title+' (time since time of detection: {} months)'.format(str(moi))
        #else:
        #    title_new=title+' (time since time of detection: {} months, threshold: below {} percentile of {})'.format(str(moi), str(pct), wvar)

        #Plot the Map
        img=latlon_plot(r_data, grid, ctype='plusminus', vmin=vmin, vmax=vmax, title=title_new, titlesize=14, figsize=(14,5), ax=ax, make_cbar=False, ymax=-70)

        bath=ax.contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
        plt.clabel(bath, [-1000])

        show_box=True
        if show_box==True:
            x=[-115, -102, -102, -115, -115]
            y=[-71.8, -71.8, -70.2,  -70.2, -71.8] 
            ax.plot(x, y, c='g', marker='o')

        return 

    
    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")
    
    anim = animation.FuncAnimation(fig, animate, frames=months[-1]-months[0], repeat=False, interval=5000)
    writergif = animation.PillowWriter(fps=2)
    anim.save('../03_output/composite/animation_composite_'+filename+'_'+today+'.gif', writer=writergif)
    return



def createCompositeMap_2(data, selected, title='Composite map', 
                       filename='map', window=24, method='mean', 
                       longwindow=12*25, 
                       unit='[m/s]',
                       members='all', larger_than=True, pct=90, startyr='1920', endyr='2013', months=[-25, 10]): 
    
    #from unixFunctionsMITgcm_final import makeCompositeAnimation, makeCompositeMap
    
    def month_diff(a, b):
        return round(((a - b)/np.timedelta64(1, 'M')))
    
    #Find the start date per ensemble member.
    start_list=pd.DataFrame()
    for ens in selected.columns:
        start=np.array([False,]*len(selected[ens]))
        for i in range(len(selected[ens])):
            if (selected[ens].iloc[i]==True) & (selected[ens].iloc[i-1]==False):
                start[i]=True
        start_list[ens]=pd.Series(start, index=selected[ens].index)
    
    
    #Load the correct data and give it a relative time name.
    columns=start_list.columns
    columns=columns.sort_values()
    
    
    months=np.arange(months[0], months[-1])
    
    data_new=None
    for month in months:
        b=None
        print('Month is:'+str(month))
        for ensnum, col in enumerate(columns):
            #print(col)
            #print(ensnum)
            true_list=start_list[col][start_list[col]==True]
            for i in range(len(true_list)):
                timestring=(true_list.index[i]+pd.DateOffset(months=month))

                if (timestring<pd.to_datetime(startyr)) | (timestring>pd.to_datetime(endyr)):
                    timestring=true_list.index[i].strftime('%Y-%m-%d')
                    a=data.sel(ens=ensnum).sel(time=timestring)*np.nan
                else:
                    timestring=timestring.strftime('%Y-%m-%d')
                    a=data.sel(ens=ensnum).sel(time=timestring)
                a['time']=[month]
                if type(b)==type(None):
                    b=a
                else:
                    b=xr.concat([b,a], dim='events')
        if type(data_new)==type(None):
            data_new=b
        else:
            data_new=xr.concat([data_new, b], dim='time')
            
            
    #Completing the title
    if members=='all':
        title+=' (full ensemble)'
    else:
        title+=' ({} members)'.format(str(members[-1]-members[0]))
        
    title=title+'\n (both {} month moving average,'.format(str(window))
    if method=='linear':
        title+=' linearly detrended)'
    elif method=='quadratic':
        title+=' quadratically detrended)'
    elif method=='mean':
        title+='  detrended with {} month rolling mean)'.format(str(longwindow))
    
    
    print('Starting with plotting')
    
    makeCompositeAnimation_2(data=data_new, months=[months[0],months[-1]], title=title, filename=filename, unit=unit, 
                           pct=pct, larger_than=larger_than)
    
    print('Finished!')
    return

def workflowCompositeMap(input_filename, var, members, deseasonalize, selected, output_filename='', title='', window=1):
    
    from functions.loading_and_processing_data import loadMITgcmData, detrend_and_average_MITgcmData, read_all_data
    #Loading the data.
    data=loadMITgcmData(filename=input_filename, members=members)

    #Composite Maps
    data_dt=detrend_and_average_MITgcmData(data, var, window=window)
    if deseasonalize==True:
             data_dt = data_dt.groupby("time.month")-data_dt.groupby("time.month").mean(dim="time")
    createCompositeMap_2(data_dt, selected, title=title, 
                       filename=output_filename,
                       window=window, members=members)
    return
