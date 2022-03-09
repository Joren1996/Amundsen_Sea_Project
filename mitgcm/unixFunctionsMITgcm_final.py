# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 17:29:53 2022

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

#%% READING AND LOADING DATA

def readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all'):
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
        if var=='THETA':
            b=b.sel(Z=slice(-200, -700)).mean(dim='Z')

        for f in fd[1:]:
            a=xr.open_dataset(f)
            a=a[[var]]
            if var=='THETA':
                a=a.sel(Z=slice(-200, -700)).mean(dim='Z')
            b=xr.concat([b,a], dim='time')

        if i==0:
            data=b
        else:
            data=xr.concat([data,b], dim='ens')
        
    #Making sure that the dataset has the dimension 'ens'
    if number==1:
        data=data.expand_dims("ens")
            
    return data

def readMITgcmData_depth_averaged(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', z0=[-700, -200], save=False):
    '''
    var: name of variable (string).
    bd: base directory, with PAS_PACE folders.
    members: members to read ('all' or list with boundaries (e.g. [0, 1] for member 1)).
    '''
    
    from mitgcm_python_master.utils import mask_2d_to_3d, apply_mask
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from mitgcm_python_master.calculus import vertical_average
    
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    mask = np.ones([grid.ny, grid.nx]).astype(bool)
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
        if var=='THETA':
            numdat=b[var].to_numpy()
            numdat = apply_mask(numdat, np.invert(mask), time_dependent=True)
            numdat=vertical_average (numdat, grid, gtype='t', time_dependent=True)
            b =xr.DataArray(numdat, dims=('time', "YC", "XC"), coords={'time': b.time, "YC": b.YC, "XC": b.XC})
            b=b.to_dataset(name='THETA')
            #b=b.sel(Z=slice(-200, -700)).mean(dim='Z')

        for f in fd[1:]:
            a=xr.open_dataset(f)
            a=a[[var]]
            if var=='THETA':
                numdat=a[var].to_numpy()
                numdat = apply_mask(numdat, np.invert(mask), time_dependent=True)
                numdat=vertical_average (numdat, grid, gtype='t', time_dependent=True)
                a = xr.DataArray(numdat, dims=('time', "YC", "XC"), coords={'time': a.time, "YC": a.YC, "XC": a.XC})
                a=a.to_dataset(name='THETA')
                #a=a.sel(Z=slice(-200, -700)).mean(dim='Z')
            b=xr.concat([b,a], dim='time')
        if save==True:
            b.to_netcdf('./data/new2_amundsen_shelf_temp_200to700_ens'+str(i+1)+'.nc')
        if i==0:
            data=b
        else:
            data=xr.concat([data,b], dim='ens')
#    for i in range(len(data.ens)):
#        
    #Making sure that the dataset has the dimension 'ens'
    if number==1:
        data=data.expand_dims("ens")
            
    return data


def loadMITgcmData(filename='zonal_winds', members='all'):
    if members=='all':
        ens_list=np.arange(20)
    else:
        ens_list=np.arange(members[0], members[-1])
        
    for i in ens_list:
        print('Loading number:'+str(i))
        a=xr.open_dataset('./data/'+filename+'_ens'+str(i+1)+'.nc')
        if i==0:
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
    
def readMITgcmTimeSeries(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', tf='timeseries_final.nc'):
    '''
    var: name of variable (string).
    bd: base directory, with PAS_PACE folders.
    members: members to read ('all' or list with boundaries (e.g. [0, 1] for member 1)).
    tf: timeseries file name (string)
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
    
#%% CREATING CORRELATION MAPS
def findOptimalLag(ts, tm):
    lags=np.arange(-24, 24)
    r=np.zeros(len(lags))

    for z, lag in enumerate(lags):
        ts_temp=ts.shift(time=lag)
        tma=tm.stack(ens_time=['ens', 'time'])
        tsa=ts_temp.stack(ens_time=['ens', 'time'])
        a=xr.corr(tma, tsa, dim='ens_time')
        r[z]=a

    plt.plot(lags, r)
    maxind=np.argmax(r)
    final_lag=lags[maxind]
    
    return final_lag

def computeCorrelation(data, ts, final_lag):
    #Correlation
    ts_temp=ts.shift(time=final_lag)
    da=data.stack(ens_time=['ens', 'time'])
    tsa=ts_temp.stack(ens_time=['ens', 'time'])

    #This is a bit tricky.... Drop the NaN values...
    da=da.drop_sel(ens_time=tsa[np.isnan(tsa)].ens_time)
    tsa=tsa.dropna(dim='ens_time')
    
    tsa=tsa.drop_sel(ens_time=da[0,0,:][np.isnan(da[0,0,:])].ens_time)
    da=da.dropna(dim='ens_time')
    
    cors=np.zeros((len(da.YC), len(da.XC)))
    p=np.zeros((len(da.YC), len(da.XC)))

    for i in range(len(da.YC)):
        for j in range(len(da.XC)):
            cors[i,j], p[i,j]=stats.pearsonr(da.values[i,j,:],tsa)
            
    print('Finished!')
    
    return cors, p

def makeCorrMap(cors, p=None, final_lag=-9, title='Correlation Map', gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/', region='amundsen_sea', ymax=None, show_p=False, show_box=True):
    #From my code to Kaitlin's code
    from mitgcm_python_master.plot_latlon import latlon_plot
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from mitgcm_python_master.utils import real_dir, daily_to_monthly, fix_lon_range, split_longitude, mask_land_ice, moving_average, index_year_start, index_year_end, index_period, mask_2d_to_3d, days_per_month, add_time_dim, z_to_xyz, select_bottom, convert_ismr, mask_except_ice, xy_to_xyz, apply_mask, var_min_max, mask_3d, average_12_months, depth_of_isoline, mask_land, axis_edges
    from mitgcm_python_master.plot_utils.latlon import shade_mask, overlay_vectors

    grid = Grid(gp)
    mask = grid.get_region_mask(region)
    mask=mask.astype(bool)
    
    #Prepare everything for the map!
    r_data = np.zeros([grid.ny, grid.nx])

    r_data[mask] = cors[mask]

    r_data = mask_land_ice(r_data, grid)
    print(('Mean r over all points and ensemble members = '+str(np.nanmean(r_data))))
    r_data[np.isnan(r_data)]=0

    vmin=np.min(r_data[mask])
    vmax=np.max(r_data[mask])      

    levels=np.arange(0, 1, 0.1)
    
    #Plot the Map
    fig, ax=plt.subplots(figsize=(14,5))

    latlon_plot(r_data, grid, ctype='plusminus', ymax=ymax, vmin=None, vmax=None, title=title+' at optimum lag of '+str(final_lag)+' months', titlesize=14, figsize=(14,5), ax=ax)
    r=ax.contour(grid.lon_2d, grid.lat_2d, r_data, levels=levels, colors=['yellow'])
    plt.clabel(r, levels)

    bath=ax.contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])

    if show_p==True:
        p_data = np.zeros([grid.ny, grid.nx])
        p_data[mask] = p[mask]
        p_data = mask_land_ice(p_data, grid)
        a=5
        ax.scatter(grid.lon_2d[::a, ::a][(p_data[::a, ::a]<=0.05) & (p_data[::a, ::a].mask==False)], grid.lat_2d[::a, ::a][(p_data[::a, ::a]<=0.05) & (p_data[::a, ::a].mask==False)], c='k', s=.02)

    if show_box==True:
        x=[-115, -102, -102, -115, -115]
        y=[-71.8, -71.8, -70.2,  -70.2, -71.8] 
        ax.plot(x, y, c='g', marker='o')
        
    return fig

def createCorrelationMap(data, ts, tm, title='Correlation map', 
                         filename='correlation_map', method='mean', 
                         window=24, long=25*12, members='all'):
    
    from mitgcm_python_master.plot_utils.windows import finished_plot, set_panels
    
    #Completing the title
    title=title+'\n (both {} month moving average,'.format(str(window))
    if method=='linear':
        title+=' linearly detrended)'
    elif method=='quadratic':
        title+=' quadratically detrended)'
    elif method=='mean':
        title+='  detrended with {} month rolling mean)'.format(str(long))
    if members=='all':
        title+=' (full ensemble)'
    else:
        title+=' ({} members)'.format(str(members[-1]-members[0]))
    
    final_lag=findOptimalLag(ts, tm)
    cors,p=computeCorrelation(data, ts, final_lag)
    fig=makeCorrMap(cors=cors, p=p, final_lag=final_lag, title=title)
    
    
    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")
    
    finished_plot(fig, fig_name='./figures/'+filename+today+'.png', dpi=300)
    
    return

#%% CREATING COMPOSITES

def read_all_data(full, total, units, longnames, kind):
    file_paths=[os.path.join(full, file) for file in os.listdir(full) if kind in file]
    
    
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

def show_selected_periods(sel, cind, indnms=['ONI'], windows=[1], xlim=[1920,2020], col=None, title=''):
    #plt.figure(figsize=(20,10))
    plt.title(title)
    for i, nm in enumerate(indnms):
        #plt.subplot(np.ceil(len(indnms)/2), 2, i+1)
        if col==None:
            cind[nm].rolling(windows[i], center=True).mean().plot()
        else:
            cind[nm][col].rolling(windows[i], center=True).mean().plot()
        try:
            sel=sel.to_series()
        except:
            sel=sel
        for j in range(len(sel)):
            if (sel.iloc[j]==True) & (sel.iloc[j-1]==False):
                c=0
                while sel.iloc[j+c]==True:
                    c+=1       
                plt.axvspan(sel.index[j], sel.index[j+c-1], color='red', alpha=0.3)
        plt.hlines(0, pd.to_datetime('1910'), pd.to_datetime('2020'), colors='k', linestyles='dashed')
        plt.grid()
        plt.xlim([pd.to_datetime(str(k)) for k in xlim])
        #plt.ylabel(nm)
        #plt.title(nm+'; {} month rolling mean (centered)'.format(windows[i]))
    return

#STEP 1: FINDING MELT EVENTS

def findHighMeltEvents(wvar='dotson_to_cosgrove_massloss', window=24,
                       pct=90, larger_than=True,
                       method='mean', longwindow=12*25, min_periods=5*12,
                       start='1920', end='2013',
                       full='/data/hpcdata/users/grejan/mitgcm/data'):
    
    total={}
    units={}
    longnames={}
    total, units, longnames=read_all_data(full, total, units, longnames, kind='isotherm')
    total, units, longnames=read_all_data(full, total, units, longnames, kind='final')
    
    
    plt.figure(figsize=(20,10))

    long=total[wvar][start:end].rolling(longwindow, min_periods=min_periods, center=True).mean()

    threshold=np.nanpercentile(total[wvar][start:end].rolling(window, center=True).mean()-long, pct)
    windows=[window]

    selected=pd.DataFrame()

    a=[pd.to_datetime(str(j)) for j in np.arange(1920, 2011, 20)]
    b=[str(j) for j in np.arange(1920, 2011, 20)]
    for i in range(len(total[wvar].columns)):    
        ax=plt.subplot(5,4,i+1)
        col=total[wvar].columns[i]
        if larger_than==True:
            sel=(total[wvar][col][start:end].rolling(window, center=True).mean()-total[wvar][col][start:end].rolling(longwindow, min_periods=min_periods, center=True).mean())>threshold
        elif larger_than==False:
            sel=(total[wvar][col][start:end].rolling(window, center=True).mean()-total[wvar][col][start:end].rolling(longwindow, min_periods=min_periods, center=True).mean())<threshold

        selected[col]=sel

        show_selected_periods(sel, total, indnms=[wvar], windows=windows, xlim=[1920,2020], col=col, title='')

        plt.ylabel(units[wvar])    
        plt.ylim([0, 250])
        plt.xticks(a,b)
        ax.grid(True)


    plt.subplots_adjust(hspace=0.4, wspace=0.4)

    plt.suptitle(longnames[wvar]+';\n {} month rolling mean (centered)'.format(windows[0])+'\n Threshold: {} percentile of: '.format(pct)+longnames[wvar])
    
    return selected, units

def makeCompositeAnimation(data, months=[-25,5], show_box=True, title='Composite map',
                           filename='map',
                           vmin_pct=0.3, vmax_pct=0.7,
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
    cors=data.sel(time=moi).mean(dim='time')
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
    r_data = mask_land_ice(r_data, grid)
    r_data[np.isnan(r_data)]=0
    if larger_than==True:
        title_new=title+'\n (time since time of detection: {} months, threshold: above {} percentile of {})'.format(str(moi), str(pct), wvar)
    else:
        title_new=title+'\n (time since time of detection: {} months, threshold: below {} percentile of {})'.format(str(moi), str(pct), wvar)

    #Plot the Map
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


        cors=data.sel(time=moi).mean(dim='time')
        try:
            cors=cors.to_array().to_numpy()[:,:]
        except:
            cors=cors.to_numpy()[:,:]

        #Prepare everything for the map!
        r_data = np.zeros([grid.ny, grid.nx])
        r_data[mask] = cors[mask]
        r_data = mask_land_ice(r_data, grid)
        r_data[np.isnan(r_data)]=0

        if larger_than==True:
            title_new=title+'\n (time since time of detection: {} months, threshold: above {} percentile of {})'.format(str(moi), str(pct), wvar)
        else:
            title_new=title+'\n (time since time of detection: {} months, threshold: below {} percentile of {})'.format(str(moi), str(pct), wvar)

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
    anim.save('./figures/animation_composite_'+filename+'_'+today+'.gif', writer=writergif)
    return

def makeCompositeMap(data, months=[-25,3], show_box=True, title='Composite map',
                     filename='map',
                     vmin=None, vmax=None, 
                     gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/',
                     unit='[m/s]',
                     region='amundsen_sea', wvar='dotson_to_cosgrove_massloss', pct=90, larger_than=True):
    
    #From my code to Kaitlin's code
    from mitgcm_python_master.plot_latlon import latlon_plot
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from mitgcm_python_master.utils import real_dir, daily_to_monthly, fix_lon_range, split_longitude, mask_land_ice, moving_average, index_year_start, index_year_end, index_period, mask_2d_to_3d, days_per_month, add_time_dim, z_to_xyz, select_bottom, convert_ismr, mask_except_ice, xy_to_xyz, apply_mask, var_min_max, mask_3d, average_12_months, depth_of_isoline, mask_land, axis_edges
    from mitgcm_python_master.plot_utils.latlon import shade_mask, overlay_vectors
    from mitgcm_python_master.plot_utils.windows import finished_plot, set_panels

    for moi in np.arange(months[0],months[-1]):
        #VISUALISATION:
        cors=data.sel(time=moi).mean(dim='time')
        try:
            cors=cors.to_array().to_numpy()[:,:]
        except:
            cors=cors.to_numpy()[:,:]
        
        grid = Grid(gp)
        region=region
        mask = grid.get_region_mask(region)
        mask=mask.astype(bool)

        #Prepare everything for the map!
        r_data = np.zeros([grid.ny, grid.nx])
        r_data[mask] = cors[mask]
        r_data = mask_land_ice(r_data, grid)
        r_data[np.isnan(r_data)]=0

        if larger_than==True:
            title_new=title+'\n (time since time of detection: {} months, threshold: above {} percentile of {})'.format(str(moi), str(pct), wvar)
        else:
            title_new=title+'\n (time since time of detection: {} months, threshold: below {} percentile of {})'.format(str(moi), str(pct), wvar)

        #Plot the Map
        fig, ax=plt.subplots(figsize=(14,5))
        img=latlon_plot(r_data, grid, ctype='plusminus', vmin=vmin, vmax=vmax, title=title_new, titlesize=14, figsize=(14,5), ax=ax, make_cbar=False)
        cbar = plt.colorbar(img)

        cbar.ax.set_ylabel(unit, rotation=270)


        bath=ax.contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
        plt.clabel(bath, [-1000])

        if show_box==True:
            x=[-115, -102, -102, -115, -115]
            y=[-71.8, -71.8, -70.2,  -70.2, -71.8] 
            ax.plot(x, y, c='g', marker='o')

    
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
    
        finished_plot(fig, fig_name='./figures/composite_'+filename+'_'+str(moi)+'_'+today+'.png', dpi=300)
    return

def createCompositeMap(data, title='Composite map', 
                       filename='map', window=24, method='mean', 
                       longwindow=12*25, 
                       unit='[m/s]',
                       members='all', larger_than=True, pct=90): 
    
    def month_diff(a, b):
        return round(((a - b)/np.timedelta64(1, 'M')))
    
    #First, we need to identify the moments of interest.
    selected, units=findHighMeltEvents(wvar='dotson_to_cosgrove_massloss', window=24,
                       pct=pct, larger_than=larger_than,
                       method='mean', longwindow=12*25, min_periods=5*12,
                       start='1920', end='2013',
                       full='/data/hpcdata/users/grejan/mitgcm/data')
    
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
    
    for i, ens in enumerate(columns):
        print('Ens is: '+str(ens))
        years=start_list[ens][start_list[ens]==True].index.year.values
        for year in years:
            a=data.sel(ens=i).sel(time=slice(str(year-2), str(year)))
            if year==years[0]:
                b=a
            else:
                b=xr.concat([b,a], dim='time')
                
        times=b.indexes['time'].to_datetimeindex()

        reltime=np.ones(len(times))*99 #relative time to the the starting month [months]

        kappa=0
        for iota, t in enumerate(times):
            if start_list[ens][start_list[ens]==True].index[kappa].year>=t.year:
                reltime[iota]=month_diff(t, start_list[ens][start_list[ens]==True].index[kappa])
            else:
                kappa+=1
                reltime[iota]=month_diff(t, start_list[ens][start_list[ens]==True].index[kappa])

        b['time']=reltime        
                
        
        if i==0:
            data_new=b
        else:
            data_new=xr.concat([data_new,b], dim='time')
            
            
    #Completing the title
    title=title+'\n (both {} month moving average,'.format(str(window))
    if method=='linear':
        title+=' linearly detrended)'
    elif method=='quadratic':
        title+=' quadratically detrended)'
    elif method=='mean':
        title+='  detrended with {} month rolling mean)'.format(str(longwindow))
    if members=='all':
        title+=' (full ensemble)'
    else:
        title+=' ({} members)'.format(str(members[-1]-members[0]))
    
    print('Starting with plotting')
    #When data is read, let's start with the plotting!
    makeCompositeMap(data=data_new, months=[-25,3], title=title, filename=filename, unit=unit, 
                     pct=pct, larger_than=larger_than)
    makeCompositeAnimation(data=data_new, months=[-25,3], title=title, filename=filename, unit=unit, 
                           pct=pct, larger_than=larger_than)
    
    print('Finished!')
    
    return
    
    
#------------------------------------------------------------------------------------------------------------

def makeCompositeAnimation2(data, months=[-25,5], show_box=True, title='Composite map or not',
                           filename='map',
                           vmin_pct=0.3, vmax_pct=0.7,
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
    r_data = mask_land_ice(r_data, grid)
    r_data[np.isnan(r_data)]=0
    
    title_new=title+' (time since time of detection: {} months)'.format(str(moi))
    #if larger_than==True:
    #    title_new=title+'\n (time since time of detection: {} months, threshold: above {} percentile of {})'.format(str(moi), str(pct), wvar)
    #else:
    #    title_new=title+'\n (time since time of detection: {} months, threshold: below {} percentile of {})'.format(str(moi), str(pct), wvar)

    #Plot the Map
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
        r_data = mask_land_ice(r_data, grid)
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
    anim.save('./figures/animation_composite_'+filename+'_'+today+'.gif', writer=writergif)
    return



def createCompositeMap2(data, selected, title='Composite map', 
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
    
    makeCompositeAnimation2(data=data_new, months=[months[0],months[-1]], title=title, filename=filename, unit=unit, 
                           pct=pct, larger_than=larger_than)
    
    print('Finished!')
    return

#---------------------------------------------------------------------------------------------------------------

def findEvents(total, units, longnames, threshold=50, wvar='dotson_to_cosgrove_massloss', window=24,
               larger_than=True, method='mean', longwindow=12*25, min_periods=5*12,
               start='1920', end='2013', filename=''):
    
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
        series_dt=total[wvar][col][start:end].rolling(window, center=True).mean()-total[wvar][col][start:end].rolling(longwindow, min_periods=min_periods, center=True).mean()
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
                plt.axvspan(sel.index[j], sel.index[j+c-1], color='red', alpha=0.3)
        plt.hlines(0, pd.to_datetime('1910'), pd.to_datetime('2020'), colors='k', linestyles='dashed')
        plt.grid()
        plt.xlim([pd.to_datetime(str(k)) for k in xlim])
        
        plt.title(col)
        plt.ylabel(units[wvar])    
        plt.ylim([-lim*2.5, lim*2.5])
        plt.xticks(a,b)
        ax.grid(True)


    plt.subplots_adjust(hspace=0.4, wspace=0.4)

    plt.suptitle(longnames[wvar]+';\n {} month rolling mean (centered)'.format(window)+'\n Threshold: {} in: '.format(str(threshold))+longnames[wvar])
    
    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")
    
    plt.savefig('./figures/selected_periods_'+filename+'_'+today+'.png')
    
    return selected, units
