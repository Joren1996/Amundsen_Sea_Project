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
import matplotlib.colors as colors

#%% READING AND LOADING DATA

def readMITgcmData(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', save=False, filename='false'):
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
        if var=='ADVx_TH':
            b=b.sel(Z=-405)
        if var=='ADVy_TH':
            b=b.sel(Z=-405)
            


        for f in fd[1:]:
            a=xr.open_dataset(f)
            a=a[[var]]
            if var=='THETA':
                a=a.sel(Z=slice(-200, -700)).mean(dim='Z')
            if var=='ADVx_TH':
                a=a.sel(Z=-405)
            if var=='ADVy_TH':
                a=a.sel(Z=-405)
            b=xr.concat([b,a], dim='time')

        if i==0:
            data=b
        else:
            data=xr.concat([data,b], dim='ens')
            
        if save==True:
            if filename=='false':
                b.to_netcdf('./data/'+var+'_ens'+str(i+1)+'.nc')
            else:
                b.to_netcdf('./data/'+filename+'_'+var+'_ens'+str(i+1)+'.nc')
        
    #Making sure that the dataset has the dimension 'ens'
    if number==1:
        data=data.expand_dims("ens")
            
    return data

def readMITgcmData_depth_averaged(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', z0=[-700, -200], save=False, bottom=False):
    '''
    var: name of variable (string).
    bd: base directory, with PAS_PACE folders.
    members: members to read ('all' or list with boundaries (e.g. [0, 1] for member 1)).
    '''
    
    from mitgcm_python_master.utils import mask_2d_to_3d, apply_mask, select_bottom
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from mitgcm_python_master.calculus import vertical_average
    
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    if bottom==False:
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
        
        if ('VEL' in var) |  ('ADV' in var):
            numdat=b[var].to_numpy()
            
            if bottom==True:
                print('Apply Bottom Mask')
                depth=b[var].Z.to_numpy()
                depth = np.stack([depth]*np.shape(grid.bathy)[0], axis=1)
                depth = np.stack([depth]*np.shape(grid.bathy)[1], axis=2)
                mask=(numdat[-1,:,:,:]==select_bottom(numdat[-1,:,:,:], grid=grid, masked=False)) | (depth<=(grid.bathy+100)) & (numdat[-1,:,:]!=0)
                numdat = apply_mask(numdat, np.invert(mask), time_dependent=True)
            numdat=vertical_average (numdat, grid, gtype='t', time_dependent=True)
            if 'U' in var:
                b =xr.DataArray(numdat, dims=('time', "YC", "XG"), coords={'time': b.time, "YC": b.YC, "XG": b.XG})
            else:
                b =xr.DataArray(numdat, dims=('time', "YG", "XC"), coords={'time': b.time, "YG": b.YG, "XC": b.XC})
            b=b.to_dataset(name=var)
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
            if ('VEL' in var) |  ('ADV' in var):
                numdat=a[var].to_numpy()
                if bottom==True:
                    numdat = apply_mask(numdat, np.invert(mask), time_dependent=True)
                numdat=vertical_average (numdat, grid, gtype='t', time_dependent=True)
                if 'U' in var:
                    a = xr.DataArray(numdat, dims=('time', "YC", "XG"), coords={'time': a.time, "YC": a.YC, "XG": a.XG})
                else:
                    a = xr.DataArray(numdat, dims=('time', "YG", "XC"), coords={'time': a.time, "YG": a.YG, "XC": a.XC})
                a=a.to_dataset(name=var)
                #a=a.sel(Z=slice(-200, -700)).mean(dim='Z')
            b=xr.concat([b,a], dim='time')
        if save==True:
            if bottom==True:
                b.to_netcdf('./data/bottom100m_averaged_'+var+'_ens'+str(i+1)+'.nc')
            else:
                b.to_netcdf('./data/depth_averaged_'+var+'_ens'+str(i+1)+'.nc')
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
        
    for j,i in enumerate(ens_list):
        print('Loading number:'+str(i))
        a=xr.open_dataset('./data/'+filename+'_ens'+str(i+1)+'.nc')
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
                           vmin_pct=0.01, vmax_pct=0.99,
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

# def findEvents(total, units, longnames, threshold=50, wvar='dotson_to_cosgrove_massloss', window=24,
#                larger_than=True, method='mean', longwindow=12*25, min_periods=5*12,
#                start='1920', end='2013', filename=''):
    
#     #Create a long period running mean for detrending. We are interested in decadal variability.
#     long=total[wvar][start:end].rolling(longwindow, min_periods=min_periods, center=True).mean()
    
    
#     #threshold=np.nanpercentile(total[wvar][start:end].rolling(window, center=True).mean()-long, pct)
#     #windows=[window]
    
#     #Creat a dataframe with the selected time periods.
#     selected=pd.DataFrame() 
    
#     #Create arrays for the time axis
#     a=[pd.to_datetime(str(j)) for j in np.arange(1920, 2011, 20)] 
#     b=[str(j) for j in np.arange(1920, 2011, 20)]
#     xlim=[1920,2020]
    
#     lim=total[wvar][start:end].rolling(window, center=True).mean().std().mean()
    
#     #Create figure
#     plt.figure(figsize=(20,10))
#     for i in range(len(total[wvar].columns)):    
#         ax=plt.subplot(5,4,i+1)
        
        
#         col=total[wvar].columns.sort_values()[i]
#         series_dt=total[wvar][col][start:end].rolling(window, center=True).mean()-total[wvar][col][start:end].rolling(longwindow, min_periods=min_periods, center=True).mean()
#         series_dt.plot()
        
#         if larger_than==True:
#             sel=series_dt>threshold
#         elif larger_than==False:
#             sel=series_dt<threshold
        
#         selected[col]=sel
#         for j in range(len(sel)):
#             if (sel.iloc[j]==True) & (sel.iloc[j-1]==False):
#                 c=0
#                 while sel.iloc[j+c]==True:
#                     c+=1       
#                 plt.axvspan(sel.index[j], sel.index[j+c-1], color='red', alpha=0.3)
#         plt.hlines(0, pd.to_datetime('1910'), pd.to_datetime('2020'), colors='k', linestyles='dashed')
#         plt.grid()
#         plt.xlim([pd.to_datetime(str(k)) for k in xlim])
        
#         plt.title(col)
#         plt.ylabel(units[wvar])    
#         plt.ylim([-lim*2.5, lim*2.5])
#         plt.xticks(a,b)
#         ax.grid(True)


#     plt.subplots_adjust(hspace=0.4, wspace=0.4)

#     plt.suptitle(longnames[wvar]+';\n {} month rolling mean (centered)'.format(window)+'\n Threshold: {} in: '.format(str(threshold))+longnames[wvar])
    
#     from datetime import date
#     today = date.today()
#     today=today.strftime("%Y%m%d")
    
#     plt.savefig('./figures/selected_periods_'+filename+'_'+today+'.png')
    
#     return selected, units



#%%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def makeLatSliceAnimation(data, var, months=[-25,5], title='Composite map or not',
                           filename='map', unit='[m/s]', wvar='dotson_to_cosgrove_massloss', pct=90, larger_than=True, x=106):
    
    #From my code to Kaitlin's code
    from mitgcm_python_master.plot_latlon import latlon_plot
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    from mitgcm_python_master.utils import real_dir, daily_to_monthly, fix_lon_range, split_longitude, mask_land_ice, moving_average, index_year_start, index_year_end, index_period, mask_2d_to_3d, days_per_month, add_time_dim, z_to_xyz, select_bottom, convert_ismr, mask_except_ice, xy_to_xyz, apply_mask, var_min_max, mask_3d, average_12_months, depth_of_isoline, mask_land, axis_edges
    from mitgcm_python_master.plot_utils.latlon import shade_mask, overlay_vectors
    from mitgcm_python_master.plot_utils.windows import finished_plot, set_panels
    from matplotlib.animation import FuncAnimation
    import matplotlib.animation as animation
    
    fig=plt.figure(figsize=(15,5))
    ax=plt.subplot(1,1,1)
    moi=months[0]
    divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=np.nanpercentile(data[var], 1), vmax=np.nanpercentile(data[var], 99))
    
    a=data[var].sel(time=moi).mean(dim='events')
    a=a.where(a != 0) 
    
    img=ax.pcolormesh( a.YC, a.Z, a, cmap='seismic', norm=divnorm)
    plt.colorbar(img, label=unit)
    ax.set_facecolor('grey')
    ax.set_xlim([-76, -70])
    ax.set_ylim([-1500, 0])
    #ax.set_title('Absolute Values', fontsize=15)
    ax.set_xlabel('Latitude [deg]', fontsize=14)
    ax.set_ylabel('Depth [m]', fontsize=14)

    title_new=title+'\n (time since time of detection: {} months)'.format(str(moi))
    ax.set_title(title_new)
    #if larger_than==True:
    #    title_new=title+'\n (time since time of detection: {} months, threshold: above {} percentile of {})'.format(str(moi), str(pct), wvar)
    #else:
    #    title_new=title+'\n (time since time of detection: {} months, threshold: below {} percentile of {})'.format(str(moi), str(pct), wvar)


    # animation function
    def animate(i):
        print(i)
        moi=i+months[0]


        a=data[var].sel(time=moi).mean(dim='events')
        a=a.where(a != 0) 

        img=ax.pcolormesh( a.YC, a.Z, a, cmap='seismic', norm=divnorm)
        #plt.colorbar(img, label=data[var].units)
        ax.set_facecolor('grey')
        ax.set_xlim([-76, -70])
        ax.set_ylim([-1500, 0])
        #ax.set_title('Absolute Values', fontsize=15)
        ax.set_xlabel('Latitude [deg]', fontsize=14)
        ax.set_ylabel('Depth [m]', fontsize=14)
        
        title_new=title+'(time since time of detection: {} months)'.format(str(moi))
        ax.set_title(title_new)
        
        return 


    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")

    anim = animation.FuncAnimation(fig, animate, frames=months[-1]-months[0], repeat=False, interval=5000)
    writergif = animation.PillowWriter(fps=2)
    anim.save('./figures/'+today+'_anim_cross'+str(x)+'N_'+var+'_'+filename+'.gif', writer=writergif)
    
    return






def createLatAnimation(data, selected, var='THETA', title='Composite map', 
                       filename='map', window=24, method='mean', 
                       longwindow=12*25, 
                       unit='[m/s]',
                       members='all', larger_than=True, pct=90, startyr='1920', endyr='2013', months=[-25, 10], x=106): 
    
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
        

    if 'events' not in data_new.dims:
        data_new=data_new.expand_dims("events") 
            
    #Completing the title
    if members=='all':
        title+=' (full ensemble)'
    else:
        title+=' ({} members)'.format(str(members[-1]-members[0]))
        
    title=title+'\n ({} month moving average'.format(str(window))
    if method=='linear':
        title+=', linearly detrended)'
    elif method=='quadratic':
        title+=', quadratically detrended)'
    elif method=='mean':
        title+=',  detrended with {} month rolling mean)'.format(str(longwindow))
    elif method==None:
        title+=')'
    
    
    print('Starting with plotting')
    
    makeLatSliceAnimation(data=data_new, var=var, months=[months[0],months[-1]], title=title, filename=filename, unit=unit, 
                           pct=pct, larger_than=larger_than, x=x)
    
    print('Finished!')
    return


def readMITgcmData_Slice(var, bd='/data/oceans_output/shelf/kaight/mitgcm', members='all', method='lon', x=106, y=-72.5):
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
        if i==0:
            data=b
        else:
            data=xr.concat([data,b], dim='ens')
            
    #Making sure that the dataset has the dimension 'ens'
    if number==1:
        data=data.expand_dims("ens")
            
    return data

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
    
    plt.savefig('./figures/defined_periods_'+filename+'_'+today+'.png')
    
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
    
    plt.savefig('./figures/selected_periods_'+filename+'_'+today+'.png')
    
    return selected, units


# def plotLonSlice(data, var, filename='', x=106):
#     fig=plt.figure(figsize=(20, 15))

#     title='Cross Sections of '+data_xr['THETA'].long_name + 'at {} degrees W'.format(str(x))
#     fig.suptitle(title, fontsize=16, x=0.45)

#     a=data[var].mean(dim='time').mean(dim='ens')
#     a=a.where(a != 0) 

#     ax=plt.subplot(2,1,1)
#     divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=np.nanpercentile(a, 1), vmax=np.nanpercentile(a, 99))
#     img=ax.pcolormesh( a.YC, a.Z, a, cmap='seismic', norm=divnorm)
#     plt.colorbar(img, label=data[var].units)
#     ax.set_facecolor('grey')
#     ax.set_xlim([-76, -70])
#     ax.set_ylim([-1500, 0])
#     ax.set_title('Absolute Values', fontsize=15)
#     ax.set_xlabel('Latitude [deg]', fontsize=14)
#     ax.set_ylabel('Depth [m]', fontsize=14)

#     from unixFunctionsMITgcm_final import detrend_and_average_MITgcmData

#     data_dt=detrend_and_average_MITgcmData(data, var, window=1)

#     a=data_dt.sel(time=slice('1925', '1930')).mean(dim='time').mean(dim='ens')
#     a=a.where(a != 0) 

#     ax=plt.subplot(2,1,2)
#     ax.set_xlim([-76, -70])
#     ax.set_ylim([-1500, 0])

#     divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=np.nanpercentile(a, 1), vmax=np.nanpercentile(a, 99))
#     img=ax.pcolormesh( a.YC, a.Z, a, cmap='seismic', norm=divnorm)
#     plt.colorbar(img, label=data[var].units)
#     ax.set_facecolor('grey')
#     ax.set_title('Anomaly', fontsize=15)
#     ax.set_xlabel('Latitude [deg]', fontsize=14)
#     ax.set_ylabel('Depth [m]', fontsize=14)

#     from datetime import date
#     today = date.today()
#     today=today.strftime("%Y%m%d")

#     fig.savefig('./figures/'+today+'_cross'+str(x)+'W_'+var+'_'+filename+'.png')
    
#     return

# def plotLatSlice(data, var, filename='', y=-72.5):
#     fig=plt.figure(figsize=(20, 15))

#     title='Cross Sections of '+data_xr['THETA'].long_name + 'at {} degrees N'.format(str(y))
#     fig.suptitle(title, fontsize=16, x=0.45)

#     a=data[var].mean(dim='time').mean(dim='ens')
#     a=a.where(a != 0) 

#     ax=plt.subplot(2,1,1)
#     divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=np.nanpercentile(a, 1), vmax=np.nanpercentile(a, 99))
#     img=ax.pcolormesh( a.XC, a.Z, a, cmap='seismic', norm=divnorm)
#     plt.colorbar(img, label=data[var].units)
#     ax.set_facecolor('grey')
#     ax.set_xlim([360-120, 360-100])
#     ax.set_ylim([-800, 0])
#     ax.set_title('Absolute Values', fontsize=15)
#     ax.set_xlabel('Longitude [deg]', fontsize=14)
#     ax.set_ylabel('Depth [m]', fontsize=14)

#     from unixFunctionsMITgcm_final import detrend_and_average_MITgcmData

#     data_dt=detrend_and_average_MITgcmData(data, var, window=1)

#     a=data_dt.sel(time=slice('1925', '1930')).mean(dim='time').mean(dim='ens')
#     a=a.where(a != 0) 

#     ax=plt.subplot(2,1,2)
#     ax.set_xlim([360-120, 360-100])
#     ax.set_ylim([-800, 0])

#     divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=np.nanpercentile(a, 1), vmax=np.nanpercentile(a, 99))
#     img=ax.pcolormesh( a.XC, a.Z, a, cmap='seismic', norm=divnorm)
#     plt.colorbar(img, label=data[var].units)
#     ax.set_facecolor('grey')
#     ax.set_title('Anomaly', fontsize=15)
#     ax.set_xlabel('Longitude [deg]', fontsize=14)
#     ax.set_ylabel('Depth [m]', fontsize=14)

#     from datetime import date
#     today = date.today()
#     today=today.strftime("%Y%m%d")

#     fig.savefig('./figures/'+today+'_cross'+str(y)+'N_'+var+'_'+filename+'.png')
    
#     return



#%-------------------------------------------------------------
#PCA

def performPCA(data, var, method='all', interest=10, filepath='./figures/', filename=''):
    import sklearn.decomposition as pca
    from scipy import signal
    from numpy import ma
    
    data=data.drop_isel(time=np.where(np.isnan(data.values[0,:,150,150]))[0])
    
    if method=='ens_mean':
        data=data.mean(dim='ens')
        tc=data.indexes['time'].to_datetimeindex()
    elif method=='all':
        data['time']=data.indexes['time'].to_datetimeindex()
        tc=data.time
        dims=[len(data.ens), len(data.time)]
        data=data.stack(ens_time=['ens', 'time'])
        data=data.transpose("ens_time", "YC", "XC")
        
        
    sla=data.values
    lat=data.YC
    lon=data.XC
    
    print(np.shape(sla))
    
    #---------------------------------------------------------------------------------------------------------------------------
    #PERFORM PCA

    #Reshape the grid to one column with position and row corresponding to time
    X = np.reshape(sla, (sla.shape[0], len(lat) * len(lon)), order='F')
    #Create Boolean array corresponding to ocean indices and mask the array where there is land
    X = ma.masked_array(X, np.isnan(X))
    land = X.sum(0).mask
    X[:,land]=0
    total_modes = 10
    EOF_recons = np.ones((total_modes, len(lat) * len(lon))) * -999.

    #Perform PCA and take the EOFs and corersponding Pcs
    skpca = pca.PCA(total_modes)
    skpca.fit(X)
    PCs = skpca.transform(X)
    EOFs = skpca.components_
    #Reconstruct the EOF fields into gridded data EOF_recons =
    np.ones((total_modes, len(lat) * len(lon))) * -999.
    for i in range(total_modes):
        EOF_recons[i,:] = EOFs[i,:]
    EOF_recons = ma.masked_values(np.reshape(EOF_recons, (total_modes, len(lat), len(lon)), order='F'), -999.)
    print(EOF_recons.shape)
    
    
    #---------------------------------------------------------------------------------------------------------------------------
    #VISUALIZE
    
    ### plot loadings and scores
    from matplotlib.colors import Normalize
    import matplotlib.colors as colors

    

    
    for j in range(0,interest):
        #norm = MidpointNormalize(midpoint=0)
        #Plot principal components
        fig, ax = plt.subplots(2,1,figsize=(12,12))
        EOF_recons[j,:,:][EOF_recons[j,:,:]==0]=np.nan
        #im = ax[0].imshow(np.flip(EOF_recons[j,:,:]*max(abs(PCs[j,:])), axis=0), norm=divnorm,
        #cmap=plt.cm.RdBu_r, interpolation='none',extent=[lon[0],lon[-
        #1],lat[0],lat[-1]])
        divnorm = colors.TwoSlopeNorm(vmin=-np.nanquantile(abs(EOF_recons[j,:,:]*max(abs(PCs[j,:]))), 0.99), vmax=np.nanquantile(abs(EOF_recons[j,:,:]*max(abs(PCs[j,:]))), 0.99), vcenter=0)
        im = ax[0].pcolor(lon, lat, EOF_recons[j,:,:]*max(abs(PCs[j,:])), norm=divnorm,
        cmap=plt.cm.RdBu_r)

        ax[0].set_facecolor('grey')
        cbar=plt.colorbar(im, fraction=0.016, pad=0.04, ax=ax[0])
        cbar.set_label('(m)')
        ax[0].set_xlabel(r'Longitude ($^\circ$E)')
        ax[0].set_ylabel(r'Latitude ($^\circ$N)')

        eof=j+1
        ax[0].set_title(r'Empirical Orthogonal Function %i' % eof)
        PCnor=PCs[:,j]/max(abs(PCs[:,j]))
        
        if method=='all':
            PCnor=PCnor.reshape(dims[0], dims[1])
            for k in range(dims[0]):
                if k==0:
                    ax[1].plot(tc,PCnor[k,:], c='grey', linewidth=0.2, label='Individual member PC amplitude')
                else:
                     ax[1].plot(tc,PCnor[k,:], c='grey', linewidth=0.2)
            ax[1].plot(tc,np.mean(PCnor, axis=0), c='black', linewidth=1.2, label='Ensemble averaged PC amplitude')
        else:
            ax[1].plot(tc, PCnor, label='PC')
        #ax[1].plot(tc, (mean-mean.mean())/max(abs(mean-mean.mean())), label='Shelf Mean: corr:{}'.format(str(np.corrcoef(mean, PCnor)[0,1])))
        ax[1].grid()
        ax[1].set_title('Principal Component %i' % eof)
        ax[1].set_xlabel('Time (year)')
        ax[1].set_ylabel('Magnitude')
        #ax[1].set_xlim([tc[0],tc[-1]])
        ax[1].set_ylim([-1,1])
        ax[1].legend()
        
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        
        plt.savefig(filepath+today+'_'+var+"_PCA_component_"+str(eof)+'_'+method+'_'+filename+".png", bbox_inches='tight')
        plt.show()
    
    variance=skpca.explained_variance_ratio_
    variancecumsum=np.cumsum(variance)
    eofnumber=np.arange(1,len(variance)+1)

    f, ax = plt.subplots(figsize=(6,4))
    ax.plot(eofnumber,variancecumsum)
    ax.plot(eofnumber[0:5],variancecumsum[0:5],'ro',markersize=4)
    ax.set_xlabel('EOF modes')
    ax.set_ylabel('Cumulative fraction')
    ax.set_title("Cumulative fraction of explained variance",
    fontsize=14)
    ax.grid()
    #ax.set_xlim([-1,313])
    #ax.set_xticks([1,50,100,150,200,250,300])
    plt.savefig(filepath+today+'_'+var+"cumulativefraction"+method+'_'+filename+".png", bbox_inches='tight')
    plt.show()
    print("Variance contained in first 10 EOFs: % .2f" %
    variancecumsum[9])
    
    np.save('./data/'+today+'_'+var+"_PCA_component"+method+'_'+filename+".npy", PCs)
    
    return PCs, tc


def plotPCACorrs(PCs, total, real_var, data_dt, var_list=['amundsen_shelf_break_uwind_avg', 'amundsen_shelf_temp_btw_200_700m', 'dotson_to_cosgrove_massloss'], interest=10, filename='', filepath='./figures/'):
    
    #print(cors)
    plt.figure(figsize=(15,10))
    
    for i, var in enumerate(var_list):
        total[var] = total[var].reindex(sorted(total[var].columns), axis=1)
        ts=total[var]['1920':'2013']
        ts=ts.rolling(24,center=True).mean()-ts.rolling(25*12, min_periods=5*12, center=True).mean()
        ts=ts.iloc[12:-11]
        mean=ts.to_numpy().T.flatten()
        
        
        #print(len(mean))
        cors=np.zeros(interest)
        for j in range(0,interest):
            PCnor=PCs[:,j]/max(abs(PCs[:,j]))
            cors[j]=np.corrcoef(mean, PCnor)[0,1]
        
        plt.subplot(len(var_list),1,i+1)
        plt.plot(range(1,interest+1), cors, 'o-', label='Pearson r')
        plt.plot(range(1,interest+1), abs(cors), 'o-', label='Absolute values')
        #plt.yticks([-1,0,1])
        plt.grid()
        plt.legend()
        plt.ylabel('Pearson R')
        plt.xlabel('Principal Component #')
        plt.title('Correlation between '+var+' and PCs from '+real_var)
    
    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")
    
    plt.subplots_adjust(hspace=0.5)
    plt.savefig(filepath+today+"_PCAcorrs_"+real_var+'_'+filename+".png", bbox_inches='tight')
    
    return



#%%-------------------------OVERVIEW_PLOT-----------------------------------------------------------------------------------


def plotLonSlice(data, var, ax, year, unit='[K]', filename='', x=106):
    import matplotlib.colors as colors 
    
    a=data.sel(time=year).mean(dim='time').mean(dim='ens')
    a=a.where(a != 0) 

    ax.set_xlim([-76, -70])
    ax.set_ylim([-1500, 0])

    divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=-0.8, vmax=0.8)
    
    img=ax.pcolormesh(a.YC, a.Z, a, cmap='seismic', norm=divnorm)
    plt.colorbar(img, label=unit, ax=ax)
    ax.set_facecolor('grey')
    ax.set_title('Cross Section at {} deg W'.format(str(x)), fontsize=15)
    ax.set_xlabel('Latitude [deg]', fontsize=14)
    ax.set_ylabel('Depth [m]', fontsize=14)

    return

def plotLatSlice(data, var, ax, year, unit='[K]', filename='', y=-72.5):
    import matplotlib.colors as colors 
    
    a=data.sel(time=year).mean(dim='time').mean(dim='ens')
    a=a.where(a != 0) 

    ax.set_xlim([360-120, 360-100])
    ax.set_ylim([-800, 0])

    divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=-0.8, vmax=0.8)
    
    img=ax.pcolormesh(a.XC, a.Z, a, cmap='seismic', norm=divnorm)
    plt.colorbar(img, label=unit, ax=ax)
    ax.set_facecolor('grey')
    ax.set_title('Cross Section at {} deg N'.format(str(y)), fontsize=15)
    ax.set_xlabel('Longitude [deg]', fontsize=14)
    ax.set_ylabel('Depth [m]', fontsize=14)

    return



def createOverviewAnimation(member=6, start_year=1935, months=12*10, window=24, deseasonalized=False, filename=''):
    
    #Load Data that we already have
    data={}
    data_dt={}
    members=[member-1, member]

#     var_list=['EXFuwind', 'THETA', 'UVEL', 'VVEL', 'THETAW106', 'THETAW120']#, 'UVELB', 'VVELB']
#     filenames={'EXFuwind':'zonal_winds', 
#                'THETA':'new2_amundsen_shelf_temp_200to700', 
#                'UVEL':'depth_averaged_UVEL', 
#                'VVEL':'depth_averaged_VVEL',
#                'curlTAU':'curl_tau', 
#                'THETAW106':'theta_W106', 
#                'THETAW120':'theta_W120', 
#                'UVELB':'bottom100m_averaged_UVEL', 
#                'VVELB':'bottom100m_averaged_VVEL'}

#     for var in var_list: 
#         print('Starting with: '+var)
#         data[var]=loadMITgcmData(filename=filenames[var], members=members)
#         if 'VELB' in var:
#             if 'U' in var:
#                 data[var]=data[var]-data['UVEL']
#                 data_dt[var]=detrend_and_average_MITgcmData(data[var], 'UVEL', window=window)
#             if 'V' in var:
#                 data[var]=data[var]-data['VVEL']
#                 data_dt[var]=detrend_and_average_MITgcmData(data[var], 'VVEL', window=window)
#         elif 'THETA' in var:
#             data_dt[var]=detrend_and_average_MITgcmData(data[var], 'THETA', window=window)
#         else:
#             data_dt[var]=detrend_and_average_MITgcmData(data[var], var, window=window)
        
#     #Load Data that we already have
#     data={}
#     data_dt={}
    
    #Reading the data.
    full='/data/hpcdata/users/grejan/mitgcm/data'
    total={}
    units={}
    longnames={}
#     total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='isotherm')
#     total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='final')
    total, units, longnames=read_all_data(full, total, units, longnames, kind='ohc')
    
    
    var_list=['EXFuwind', 'THETA', 'UVEL', 'VVEL', 'THETAW106', 'THETAS72-5', 'UVELB', 'VVELB']
    filenames={'EXFuwind':'zonal_winds', 
               'THETA':'new2_amundsen_shelf_temp_200to700', 
               'UVEL':'depth_averaged_UVEL', 
               'VVEL':'depth_averaged_VVEL',
               'curlTAU':'curl_tau', 
               'THETAW106':'theta_W106', 
               'THETAW120':'theta_W120',
               'THETAW103-5':'theta_W103-5',
               'THETAW100':'theta_W100',
               'THETAS72-5':'theta_S72-5',
               'UVELB':'bottom100m_averaged_UVEL', 
               'VVELB':'bottom100m_averaged_VVEL'}

    for var in var_list: 
        print('Starting with: '+var)
        data[var]=fc.loadMITgcmData(filename=filenames[var], members=members)
        if 'VELB' in var:
            if 'U' in var:
                data[var]=data[var]-data['UVEL']
                data_dt[var]=fc.detrend_and_average_MITgcmData(data[var], 'UVEL', window=window)
            elif 'VVEL' in var:
                data[var]=data[var]-data['VVEL']
                data_dt[var]=fc.detrend_and_average_MITgcmData(data[var], 'VVEL', window=window)
        elif 'THETA' in var:
            data_dt[var]=fc.detrend_and_average_MITgcmData(data[var], 'THETA', window=window)
        else:
            data_dt[var]=fc.detrend_and_average_MITgcmData(data[var], var, window=window)
        
        if deseasonalized==True:
            data_dt[var] = data_dt[var].groupby("time.month")-data_dt[var].groupby("time.month").mean(dim="time")
            

    def showOverview(i):
        import matplotlib.colors as colors
        plt.rcParams.update({'font.size': 14})

        #start_year=1935
        #import pandas as pd
        a=pd.to_datetime(str(start_year))+pd.DateOffset(months=i)
        year=a.strftime('%Y-%m')

        import matplotlib.gridspec as gridspec
        #plt.figure(figsize=(30,10))
        gs = gridspec.GridSpec(4, 4)
        gs.update(wspace=0.5)
        ax1 = plt.subplot(gs[0, :2])
        ax2 = plt.subplot(gs[0, 2:])
        ax3 = plt.subplot(gs[1, :2])
        ax4 = plt.subplot(gs[1, 2:])
        ax5 = plt.subplot(gs[2, 0])
        ax6 = plt.subplot(gs[2, 1])
        ax7 = plt.subplot(gs[2, 2])
        ax8 = plt.subplot(gs[2, 3])
        ax9 = plt.subplot(gs[3, 1:3])
        ax10 = plt.subplot(gs[3, 3])

        from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
        gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
        grid = Grid(gp)


        var='EXFuwind'

        xcor=[m for m in data_dt[var].coords if 'X' in m][0]
        ycor=[m for m in data_dt[var].coords if 'Y' in m][0]

        divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=np.nanpercentile(data_dt[var], 1), vmax=np.nanpercentile(data_dt[var], 99))

        data_dt[var]=data_dt[var].where(data_dt[var]!=0)
        img=ax1.pcolor(data_dt[var][xcor], data_dt[var][ycor], data_dt[var].sel(time=year).mean(dim='time').mean(dim='ens').values, cmap='seismic', norm=divnorm)
        ax1.set_facecolor('grey')
        bath=ax1.contour(data_dt[var][xcor], data_dt[var][ycor], grid.bathy, levels=[-1000], colors=['black'])
        plt.clabel(bath, [-1000])
        plt.colorbar(img, label='[m/s]', ax=ax1)
        ax1.set_title('Zonal Wind')

        x=[360-120, 360-120]
        y=[-76, -70]
        ax1.plot(x, y, 'o--g')

        x=[360-106, 360-106]
        y=[-76, -70]
        ax1.plot(x, y, 'o--g')

        x=[360-103.5, 360-103.5]
        y=[-76, -70]
        ax1.plot(x, y, 'o--g')

        x=[360-100, 360-100]
        y=[-76, -70]
        ax1.plot(x, y, 'o--g')

        x=[240, 260]
        y=[-72.5, -72.5]
        ax1.plot(x, y, 'o--c')


        var='THETA'
        divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=np.nanpercentile(data_dt[var], 1), vmax=np.nanpercentile(data_dt[var], 99))
        img=ax2.pcolor(data_dt[var][xcor], data_dt[var][ycor], data_dt[var].sel(time=year).mean(dim='time').mean(dim='ens').values, cmap='seismic', norm=divnorm)
        plt.colorbar(img, label='[m/s]', ax=ax2)
        ax2.set_facecolor('grey')
        bath=ax2.contour(data_dt[var][xcor], data_dt[var][ycor], grid.bathy, levels=[-1000], colors=['black'])
        plt.clabel(bath, [-1000])
        ax2.set_title('Average Temperature [200-700m]')
        ax2.set_ylim([-76, -70])


        var1='UVEL'
        var2='VVEL'
        vfield=np.sqrt(data_dt[var1].sel(time=year).mean(dim='time').mean(dim='ens').values**2
                       +data_dt[var2].sel(time=year).mean(dim='time').mean(dim='ens').values**2)

        img=ax3.pcolor(data_dt[var][xcor], data_dt[var][ycor], vfield, cmap='rainbow', vmax=0.02)
        plt.colorbar(img, label='[m/s]', ax=ax3)
        ax3.set_facecolor('grey')
        step=10
        ax3.quiver(data_dt[var][xcor][::step], data_dt[var][ycor][::step], 
                       data_dt['UVEL'].sel(time=year).mean(dim='time').mean(dim='ens').values[::step, ::step], 
                       data_dt['VVEL'].sel(time=year).mean(dim='time').mean(dim='ens').values[::step, ::step], 
                       color='black')
        bath=ax3.contour(data_dt[var][xcor], data_dt[var][ycor], grid.bathy, levels=[-1000], colors=['white'])
        plt.clabel(bath, [-1000])
        ax3.set_ylim([-76, -70])
        ax3.set_title('Depth Averaged Ocean Velocity')


        var1='UVELB'
        var2='VVELB'
        vfield=np.sqrt(data_dt[var1].sel(time=year).mean(dim='time').mean(dim='ens').values**2
                       +data_dt[var2].sel(time=year).mean(dim='time').mean(dim='ens').values**2)

        img=ax4.pcolor(data_dt[var][xcor], data_dt[var][ycor], vfield, cmap='rainbow', vmax=0.02)
        plt.colorbar(img, label='[m/s]', ax=ax4)
        ax4.set_facecolor('grey')
        step=10
        ax4.quiver(data_dt[var][xcor][::step], data_dt[var][ycor][::step], 
                       data_dt['UVELB'].sel(time=year).mean(dim='time').mean(dim='ens').values[::step, ::step], 
                       data_dt['VVELB'].sel(time=year).mean(dim='time').mean(dim='ens').values[::step, ::step], 
                       color='black')
        bath=ax4.contour(data_dt[var][xcor], data_dt[var][ycor], grid.bathy, levels=[-1000], colors=['white'])
        plt.clabel(bath, [-1000])
        ax4.set_ylim([-76, -70])
        ax4.set_title('Bottom 100m Ocean Velocity - Depth Averaged')





        plotLonSlice(data_dt['THETAW120'], var='THETA', ax=ax6, year=year, filename='', x=120)
        plotLonSlice(data_dt['THETAW106'], var='THETA', ax=ax6, year=year, filename='', x=106)
        plotLonSlice(data_dt['THETAW103-5'], var='THETA', ax=ax7, year=year, filename='', x=103.5)
        plotLonSlice(data_dt['THETAW100'], var='THETA', ax=ax8, year=year, filename='', x=100)

        plotLatSlice(data_dt['THETAS72-5'], var='THETA', ax=ax9, year=year, filename='', y=-72.5)

        plt.subplots_adjust(wspace=0.4, hspace=0.4)
        plt.suptitle('Averaged Anomalies in Member '+str(members[-1])+' over ' +year)

        #years=[1935, 1945]
        total['amundsen_shelf_ohc_below_0m']['ens%02d' % (members[-1],)].rolling(24).mean()[str(start_year):str(start_year+10)].plot(grid=True, ax=ax10, rot=90)
        ax10.set_ylabel('OHC')
        ax10.set_title('Ocean Heat Content in Event of Interest', pad=20)
        ax10.axvline(x=total['amundsen_shelf_ohc_below_0m'].loc[year].index[0], color='r')


        return

    import matplotlib.animation as animation

    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")

    fig=plt.figure(figsize=(30,10))
    anim = animation.FuncAnimation(fig, showOverview, frames=months, repeat=False, interval=5000)
    writergif = animation.PillowWriter(fps=2)
    anim.save('./figures/animation_overview_'+filename+'_'+today+'.gif', writer=writergif)
    
    return














#%%%%%%%%%%%%%%%%%%_-------------------------------------------------------------------------------------------------------------------------------------------

def createVelocityCompositeMap(dataU, dataV, selected, title='Composite map', 
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
                           pct=pct, larger_than=larger_than)
    
    print('Finished!')
    return

def makeVelocityCompositeAnimation(dataU, dataV, months=[-25,5], show_box=True, title='Composite map or not',
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
    anim.save('./figures/animation_composite_'+filename+'_'+today+'.gif', writer=writergif)
    return
