# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain test!.

@author: joren
"""

from functions.loading_and_processing_data import loadMITgcmData, loadPACEData, detrend_and_average_MITgcmData, read_all_data
from functions.composite_maps import createCompositeMap, createVelocityCompositeMap, defineEvents, findEvents, createCompositeMap_2, workflowCompositeMap
import numpy as np
import xarray as xr


#00-INITIALISATION------------------------------------------------------------

members='all'#[0,1]
months=[-25, 15]

# members=[0,2]
# months=[-3, 1]

deseasonalize=True
window=1
output_filename='new_warm'
suptitle='Composite of warm events'

#First, read time series so that we can define events.
full='/data/hpcdata/users/grejan/mitgcm/data'
total={}
units={}
longnames={}
total, units, longnames=read_all_data(full, total, units, longnames, kind='isotherm')
total, units, longnames=read_all_data(full, total, units, longnames, kind='final')
total, units, longnames=read_all_data(full, total, units, longnames, kind='ohc')

#Select events that we are interested in.
threshold=defineEvents(total=total, longnames=longnames, units=units, wvar='amundsen_shelf_ohc_below_0m', larger_than=True, pct=75, detrend=True, window=12*2, filename='warm75_colored_in')
print(threshold)
selected,_=findEvents(total=total, longnames=longnames, units=units, threshold=threshold, wvar='amundsen_shelf_ohc_below_0m', larger_than=True, detrend=True, window=12*2, filename='warm75_colored_in')


import pandas as pd
import cftime


startyr='1920'
endyr='2013'

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

if members!='all':
    columns=columns[members[0]:members[-1]]

months=np.arange(months[0], months[-1])

#----------------------------------------------------------
#02-Overview-Composites------------------------------------------------

#Load Data that we already have


var_list=['UVELW115', 'EXFuwind', 'THETA', 'THETAPITE', 'UVEL', 'UVELB', 'PSL']

filenames={'EXFuwind':'zonal_winds',
           'THETA':'THETA_depth_integrated',
           'UVEL':'depth_averaged_UVEL', 
           'VVEL':'depth_averaged_VVEL',
           'curlTAU':'curl_tau', 
           'THETAW106':'theta_W106', 
           'THETAW120':'theta_W120',
           'UVELW115':'UVEL_W115',
           'THETAW103-5':'theta_W103-5',
           'THETAW100':'theta_W100',
           'THETAS72-5':'theta_S72-5',
           'THETAPITE':['THETA_transect_PITE_full_0', 'THETA_transect_PITE_full_1', 
                        'THETA_transect_PITE_full_2', 'THETA_transect_PITE_full_3', 
                        'THETA_transect_PITE_full_4'],
           'UVELB':'bottom100m_averaged_UVEL', 
           'VVELB':'bottom100m_averaged_VVEL'}
#'THETA':'new2_amundsen_shelf_temp_200to700',

data_new={}
data={}
for var in var_list:
    data_dt={}
    kind='maps'
    print('Starting with: '+var)
    if var=='THETAPITE':
        kind='transects'
        for fn in filenames[var]:
            data[var+'_'+fn[-1]]=loadMITgcmData(filename=fn, members=members, kind=kind)
    elif var=='PSL':
        data[var]=loadPACEData(var=var, members=members)    
    else:
        if var in ['EXFuwind', 'UVEL', 'UVELB']:
            kind='old'
        elif ('W1' in var) or ('S7' in var):
            kind='slices'
        print('Continue')
        if var in ['UVEL', 'UVELB']:
            data[var]=loadMITgcmData(filename=filenames[var], members=members, kind=kind)
            vary=var.replace('U', 'V')
            data[vary]=loadMITgcmData(filename=filenames[vary], members=members, kind=kind)    
        else:
            data[var]=loadMITgcmData(filename=filenames[var], members=members, kind=kind)    
        
    print('Detrending: '+var)    
    
    if 'VELB' in var:
        data[var]=data[var]-data['UVEL']
        data_dt[var]=detrend_and_average_MITgcmData(data[var], 'UVEL', window=window)
        vary=var.replace('U', 'V')
        data[vary]=data[vary]-data['VVEL']
        data_dt[vary]=detrend_and_average_MITgcmData(data[vary], 'VVEL', window=window)
        del data['UVEL']
        del data['VVEL']
    elif ('W1' in var) or ('S7' in var):
        new = ''.join([i for i in var if not i.isdigit()])[:-1]
        data_dt[var]=detrend_and_average_MITgcmData(data[var], new, window=window)        
    elif 'THETAPITE' in var:
        for fn in filenames[var]:
            data_dt[var+'_'+fn[-1]]=detrend_and_average_MITgcmData(data[var+'_'+fn[-1]], 'THETA', window=window)
            #Introduce shitty absolute data
            print('Adding absolute '+var)
            data_dt[var+'_abs_'+fn[-1]]=data[var+'_'+fn[-1]]['THETA']
    else:
        data_dt[var]=detrend_and_average_MITgcmData(data[var], var, window=window)
    
    if var in ['UVEL', 'UVELB']:
        vary=var.replace('U', 'V')
        varxy=[var, vary]
        x=data[varxy[0]]['UVEL'].values
        y=data[varxy[1]]['VVEL'].values
        vfield=np.sqrt(x**2+y**2)
        vel=var.replace('U','')
        data[vel+'field'] = xr.DataArray(vfield, coords=[data[varxy[0]].ens, data[varxy[0]].time, data[varxy[0]].YC, data[varxy[1]].XC], dims=[ "ens", "time", "YC", "XC"])
        data[vel+'field']=data[vel+'field'].to_dataset(name=vel+'field')
        data_dt[vel+'field']=detrend_and_average_MITgcmData(data[vel+'field'], vel+'field', window=window)
        
    if var=='UVEL':
        vary=var.replace('U', 'V')
        data_dt[vary]=detrend_and_average_MITgcmData(data[vary], vary, window=window)
    else:
        if 'THETAPITE' in var:
            for fn in filenames[var]:
                del data[var+'_'+fn[-1]]
        else:
            del data[var]
    
    if deseasonalize==True:
        #print('Deseasonalize it')
        #for var in list(data_dt):
        print('Deseasonalize: '+var)
        if 'THETAPITE' in var:
            for fn in filenames[var]:
                data_dt[var+'_abs_'+fn[-1]] = data_dt[var+'_abs_'+fn[-1]].groupby("time.month")-data_dt[var+'_'+fn[-1]].groupby("time.month").mean(dim="time")
                data_dt[var+'_'+fn[-1]] = data_dt[var+'_'+fn[-1]].groupby("time.month")-data_dt[var+'_'+fn[-1]].groupby("time.month").mean(dim="time")
        else:
            data_dt[var] = data_dt[var].groupby("time.month")-data_dt[var].groupby("time.month").mean(dim="time")
            if ('VEL' in var) & ('W1' not in var):
                vary=var.replace('U', 'V')
                vel=var.replace('U', '')+'field'
                data_dt[vary] = data_dt[vary].groupby("time.month")-data_dt[vary].groupby("time.month").mean(dim="time")
                data_dt[vel] = data_dt[vel].groupby("time.month")-data_dt[vel].groupby("time.month").mean(dim="time")
        
    #-------------------------------------------------------------------------------------
    print('Start selected months with with: '+var)
    
    def groupByMonth(data_new, var, data_dt):
        data_new[var]=None
        try:
            data_dt[var]['time']=data_dt[var].indexes['time'].to_datetimeindex()
        except:
            print(var+': already datetime index')
        for month in months:
            b=None
            print('Month is:'+str(month))
            for ensnum, col in enumerate(columns):
                print(col)
                print(ensnum)
                true_list=start_list[col][start_list[col]==True]
                tlist=true_list.index+pd.DateOffset(months=month)
                tlist=[t for t in tlist if (t>pd.to_datetime(startyr)) & (t<pd.to_datetime(endyr))] #Check that all values are within the correct time window
                a=data_dt[var].sel(ens=ensnum).sel(time=tlist)
                a=a.rename({'time':'events'})
                a['events']=np.arange(len(a.events.values))+ensnum*100
                a=a.expand_dims("time")
                a['time']=[month]
                if type(b)==type(None):
                    b=a
                else:
                    b=xr.concat([b,a], dim='events')
            if type(data_new[var])==type(None):
                data_new[var]=b
            else:
                data_new[var]=xr.concat([data_new[var], b], dim='time')
        data_new[var]=data_new[var].mean(dim='events')
        return data_new
    
    if 'THETAPITE' in var:
        for fn in filenames[var]:
            data_new=groupByMonth(data_new, var+'_'+fn[-1], data_dt)
            data_new=groupByMonth(data_new, var+'_abs_'+fn[-1], data_dt)
    else:
        data_new=groupByMonth(data_new, var, data_dt)
        if ('VEL' in var) & ('W1' not in var):
            vary=var.replace('U', 'V')
            vel=var.replace('U', '')+'field'
            data_new=groupByMonth(data_new, vary, data_dt)
            data_new=groupByMonth(data_new, vel, data_dt)
    del data_dt

total_dt=total['amundsen_shelf_ohc_below_0m']-total['amundsen_shelf_ohc_below_0m'].rolling(25*12, center=True, min_periods=5*12).mean()
events=pd.DataFrame()
for col in columns:
    for i in start_list[col][start_list[col]==True].index:
        events[col+'_'+str(i.year)+'-'+str(i.month)]=total_dt[col].rolling(24, center=True).mean()[i+pd.DateOffset(months=months[0]):i+pd.DateOffset(months=months[-1])].values


#ANIMATION
def OverviewAnimation_new(data_dt, suptitle='', filename='test'):
    from matplotlib.animation import FuncAnimation
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

    plt.rcParams.update({'font.size': 14})
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)

    year=months[0]

    def plotSimpleMap(var, ax, year, title='Zonal Wind', unit='[m/s]', ymax=None, cbar=False):
        xcor=[m for m in data_dt[var].coords if 'X' in m][0]
        ycor=[m for m in data_dt[var].coords if 'Y' in m][0]
        data_dt[var]=data_dt[var].where(data_dt[var]!=0)
        divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=np.nanpercentile(data_dt[var], 1), vmax=np.nanpercentile(data_dt[var], 99))
        img=ax.pcolor(data_dt[var][xcor], data_dt[var][ycor], data_dt[var].sel(time=year).values, cmap='seismic', norm=divnorm)
        ax.set_facecolor('grey')
        bath=ax.contour(data_dt[var][xcor], data_dt[var][ycor], grid.bathy, levels=[-1000], colors=['black'])
        plt.clabel(bath, [-1000])
        if type(ymax)!=None:
            ax.set_ylim([-75.5, ymax])
        if cbar:
            plt.colorbar(img, label=unit, ax=ax)
        ax.set_title(title)
        return

    def plotGlobalMap(var,  year, ax_nr=331, title='Sea Level Pressure', unit='[Pa]', cbar=False):
        ax = plt.subplot(ax_nr, projection=ccrs.PlateCarree(central_longitude=180.0))
        data_dt[var]=data_dt[var].where(data_dt[var]!=0)
        divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=np.nanpercentile(data_dt[var], 1), vmax=np.nanpercentile(data_dt[var], 99))
        img=ax.pcolor(data_dt[var].lon, data_dt[var].lat, 
                      data_dt[var].sel(time=year).values, 
                      cmap='seismic', norm=divnorm, transform=ccrs.PlateCarree())
        ax.set_facecolor('grey')
        if cbar:
            plt.colorbar(img, label=unit, ax=ax, pad=0.15)
        ax.coastlines()
        ax.set_title(title)
        ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(central_longitude=180.0))

    def plotLonSlice(var, ax,  year, unit='[K]', title='Zonal Velocity', x=106, 
                     xlim=[-71.8, -71.1], ylim=[-2500, 0], cbar=False):
        import matplotlib.colors as colors 

        a=data_dt[var].sel(time=year)
        a=a.where(a != 0) 

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=np.nanpercentile(data_dt[var], 1), vmax=np.nanpercentile(data_dt[var], 99))
        img=ax.pcolormesh(a.YC, a.Z, a, cmap='seismic', norm=divnorm)
        if cbar:
            plt.colorbar(img, label=unit, ax=ax)
        ax.set_facecolor('grey')
        ax.set_title(title+' at {} deg W'.format(str(x)), fontsize=15)
        ax.set_xlabel('Latitude [deg]', fontsize=14)
        ax.set_ylabel('Depth [m]', fontsize=14)

        return

    def plotTransect(var, ax,  year, unit='[K]', title='Theta through PITE trough', cbar=False):
        if 'abs' in var:
            sublist=[k for k in data_dt if var in k]
        else:
            sublist=[k for k in data_dt if (var in k) & ('abs' not in k)]
        for i, sub in enumerate(sublist):
            if i==0:
                dat=data_dt[sub].sel(time=year)
            else:
                dat1=data_dt[sub].sel(time=year)
                dat1['Dist']=dat1.Dist+np.max(dat.Dist)
                dat1['Z']=dat.Z
                dat=xr.concat([dat, dat1], dim='Dist')
        dat=dat.where(dat!=0)
        ax.set_ylim([-2000, 0])
        img=ax.pcolor(dat.Dist, dat.Z, dat, cmap='seismic')
        if cbar:
            plt.colorbar(img, label=unit, ax=ax)
        ax.set_facecolor('grey')
        ax.set_title(title)
        return

    def plotVelMap(var, varxy,  year, ax, title='Depth integrated velocity', unit='[m/s]', ymax=None, cbar=False):
        xcor=[m for m in data_dt[var].coords if 'X' in m][0]
        ycor=[m for m in data_dt[var].coords if 'Y' in m][0]
        data_dt[var]=data_dt[var].where(data_dt[var]!=0)
        divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=np.nanpercentile(data_dt[var], 1), vmax=np.nanpercentile(data_dt[var], 99))
        z=data_dt[var].sel(time=year).values.copy()
        z[np.isfinite(z)]=0
        img=ax.pcolor(data_dt[var][xcor], data_dt[var][ycor], z, cmap='seismic', norm=divnorm)
        ax.set_facecolor('grey')
        #bath=ax.contour(data_dt[var][xcor], data_dt[var][ycor], grid.bathy, levels=[-1000], colors=['black'])
        #plt.clabel(bath, [-1000])
        if type(ymax)!=None:
            ax.set_ylim([-75.5, ymax])

        ax.set_title(title)

        #     ax.quiver(data_dt[var][xcor][::step], data_dt[var][ycor][::step], 
        #                data_dt[varxy[0]].sel(time=year).values[::step, ::step], 
        #                data_dt[varxy[1]].sel(time=year).values[::step, ::step], edgecolor='white',width=0.02, headwidth=2, minshaft=5, headlength=5, linewidth = 0.5,
        #                color='black')
        step=4
        step2=10
        img=ax.quiver(data_dt[var][xcor][::step2], data_dt[var][ycor][::step], 
                   data_dt[varxy[0]].sel(time=year).values[::step, ::step2], 
                   data_dt[varxy[1]].sel(time=year).values[::step, ::step2], data_dt[var].sel(time=year).values[::step, ::step2], scale=0.2, width=0.01, edgecolor='black', linewidth = 0.5, cmap='seismic', norm=divnorm)
        if cbar:
            plt.colorbar(img, label=unit, ax=ax)
        return


    fig, ax=plt.subplots(3,3,figsize=(20,10))
    plotGlobalMap(var='PSL', ax_nr=331, year=year, title='Sea Level Pressure', unit='[Pa]', cbar=True)
    plotSimpleMap(var='EXFuwind', ax=ax[0,1], year=year, title='Zonal Wind', unit='[m/s]', cbar=True)
    plotLonSlice(var='UVELW115', ax=ax[0,2], year=year, unit='[m/s]', title='Zonal Velocity', x=115, cbar=True)
    plotSimpleMap(var='THETA', ax=ax[1,0], year=year, title='Depth Integrated THETA', unit='[K]', ymax=-70, cbar=True)

    plotVelMap(var='VELfield', varxy=['UVEL', 'VVEL'], year=year, ax=ax[1,1], title='Depth averaged velocity', 
               unit='[m/s]', ymax=-70, cbar=True)
    plotVelMap(var='VELBfield', varxy=['UVELB', 'VVELB'], year=year, ax=ax[1,2], 
               title='Bottom 100m velocity -depth averaged', unit='[m/s]', ymax=-70, cbar=True)

    plotTransect(var='THETAPITE_abs', ax=ax[2,0], year=year, unit='[K]', title='Abs Theta through PITE trough', cbar=True)
    plotTransect(var='THETAPITE', year=year, ax=ax[2,1], unit='[K]', title='Delta Theta through PITE trough', cbar=True)

    for ev in list(events):
        ax[2,2].plot(months, events[ev] ,linewidth=0.1, c='grey')
    ax[2,2].plot(months, events.mean(axis=1) ,linewidth=1.2, c='k')
    ax[2,2].fill_between(months, events.mean(axis=1)-events.std(axis=1), events.mean(axis=1)+events.std(axis=1), color='b', alpha=0.2)
    ax[2,2].set_ylabel('OHC')
    ax[2,2].set_title('Ocean Heat Content during Events', pad=20)
    ax[2,2].set_xlabel('Months from event')
    ax[2,2].axvline(x=year, color='r')

    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    fig.suptitle(suptitle+'; Month wrt event: '+str(year))
    
    def animate(i):
        print(i)
        year=i+months[0]
        #fig, ax=plt.subplots(3,3,figsize=(20,10))
        plotGlobalMap(var='PSL', ax_nr=331, year=year, title='Sea Level Pressure', unit='[Pa]', cbar=True)
        plotSimpleMap(var='EXFuwind', ax=ax[0,1], year=year, title='Zonal Wind', unit='[m/s]')
        plotLonSlice(var='UVELW115', ax=ax[0,2], year=year, unit='[m/s]', title='Zonal Velocity', x=115)
        plotSimpleMap(var='THETA', ax=ax[1,0], year=year, title='Depth Integrated THETA', unit='[K]', ymax=-70)

        plotVelMap(var='VELfield', varxy=['UVEL', 'VVEL'], year=year, ax=ax[1,1], title='Depth averaged velocity', unit='[m/s]', ymax=-70)
        plotVelMap(var='VELBfield', varxy=['UVELB', 'VVELB'], year=year, ax=ax[1,2], title='Bottom 100m velocity -depth averaged', unit='[m/s]', ymax=-70)

        plotTransect(var='THETAPITE_abs', ax=ax[2,0], year=year, unit='[K]', title='Abs Theta through PITE trough')
        plotTransect(var='THETAPITE', ax=ax[2,1], year=year, unit='[K]', title='Delta Theta through PITE trough')

#         for ev in list(events):
#             ax[2,2].plot(months, events[ev] ,linewidth=0.1, c='grey')
#         ax[2,2].plot(months, events.mean(axis=1) ,linewidth=1.2, c='k')
#         ax[2,2].fill_between(months, events.mean(axis=1)-events.std(axis=1), events.mean(axis=1)+events.std(axis=1), color='b', alpha=0.2)
#         ax[2,2].set_ylabel('OHC')
#         ax[2,2].set_title('Ocean Heat Content during Events', pad=20)
#         ax[2,2].set_xlabel('Months from event')
        ax[2,2].axvline(x=year, color='r')
        
        fig.suptitle(suptitle+'; Month wrt event: '+str(year))
        return


    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")
    anim = animation.FuncAnimation(fig, animate, frames=months[-1]-months[0], repeat=False, interval=5000)
    writergif = animation.PillowWriter(fps=1)
    anim.save('../03_output/composite/'+today+'_animation_overview_composite_'+filename+'_.gif', writer=writergif)
    return

print(data_new)
OverviewAnimation_new(data_dt=data_new, suptitle=suptitle, filename=output_filename)