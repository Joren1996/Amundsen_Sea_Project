# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 16:53:42 2022

@author: joren
"""

#%% IMPORTING
import unixFunctionsMITgcm_final as fc
import xarray as xr

#%% Composite Maps

#Reading the data.
full='/data/hpcdata/users/grejan/mitgcm/data'
total={}
units={}
longnames={}
total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='isotherm')
total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='final')

#Test
selected,_=fc.findEvents(total, units, longnames, threshold=0.15, wvar='amundsen_shelf_temp_btw_200_700m', larger_than=True, filename='test')


#----------------------------------------------------------------------------------
#%% INITIALISATION
members='all' #[0,1]
var='THETA'
# #Loading the data.
data=fc.loadMITgcmData(filename='new2_amundsen_shelf_temp_200to700', members=members)

#Composite Maps
window=24
data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)

#Warm Events
selected,_=fc.findEvents(total, units, longnames, threshold=0.15, wvar='amundsen_shelf_temp_btw_200_700m', larger_than=True, filename='warm')
fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during warm events (>0.15K)', 
                   filename='24month_mean_amundstemp_warm',
                   window=window)

#High Melt Events
selected,_=fc.findEvents(total, units, longnames, threshold=30, wvar='dotson_to_cosgrove_massloss', larger_than=True, filename='highmelt')
fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during high melt events (>30Gt/yr)', 
                   filename='24month_mean_amundstemp_high_melt',
                   window=window)

#Cold Events
selected,_=fc.findEvents(total, units, longnames, threshold=-0.15, wvar='amundsen_shelf_temp_btw_200_700m', larger_than=False, filename='cold')
fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during cold events (<-0.15K)', 
                   filename='24month_mean_amundstemp_cold',
                   window=window)

#Low Melt Events
selected,_=fc.findEvents(total, units, longnames, threshold=-30, wvar='dotson_to_cosgrove_massloss', larger_than=False, filename='lowmelt')
fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during low melt events (<-30Gt/yr)', 
                   filename='24month_mean_amundstemp_low_melt',
                   window=window)


#1940s event!
selected[selected==True]=False
selected.loc['1941-01']=True
fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during 1940s (wrt 1941-01)', 
                   filename='1940s_mean_amundstemp',
                   window=window, months=[-25, 24])


#Monthly picture
window=1
data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)

deseasonalize=True
if deseasonalize==True:
     data_dt = data_dt.groupby("time.month")-data_dt.groupby("time.month").mean(dim="time")

#Warm Events
selected,_=fc.findEvents(total, units, longnames, threshold=0.15, wvar='amundsen_shelf_temp_btw_200_700m', larger_than=True, filename='warm')
fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during warm events (>0.15K)', 
                   filename='1month_mean_amundstemp_warm',
                   window=window)








#%% Composite Maps

# window=24
# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)

# fc.createCompositeMap(data_dt, title='Composite map of local surface stress curl during high melt events', 
#                    filename='24monthmean_stress_curl',
#                    unit='[N/m^3]',
#                    window=window)


# var='THETA'
# tmvar='amundsen_shelf_temp_btw_200_700m'


#var='oceTAUX'
#dataX=fc.readMITgcmData(var=var, members=members)


#dataX=dataX.rename({'oceTAUX':'curlTAU', 'XG':'XC'})
#dataX=dataX.differentiate('YC')

#var='oceTAUY'
#data=fc.readMITgcmData(var=var, members=members)
#data=data.rename({'oceTAUY':'curlTAU', 'YG':'YC'})
#data=data.differentiate('XC')

#data=xr.Dataset(data_vars=dict(curlTAU=(['ens', 'time', 'YC', 'XC'], data['curlTAU'].values-dataX['curlTAU'].values)), coords=dict(ens=dataX.ens, time=dataX.time, XC=dataX.XC, YC=dataX.YC))

#var='curlTAU'

#save=True
#deseasonalize=True

#%%Reading the data.
# =============================================================================
# data=fc.readMITgcmData(var=var, members=members)
# 
#if save==True:
#    for i in range(len(data.ens)):
#        data.sel(ens=i).to_netcdf('./data/curl_tau_ens'+str(i+1)+'.nc')
# =============================================================================


#data=fc.loadMITgcmData(filename='amundsen_shelf_temp_200to700', members=members)
#%% Correlation Maps

# =============================================================================
# tm=fc.readMITgcmTimeSeries(var=tmvar, members=members)
# 
# window=24
# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)
# tm_dt=fc.detrend_and_average_MITgcmData(tm, tmvar, window=window)
# 
# tsvar='amundsen_shelf_break_uwind_avg'
# ts=fc.readMITgcmTimeSeries(var=tsvar, members=members)
# ts=fc.detrend_and_average_MITgcmData(ts, tsvar, window=window)
# fc.createCorrelationMap(data=data_dt, ts=ts, tm=tm_dt, 
#                         title='Correlation map of local Amundsen sea temperature (200-700m) with zonal shelf winds over PITT', 
#                         filename='correlation_map_amundstemp_zonalwind', members=members, window=window)
# 
# ts=ts.cumsum(dim='time')
# fc.createCorrelationMap(data=data_dt, ts=ts, tm=tm_dt, 
#                         title='Correlation map of local Amundsen sea temperature (200-700m) with cumulative zonal shelf winds over PITT', 
#                         filename='correlation_map_amundstemp_cumul_zonalwind', members=members, window=window)
# =============================================================================



# window=1
# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)

# if deseasonalize==True:
#     data_dt = data_dt.groupby("time.month")-data_dt.groupby("time.month").mean(dim="time")

# fc.createCompositeMap(data_dt, title='Composite map of local winds stress curl during high melt events', 
#                    filename='1monthmean_curlTAU_deseasonalized',
#                    unit='[N/m^2]',
#                    window=window)



# window=1
# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)

# if deseasonalize==True:
#    data_dt = data_dt.groupby("time.month")-data_dt.groupby("time.month").mean(dim="time")

# fc.createCompositeMap(data_dt, title='Composite map of local Amundsen sea temperature (200-700m) during high melt events', 
#                    filename='1monthmean_amundstemp_deseasonalized',
#                    unit='[K]',
#                    window=window)

# fc.createCompositeMap(data_dt, title='Composite map of local Amundsen sea temperature (200-700m) during low melt events', 
#                    filename='1monthmean_amundstemp_lowmelt_deseasonalized',
#                    unit='[K]',
#                    window=window,
#                    larger_than=False, pct=10)

# =============================================================================
# fc.createCompositeMap(data_dt, title='Composite map of local Amundsen sea temperature (200-700m) during high melt events', 
#                    filename='24monthmean_amundstemp',
#                    tmvar=tmvar,
#                    window=window)
# 
# fc.createCompositeMap(data_dt, title='Composite map of local Amundsen sea temperature (200-700m) during low melt events', 
#                    filename='24monthmean_amundstemp_lowmelt',
#                    tmvar=tmvar,
#                    window=window,
#                    larger_than=False, pct=10)
# =============================================================================



