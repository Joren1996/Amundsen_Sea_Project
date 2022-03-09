# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 16:53:42 2022

@author: joren
"""

#%% IMPORTING
import unixFunctionsMITgcm_final as fc
import xarray as xr


#%% INITIALISATION
members='all' #[0,1]
var='THETA'

# #Loading the data.
data=fc.readMITgcmData_depth_averaged(var=var, members=members, save=True)


#data=fc.loadMITgcmData(filename='new_amundsen_shelf_temp_200to700', members='all')


#%% Composite Maps
window=24
data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)

fc.createCompositeMap(data_dt, title='Composite map of Amundsen Shelf Temperature between 200-700m during high melt events', 
                   filename='24month_mean_amundtemp_new',
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



