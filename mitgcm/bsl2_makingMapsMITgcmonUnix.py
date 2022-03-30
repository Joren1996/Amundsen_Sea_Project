# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 16:53:42 2022

@author: joren
"""

#%% IMPORTING
import xarray as xr
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import os
import numpy as np
import unixFunctionsMITgcm_final as fc

#%% Correlation Maps

filename='zonal_winds'
var='EXFuwind'
tmvar='amundsen_shelf_break_uwind_avg'



members='all'

# =============================================================================
data=fc.loadMITgcmData(filename=filename, members=members)
tm=fc.readMITgcmTimeSeries(var=tmvar, members=members)

window=24
data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)
tm_dt=fc.detrend_and_average_MITgcmData(tm, tmvar, window=window)

#tsvar='amundsen_shelf_break_uwind_avg'
#ts=fc.readMITgcmTimeSeries(var=tsvar, members=members)
#ts=fc.detrend_and_average_MITgcmData(ts, tsvar, window=window)

ts1=xr.open_dataarray('./data/PITE_depth_averaged_vflow.nc')
ts1=ts1.drop('YG')

ts2=xr.open_dataarray('./data/PITE_bottom100m_vflow.nc')
ts2=ts2.drop('YG')

ts=ts2-ts1

#ts=ts.sel(ens=slice(0,2))

fc.createCorrelationMap(data=data_dt, ts=ts, tm=tm_dt, 
                        title='Correlation map of local zonal wind with PITE baroclinic flow', 
                        filename='correlation_map_zonal_wind_PITE_bottom_anomaly', members=members, window=window)




# ts=xr.open_dataarray('./data/BRE_depth_averaged_vflow.nc')
# ts=ts.drop('YG')
# #ts=ts.sel(ens=slice(0,2))

# fc.createCorrelationMap(data=data_dt, ts=ts, tm=tm_dt, 
#                         title='Correlation map of local zonal wind with BRE depth averaged meridional flow', 
#                         filename='correlation_map_zonal_wind_BRE_flow', members=members, window=window)

# ts=xr.open_dataarray('./data/BRW_depth_averaged_vflow.nc')
# ts=ts.drop('YG')
# #ts=ts.sel(ens=slice(0,2))

# fc.createCorrelationMap(data=data_dt, ts=ts, tm=tm_dt, 
#                         title='Correlation map of local zonal wind with BRW depth averaged meridional flow', 
#                         filename='correlation_map_zonal_wind_BRW_flow', members=members, window=window)


# data_dt=data_dt.cumsum(dim='time')



# ts=xr.open_dataarray('./data/PITE_depth_averaged_vflow.nc')
# ts=ts.drop('YG')
# #ts=ts.sel(ens=slice(0,2))

# fc.createCorrelationMap(data=data_dt, ts=ts, tm=tm_dt, 
#                         title='Correlation map of local cumulative zonal wind with PITE depth averaged meridional flow', 
#                         filename='correlation_map_cum_zonal_wind_PITE_flow', members=members, window=window)


# ts=xr.open_dataarray('./data/BRE_depth_averaged_vflow.nc')
# ts=ts.drop('YG')
# #ts=ts.sel(ens=slice(0,2))

# fc.createCorrelationMap(data=data_dt, ts=ts, tm=tm_dt, 
#                         title='Correlation map of local cumulative zonal wind with BRE depth averaged meridional flow', 
#                         filename='correlation_map_cum_zonal_wind_BRE_flow', members=members, window=window)

# ts=xr.open_dataarray('./data/BRW_depth_averaged_vflow.nc')
# ts=ts.drop('YG')
# #ts=ts.sel(ens=slice(0,2))

# fc.createCorrelationMap(data=data_dt, ts=ts, tm=tm_dt, 
#                         title='Correlation map of local cumulative zonal wind with BRW depth averaged meridional flow', 
#                         filename='correlation_map_cum_zonal_wind_BRW_flow', members=members, window=window)





#ts=ts.cumsum(dim='time')
#fc.createCorrelationMap(data=data_dt, ts=ts, tm=tm_dt, 
#                        title='Correlation map of local Amundsen sea temperature (200-700m) with cumulative zonal shelf winds over PITT', 
#                        filename='correlation_map_amundstemp_cumul_zonalwind', members=members, window=window)
# =============================================================================



#data_dt=data_dt.cumsum(dim='time')
























# #Reading the data.
# full='/data/hpcdata/users/grejan/mitgcm/data'
# total={}
# units={}
# longnames={}
# total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='isotherm')
# total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='final')
# total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='ohc')


# filename='zonal_winds'
# var='EXFuwind'

# members='all'
# window=24


# data=fc.loadMITgcmData(filename=filename, members=members)
# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)
# data_dt=data_dt.cumsum(dim='time')

# PCs, tc=fc.performPCA(data=data_dt, var=var, method='all', interest=10, filepath='./figures/', filename='full_ensemble_cumsum') 

# PCs, tc=fc.performPCA(data=data_dt, var=var, method='ens_mean', interest=10, filepath='./figures/', filename='ensemble_mean_cumsum') 


# data_iv=data_dt-data_dt.mean(dim='ens')
# PCs, tc=fc.performPCA(data=data_iv, var=var, method='all', interest=10, filepath='./figures/', filename='internal_variability_cumsum') 




#fc.plotPCACorrs(PCs=PCs, total=total, real_var=var, data_dt=data_dt, filename='full_ensemble')





































# #%% Composite Maps


# #Reading the data.
# full='/data/hpcdata/users/grejan/mitgcm/data'
# total={}
# units={}
# longnames={}
# total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='isotherm')
# total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='final')
# total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='ohc')



# threshold=fc.defineEvents(total=total, longnames=longnames, units=units, wvar='amundsen_shelf_ohc_below_0m', larger_than=True, pct=75, detrend=True, window=12*2, filename='warm75_colored_in')
# print(threshold)
# selected,_=fc.findEvents(total=total, longnames=longnames, units=units, threshold=threshold, wvar='amundsen_shelf_ohc_below_0m', larger_than=True, detrend=True, window=12*2, filename='warm75_colored_in')


# # save=True

# # data=fc.readMITgcmData_Slice(var='THETA', members='all', x=106)

# # if save==True:
# #     for i in range(len(data.ens)):
# #         data.sel(ens=i).to_netcdf('./data/theta_W106_ens'+str(i+1)+'.nc')

# members='all'
# data=fc.loadMITgcmData(filename='theta_W106', members=members)

        
# var='THETA'
# x=106

# window=1
# fc.createLatAnimation(data, selected, var=var, title='Cross Sections of '+data[var].long_name + 'at {} degrees W (Warm Events)'.format(str(x)), 
#                    filename='1month_absolute_warm', window=window, method=None, 
#                    longwindow=12*25, 
#                    unit='[degC]',
#                    members='all', larger_than=True, pct=90, startyr='1920', endyr='2013', months=[-36, 24], x=106)

# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)
# fc.createLatAnimation(data_dt.to_dataset(), selected, var=var, title='Cross Sections of Anomaly of '+data[var].long_name + 'at {} degrees W (Warm Events)'.format(str(x)), 
#                    filename='1month_anomaly_warm', window=1, method='mean', 
#                    longwindow=12*25, 
#                    unit='[degC]',
#                    members='all', larger_than=True, pct=90, startyr='1920', endyr='2013', months=[-36, 24], x=106)

# window=24
# fc.createLatAnimation(data, selected, var=var, title='Cross Sections of '+data[var].long_name + 'at {} degrees W (Warm Events)'.format(str(x)), 
#                    filename='24month_absolute_warm', window=window, method=None, 
#                    longwindow=12*25, 
#                    unit='[degC]',
#                    members='all', larger_than=True, pct=90, startyr='1920', endyr='2013', months=[-36, 24], x=106)

# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)
# fc.createLatAnimation(data_dt.to_dataset(), selected, var=var, title='Cross Sections of Anomaly of '+data[var].long_name + 'at {} degrees W (Warm Events)'.format(str(x)), 
#                    filename='24month_anomaly_warm', window=1, method='mean', 
#                    longwindow=12*25, 
#                    unit='[degC]',
#                    members='all', larger_than=True, pct=90, startyr='1920', endyr='2013', months=[-36, 24], x=106)







# threshold=fc.defineEvents(total=total, longnames=longnames, units=units, wvar='amundsen_shelf_ohc_below_0m', larger_than=False, pct=25, detrend=True, window=12*2, filename='cold25_colored_in')
# print(threshold)
# selected,_=fc.findEvents(total=total, longnames=longnames, units=units, threshold=threshold, wvar='amundsen_shelf_ohc_below_0m', larger_than=False, detrend=True, window=12*2, filename='cold25_colored_in')

# window=1
# fc.createLatAnimation(data, selected, var=var, title='Cross Sections of '+data[var].long_name + 'at {} degrees W (Cold Events)'.format(str(x)), 
#                    filename='1month_absolute_cold', window=window, method=None, 
#                    longwindow=12*25, 
#                    unit='[degC]',
#                    members='all', larger_than=True, pct=90, startyr='1920', endyr='2013', months=[-36, 24], x=106)

# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)
# fc.createLatAnimation(data_dt.to_dataset(), selected, var=var, title='Cross Sections of Anomaly of '+data[var].long_name + 'at {} degrees W (Cold Events)'.format(str(x)), 
#                    filename='1month_anomaly_cold', window=1, method='mean', 
#                    longwindow=12*25, 
#                    unit='[degC]',
#                    members='all', larger_than=True, pct=90, startyr='1920', endyr='2013', months=[-36, 24], x=106)

# window=24
# fc.createLatAnimation(data, selected, var=var, title='Cross Sections of '+data[var].long_name + 'at {} degrees W (Cold Events)'.format(str(x)), 
#                    filename='24month_absolute_cold', window=window, method=None, 
#                    longwindow=12*25, 
#                    unit='[degC]',
#                    members='all', larger_than=True, pct=90, startyr='1920', endyr='2013', months=[-36, 24], x=106)

# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)
# fc.createLatAnimation(data_dt.to_dataset(), selected, var=var, title='Cross Sections of Anomaly of '+data[var].long_name + 'at {} degrees W (Cold Events)'.format(str(x)), 
#                    filename='24month_anomaly_cold', window=1, method='mean', 
#                    longwindow=12*25, 
#                    unit='[degC]',
#                    members='all', larger_than=True, pct=90, startyr='1920', endyr='2013', months=[-36, 24], x=106)















# #----------------------------------------------------------------------------------
# #%% INITIALISATION
# members='all' #[0,1]
# var='THETA'
# # #Loading the data.
# data=fc.loadMITgcmData(filename='new2_amundsen_shelf_temp_200to700', members=members)

# #Composite Maps
# window=24
# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)

# #Warm Events
# selected,_=fc.findEvents(total, units, longnames, threshold=0.15, wvar='amundsen_shelf_temp_btw_200_700m', larger_than=True, filename='warm')
# fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during warm events (>0.15K)', 
#                    filename='24month_mean_amundstemp_warm',
#                    window=window)

# #High Melt Events
# selected,_=fc.findEvents(total, units, longnames, threshold=30, wvar='dotson_to_cosgrove_massloss', larger_than=True, filename='highmelt')
# fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during high melt events (>30Gt/yr)', 
#                    filename='24month_mean_amundstemp_high_melt',
#                    window=window)

# #Cold Events
# selected,_=fc.findEvents(total, units, longnames, threshold=-0.15, wvar='amundsen_shelf_temp_btw_200_700m', larger_than=False, filename='cold')
# fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during cold events (<-0.15K)', 
#                    filename='24month_mean_amundstemp_cold',
#                    window=window)

# #Low Melt Events
# selected,_=fc.findEvents(total, units, longnames, threshold=-30, wvar='dotson_to_cosgrove_massloss', larger_than=False, filename='lowmelt')
# fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during low melt events (<-30Gt/yr)', 
#                    filename='24month_mean_amundstemp_low_melt',
#                    window=window)


# #1940s event!
# selected[selected==True]=False
# selected.loc['1941-01']=True
# fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during 1940s (wrt 1941-01)', 
#                    filename='1940s_mean_amundstemp',
#                    window=window, months=[-25, 24])


# #Monthly picture
# window=1
# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)

# deseasonalize=True
# if deseasonalize==True:
#      data_dt = data_dt.groupby("time.month")-data_dt.groupby("time.month").mean(dim="time")

# #Warm Events
# selected,_=fc.findEvents(total, units, longnames, threshold=0.15, wvar='amundsen_shelf_temp_btw_200_700m', larger_than=True, filename='warm')
# fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during warm events (>0.15K)', 
#                    filename='1month_mean_amundstemp_warm',
#                    window=window)


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



