# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 16:53:42 2022

@author: joren
"""

#%% IMPORTING
import unixFunctionsMITgcm_final as fc
import xarray as xr

#Reading the data.
full='/data/hpcdata/users/grejan/mitgcm/data'
total={}
units={}
longnames={}
total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='isotherm')
total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='final')
total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='ohc')
total, units, longnames=fc.read_all_data(full, total, units, longnames, kind='trough')
        
threshold=fc.defineEvents(total=total, longnames=longnames, units=units, wvar='amundsen_shelf_ohc_below_0m', larger_than=True, pct=75, detrend=True, window=12*2, filename='warm75_colored_in')
selected,_=fc.findEvents(total=total, longnames=longnames, units=units, threshold=threshold, wvar='amundsen_shelf_ohc_below_0m', larger_than=True, detrend=True, window=12*2, filename='warm75_colored_in')



#READ SEA ICE FLUX

save=True

#%% INITIALISATION
members='all' #[0,1]
var='SHIfwFlx'


data=fc.loadMITgcmData(filename='seaiceflux', members=members)

# #Loading the data.
# data=fc.readMITgcmData(var=var, members=members, save=True)
# #data=fc.readMITgcmData_depth_averaged(var=var, members=members, save=True)

# if save==True:
#     for i in range(len(data.ens)):
#         data.sel(ens=i).to_netcdf('./data/seaiceflux_ens'+str(i+1)+'.nc')


#Covective Maps
window=24
data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)

choose='only'
if choose=='full':
    title='\n Selecting convection based on Definition Naughten'
    convection=(total['amundsen_shelf_isotherm_0.5C_below_100m'].rolling(24, center=True).mean()<-440) & (total['pine_island_bay_isotherm_0C_below_100m'].rolling(24, center=True).mean()<-470) & (total['dotson_bay_isotherm_-1C_below_100m'].rolling(24, center=True).mean()<-430)
elif choose=='only':
    title='\n Selecting convection based on Isotherms in PIB'
    convection=(
                    #total[ 'amundsen_shelf_isotherm_0.5C_below_100m'].rolling(24).mean()<-430 #| 
                    total['pine_island_bay_isotherm_0C_below_100m'].rolling(24, center=True).mean()<-470#/600 
                    #total['dotson_bay_isotherm_-1C_below_100m'].rolling(24, center=True).mean()<-430
                    )
if choose=='two':
    title='\n Selecting convection based on Isotherms in PIB and Amundsen Shelf'
    convection=(total['amundsen_shelf_isotherm_0.5C_below_100m'].rolling(24, center=True).mean()<-440) & (total['pine_island_bay_isotherm_0C_below_100m'].rolling(24, center=True).mean()<-470)# & (total['dotson_bay_isotherm_-1C_below_100m'].rolling(24, center=True).mean()<-430)
    
selected=convection['1920':'2013']        
fc.createCompositeMap2(data_dt, selected, title='Composite map of Sea Ice Flux during PIB convection', 
                   filename='24month_mean_seaice_PIB_convection',
                   window=window)

choose='full'
if choose=='full':
    title='\n Selecting convection based on Definition Naughten'
    convection=(total['amundsen_shelf_isotherm_0.5C_below_100m'].rolling(24, center=True).mean()<-440) & (total['pine_island_bay_isotherm_0C_below_100m'].rolling(24, center=True).mean()<-470) & (total['dotson_bay_isotherm_-1C_below_100m'].rolling(24, center=True).mean()<-430)
elif choose=='only':
    title='\n Selecting convection based on Isotherms in PIB'
    convection=(
                    #total[ 'amundsen_shelf_isotherm_0.5C_below_100m'].rolling(24).mean()<-430 #| 
                    total['pine_island_bay_isotherm_0C_below_100m'].rolling(24, center=True).mean()<-470#/600 
                    #total['dotson_bay_isotherm_-1C_below_100m'].rolling(24, center=True).mean()<-430
                    )
if choose=='two':
    title='\n Selecting convection based on Isotherms in PIB and Amundsen Shelf'
    convection=(total['amundsen_shelf_isotherm_0.5C_below_100m'].rolling(24, center=True).mean()<-440) & (total['pine_island_bay_isotherm_0C_below_100m'].rolling(24, center=True).mean()<-470)# & (total['dotson_bay_isotherm_-1C_below_100m'].rolling(24, center=True).mean()<-430)
    
selected=convection['1920':'2013']      
fc.createCompositeMap2(data_dt, selected, title='Composite map of Sea Ice Flux during full convection', 
                   filename='24month_mean_seaice_full_convection',
                   window=window)

choose='only'
if choose=='full':
    title='\n Selecting convection based on Definition Naughten'
    convection=(total['amundsen_shelf_isotherm_0.5C_below_100m'].rolling(24, center=True).mean()<-440) & (total['pine_island_bay_isotherm_0C_below_100m'].rolling(24, center=True).mean()<-470) & (total['dotson_bay_isotherm_-1C_below_100m'].rolling(24, center=True).mean()<-430)
elif choose=='only':
    title='\n Selecting convection based on Isotherms in PIB'
    convection=(
                    total[ 'amundsen_shelf_isotherm_0.5C_below_100m'].rolling(24).mean()<-430 #| 
                    #total['pine_island_bay_isotherm_0C_below_100m'].rolling(24, center=True).mean()<-470#/600 
                    #total['dotson_bay_isotherm_-1C_below_100m'].rolling(24, center=True).mean()<-430
                    )
if choose=='two':
    title='\n Selecting convection based on Isotherms in PIB and Amundsen Shelf'
    convection=(total['amundsen_shelf_isotherm_0.5C_below_100m'].rolling(24, center=True).mean()<-440) & (total['pine_island_bay_isotherm_0C_below_100m'].rolling(24, center=True).mean()<-470)# & (total['dotson_bay_isotherm_-1C_below_100m'].rolling(24, center=True).mean()<-430)
    
selected=convection['1920':'2013']    
fc.createCompositeMap2(data_dt, selected, title='Composite map of Sea Ice Flux during Amundsen Shelf convection', 
                   filename='24month_mean_seaice_amundsen_convection',
                   window=window)
        

tmvar='amundsen_shelf_break_uwind_avg'
tm=fc.readMITgcmTimeSeries(var=tmvar, members=members)
tm_dt=fc.detrend_and_average_MITgcmData(tm, tmvar, window=window)

tsvar='amundsen_shelf_temp_btw_200_700m'
ts=fc.readMITgcmTimeSeries(var=tsvar, members=members)
ts_dt=fc.detrend_and_average_MITgcmData(ts, tsvar, window=window)


selected=convection['1920':'2013']
fc.createCorrelationMap(data=data_dt, ts=ts, tm=tm_dt, 
                        title='Correlation map of local ice flux and Amundsen Shelf temperature', 
                        filename='correlation_map_seaice_amundstemp', members=members, window=window)
    
# fc.createCorrelationMap(data=data_dt, ts=ts, tm=tm_dt, 
#                          title='Correlation map of local zonal wind with PITE depth averaged meridional flow', 
#                          filename='correlation_map_zonal_wind_PITE_flow', members=members, window=window) 
    
    
    
    
        
# data=fc.readMITgcmData_Slice(var='THETA', members='all', method='lon', x=100)

# if save==True:
#     for i in range(len(data.ens)):
#         data.sel(ens=i).to_netcdf('./data/theta_W100_ens'+str(i+1)+'.nc')


# data=fc.readMITgcmData_Slice(var='THETA', members='all', method='lon', x=103.5)

# if save==True:
#     for i in range(len(data.ens)):
#         data.sel(ens=i).to_netcdf('./data/theta_W103-5_ens'+str(i+1)+'.nc')


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


# #%% INITIALISATION
# members='all' #[0,1]
# var='THETA'
# # #Loading the data.
# data=fc.loadMITgcmData(filename='new2_amundsen_shelf_temp_200to700', members=members)

# #Composite Maps
# window=24
# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)



# fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during warm events', 
#                    filename='24month_mean_amundstemp_warm',
#                    window=window)



# threshold=fc.defineEvents(total=total, longnames=longnames, units=units, wvar='amundsen_shelf_ohc_below_0m', larger_than=False, pct=25, detrend=True, window=12*2, filename='cold25_colored_in')
# print(threshold)
# selected,_=fc.findEvents(total=total, longnames=longnames, units=units, threshold=threshold, wvar='amundsen_shelf_ohc_below_0m', larger_than=False, detrend=True, window=12*2, filename='cold25_colored_in')


# #Cold Events
# fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during cold events', 
#                    filename='24month_mean_amundstemp_cold',
#                    window=window)




# #1940s event!
# selected[selected==True]=False
# selected.loc['1941-01']=True
# fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during 1940s (wrt 1941-01)', 
#                    filename='1940s_24month_mean_amundstemp',
#                    window=window, months=[-25, 24])




# #Monthly picture
# window=1
# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)

# deseasonalize=True
# if deseasonalize==True:
#      data_dt = data_dt.groupby("time.month")-data_dt.groupby("time.month").mean(dim="time")

        

# threshold=fc.defineEvents(total=total, longnames=longnames, units=units, wvar='amundsen_shelf_ohc_below_0m', larger_than=True, pct=75, detrend=True, window=12*2, filename='warm75_colored_in')
# print(threshold)
# selected,_=fc.findEvents(total=total, longnames=longnames, units=units, threshold=threshold, wvar='amundsen_shelf_ohc_below_0m', larger_than=True, detrend=True, window=12*2, filename='warm75_colored_in')


# fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during warm events', 
#                    filename='1month_mean_amundstemp_warm',
#                    window=window)


# threshold=fc.defineEvents(total=total, longnames=longnames, units=units, wvar='amundsen_shelf_ohc_below_0m', larger_than=False, pct=25, detrend=True, window=12*2, filename='cold25_colored_in')
# print(threshold)
# selected,_=fc.findEvents(total=total, longnames=longnames, units=units, threshold=threshold, wvar='amundsen_shelf_ohc_below_0m', larger_than=False, detrend=True, window=12*2, filename='cold25_colored_in')


# #Cold Events
# fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during cold events', 
#                    filename='1month_mean_amundstemp_cold',
#                    window=window)


# #1940s event!
# selected[selected==True]=False
# selected.loc['1941-01']=True
# fc.createCompositeMap2(data_dt, selected, title='Composite map of local Amundsen Shelf Temperature (200-700m) during 1940s (wrt 1941-01)', 
#                    filename='1940s_1month_mean_amundstemp',
#                    window=window, months=[-25, 24])






# #%% INITIALISATION
# members='all' #[0,1]
# var='THETA'

# # #Loading the data.
# #data=fc.readMITgcmData_depth_averaged(var=var, members=members, save=True)


# data=fc.loadMITgcmData(filename='new_amundsen_shelf_temp_200to700', members='all')


# #%% Composite Maps
# window=24
# data_dt=fc.detrend_and_average_MITgcmData(data, var, window=window)
# fc.createCompositeMap(data_dt, title='Composite map of Amundsen Shelf Temperature between 200-700m during high melt events', 
#                    filename='24month_mean_amundtemp_new',
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



