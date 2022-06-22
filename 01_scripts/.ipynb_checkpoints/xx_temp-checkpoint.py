# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain test!.

@author: joren
"""
from functions.loading_and_processing_data import loadMITgcmData, loadPACEData, detrend_and_average_MITgcmData, read_all_data
from functions.correlation_maps import simpleCorrelationMap, correlationMap_1D_v3, correlationMap_2D_v3,  showCombinedCorrelationMap_2D_v3#, correlationMap_2D, correlationMap_1D,
from functions.loading_and_processing_data import loadMITgcmData, loadTimeSeriesData
from functions.other_functions import computeBuoyancy
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


wind=loadTimeSeriesData(filename='timeseries_final', members='all', kind='old')['amundsen_shelf_break_uwind_avg']
wind=wind.mean(dim='ens')



# window=1

# data=xr.open_dataset('/data/hpcdata/users/grejan/mitgcm/data/new2_amundsen_shelf_temp_200to700_ensmean.nc')
# correlationMap_1D_v3(members=[20,21],
#                   data=data.expand_dims('ens'), datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
#                   ind=wind.to_dataset(name='amundsen_shelf_break_uwind_avg').expand_dims('ens'), indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation Map of Ensemble Mean Amundsen Shelf Temperature (200 to 700m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='ensmean_amundsen_shelf_temp_200to700_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, save_results=True)

u_ensmean=xr.open_dataset('/data/hpcdata/users/grejan/mitgcm/02_data/maps/bottom100m_averaged_UVEL_ensmean.nc')
v_ensmean=xr.open_dataset('/data/hpcdata/users/grejan/mitgcm/02_data/maps/bottom100m_averaged_VVEL_ensmean.nc')

# correlationMap_2D_v3(members=[20,21],
#                   data=[u_ensmean.expand_dims('ens'), v_ensmean.expand_dims('ens')], datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=wind.to_dataset(name='amundsen_shelf_break_uwind_avg').expand_dims('ens'), indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation Map of Ensemble Mean Average Bottom Current (lowest 100m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='ensmean_bottom100m_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, save_results=True)

window=60

# correlationMap_1D_v3(members=[20,21],
#                   data=data.expand_dims('ens'), datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
#                   ind=wind.to_dataset(name='amundsen_shelf_break_uwind_avg').expand_dims('ens'), indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation Map of Ensemble Mean Amundsen Shelf Temperature (200 to 700m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='ensmean_amundsen_shelf_temp_200to700_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, save_results=True)

correlationMap_2D_v3(members=[20,21],
                  data=[u_ensmean.expand_dims('ens'), v_ensmean.expand_dims('ens')], datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
                  ind=wind.to_dataset(name='amundsen_shelf_break_uwind_avg').expand_dims('ens'), indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
                  detrend=True, window=window, title='Correlation Map of Ensemble Mean Average Bottom Current (lowest 100m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='ensmean_bottom100m_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, save_results=True)


melt=loadTimeSeriesData(filename='timeseries_final', members='all', kind='old')['dotson_to_cosgrove_massloss']
melt=melt.mean(dim='ens')
data=xr.open_dataset('/data/hpcdata/users/grejan/mitgcm/02_data/maps/EXFuwind_ensmean.nc')

window=60
correlationMap_1D_v3(members=[20,21],
                  data=data.expand_dims('ens'), datafn='EXFuwind', datavar='EXFuwind', datakind='old',
                  ind=melt.to_dataset(name='dotson_to_cosgrove_massloss').expand_dims('ens'), indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  detrend=True, window=window, title='Ensemble Mean EXFuwind vs Dotson to Cosgrove Mass Loss', ymax=None, save=True, savefn='ensmean_EXFuwind_dotson_to_cosgrove_massloss_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, save_results=True)


# window=60
# members='all'

# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Temperature (200-700m) vs Dotson to Cosgrove Massloss', ymax=-70, save=True, savefn='amundsentemp_200-700_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='Temperature (200-700m) vs Amundsen Shelf Break Uwind ('+str(window)+' month rolling mean)', ymax=-70, save=True, savefn='amundsentemp_200-700_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)




# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Average Bottom Current (lowest 100m) vs Dotson to Cosgrove Melt \n ('+str(window)+' month rolling mean)', ymax=-70, save=True, savefn='bottom100m_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Depth Averaged Current vs Dotson to Cosgrove Melt ('+str(window)+' month rolling mean)', ymax=-70, save=True, savefn='depth_averaged_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='Average Bottom Current (lowest 100m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='bottom100m_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='Depth Averaged Current vs Amundsen Shelf Break Uwind ('+str(window)+' month rolling mean)', ymax=-70, save=True, savefn='depth_averaged_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)


# showCombinedCorrelationMap_2D_v3(members=members,
#                       data=None, datafn='baroclinic_UVEL', datavar='UVEL', datakind='old',
#                       ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                       detrend=True, window=window, title='Baroclinic Current (lowest 100m-depth averaged) vs Dotson to Cosgrove Melt ('+str(window)+' month rolling mean)', ymax=-70, save=True, 
#                       savefn='baroclinic_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False)

# showCombinedCorrelationMap_2D_v3(members=members,
#                       data=None,datafn='baroclinic_UVEL', datavar='UVEL', datakind='old',
#                       ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                       detrend=True, window=window, title='Baroclinic Current (lowest 100m-depth averaged) vs Amundsen Shelf Break Uwind ('+str(window)+' month rolling mean)', ymax=-70, save=True, savefn='baroclinic_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False)

# showCombinedCorrelationMap_2D_v3(members=members,
#                   data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Zonal Wind vs Dotson to Cosgrove Massloss ('+str(window)+' month rolling mean)', ymax=-70, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)




# window=24
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='curlTAU', datavar='curlTAU', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Surface Stress Curl vs Dotson to Cosgrove Massloss ('+str(window)+' month rolling mean)', ymax=-70, save=True, savefn='curlTAU_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)