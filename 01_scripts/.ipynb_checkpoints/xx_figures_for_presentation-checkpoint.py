# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain test!.

@author: joren
"""

from functions.correlation_maps import simpleCorrelationMap, correlationMap_1D_v3, correlationMap_2D_v3,  showCombinedCorrelationMap_2D_v3#, correlationMap_2D, correlationMap_1D,
from functions.loading_and_processing_data import loadMITgcmData
from functions.other_functions import computeBuoyancy
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt



for window in [1, 60]:
    showCombinedCorrelationMap_2D_v3(members='all',
                      data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      detrend=True, window=window, title='Average Bottom Current (lowest 100m) vs Dotson to Cosgrove Ice Shelf Basal Mass Loss', ymax=-70, save=True, savefn='bottom100m_VEL_dotson_to_cosgrove_'+str(window)+'month_dt_dS',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False, deseasonalize=True)


    showCombinedCorrelationMap_2D_v3(members='all',
                      data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      detrend=True, window=window, title='Depth Averaged Current vs Dotson to Cosgrove Ice Shelf Basal Mass Loss', ymax=-70, save=True, savefn='depth_averaged_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False)

    showCombinedCorrelationMap_2D_v3(members='all',
                      data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
                      detrend=True, window=window, title='Average Bottom Current (lowest 100m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='bottom100m_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False)

    showCombinedCorrelationMap_2D_v3(members='all',
                      data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
                      detrend=True, window=window, title='Depth Averaged Current vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='depth_averaged_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False)


    showCombinedCorrelationMap_2D_v3(members='all',
                          data=None, datafn='baroclinic_UVEL', datavar='UVEL', datakind='old',
                          ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                          detrend=True, window=window, title='Baroclinic Current (lowest 100m-depth averaged) vs Dotson to Cosgrove Ice Shelf Basal Mass Loss', ymax=-70, save=True, 
                          savefn='baroclinic_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
                          draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                          inverse=False,
                          give_back=False)

    showCombinedCorrelationMap_2D_v3(members='all',
                          data=None,datafn='baroclinic_UVEL', datavar='UVEL', datakind='old',
                          ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
                          detrend=True, window=window, title='Baroclinic Current (lowest 100m-depth averaged) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='baroclinic_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                          draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                          inverse=False,
                          give_back=False)

    
    

window=60
showCombinedCorrelationMap_2D_v3(members='all',
                  data=None, datafn='curlTAU', datavar='curlTAU', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  detrend=True, window=window, title='Surface Stress Curl vs Dotson to Cosgrove Ice Shelf Basal Mass Loss', ymax=-70, save=True, savefn='curlTAU_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False)

showCombinedCorrelationMap_2D_v3(members='all',
                  data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  detrend=True, window=window, title='Zonal Wind vs Dotson to Cosgrove Ice Shelf Basal Mass Loss', ymax=-70, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False)

    
window=60

correlationMap_1D_v3(members=[20,21],
                  data=data.expand_dims('ens'), datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
                  ind=wind.to_dataset(name='amundsen_shelf_break_uwind_avg').expand_dims('ens'), indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
                  detrend=True, window=window, title='Ensemble Mean Amundsen Shelf Temperature (200 to 700m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='ensmean_amundsen_shelf_temp_200to700_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, save_results=True)

correlationMap_2D_v3(members=[20,21],
                  data=[u_ensmean.expand_dims('ens'), v_ensmean.expand_dims('ens')], datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
                  ind=wind.to_dataset(name='amundsen_shelf_break_uwind_avg').expand_dims('ens'), indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
                  detrend=True, window=window, title='Ensemble Mean Average Bottom Current (lowest 100m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='ensmean_bottom100m_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False,
                  give_back=False, save_results=True)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
# window=1

# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Average Bottom Current (lowest 100m) vs Dotson to Cosgrove Melt', ymax=-70, save=True, savefn='bottom100m_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Depth Averaged Current vs Dotson to Cosgrove Melt', ymax=-70, save=True, savefn='depth_averaged_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Average Bottom Current (lowest 100m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='bottom100m_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Depth Averaged Current vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='depth_averaged_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)


# window=60
# members='all'

# showCombinedCorrelationMap_2D_v3(members=members,
#                       data=None, datafn='oceBUOY', datavar='oceBUOY', datakind='maps',
#                       ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                       detrend=True, window=window, title='Correlation map of Buoyancy Forcing vs Dotson to Cosgrove Massloss ('+str(window)+' month rolling mean)', ymax=-70, save=True, savefn='oceBUOY_dotson_to_cosgrove_'+str(window)+'month_dt',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False)
    
# showCombinedCorrelationMap_2D_v3(members=members,
#                       data=None, datafn='oceBUOY', datavar='oceBUOY', datakind='maps',
#                       ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                       detrend=True, window=window, title='Correlation map of Buoyancy Forcing vs Amundsen Shelf Break Uwind ('+str(window)+' month rolling mean)', ymax=-70, save=True, savefn='oceBUOY_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False)


# window=24
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='curlTAU', datavar='curlTAU', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Surface Stress Curl vs Dotson to Cosgrove Massloss ('+str(window)+' month rolling mean)', ymax=-70, save=True, savefn='curlTAU_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)