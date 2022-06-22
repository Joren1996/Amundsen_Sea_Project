# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain all the code to create the correlation maps!

@author: joren
"""

#from functions.loading_and_processing_data import loadMITgcmData, loadPACEData, detrend_and_average_MITgcmData, read_all_data
from functions.correlation_maps import simpleCorrelationMap, correlationMap_1D_v3, correlationMap_2D_v3#, correlationMap_2D, correlationMap_1D,
from functions.loading_and_processing_data import loadMITgcmData
from functions.other_functions import computeBuoyancy
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt



# for ens in range(13, 14):
#     members=[ens-1,ens]
#     window=60
#     data=computeBuoyancy(member=ens)

#     correlationMap_1D_v3(members=members,
#                       data=data, datafn='oceBUOY', datavar='oceBUOY', datakind='maps',
#                       ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                       detrend=True, window=window, title='Correlation map of Buoyancy Forcing vs Dotson to Cosgrove Massloss', ymax=-70, save=True, savefn='oceBUOY_dotson_to_cosgrove_'+str(window)+'month_dt',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False)
    
#     correlationMap_1D_v3(members=members,
#                       data=data, datafn='oceBUOY', datavar='oceBUOY', datakind='maps',
#                       ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                       detrend=True, window=window, title='Correlation map of Buoyancy Forcing vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='oceBUOY_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False)

#window=60
# correlationMap_1D_v3(members='all',
#                   data=None, datafn='curlTAU', datavar='curlTAU', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Surface Stress Curl vs Dotson to Cosgrove Massloss', ymax=-70, save=True, savefn='curlTAU_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)


# correlationMap_1D_v2(members=[1,2],
#                   data=None, datafn='curlWIND', datavar='curlWIND', datakind='maps',
#                     ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=24, title='Correlation map of Wind Curl vs Dotson to Cosgrove Massloss', ymax=-70, save=True, savefn='curlWIND_dotson_to_cosgrove_24month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# window=60
# correlationMap_1D_v3(members='all',
#                   data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(0, 24),
#                   detrend=True, window=window, title='Correlation map of Zonal Wind vs Dotson to Cosgrove Massloss', ymax=None, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True)




# correlationMap_1D_v2(members=[1,2],
#                   data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=24, title='Correlation map of Meridional Wind vs Dotson to Cosgrove Massloss', ymax=-70, save=True, savefn='EXFvwind_dotson_to_cosgrove_24month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)


# window=1
# correlationMap_2D_v3(members='all',
#                   data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Average Bottom Current (lowest 100m) vs Dotson to Cosgrove Melt', ymax=-70, save=True, savefn='bottom100m_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True)

# correlationMap_2D_v3(members='all',
#                   data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Depth Averaged Current vs Dotson to Cosgrove Melt', ymax=-70, save=True, savefn='depth_averaged_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True)

# correlationMap_2D_v3(members='all',
#                   data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Average Bottom Current (lowest 100m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='bottom100m_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True)

# correlationMap_2D_v3(members='all',
#                   data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Depth Averaged Current vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='depth_averaged_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True)

window=1
for ens in range(16, 21):
    members=[ens-1,ens]
    data1=loadMITgcmData(filename='bottom100m_averaged_UVEL', members=members, kind='old')
    data2=loadMITgcmData(filename='depth_averaged_UVEL', members=members, kind='old')
    dataU=data1-data2
    data1=loadMITgcmData(filename='bottom100m_averaged_VVEL', members=members, kind='old')
    data2=loadMITgcmData(filename='depth_averaged_VVEL', members=members, kind='old')
    dataV=data1-data2
    data=[dataU,dataV]
    
    correlationMap_2D_v3(members=members,
                      data=data, datafn='baroclinic_UVEL', datavar='UVEL', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      detrend=True, window=window, title='2D Correlation Map of Baroclinic Current (lowest 100m-depth averaged) vs Dotson to Cosgrove Melt', ymax=-70, save=True, 
                      savefn='baroclinic_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False,
                      give_back=False, deseasonalize=True)

    correlationMap_2D_v3(members=members,
                          data=data,datafn='baroclinic_UVEL', datavar='UVEL', datakind='old',
                          ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
                          detrend=True, window=window, title='2D Correlation Map of Baroclinic Current (lowest 100m-depth averaged) vs Amundsen Shelf Break Uwind', ymax=-70, save=True,
                          savefn='baroclinic_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                          draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                          inverse=False,
                          give_back=False, deseasonalize=True)


# correlationMap_2D_v2(members=[1,2],
#                   data=None, datafn='ADVx_TH_corrected_bottom100m_averaged', datavar='ADVx_TH', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=24, title='2D Correlation Map of Average Bottom Heat Advection (lowest 100m) vs Dotson to Cosgrove Melt', ymax=-70, save=True, savefn='bottom100m_ADV_TH_dotson_to_cosgrove_24month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# correlationMap_2D_v2(members=[1,2],
#                   data=None,datafn='ADVx_TH_corrected_depth_averaged', datavar='ADVx_TH', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=24, title='2D Correlation Map of Depth Averaged Heat Advection vs Dotson to Cosgrove Melt', ymax=-70, save=True, savefn='depth_averaged_ADV_TH_dotson_to_cosgrove_24month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# correlationMap_2D_v2(members=[1,2],
#                   data=None, datafn='ADVx_TH_corrected_bottom100m_averaged', datavar='ADVx_TH', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=24, title='2D Correlation Map of Average Heat Advection (lowest 100m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='bottom100m_ADV_TH_amundsen_shelf_break_uwind_24month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# correlationMap_2D_v2(members=[1,2],
#                   data=None, datafn='ADVx_TH_corrected_depth_averaged', datavar='ADVx_TH', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=24, title='2D Correlation Map of Depth Averaged Current vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='depth_averaged_ADV_TH_amundsen_shelf_break_uwind_24month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# members=[1,2]
# data1=loadMITgcmData(filename='ADVx_TH_corrected_bottom100m_averaged', members=members, kind='maps')
# data2=loadMITgcmData(filename='ADVx_TH_corrected_depth_averaged', members=members, kind='maps')
# dataU=data1-data2
# data1=loadMITgcmData(filename='ADVy_TH_corrected_bottom100m_averaged', members=members, kind='maps')
# data2=loadMITgcmData(filename='ADVy_TH_corrected_depth_averaged', members=members, kind='maps')
# dataV=data1-data2
# data=[dataU,dataV]

# correlationMap_2D_v2(members=members,
#                   data=data, datafn='baroclinic_ADVx_TH', datavar='ADVx_TH', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=24, title='2D Correlation Map of Baroclinic Heat Advection (lowest 100m-depth averaged) vs Dotson to Cosgrove Melt', ymax=-70, save=True, 
#                   savefn='baroclinic_ADV_TH_dotson_to_cosgrove_24month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# correlationMap_2D_v2(members=members,
#                   data=data,datafn='baroclinic_ADVx_TH', datavar='ADVx_TH', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=24, title='2D Correlation Map of Baroclinic Heat Advection (lowest 100m-depth averaged) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='baroclinic_ADV_TH_amundsen_shelf_break_uwind_24month_dt',
#                   draw_box=False, box_x=[360-115, 360-102], box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)






#More complicated correlation maps

# window=60

# correlationMap_1D_v3(members='all',
#                   data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Zonal Wind vs Dotson to Cosgrove Massloss', ymax=-70, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)


# correlationMap_1D_v3(members='all',
#                   data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Meridional Wind vs Dotson to Cosgrove Massloss', ymax=-70, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)



# correlationMap_1D_v3(members='all',
#                   data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Temperature (200-700m) vs Dotson to Cosgrove Massloss', ymax=-70, save=True, savefn='amundsentemp_200-700_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# correlationMap_1D_v3(members='all',
#                   data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Temperature (200-700m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='amundsentemp_200-700_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)


# correlationMap_1D_v2(members=[1,2],
#                   data=None, datafn='THETA_depth_integrated', datavar='THETA', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=24, title='Correlation map of Depth Integrated Temperature vs Dotson to Cosgrove Massloss', ymax=-70, save=True, savefn='amundsentemp_depthintegrated_dotson_to_cosgrove_24month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)


# correlationMap_1D_v2(members=[1,2],
#                   data=None, datafn='THETA_depth_integrated', datavar='THETA', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=24, title='Correlation map of Depth Integrated Temperature vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='amundsentemp_depthintegrated_amundsen_shelf_break_uwind_24month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)





#SIMPLE CORRELATION MAPS



# simpleCorrelationMap(members=[1,2],
#                      data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
#                      ind=None, indfn='undercurrent_W115', indvar='UVEL',
#                      detrend=True, window=24, 
#                      title='Correlation between local zonal wind and shelf break undercurrent at W115', 
#                      ymax=-70, save=True, savefn='EXFuwind_UW115_detrend_24month',
#                      give_back=False)

# simpleCorrelationMap(members=[1,2],
#                      data=None, datafn='oceQnet', datavar='oceQnet', datakind='maps',
#                      ind=None, indfn='undercurrent_W115', indvar='UVEL',
#                      detrend=True, window=24, 
#                      title='Correlation between local net heat uptake and shelf break undercurrent at W115', 
#                      ymax=-70, save=True, savefn='oceQnet_UW115_detrend_24month',
#                      give_back=False)

# simpleCorrelationMap(members=[1,2],
#                      data=None, datafn='oceTAUX', datavar='oceTAUX', datakind='maps',
#                      ind=None, indfn='undercurrent_W115', indvar='UVEL',
#                      detrend=True, window=24, 
#                      title='Correlation between local zonal wind stress and shelf break undercurrent at W115', 
#                      ymax=-70, save=True, savefn='oceTAUX_UW115_detrend_24month',
#                      give_back=False)

# simpleCorrelationMap(members=[1,2],
#                      data=None, datafn='THETA_depth_integrated', datavar='THETA', datakind='maps',
#                      ind=None, indfn='undercurrent_W115', indvar='UVEL',
#                      detrend=True, window=24, 
#                      title='Correlation between local depth integrated temperature and shelf break undercurrent at W115', 
#                      ymax=-70, save=True, savefn='THETA_depthintegrated_UW115_detrend_24month',
#                      give_back=False)

# data=xr.open_dataset('/data/hpcdata/users/grejan/mitgcm/02_data/maps/curlWIND_ensmean.nc')
# data=data.expand_dims('ens')
# data1, ind1, res1= simpleCorrelationMap(members='all',
#                      data=data, datafn='curlWIND', datavar='curlWIND', datakind='maps',
#                      ind=None, indfn='undercurrent_W115', indvar='UVEL',
#                      detrend=True, window=24, 
#                      title='Correlation between local wind curl and shelf break undercurrent at W115', 
#                      ymax=-70, save=True, savefn='curlWIND_UW115_detrend_24month',
#                      give_back=True)

# data1, ind1, res1= simpleCorrelationMap(members='all',
#                      data=data, datafn='curlWIND', datavar='curlWIND', datakind='maps',
#                      ind=None, indfn='ADVy_TH_PITE', indvar='ADVy_TH',
#                      detrend=True, window=24, 
#                      title='Correlation between local wind curl and PITE advection', 
#                      ymax=-70, save=True, savefn='curlWIND_UW115_detrend_24month',
#                      give_back=True)


# data=xr.open_dataset('/data/hpcdata/users/grejan/mitgcm/02_data/maps/curlWIND_ensmean.nc')
# data=data.expand_dims('ens')
# data1, ind1, res1= simpleCorrelationMap(members='all',
#                      data=data, datafn='curlWIND', datavar='curlWIND', datakind='maps',
#                      ind=None, indfn='undercurrent_W115', indvar='UVEL',
#                      detrend=False, window=24, 
#                      title='Correlation between local wind curl and shelf break undercurrent at W115', 
#                      ymax=-70, save=True, savefn='curlWIND_UW115_24month',
#                      give_back=True)

# data1, ind1, res1= simpleCorrelationMap(members='all',
#                      data=data, datafn='curlWIND', datavar='curlWIND', datakind='maps',
#                      ind=None, indfn='ADVy_TH_PITE', indvar='ADVy_TH',
#                      detrend=False, window=24, 
#                      title='Correlation between local wind curl and PITE advection', 
#                      ymax=-70, save=True, savefn='curlWIND_UW115_24month',
#                      give_back=True)



# data, ind, res= simpleCorrelationMap(members='all',
#                      data=None, datafn='curlTAU', datavar='curlTAU', datakind='maps',
#                      ind=None, indfn='undercurrent_W115', indvar='UVEL',
#                      detrend=True, window=24, 
#                      title='Correlation between local stress curl and shelf break undercurrent at W115', 
#                      ymax=-70, save=True, savefn='new_curlTAU_UW115_detrend_24month',
#                      give_back=True)