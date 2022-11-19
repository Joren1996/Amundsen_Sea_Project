# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain test!.

@author: joren
"""

from functions.correlation_maps import correlationMap_1D_v3, correlationMap_2D_v3,  showCombinedCorrelationMap_2D_v3, showMapPACE, createPACECorrelationMaps, createWindCorrelationMaps, createFlowCorrelationMaps, createMapsOfProcess#, correlationMap_2D, correlationMap_1D, showCombinedCorMap_Barotropic_Baroclinic,
from functions.loading_and_processing_data import loadMITgcmData
from functions.other_functions import computeBuoyancy
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

#Figure 5
createWindCorrelationMaps(save=True)

#Figure 7
createFlowCorrelationMaps(save=True)

#Figure 2
window=60
showCombinedCorrelationMap_2D_v3(members='all',
                      data=None, datafn='ADVx_TH_corrected_bottom100m_averaged', datavar='ADVx_TH', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      detrend=True, window=window, title='Correlation between the Average Bottom Heat Advection (lowest 100m) \n and Dotson to Cosgrove Mass Loss', ymax=-70, save=True, savefn='bottom100m_ADV_TH_dotson_to_cosgrove_'+str(window)+'month_dt',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      inverse=False, deseasonalize=False,
                      give_back=False)


#Figure 4
window=60
showCombinedCorrelationMap_2D_v3(members='all',
                  data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
                  ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
                  detrend=True, window=window, title='Correlation between the Average Bottom Current (lowest 100m) \n and Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='bottom100m_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, 
                  give_back=False, deseasonalize=False)



#Figure 8
window=60
showCombinedCorrelationMap_2D_v3(members='all',
                  data=None, datafn='SALT_averaged_-200to0', datavar='SALT', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  detrend=True, window=window, title='Correlation between Average Salinity (-200 to 0m) \n and Dotson to Cosgrove Massloss', ymax=-70, save=True, 
                  savefn='SALT-200to0_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=False)


#NOTICE UP TO 12 month LEAD...
window=60
showCombinedCorrelationMap_2D_v3(members='all',
                  data=None, datafn='oceFWflx', datavar='oceFWflx', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 12),
                  detrend=True, window=window, title='Correlation between Surface Freshwater Flux \n and Dotson to Cosgrove Mass Loss', ymax=-70, save=True, savefn='oceFWflx_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=False)

#Figure 11
createPACECorrelationMaps(save=True)


#Appendix
createMapsOfProcess()


#Extra statistics:
window=60
showCombinedCorrelationMap_2D_v3(members='all',
                  data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Correlation of Zonal Wind with Ice Shelf Basal Mass Loss (Dotson to Cosgrove)', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=False )

showCombinedCorrelationMap_2D_v3(members='all',
                  data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Meridional Wind', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=False )


showCombinedCorrelationMap_2D_v3(members='all',
                  data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Cumulative Zonal Wind', ymax=None, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=True)
    
showCombinedCorrelationMap_2D_v3(members='all',
                  data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
                  detrend=True, window=window, title='Cumulative Meridional Wind', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, mask='land_and_ice_mask',
                  give_back=False, deseasonalize=True, cumsum_map=True)


#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS
#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS
#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS
#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS
#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUS#%%----BONUSv

# window=60
# showCombinedCorrelationMap_2D_v3(members=[0,10],
#                   data=None, datafn='UVEL_averaged_-50to0', datavar='UVEL', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Average Surface Current (top 50m) vs Dotson to Cosgrove Melt', ymax=-70, save=True, savefn='surface_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True)


# showCombinedCorrelationMap_2D_v3(members=[0,10],
#                   data=None, datafn='UVEL_averaged_-50to0', datavar='UVEL', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Average Surface Current (top 50m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='surface_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True)




# window=60
# showMapPACE(members='all',
#                    data=None, datavar='PSL',
#                    ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 0),
#                    years=['1920', '2013'],
#                    detrend=True, deseasonalize=True, window=window, title='Correlation Map of Surface Air Pressure vs Amundsen Shelf Break Zonal Wind', ymax=-70, 
#                    save=True, savefn='PACE_PSL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                    draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                    inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False)
# window=60
# showMapPACE(members='all',
#                    data=None, datavar='SST',
#                    ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 0),
#                    years=['1920', '2013'],
#                    detrend=True, deseasonalize=True, window=window, title='Correlation Map of Sea Surface Temperature vs Amundsen Shelf Break Zonal Wind', ymax=-70, 
#                    save=True, savefn='PACE_SST_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                    draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                    inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False)

# window=60
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_forcing', indvar='amundsen_shelf_seaice_freeze', lags=range(-24, 24),
#                   detrend=True, window=window, title='Average Bottom Current (lowest 100m) vs Amundsen Shelf Sea Ice Freeze', ymax=-70, save=True, savefn='bottom100m_VEL_amundsen_seaice_freeze_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True)

# #%%BSL4
# window=60
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_forcing', indvar='amundsen_shelf_seaice_melt', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Average Bottom Current (lowest 100m) vs Amundsen Shelf Sea Ice Melt', ymax=-70, save=True, savefn='bottom100m_VEL_amundsen_seaice_melt_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True)

# window=60
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='SALT_averaged_-200to0', datavar='SALT', datakind='maps',
#                   ind=None, indfn='timeseries_forcing', indvar='amundsen_shelf_seaice_freeze', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Average Salinity (-200 to 0m) vs Amundsen Sea Ice Freeze', ymax=-70, save=True, 
#                   savefn='SALT-200to0_amundsen_shelf_sea_ice_freeze_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=False)


##--------------------------------------------------------------------------------------

# window=1
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='SALT_averaged_-200to0', datavar='SALT', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Average Salinity (-200 to 0m) vs Dotson to Cosgrove Massloss', ymax=-70, save=True, 
#                   savefn='SALT-200to0_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=False)


# window=1
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_forcing', indvar='amundsen_shelf_seaice_freeze', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Average Bottom Current (lowest 100m) vs Amundsen Shelf Sea Ice Freeze', ymax=-70, save=True, savefn='bottom100m_VEL_amundsen_seaice_freeze_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True)

# window=1
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='SALT_averaged_-200to0', datavar='SALT', datakind='maps',
#                   ind=None, indfn='timeseries_forcing', indvar='amundsen_shelf_seaice_freeze', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Average Salinity (-200 to 0m) vs Amundsen Sea Ice Freeze', ymax=-70, save=True, 
#                   savefn='SALT-200to0_amundsen_shelf_seaice_freeze_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=False)


# print('Start Map')
# window=60
# showCombinedCorrelationMap_2D_v3(members=[20,21],
#                   data=None, datafn='baroclinic_VEL', datavar='UVEL', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   mask='land', years=['1920', '2013'],
#                   detrend=True, deseasonalize=True, window=window, title='1D Correlation Map Ensemble Mean Baroclinic Flow and Dotson to Cosgrove Melt', ymax=-70, 
#                   save=True, savefn='baroclinic_VEL_dotson_to_cosgrove_massloss_ensmean_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False, save_results=True, give_back=False, cumsum_map=False)
# print('Finish')

# window=1
# showCombinedCorrelationMap_2D_v3(members=[20,21],
#                   data=None, datafn='baroclinic_VEL', datavar='UVEL', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   mask='land', years=['1920', '2013'],
#                   detrend=True, deseasonalize=True, window=window, title='1D Correlation Map Ensemble Mean Baroclinic Flow and Dotson to Cosgrove Melt', ymax=-70, 
#                   save=True, savefn='baroclinic_VEL_dotson_to_cosgrove_massloss_ensmean_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False, save_results=True, give_back=False, cumsum_map=False)

# window=60


# showCombinedCorrelationMap_2D_v3(members=[20,21],
#                   data=None, datafn='depth_averaged_VEL', datavar='UVEL', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   mask='land', years=['1920', '2013'],
#                   detrend=True, deseasonalize=True, window=window, title='1D Correlation Map Ensemble Mean Barotropic Flow and Dotson to Cosgrove Melt', ymax=-70, 
#                 save=True, savefn='barotropic_VEL_dotson_to_cosgrove_massloss_ensmean_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False, save_results=True, give_back=False, cumsum_map=False)
# window=1
# showCombinedCorrelationMap_2D_v3(members=[20,21],
#                   data=None, datafn='depth_averaged_VEL', datavar='UVEL', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   mask='land', years=['1920', '2013'],
#                   detrend=True, deseasonalize=True, window=window, title='1D Correlation Map Ensemble Mean Barotropic Flow and Dotson to Cosgrove Melt', ymax=-70, 
#                 save=True, savefn='barotropic_VEL_dotson_to_cosgrove_massloss_ensmean_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False, save_results=True, give_back=False, cumsum_map=False, cumsum_ind=False)


# window=60
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
#                   detrend=True, window=window, title='Correlation map of Meridional Wind vs Dotson to Cosgrove Mass Loss', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=False)

# window=60
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
#                   detrend=True, window=window, title='Correlation map of Cumulative Meridional Wind vs Dotson to Cosgrove Mass Loss', ymax=None, save=True, savefn='cumEXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=True)


# window=60
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
#                   detrend=True, window=window, title='Correlation map of Zonal Wind vs Dotson to Cosgrove Mass Loss', ymax=None, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=False)

# window=60
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
#                   detrend=True, window=window, title='Correlation map of Cumulative Zonal Wind vs Dotson to Cosgrove Mass Loss', ymax=None, save=True, savefn='cumEXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=True)


# window=60

# showCombinedCorrelationMap_2D_v3(members='all',
#                       data=None, datafn='ADVx_TH_corrected_bottom100m_averaged', datavar='ADVx_TH', datakind='old',
#                       ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                       detrend=True, window=window, title='Average Bottom Heat Advection (lowest 100m) vs Dotson to Cosgrove Mass Loss', ymax=-70, save=True, savefn='bottom100m_ADV_TH_dotson_to_cosgrove_'+str(window)+'month_dt',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False, deseasonalize=False,
#                       give_back=False)

# # showCombinedCorrelationMap_2D_v3(members='all',
# #                       data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
# #                       ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
# #                       detrend=True, window=window, title='Average Bottom Current (lowest 100m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='bottom100m_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
# #                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
# #                       inverse=False, deseasonalize=False,
# #                       give_back=False)

# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Temperature (200-700m) vs Dotson to Cosgrove Mass Loss', ymax=-70, save=True, savefn='amundsentemp_200-700_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False, deseasonalize=False,
#                   give_back=False)

# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='SALT_averaged_-200to0', datavar='SALT', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Average Salinity (-200 to 0m) vs Dotson to Cosgrove Mass Loss', ymax=-70, save=True, 
#                   savefn='SALT-200to0_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=False)

# #NOTICE: UP TO 12 months lead!
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='oceFWflx', datavar='oceFWflx', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 12),
#                   detrend=True, window=window, title='Correlation map of Surface Freshwater Flux vs Dotson to Cosgrove Mass Loss', ymax=-70, save=True, savefn='oceFWflx_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=False)



# showCombinedCorrelationMap_2D_v3(members=[0,5],
#                   data=None, datafn='UVEL_averaged_-50to0', datavar='UVEL', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Average Surface Current (top 50m) vs Dotson to Cosgrove Melt', ymax=-70, save=True, savefn='surface_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True)




# for window in [1, 60]:
#     showCombinedCorrelationMap_2D_v3(members='all',
#                       data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                       ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                       detrend=True, window=window, title='Average Bottom Current (lowest 100m) vs Dotson to Cosgrove Ice Shelf Basal Mass Loss', ymax=-70, save=True, savefn='bottom100m_VEL_dotson_to_cosgrove_'+str(window)+'month_dt_dS',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False, deseasonalize=True)


#     showCombinedCorrelationMap_2D_v3(members='all',
#                       data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
#                       ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                       detrend=True, window=window, title='Depth Averaged Current vs Dotson to Cosgrove Ice Shelf Basal Mass Loss', ymax=-70, save=True, savefn='depth_averaged_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False)

#     showCombinedCorrelationMap_2D_v3(members='all',
#                       data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                       ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                       detrend=True, window=window, title='Average Bottom Current (lowest 100m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='bottom100m_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False)

#     showCombinedCorrelationMap_2D_v3(members='all',
#                       data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
#                       ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                       detrend=True, window=window, title='Depth Averaged Current vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='depth_averaged_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                       draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                       inverse=False,
#                       give_back=False)


#     showCombinedCorrelationMap_2D_v3(members='all',
#                           data=None, datafn='baroclinic_UVEL', datavar='UVEL', datakind='old',
#                           ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                           detrend=True, window=window, title='Baroclinic Current (lowest 100m-depth averaged) vs Dotson to Cosgrove Ice Shelf Basal Mass Loss', ymax=-70, save=True, 
#                           savefn='baroclinic_VEL_dotson_to_cosgrove_'+str(window)+'month_dt',
#                           draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                           inverse=False,
#                           give_back=False)

#     showCombinedCorrelationMap_2D_v3(members='all',
#                           data=None,datafn='baroclinic_UVEL', datavar='UVEL', datakind='old',
#                           ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                           detrend=True, window=window, title='Baroclinic Current (lowest 100m-depth averaged) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='baroclinic_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                           draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                           inverse=False,
#                           give_back=False)

    
    

# window=60
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='curlTAU', datavar='curlTAU', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Surface Stress Curl vs Dotson to Cosgrove Ice Shelf Basal Mass Loss', ymax=-70, save=True, savefn='curlTAU_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Zonal Wind vs Dotson to Cosgrove Ice Shelf Basal Mass Loss', ymax=-70, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False)

    
# window=60

# showCombinedCorrelationMap_2D_v3(members=[20,21],
#                   data=data.expand_dims('ens'), datafn='new2_amundsen_shelf_temp_200to700', datavar='THETA', datakind='old',
#                   ind=wind.to_dataset(name='amundsen_shelf_break_uwind_avg').expand_dims('ens'), indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='Ensemble Mean Amundsen Shelf Temperature (200 to 700m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='ensmean_amundsen_shelf_temp_200to700_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, save_results=True)

# showCombinedCorrelationMap_2D_v3(members=[20,21],
#                   data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='amundsen_shelf_break_uwind_avg', lags=range(-24, 24),
#                   detrend=True, window=window, title='Ensemble Mean Average Bottom Current (lowest 100m) vs Amundsen Shelf Break Uwind', ymax=-70, save=True, savefn='ensmean_bottom100m_VEL_amundsen_shelf_break_uwind_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, save_results=True)    
    
    
# showCombinedCorrelationMap_2D_v3(members=[20,21],
#                   data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
#                   detrend=True, window=window, title='Ensemble Mean EXFuwind vs Dotson to Cosgrove Mass Loss', ymax=None, save=True, savefn='ensmean_EXFuwind_dotson_to_cosgrove_massloss_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, save_results=True)
    
    
    
# window=60
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Cumulative Zonal Wind vs Dotson to Cosgrove Massloss', ymax=None, save=True, savefn='cumEXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=True)
    
    
    
# window=60
# showCombinedCorrelationMap_2D_v3(members='all',
#                   data=None, datafn='SALT_averaged_-200to0', datavar='SALT', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Average Salinity (-200 to 0m) vs Dotson to Cosgrove Massloss', ymax=None, save=True, 
#                   savefn='SALT-200to0_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=False)

    
    
    
    
    
    
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