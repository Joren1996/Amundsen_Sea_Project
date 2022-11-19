'''
Perform some correlations for checking if convection influenced the results.



'''

from functions.convection import convectionCorMap, showConvectionCorMap, createConvectionMaps

#Perform correlations without convection.

#-----SALT 0-200m----------------------------------------------------
# window=60
# showConvectionCorMap(members='all',
#                   data=None, datafn='SALT_averaged_-200to0', datavar='SALT', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', 
#                   lags=range(-24, 24),
#                   detrend=True, window=window, 
#                   title='Correlation map of Average Salinity (-200 to 0m) vs \
#                   Dotson to Cosgrove Massloss', ymax=-70, save=True, 
#                   savefn='SALT-200to0_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=False)

#-----oceFWflx----------------------------------------------------
# window=60
# showConvectionCorMap(members='all',
#                   data=None, datafn='oceFWflx', datavar='oceFWflx', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='Correlation map of Surface Freshwater Flux vs Dotson to Cosgrove Massloss', ymax=None, save=True, savefn='oceFWflx_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=False)

#-----EXFuwind----------------------------------------------------
# window=60
# showConvectionCorMap(members='all',
#                   data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
#                   detrend=True, window=window, title='Correlation map of Zonal Wind vs Dotson to Cosgrove Mass Loss', ymax=None, save=True, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=False)

#-----EXFvwind----------------------------------------------------
# window=60
# showConvectionCorMap(members='all',
#                   data=None, datafn='EXFvwind', datavar='EXFvwind', datakind='maps',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 0),
#                   detrend=True, window=window, title='Correlation map of Meridional Wind vs Dotson to Cosgrove Mass Loss', ymax=None, save=True, savefn='EXFvwind_dotson_to_cosgrove_'+str(window)+'month_dt',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True, cumsum_map=False)

#--------------------------Create Final Map for in the report.--------------------------
#createConvectionMaps(save=True)