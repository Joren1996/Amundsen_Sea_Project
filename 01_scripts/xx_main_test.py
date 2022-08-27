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



def showCombinedCorrelationMap_2D_v3(members='all',
                      data=None, datafn='depth_averaged_UVEL', datavar='UVEL', datakind='old',
                      ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                      detrend=True, deseasonalize=True, window=24, title='2D Correlation Map', ymax=-70, save=False, savefn='test',
                      draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                      mask='land', years=['1920', '2013'],
                      inverse=False, save_results=True,
                      give_back=False, cumsum_map=False):
    '''
    Input: Always give U component as input!
    
    '''
    from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    print('Start preparing for correlation map!')
    print('NOTE: this is version 3!')
    plt.rcParams.update({'font.size': 22})
    
    if members=='all':
        members=[0,20]
    
    for i, ens in enumerate(range(members[0]+1,members[-1]+1)):
        if deseasonalize==True:
            if cumsum_map==True:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxp=xr.open_dataarray('../02_data/correlation/maxp_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_cum'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            else:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxp=xr.open_dataarray('../02_data/correlation/maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
        else:
            try:
                temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
            except:
                try:
                    temp_opt_lag=np.load('../02_data/correlation/opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxp=xr.open_dataarray('../02_data/correlation/maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                except:
                    temp_opt_lag=np.load('../02_data/correlation/2D_opt_lag_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.npy')
                    temp_maxp=xr.open_dataarray('../02_data/correlation/2D_maxp_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    temp_maxcor=xr.open_dataarray('../02_data/correlation/2D_maxcor_'+datafn+'_'+indvar+'_dt'+str(detrend)+'_dS'+str(deseasonalize)+'_w'+str(window)+'_ens'+str(ens)+'.nc')
                    
        if i==0:
            opt_lag=temp_opt_lag
            maxp=temp_maxp
            maxcor=temp_maxcor
        else:
            opt_lag=np.dstack((opt_lag, temp_opt_lag))
            maxp=xr.concat([maxp, temp_maxp], dim='ens')
            maxcor=xr.concat([maxcor, temp_maxcor], dim='ens')
    if i==0:
        maxp=maxp.expand_dims('ens')
        maxcor=maxcor.expand_dims('ens')
    if members[-1]!=21:
        std_lag=np.nanstd(opt_lag, axis=2)
        opt_lag=np.nanmean(opt_lag, axis=2)
    else:
        std_lag=np.zeros(np.shape(opt_lag))
    
    stdp=maxp.mean(dim='ens')
    maxp=maxp.mean(dim='ens')
    
    stdcor=maxcor.std(dim='ens')
    maxcor=maxcor.mean(dim='ens')
    
    #Reading Grid
    gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
    grid = Grid(gp)
    
    print('...Plotting')
    fig, ax=plt.subplots(2,1, figsize=(20,14))
    
    if members[-1]==21:
        fig.suptitle(title+'\n Ensemble Mean; '+str(window)+' month rolling mean; p<0.05')
    else:
        fig.suptitle(title+'\n Mean over members '+ str(members[0]+1)+' to '+ str(members[-1])+'; '+str(window)+' month rolling mean; p<0.05')


    #PLOT CORRELATION COEFFICIENTS
    img=ax[0].pcolor(maxcor.XC, maxcor.YC, maxcor, cmap='seismic', vmin=-1, vmax=1)
    plt.colorbar(img, ax=ax[0], label='Max Pearson Correlation Coefficient')
    CS=ax[0].contour(maxcor.XC, maxcor.YC, maxcor, [-0.5, 0.5], colors='yellow')
    ax[0].clabel(CS, CS.levels, inline=True,fontsize=10)

    bath=ax[0].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax[0].set_ylim([-75.5, ymax])
    if 'ice' in mask:
        ax[0].pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
    ax[0].set_facecolor('grey')
    ax[0].set_xlabel('Longitude [deg]')
    ax[0].set_ylabel('Latitude [deg]')
    ax[0].set_title('Optimal Correlation Coefficient')

    if draw_box==True:
        ax[0].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')

    #PLOT OPTIMAL LAG
    img=ax[1].pcolor(maxcor.XC, maxcor.YC, opt_lag, cmap='seismic', vmin=lags[0], vmax=lags[-1])
    plt.colorbar(img, ax=ax[1], label='Optimal Lag [months]')
    img=ax[1].pcolor(maxp.where(maxp>0.05).XC, maxp.where(maxp>0.05).YC, maxp.where(maxp>0.05)*0, cmap='ocean')

    bath=ax[1].contour(grid.lon_2d+360, grid.lat_2d, grid.bathy, levels=[-1000], colors=['black'])
    plt.clabel(bath, [-1000])
    if type(ymax)!=type(None):
        ax[1].set_ylim([-75.5, ymax])

    if 'ice' in mask:
        ax[1].pcolor(grid.lon_2d+360, grid.lat_2d, np.where(grid.ice_mask, np.ones(np.shape(grid.ice_mask)), np.nan*np.ones(np.shape(grid.ice_mask))), cmap='cool')
    ax[1].set_facecolor('grey')
    ax[1].set_xlabel('Longitude [deg]')
    ax[1].set_ylabel('Latitude [deg]')

    if draw_box==True:
        ax[1].plot([box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], 'o-k')
    ax[1].set_title('Optimal Lag \n (positive = timeseries leading)')

    fig.subplots_adjust(hspace=0.5, top=0.85)

    if save==True:
        print('...Saving figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_2Dcormap_'+savefn+'_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
    fig=plt.figure(figsize=(20,10))
    plt.suptitle('Extra statistics over ensemble members for : \n'+title)
    plt.subplot(2,2,1)
    img=plt.pcolor(maxp.XC, maxp.YC, maxp, vmax=0.05, cmap='seismic')
    plt.colorbar(img, label='Pvalue')
    #CS=plt.contour(maxp.XC, maxp.YC, maxp, [1e-4], colors='yellow')
    #plt.clabel(CS, CS.levels, inline=True,fontsize=10)
    plt.title('Mean P value with highest lag')# \n Yellow line indicates 1e-4')
    
    plt.subplot(2,2,2)
    img=plt.pcolor(stdp.XC, stdp.YC, stdp, vmax=0.05, cmap='seismic')
    plt.colorbar(img, label='Std of Pvalue')
    #CS=plt.contour(stdp.XC, stdp.YC, stdp, [1e-4], colors='yellow')
    #plt.clabel(CS, CS.levels, inline=True,fontsize=10)
    plt.title('Standard deviation of P value')# \n Yellow line indicates 1e-4')
    
    plt.subplot(2,2,3)
    img=plt.pcolor(stdcor.XC, stdcor.YC, stdcor, cmap='seismic', vmax=0.2)
    plt.colorbar(img, label='Std of cor')
    #CS=plt.contour(stdcor.XC, stdcor.YC, stdcor, [0.1], colors='yellow')
    #plt.clabel(CS, CS.levels, inline=True,fontsize=10)
    plt.title('Standard deviation of correlation coefficient')#  \n Yellow line indicates 0.1')
    
    plt.subplot(2,2,4)
    img=plt.pcolor(stdcor.XC, stdcor.YC, std_lag, cmap='seismic', vmax=10)
    plt.colorbar(img, label='Std of lag')
    #CS=plt.contour(stdcor.XC, stdcor.YC, std_lag, [4], colors='yellow')
    #plt.clabel(CS, CS.levels, inline=True,fontsize=10)
    plt.title('Standard deviation of lag')# \n Yellow line indicates 4 months')
    
    fig.subplots_adjust(hspace=0.5, wspace=0.3, top=0.8)
    
    if save==True:
        print('...Saving second figure')
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        fig.savefig('../03_output/correlation/'+today+'_BONUS_2Dcormap_'+savefn+'_members'+str(members[0])+'-'+str(members[-1])+'.png')
    
    return


window=60
showCombinedCorrelationMap_2D_v3(members='all',
                  data=None, datafn='EXFuwind', datavar='EXFuwind', datakind='maps',
                  ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
                  detrend=True, window=window, title='Zonal Wind vs Dotson to Cosgrove Ice Shelf Basal Mass Loss', ymax=-70, save=False, savefn='EXFuwind_dotson_to_cosgrove_'+str(window)+'month_dt',
                  draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
                  inverse=False, cumsum_map=True,
                  give_back=False)






# window=1
# showcorrelationMap_2D_v3(members='all',
#                   data=None, datafn='bottom100m_averaged_UVEL', datavar='UVEL', datakind='old',
#                   ind=None, indfn='timeseries_final', indvar='dotson_to_cosgrove_massloss', lags=range(-24, 24),
#                   detrend=True, window=window, title='2D Correlation Map of Average Bottom Current (lowest 100m) vs Dotson to Cosgrove Melt (1 month rolling mean)', ymax=-70, save=True, savefn='bottom100m_VEL_dotson_to_cosgrove_'+str(window)+'month_dt_dS',
#                   draw_box=False, box_x=[360-115, 360-102],box_y=[-71.8, -70.2],
#                   inverse=False,
#                   give_back=False, deseasonalize=True)

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