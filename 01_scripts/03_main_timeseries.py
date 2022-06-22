# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain test!.

@author: joren
"""

from functions.loading_and_processing_data import loadMITgcmData, loadPACEData, detrend_and_average_MITgcmData, read_all_data
from functions.composite_maps import createCompositeMap, createVelocityCompositeMap, defineEvents, findEvents, createCompositeMap_2, workflowCompositeMap
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

fn_map='EXFuwind'
var='EXFuwind'
data=loadMITgcmData(filename=fn_map, members='all', kind='maps')

def simpleCorrelationMap(fn_map='EXFuwind', mapvar='EXFuwind', fn_ind='undercurrent_W115'):
    #data=loadMITgcmData(filename=fn_map, members='all', kind='maps')

    UW115=xr.open_dataarray('/data/hpcdata/users/grejan/mitgcm/02_data/indices/'+fn_ind+'.nc')
    #UW115['time']=UW115.indexes['time'].to_datetimeindex()
    res=xr.corr(data[mapvar].mean(dim='ens'), UW115.mean(dim='ens'), dim='time')
    
    plt.figure(figsize=(20,10))
    plt.title('ENS MEAN')
    plt.pcolor( res.XC, res.YC, res, cmap='seismic', vmin=-1, vmax=1)
    plt.colorbar()
    CS=plt.contour( res.XC, res.YC, res, np.arange(-1, 1, 0.1), labels=np.arange(-1, 1, 0.1))
    plt.clabel(CS, CS.levels, inline=True,fontsize=10)
    plt.ylim([-75, -70])
    
    return

simpleCorrelationMap(fn_map='EXFuwind', mapvar='EXFuwind', fn_ind='undercurrent_W115')

dataU=xr.open_dataset('../02_data/maps/EXFuwind_ensmean.nc')
dataV=xr.open_dataset('../02_data/maps/EXFvwind_ensmean.nc')