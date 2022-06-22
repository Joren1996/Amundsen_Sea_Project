# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain all functions to read data from the output of MITgcm to smaller files that are useful for the analysis.

@author: joren
"""

import xarray as xr
import numpy as np
import pandas as pd

def computeBuoyancy(member=1):
    print('...computing Buoyancy Forcing using EOS for member '+str(member))
    g = 9.81 #m/s^2 is gravity
    tAlpha = 3.733e-5 #K^-1 is thermal expansion coefficient
    Cp = 4180 #J/K/kg is specific heat of seawater
    sBeta = 7.843e-4 #psu^-1 is haline contraction coefficient
    #SSS is surface salinity
    
    oceFWflx=xr.open_dataset('../02_data/maps/oceFWflx_ens'+str(member)+'.nc')
    oceQnet=xr.open_dataset('../02_data/maps/oceQnet_ens'+str(member)+'.nc')
    SSS=xr.open_dataset('../02_data/maps/Z5_SALT_ens'+str(member)+'.nc')
    
    oceBUOY= g * (tAlpha/Cp * oceQnet['oceQnet'] + sBeta * SSS['SALT'] * oceFWflx['oceFWflx'])

    return oceBUOY.to_dataset(name='oceBUOY').expand_dims('ens')