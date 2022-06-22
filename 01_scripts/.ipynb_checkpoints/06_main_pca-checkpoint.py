# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain all important figures that were created with PCA.

@author: joren
"""
import xarray as xr
import numpy as np
from functions.loading_and_processing_data import loadMITgcmData, detrend_and_average_MITgcmData
import functions.pca as fc


#%%CREATE ANIMATIONS OF DEPTH INTEGRATED THETA
# input_filename='THETA_depth_integrated'
# var='THETA'
# members='all'
# window=1
# deseasonalize=True
# method='ens_mean'
# pcnumbers=10
# start_year=[1936, 1966, 1991]
# output_filename='integrated_THETA' 
# frames=48

#%%CREATE ANIMATIONS OF Curl TAU
input_filename='curlTAU'
var='curlTAU'
members=[0,1]
window=24
deseasonalize=True
method='ens_mean'
pcnumbers=10
start_year=[1936, 1966, 1991]
output_filename='curlTAU' 
frames=48


def workflowPCA(input_filename, var, members, window, deseasonalize, method, pcnumbers, start_year=[1938, 1968, 1993], output_filename='', frames=36, maxy=-70, animate=True):

    if method=='ens_mean':
        try:
            print('...Directly loading '+input_filename+' ens_mean' )
            data=xr.open_dataset('../02_data/maps/'+input_filename+'_ensmean.nc')
            if 'ens' not in data.dims:
                data=data.expand_dims('ens')
        except:
            print('...or not :p' )
            data=loadMITgcmData(filename=input_filename, members=members)
    else:
        print('....Normally loading data')
        data=loadMITgcmData(filename=input_filename, members=members)
     
    print('...Detrending data' )
    data_dt=detrend_and_average_MITgcmData(data, var, window=window)
    if deseasonalize==True:
         data_dt = data_dt.groupby("time.month")-data_dt.groupby("time.month").mean(dim="time")
    
    if 'YC' in data_dt.dims:
        data_dt=data_dt.where(data_dt.YC<maxy, drop=True)
    else:
        data_dt=data_dt.where(data_dt.YG<maxy, drop=True)
    
    print('...Perform PCA')
    EOF_recons, PCs, tc=fc.performPCA(data=data_dt, var=var, method=method, interest=pcnumbers, filepath='../03_output/pca/', filename=output_filename)
    if animate==True:
        for year in start_year:
            print('...Create Animation for '+str(year))
            fc.makePCAAnimation(data_dt=data_dt, EOF_recons=EOF_recons, PCs=PCs, tc=tc, var=var, pcnumber=pcnumbers, frames=frames, start_year=year, filename=output_filename+'_'+str(year), filepath='../03_output/pca/')
    print('Finished!')
    return

workflowPCA(input_filename, var, members, window, deseasonalize, method, pcnumbers, start_year, output_filename, frames, animate=False)
