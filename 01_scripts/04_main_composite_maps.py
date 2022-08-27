# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should eventually contain the code to create all composite maps for the final report.

@author: joren
"""

from functions.loading_and_processing_data import loadMITgcmData, detrend_and_average_MITgcmData, read_all_data
from functions.composite_maps import createCompositeMap, createVelocityCompositeMap, defineEvents, findEvents, createCompositeMap_2, workflowCompositeMap
import numpy as np
import xarray as xr

#00-INITIALISATION------------------------------------------------------------

#First, read time series so that we can define events.
full='/data/hpcdata/users/grejan/mitgcm/data'
total={}
units={}
longnames={}
total, units, longnames=read_all_data(full, total, units, longnames, kind='isotherm')
total, units, longnames=read_all_data(full, total, units, longnames, kind='final')
total, units, longnames=read_all_data(full, total, units, longnames, kind='ohc')

#Select events that we are interested in.
threshold=defineEvents(total=total, longnames=longnames, units=units, wvar='dotson_to_cosgrove_massloss', larger_than=True, pct=75, detrend=True, window=12*2, filename='melt75_colored_in')
print(threshold)
selected,_=findEvents(total=total, longnames=longnames, units=units, threshold=threshold, wvar='dotson_to_cosgrove_massloss', larger_than=True, detrend=True, window=12*2, filename='melt75_colored_in')

# #Select events that we are interested in.
# threshold=defineEvents(total=total, longnames=longnames, units=units, wvar='amundsen_shelf_ohc_below_0m', larger_than=False, pct=25, detrend=True, window=12*2, filename='cold25_colored_in')
# print(threshold)
# selected2,_=findEvents(total=total, longnames=longnames, units=units, threshold=threshold, wvar='amundsen_shelf_ohc_below_0m', larger_than=False, detrend=True, window=12*2, filename='cold25_colored_in')

#01-Normal-Composite-Map-------------------------------------------------------------------

#selected=selected[['ens01', 'ens02', 'ens03']]

#%% INITIALISATION
members='all' #[0,1]
var='SALT'
input_filename='SALT_averaged_-700to-200'
output_filename='test_SALT_-700to-200_highmelt'
title='SALT (averaged between 700 and 200m) during high melt event'
window=1
larger_than=True
deseasonalize=True
pct=75
workflowCompositeMap(input_filename=input_filename, var=var, 
                     members=members, deseasonalize=deseasonalize, output_filename=output_filename, 
                     title=title, window=window, selected=selected)



members=[0,3] #[0,1]
var='oceFWflx'
input_filename='oceFWflx'
output_filename='test_oceFWflx_highmelt'
title='oceFWflx during high melt event'
window=12
larger_than=True
deseasonalize=True
pct=75
workflowCompositeMap(input_filename=input_filename, var=var, 
                     members=members, deseasonalize=deseasonalize, output_filename=output_filename, 
                     title=title, window=window, selected=selected)

#02-Overview-Composites------------------------------------------------

















































# def workflowOverview(members, deseasonalize, output_filename='', title='', window=1):
#     #Load Data that we already have
#     data={}
#     data_dt={}

#     #var_list=['EXFuwind', 'THETA', 'UVEL', 'VVEL', 'THETAW106', 'THETAS72-5', 'THETAW120', 
#     #          'THETAW100', 'THETAW103-5', 'THETAPITE', 'UVELB', 'VVELB']
#     var_list=['EXFuwind', 'THETA', 'THETAPITE', 'UVELW115', 'UVEL', 'VVEL', 'UVELB', 'VVELB']

#     filenames={'EXFuwind':'zonal_winds', 
#                'THETA':'new2_amundsen_shelf_temp_200to700', 
#                'UVEL':'depth_averaged_UVEL', 
#                'VVEL':'depth_averaged_VVEL',
#                'curlTAU':'curl_tau', 
#                'THETAW106':'theta_W106', 
#                'THETAW120':'theta_W120',
#                'UVELW115':'UVEL_W115',
#                'THETAW103-5':'theta_W103-5',
#                'THETAW100':'theta_W100',
#                'THETAS72-5':'theta_S72-5',
#                'THETAPITE':['THETA_transect_PITE_full_0', 'THETA_transect_PITE_full_1', 
#                             'THETA_transect_PITE_full_2', 'THETA_transect_PITE_full_3', 
#                             'THETA_transect_PITE_full_4'],
#                'UVELB':'bottom100m_averaged_UVEL', 
#                'VVELB':'bottom100m_averaged_VVEL'}

#     for var in var_list: 
#         print('Starting with: '+var)
#         if var=='THETAPITE':
#             for fn in filenames[var]:
#                 data[var+'_'+fn[-1]]=fc.loadMITgcmData1(filename=fn, members=members)
#         else:
#             data[var]=fc.loadMITgcmData(filename=filenames[var], members=members)
#         if 'VELB' in var:
#             if 'U' in var:
#                 data[var]=data[var]-data['UVEL']
#                 data_dt[var]=fc.detrend_and_average_MITgcmData(data[var], 'UVEL', window=window)
#             elif 'VVEL' in var:
#                 data[var]=data[var]-data['VVEL']
#                 data_dt[var]=fc.detrend_and_average_MITgcmData(data[var], 'VVEL', window=window)

#         elif 'THETA' in var:
#             if 'PITE' in var:
#                 for fn in filenames[var]:
#                     data_dt[var+'_'+fn[-1]]=fc.detrend_and_average_MITgcmData(data[var+'_'+fn[-1]], 'THETA', window=window)
#             else:
#                 data_dt[var]=fc.detrend_and_average_MITgcmData(data[var], 'THETA', window=window)
#         else:
#             data_dt[var]=fc.detrend_and_average_MITgcmData(data[var], var, window=window)
    
    
    
    
    
    
# #     data=loadMITgcmData(filename=input_filename, members=members)
    
# #     #Composite Maps
# #     data_dt=detrend_and_average_MITgcmData(data, var, window=window)
# #     if deseasonalize==True:
# #              data_dt = data_dt.groupby("time.month")-data_dt.groupby("time.month").mean(dim="time")
# #     createCompositeMap(data_dt, selected, title=title, 
# #                        filename=output_filename,
# #                        window=window)
#     return


# members=[0,1]
# deseasonalize=True
# output_filename='test_overview'


# def workflowCompositeMap_vectorfield(input_filename_list, var_list, members, deseasonalize, output_filename='', title='', window=1, unit='[m/s]', vmax=None):
#     #Loading the data.
#     dataU=loadMITgcmData(filename=input_filename_list[0], members=members)
#     dataV=loadMITgcmData(filename=input_filename_list[1], members=members)
#     #Composite Maps
#     print(dataU)
#     data_dtU=detrend_and_average_MITgcmData(dataU, var_list[0], window=window)
#     data_dtV=detrend_and_average_MITgcmData(dataV, var_list[1], window=window)
    
#     if deseasonalize==True:
#             data_dtU = data_dtU.groupby("time.month")-data_dtU.groupby("time.month").mean(dim="time")
#             data_dtV = data_dtV.groupby("time.month")-data_dtV.groupby("time.month").mean(dim="time")
            
#     createVelocityCompositeMap(dataU=data_dtU, dataV=data_dtV, selected=selected, title=title, 
#                        filename=output_filename, window=window, 
#                        unit=unit, members=members, vmax=vmax)
#     return

# def workflowCompositeMap_vectorfield_2(input_filename_list, var_list, members, deseasonalize, output_filename='', output_filename2='', title='', title2='', window=1, unit='[K*m/s]', vmax=None):
#     #Loading the data.
#     dataU=loadMITgcmData(filename=input_filename_list[0], members=members)
#     dataV=loadMITgcmData(filename=input_filename_list[1], members=members)
#     #Composite Maps
    
    
#     vfield=np.sqrt(dataU[var_list[0]].values**2+dataV[var_list[-1]].values**2)
#     data = xr.DataArray(vfield, coords=[dataU[var_list[0]].ens, dataU[var_list[0]].time, dataU[var_list[0]].YC, dataV[var_list[-1]].XC], dims=[ "ens", "time", "YC", "XC"])
#     data=data.to_dataset(name=var_list[0])
    
#     #Composite Maps
#     data_dt=detrend_and_average_MITgcmData(data, var=var_list[0], window=window)
#     if deseasonalize==True:
#              data_dt = data_dt.groupby("time.month")-data_dt.groupby("time.month").mean(dim="time")
#     createCompositeMap_2(data_dt, selected, title=title, 
#                        filename=output_filename,
#                        window=window)
    
    
#     createCompositeMap_2(data_dt, selected2, title=title2, 
#                        filename=output_filename2,
#                        window=window)
    
# #     #print(dataU)
# #     #data_dtU=detrend_and_average_MITgcmData(dataU, var_list[0], window=window)
# #     #data_dtV=detrend_and_average_MITgcmData(dataV, var_list[1], window=window)
    
# #     if deseasonalize==True:
# #             data_dtU = data_dtU.groupby("time.month")-data_dtU.groupby("time.month").mean(dim="time")
# #             data_dtV = data_dtV.groupby("time.month")-data_dtV.groupby("time.month").mean(dim="time")
            
# #     createVelocityCompositeMap(dataU=data_dtU, dataV=data_dtV, selected=selected, title=title, 
# #                        filename=output_filename, window=window, 
# #                        unit=unit, members=members, vmax=vmax)
#     return









# # for i in range(2, 21):
# #     selected['ens'+str(i).zfill(2)]=False

# #%% HEAT ADVECTION PLOT
# members='all'
# var_list=['ADVx_TH', 'ADVy_TH']
# #input_filename_list=['ADVx_TH_corrected_depth_averaged', 'ADVy_TH_corrected_depth_averaged']
# input_filename_list=['ADVx_TH_corrected_bottom100m_averaged', 'ADVy_TH_corrected_bottom100m_averaged']
# output_filename='corrected_bottom100m_averaged_advection_warm_test'
# title='Composite Map of Bottom 100m Averaged Heat Advection during Warm Episodes'


# output_filename2='corrected_bottom100m_averaged_advection_cold_test'
# title2='Composite Map of Bottom 100m Averaged Heat Advection during Cold Episodes'

# window=1
# deseasonalize=True
# workflowCompositeMap_vectorfield_2(input_filename_list=input_filename_list, var_list=var_list, 
#                                  members=members,  deseasonalize=deseasonalize, output_filename=output_filename, output_filename2=output_filename2,
#                                  title=title, title2=title2, window=window, unit='[K*m/s]', vmax=10000)




# members='all'
# var_list=['ADVx_TH', 'ADVy_TH']
# input_filename_list=['ADVx_TH_corrected_bottom100m_averaged', 'ADVy_TH_corrected_bottom100m_averaged']
# output_filename='corrected_bottom100m_averaged_advection_warm'
# title='Composite Map of Bottom 100m Averaged Heat Advection during Warm Episodes'
# window=1
# deseasonalize=True
# workflowCompositeMap_vectorfield(input_filename_list=input_filename_list, var_list=var_list, 
#                                  members=members,  deseasonalize=deseasonalize, output_filename=output_filename, 
#                                  title=title, window=window, unit='[K*m/s]', vmax=10000)







