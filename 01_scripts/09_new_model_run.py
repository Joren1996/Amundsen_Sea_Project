# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain test!.

@author: joren
"""

from functions.loading_and_processing_data import loadMITgcmData, loadPACEData, detrend_and_average_MITgcmData
from functions.composite_maps import createCompositeMap, createVelocityCompositeMap, defineEvents, findEvents, createCompositeMap_2, workflowCompositeMap
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
import pandas as pd
import datetime

import sys
sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
from mitgcm_python_master.file_io import netcdf_time

def read_all_data(full, total, units, longnames, kind):
    '''
    Old function to read data, this file still needs to be updated...
    INPUT:
    full (string): path where to find the time series
    total (dictionary): dictionary to contain the data.
    units (dictionary): dictionary to contain the units.
    longnames (dictionary): dictionary to contain the extended names of the variables.
    kind (string): suffix of timeseries to read
    
    OUTPUT:
    total (dictionary): dictionary with the data.
    units (dictionary): dictionary with the units.
    longnames (dictionary): dictionary with the extended names of the variables.
    '''
    
    file_paths=[os.path.join(full, file) for file in os.listdir(full) if kind in file]
    
    if kind!='trough':    
        names=['ens'+str(file[-5:-3]) for file in file_paths]
        
        varlist=list(xr.open_dataset(file_paths[1]).variables)

        for var in varlist:
            total[var]=pd.DataFrame()    
        for i, file in enumerate(file_paths):
            data=xr.open_dataset(file)

            timeind=netcdf_time(file, monthly=False)

            for var in varlist:
                total[var][names[i]]=pd.Series(data[var], index=timeind)
                if (i==0) & (var!='time'):
                    units[var]=data[var].units
                    longnames[var]=data[var].long_name
    
#     else:
#         varlist=['PITE', 'Burke', 'Dotson', 'BRE', 'BRW']
#         for var in varlist:
#             ts=xr.open_dataarray('./data/'+var+'_depth_averaged_vflow.nc')
#             ts=ts.drop('YG')
#             names=['ens'+str(i+1).zfill(2) for i in ts.ens.values]
#             total[var]=pd.DataFrame(data=ts.values.T, columns=names, index=ts.indexes['time'].to_datetimeindex())
#             units[var]='[m/s]'
#             longnames[var]='Depth Averaged Meridional Flow through Cross Section of '+var
            
#             varb=var+'B'
#             ts=xr.open_dataarray('./data/'+var+'_bottom100m_vflow.nc')
#             ts=ts.drop('YG')
#             names=['ens'+str(i+1).zfill(2) for i in ts.ens.values]
#             total[varb]=pd.DataFrame(data=ts.values.T, columns=names, index=ts.indexes['time'].to_datetimeindex())
#             units[varb]='[m/s]'
#             longnames[varb]='Bottom 100m Meridional Flow (depth averaged flow subtracted) through Cross Section of '+var
            
#             varb=var+'_full'
#             ts=xr.open_dataarray('./data/'+var+'_depth_averaged_vflow_full.nc')
#             ts=ts.drop('YG')
#             names=['ens'+str(i+1).zfill(2) for i in ts.ens.values]
#             total[varb]=pd.DataFrame(data=ts.values.T, columns=names, index=ts.indexes['time'].to_datetimeindex())
#             units[varb]='[m/s]'
#             longnames[varb]='Depth Averaged Meridional Flow through Cross Section of '+var+' (Not Detrended)'
            
#             varb=var+'B_full'
#             ts=xr.open_dataarray('./data/'+var+'_bottom100m_vflow_full.nc')
#             ts=ts.drop('YG')
#             names=['ens'+str(i+1).zfill(2) for i in ts.ens.values]
#             total[varb]=pd.DataFrame(data=ts.values.T, columns=names, index=ts.indexes['time'].to_datetimeindex())
#             units[varb]='[m/s]'
#             longnames[varb]='Bottom 100m Meridional Flow (depth averaged flow subtracted) through Cross Section of '+var
    return total, units, longnames


# def netcdf_time (file_path, var_name='time', t_start=None, t_end=None, return_date=True, monthly=True, return_units=False):

#     import netCDF4 as nc

#     # Open the file and get the length of the record
#     id = nc.Dataset(file_path, 'r')
#     time_id = id.variables[var_name]
#     units = time_id.units
#     try:
#         calendar = time_id.calendar
#     except(AttributeError):
#         calendar = 'standard'
#     num_time = time_id.size

#     # Choose range of time values to consider
#     # If t_start and/or t_end are already set, use those bounds
#     # Otherwise, start at the first time_index and/or end at the last time_index in the file
#     if t_start is None:
#         t_start = 0
#     if t_end is None:
#         t_end = num_time

#     # Read the variable
#     if return_date:
#         # Return as handy Date objects
#         time = nc.num2date(time_id[t_start:t_end], units=units, calendar=calendar)
#     else:
#         # Return just as scalar values
#         time = time_id[t_start:t_end]
#     id.close()

#     if return_date:
#         # Want to convert to a datetime object
#         if monthly:
#             # Back up to previous month
#             for t in range(time.size):
#                 month = time[t].month-1
#                 year = time[t].year
#                 if month < 1:
#                     month += 12
#                     year -= 1
#                 time[t] = datetime.datetime(year, month, 1)
#         else:
#             for t in range(time.size):
#                 time[t] = datetime.datetime(time[t].year, time[t].month, time[t].day)             

#     if return_units:
#         return time, units, calendar
#     else:
#         return time
    


#First, show which members were selected for the new run.
#Reading the data.
full='/data/hpcdata/users/grejan/mitgcm/data'
total={}
units={}
longnames={}
total, units, longnames=read_all_data(full, total, units, longnames, kind='final')

#Members that we want to show
wishlist=[11, 12, 6, 10, 4]

#Selection Criteria
print((total['dotson_to_cosgrove_massloss']['1940':'1941'].mean()-total['dotson_to_cosgrove_massloss']['1938':'1939'].mean()).sort_values())
print(total['dotson_to_cosgrove_massloss']['1970':'1990'].mean().sort_values())



#START PLOTTING THE SELECTED MEMBERS
plt.rcParams.update({'font.size': 18})
fig=plt.figure(figsize=(15,15))

#First, show high melt members
wishlist=[6, 10]#, 11, 12,  4]
wl=['ens'+str(i).zfill(2) for i in wishlist]
ax=plt.subplot(2,2,1)
#Plot ensemble mean
total['dotson_to_cosgrove_massloss']['1920':'2013'].mean(axis=1).rolling(60, center=True).mean().plot(c='k', linewidth=2, label='')
#Plot normal members
for i in range(20):
        total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i+1).zfill(2)].rolling(60, center=True).mean().plot(c='grey', linewidth=0.5, label='')
#Plot High melt members
cs=['red', 'orange']
for j, i in enumerate(wishlist):
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i).zfill(2)].rolling(60, center=True).mean().plot(linewidth=2, label='Member '+str(i), c=cs[j])
#Indicate time interval of interest.
data=total['dotson_to_cosgrove_massloss']
ax.axvspan(data['1940'].index[0], 
                   data['1941'].index[-1], color='red', alpha=0.3)
ax.axvspan(data['1938'].index[0], 
                   data['1939'].index[-1], color='blue', alpha=0.3)  
plt.grid(True)
plt.ylabel('Dotson to Cosgrove \n Basal Mass Loss [Gt/yr]')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.title('Two members with strong 1940')



#Second, show low melt members
wishlist=[11, 12]
wl=['ens'+str(i).zfill(2) for i in wishlist]
ax=plt.subplot(2,2,2)
#Plot ensemble mean
total['dotson_to_cosgrove_massloss']['1920':'2013'].mean(axis=1).rolling(60, center=True).mean().plot(c='k', linewidth=2, label='')
#Plot normal members
for i in range(20):
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i+1).zfill(2)].rolling(60, center=True).mean().plot(c='grey', linewidth=0.5, label='')
#Plot low melt members
cs=['green', 'blue']
for j, i in enumerate(wishlist):
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i).zfill(2)].rolling(60, center=True).mean().plot(linewidth=2, label='Member '+str(i), c=cs[j])
#Indicate time interval of interest.
data=total['dotson_to_cosgrove_massloss']
ax.axvspan(data['1940'].index[0], 
                   data['1941'].index[-1], color='red', alpha=0.3)
ax.axvspan(data['1938'].index[0], 
                   data['1939'].index[-1], color='blue', alpha=0.3)
plt.grid(True)
plt.title('Two members with weak 1940')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel('Dotson to Cosgrove \n Basal Mass Loss [Gt/yr]')


#Third, show balancing member
wishlist=[11, 12, 6, 10]
wl=['ens'+str(i).zfill(2) for i in wishlist]
ax=plt.subplot(2,2,3)
plt.title('Balancing member (least melt 1970-1990)')
#Plot ensemble mean
total['dotson_to_cosgrove_massloss']['1920':'2013'].mean(axis=1).rolling(60, center=True).mean().plot(c='k', linewidth=2, label='')
#Plot normal members
for i in range(20):
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i+1).zfill(2)].rolling(60, center=True).mean().plot(c='grey', linewidth=0.5, label='')
#Plot balancing member
for i in [4]:
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i).zfill(2)].rolling(60, center=True).mean().plot(linewidth=2, label='Member '+str(i), c='brown')
#Plot selected mean - balancing member
total['dotson_to_cosgrove_massloss']['1920':'2013'][wl].mean(axis=1).rolling(60, center=True).mean().plot(c='m', linewidth=2, linestyle='--', label='Mean of selected (w/o 4)')
plt.grid(True)
ax.axvspan(data['1970'].index[0], 
                   data['1990'].index[-1], color='green', alpha=0.3)  
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel('Dotson to Cosgrove \n Basal Mass Loss [Gt/yr]')


#Fourth, show selected members
wishlist=[11, 12, 6, 10, 4]
wl=['ens'+str(i).zfill(2) for i in wishlist]
ax=plt.subplot(2,2,4)
plt.title('Mean')
#Plot ensemble mean
total['dotson_to_cosgrove_massloss']['1920':'2013'].mean(axis=1).rolling(60, center=True).mean().plot(c='k', linewidth=2, label='')
#Plot normal members
for i in range(20):
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i+1).zfill(2)].rolling(60, center=True).mean().plot(c='grey', linewidth=0.5, label='')
#Plot mean of selected members
total['dotson_to_cosgrove_massloss']['1920':'2013'][wl].mean(axis=1).rolling(60, center=True).mean().plot(c='m', linewidth=2,label='Mean of selected')
plt.grid(True)
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel('Dotson to Cosgrove \n Basal Mass Loss [Gt/yr]')

#Create the legend manually.....
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
black_line=Line2D([0], [0], color='black', lw=2, label='Ensemble Mean')
grey_line=Line2D([0], [0], color='grey', lw=0.5, label='Member (not selected)')
red_patch = mpatches.Patch(facecolor='red', edgecolor='black', alpha=0.3, label='1940 and 1941')
blue_patch = mpatches.Patch(facecolor='blue', edgecolor='black', alpha=0.3, label='1938 and 1939')
green_patch = mpatches.Patch(facecolor='green', edgecolor='black', alpha=0.3, label='1970 up to and including 1990')

#Add letters to subplot
import string
for n, ax in enumerate(fig.axes):
    ax.text(-0.25, 0.9, string.ascii_uppercase[n]+'.', transform=ax.transAxes, 
            size=25, weight='bold')
#Create legend
lgd=fig.legend(handles=[black_line, grey_line, red_patch, blue_patch, green_patch], ncol=3,
           loc='upper center', bbox_to_anchor=(0.5, 0.15))

#Add suptitle
st=plt.suptitle('Selected Members for Model Runs with Advanced Grounding Line')

#Save
plt.subplots_adjust(hspace=0.3, wspace=0.5, bottom=0.2)
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/new_run/selection_60months.png', bbox_extra_artists=(lgd, st), bbox_inches='tight')
