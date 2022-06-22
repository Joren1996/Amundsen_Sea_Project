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


#First, show which members were selected for the new run.
#Reading the data.
full='/data/hpcdata/users/grejan/mitgcm/data'
total={}
units={}
longnames={}
total, units, longnames=read_all_data(full, total, units, longnames, kind='final')

wishlist=[11, 12, 6, 10, 4]

#Selection Criteria
print((total['dotson_to_cosgrove_massloss']['1940':'1941'].mean()-total['dotson_to_cosgrove_massloss']['1938':'1939'].mean()).sort_values())
print(total['dotson_to_cosgrove_massloss']['1970':'1990'].mean().sort_values())


plt.figure(figsize=(15,15))
plt.suptitle('Selected Members for New Model Run')
wishlist=[6, 10]#, 11, 12,  4]

plt.subplot(4,2,1)
plt.title('Two members with strong 1940')
total['amundsen_shelf_break_uwind_avg']['1920':'2013'].mean(axis=1).rolling(24, center=True).mean().plot(c='k', linewidth=1.2, label='ensmean')
for i in range(20):
    total['amundsen_shelf_break_uwind_avg']['1920':'2013']['ens'+str(i+1).zfill(2)].rolling(24, center=True).mean().plot(c='grey', linewidth=0.1, label='')
for i in wishlist:
    total['amundsen_shelf_break_uwind_avg']['1920':'2013']['ens'+str(i).zfill(2)].rolling(24, center=True).mean().plot(linewidth=1, label=str(i))

# wl=['ens'+str(i).zfill(2) for i in wishlist]
# total['amundsen_shelf_break_uwind_avg']['1920':'2013'][wl].mean(axis=1).rolling(24, center=True).mean().plot(c='b', linewidth=1.2, label='Mean of selected')   

plt.grid(True)
#plt.legend()
plt.ylabel('amundsen_shelf_break_uwind_avg')
plt.subplot(4,2,2)
total['dotson_to_cosgrove_massloss']['1920':'2013'].mean(axis=1).rolling(24, center=True).mean().plot(c='k', linewidth=1.2, label='ensmean')
for i in range(20):
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i+1).zfill(2)].rolling(24, center=True).mean().plot(c='grey', linewidth=0.1, label='')
for i in wishlist:
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i).zfill(2)].rolling(24, center=True).mean().plot(linewidth=1, label=str(i))
# total['dotson_to_cosgrove_massloss']['1920':'2013'][wl].mean(axis=1).rolling(24, center=True).mean().plot(c='b', linewidth=1.2,label='Mean of selected')
plt.grid(True)
plt.ylabel('dotson_to_cosgrove_massloss')
plt.legend()


plt.subplot(4,2,3)
plt.title('Two members with weak 1940')
wishlist=[11, 12]
total['amundsen_shelf_break_uwind_avg']['1920':'2013'].mean(axis=1).rolling(24, center=True).mean().plot(c='k', linewidth=1.2, label='ensmean')
for i in range(20):
    total['amundsen_shelf_break_uwind_avg']['1920':'2013']['ens'+str(i+1).zfill(2)].rolling(24, center=True).mean().plot(c='grey', linewidth=0.1, label='')
for i in wishlist:
    total['amundsen_shelf_break_uwind_avg']['1920':'2013']['ens'+str(i).zfill(2)].rolling(24, center=True).mean().plot(linewidth=1, label=str(i))

# wl=['ens'+str(i).zfill(2) for i in wishlist]
# total['amundsen_shelf_break_uwind_avg']['1920':'2013'][wl].mean(axis=1).rolling(24, center=True).mean().plot(c='b', linewidth=1.2, label='Mean of selected')   

plt.grid(True)
#plt.legend()
#plt.ylabel('amundsen_shelf_break_uwind_avg')
plt.subplot(4,2,4)
total['dotson_to_cosgrove_massloss']['1920':'2013'].mean(axis=1).rolling(24, center=True).mean().plot(c='k', linewidth=1.2, label='ensmean')
for i in range(20):
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i+1).zfill(2)].rolling(24, center=True).mean().plot(c='grey', linewidth=0.1, label='')
for i in wishlist:
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i).zfill(2)].rolling(24, center=True).mean().plot(linewidth=1, label=str(i))
# total['dotson_to_cosgrove_massloss']['1920':'2013'][wl].mean(axis=1).rolling(24, center=True).mean().plot(c='b', linewidth=1.2,label='Mean of selected')
plt.grid(True)
#plt.ylabel('dotson_to_cosgrove_massloss')
plt.legend()


plt.subplot(4,2,5)
plt.title('Balancing Member (least melt 1970-1990)')
wishlist=[11, 12, 6, 10]
total['amundsen_shelf_break_uwind_avg']['1920':'2013'].mean(axis=1).rolling(24, center=True).mean().plot(c='k', linewidth=1.2, label='ensmean')
for i in range(20):
    total['amundsen_shelf_break_uwind_avg']['1920':'2013']['ens'+str(i+1).zfill(2)].rolling(24, center=True).mean().plot(c='grey', linewidth=0.1, label='')
for i in [4]:
    total['amundsen_shelf_break_uwind_avg']['1920':'2013']['ens'+str(i).zfill(2)].rolling(24, center=True).mean().plot(linewidth=1, label=str(i))

wl=['ens'+str(i).zfill(2) for i in wishlist]
total['amundsen_shelf_break_uwind_avg']['1920':'2013'][wl].mean(axis=1).rolling(24, center=True).mean().plot(c='b', linewidth=1.2, label='Mean of selected (w/o 4)')   

plt.grid(True)
#plt.legend()
#plt.ylabel('amundsen_shelf_break_uwind_avg')
plt.subplot(4,2,6)
total['dotson_to_cosgrove_massloss']['1920':'2013'].mean(axis=1).rolling(24, center=True).mean().plot(c='k', linewidth=1.2, label='ensmean')
for i in range(20):
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i+1).zfill(2)].rolling(24, center=True).mean().plot(c='grey', linewidth=0.1, label='')
for i in [4]:
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i).zfill(2)].rolling(24, center=True).mean().plot(linewidth=1, label=str(i))
total['dotson_to_cosgrove_massloss']['1920':'2013'][wl].mean(axis=1).rolling(24, center=True).mean().plot(c='b', linewidth=1.2,label='Mean of selected (w/o 4)')
plt.grid(True)
#plt.ylabel('dotson_to_cosgrove_massloss')
plt.legend()





plt.subplot(4,2,7)
plt.title('Mean')
wishlist=[11, 12, 6, 10, 4]
total['amundsen_shelf_break_uwind_avg']['1920':'2013'].mean(axis=1).rolling(24, center=True).mean().plot(c='k', linewidth=1.2, label='ensmean')
for i in range(20):
    total['amundsen_shelf_break_uwind_avg']['1920':'2013']['ens'+str(i+1).zfill(2)].rolling(24, center=True).mean().plot(c='grey', linewidth=0.1, label='')
# for i in wishlist:
#     total['amundsen_shelf_break_uwind_avg']['1920':'2013']['ens'+str(i).zfill(2)].rolling(24, center=True).mean().plot(linewidth=1, label=str(i))

wl=['ens'+str(i).zfill(2) for i in wishlist]
total['amundsen_shelf_break_uwind_avg']['1920':'2013'][wl].mean(axis=1).rolling(24, center=True).mean().plot(c='b', linewidth=1.2, label='Mean of selected')   

plt.grid(True)
#plt.legend()
#plt.ylabel('amundsen_shelf_break_uwind_avg')
plt.subplot(4,2,8)
total['dotson_to_cosgrove_massloss']['1920':'2013'].mean(axis=1).rolling(24, center=True).mean().plot(c='k', linewidth=1.2, label='ensmean')
for i in range(20):
    total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i+1).zfill(2)].rolling(24, center=True).mean().plot(c='grey', linewidth=0.1, label='')
# for i in wishlist:
#     total['dotson_to_cosgrove_massloss']['1920':'2013']['ens'+str(i).zfill(2)].rolling(24, center=True).mean().plot(linewidth=1, label=str(i))
total['dotson_to_cosgrove_massloss']['1920':'2013'][wl].mean(axis=1).rolling(24, center=True).mean().plot(c='b', linewidth=1.2,label='Mean of selected')
plt.grid(True)
#plt.ylabel('dotson_to_cosgrove_massloss')
plt.legend()


plt.subplots_adjust(hspace=0.3)
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/new_run/selection.png')