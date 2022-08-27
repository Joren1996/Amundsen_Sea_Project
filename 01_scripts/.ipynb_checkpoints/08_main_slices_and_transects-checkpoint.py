# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain test!.

@author: joren
"""
#%%IMPORTING
import xarray as xr
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import pandas as pd
import datetime
import matplotlib.colors as colors

divnorm = colors.TwoSlopeNorm(vcenter=0)

plt.rcParams.update({'font.size': 18})

from functions.loading_and_processing_data import loadMITgcmData, detrend_and_average_MITgcmData, read_all_data
from functions.composite_maps import createCompositeMap, createVelocityCompositeMap, defineEvents, findEvents

#%--------COMPARE DIFFERENT TIMES-------------------------------------------------------------------------------

members=[18,19]
#01 MAKE CROSS SECTIONS A LA ALESSANDRO
lon=123

dataS=loadMITgcmData(filename='SALT_W'+str(lon), members=members, kind='slices')
dataT=loadMITgcmData(filename='THETA_W'+str(lon), members=members, kind='slices')

#Convert everything to density.
import sys
sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
from mitgcm_python_master.diagnostics import density
ref_depth=0
eosType='MDJWF'
dataS=dataS.where(dataS!=0, drop=True)
dataT=dataT.where(dataT!=0, drop=True)

salt=dataS['SALT'].to_numpy()
temp=dataT['THETA'].to_numpy()
q=density (eosType, salt, temp, press=ref_depth)
q=xr.DataArray(q, dims=('ens', 'time', "Z", "YC"), coords={'ens': dataS['SALT'].ens, 'time': dataS['SALT'].time, "YC": dataS['SALT'].YC, "Z": dataS['SALT'].Z})


fig=plt.figure(figsize=(20,10))
fig.suptitle('Longitudinal cross sections of density in member '+str(members[-1])+' across shelf break at '+str(lon)+' degrees W')
ax=plt.subplot(2,2,1)
ax.set_title('Mean over 1920-1950')
img=ax.pcolor(q.YC, q.Z, q.mean(dim='ens').sel(time=slice('1920', '1950')).mean(dim='time')-1027, cmap='seismic', vmin=0.4)
cn=ax.contour(q.YC, q.Z, q.mean(dim='ens').sel(time=slice('1920', '1950')).mean(dim='time'), [1027.5, 1027.6, 1027.7, 1027.8], colors='k')
#ax.clabel(cn, [1027.5, 1027.6, 1027.7, 1027.8])
plt.colorbar(img, label='Density - 1027 [kg/m3]', extend='min')
if lon==115:
    ax.set_xlim([-71.8, -71.1])
elif lon==123:
    ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')

ax=plt.subplot(2,2,2)
ax.set_title('Mean over 1980-2010')
img=ax.pcolor(q.YC, q.Z, q.mean(dim='ens').sel(time=slice('1980', '2010')).mean(dim='time')-1027, cmap='seismic', vmin=0.4)
cn=ax.contour(q.YC, q.Z, q.mean(dim='ens').sel(time=slice('1980', '2010')).mean(dim='time'), [1027.5, 1027.6, 1027.7, 1027.8], colors='k')
#ax.clabel(cn, [1027.5, 1027.6, 1027.7, 1027.8])
plt.colorbar(img, label='Density - 1027 [kg/m3]', extend='min')
if lon==115:
    ax.set_xlim([-71.8, -71.1])
elif lon==123:
    ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')


ax=plt.subplot(2,2,3)
ax.set_title('Difference Mean1980to2010 - Mean1920to1950')
img=ax.pcolor(q.YC, q.Z, q.mean(dim='ens').sel(time=slice('1980', '2010')).mean(dim='time')-q.mean(dim='ens').sel(time=slice('1920', '1950')).mean(dim='time'), cmap='seismic', vmin=-0.08, vmax=0.08)
#cn=ax.contour(q.YC, q.Z, q.mean(dim='ens').sel(time=slice('1980', '2010')).mean(dim='time'), [1027.5, 1027.6, 1027.7, 1027.8], colors='k')
#ax.clabel(cn, [1027.5, 1027.6, 1027.7, 1027.8])
plt.colorbar(img)
if lon==115:
    ax.set_xlim([-71.8, -71.1])
elif lon==123:
    ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')

ax=plt.subplot(2,2,4)
ax.set_title('Trend 1920-2013')
img=ax.pcolor(q.YC, q.Z, q.sel(time=slice('1920', '2010')).polyfit(dim='time', deg=1).sel(degree=1, ens=0)['polyfit_coefficients'], cmap='seismic', vmin=-1e-20, vmax=1e-20)
#cn=ax.contour(q.YC, q.Z, q.mean(dim='ens').sel(time=slice('1980', '2010')).mean(dim='time'), [1027.5, 1027.6, 1027.7, 1027.8], colors='k')
#ax.clabel(cn, [1027.5, 1027.6, 1027.7, 1027.8])
plt.colorbar(img)
if lon==115:
    ax.set_xlim([-71.8, -71.1])
elif lon==123:
    ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')


fig.subplots_adjust(hspace=0.4, wspace=0.4)

fig.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/RHO_slicesintime_member'+str(members[-1])+'_at_lon'+str(lon)+'.png')




#%--------COMPARE DIFFERENT MEMBERS-------------------------------------------------------------------------------

#01 MAKE CROSS SECTIONS A LA ALESSANDRO
lon=115

dataS=loadMITgcmData(filename='SALT_W'+str(lon), members='all', kind='slices')
dataT=loadMITgcmData(filename='THETA_W'+str(lon), members='all', kind='slices')

#Convert everything to density.
import sys
sys.path.append('/data/hpcdata/users/grejan/mitgcm/') #Make sure we can also import Kaitlins code.
from mitgcm_python_master.diagnostics import density
ref_depth=0
eosType='MDJWF'
dataS=dataS.where(dataS!=0, drop=True)
dataT=dataT.where(dataT!=0, drop=True)

salt=dataS['SALT'].to_numpy()
temp=dataT['THETA'].to_numpy()
q=density (eosType, salt, temp, press=ref_depth)
q=xr.DataArray(q, dims=('ens', 'time', "Z", "YC"), coords={'ens': dataS['SALT'].ens, 'time': dataS['SALT'].time, "YC": dataS['SALT'].YC, "Z": dataS['SALT'].Z})


dataU=loadMITgcmData(filename='UVEL_W'+str(lon), members='all', kind='slices')
dataU=dataU.where(dataU!=0, drop=True)



#------------------------------------------------------------------------VISUALISATION
#ENSEMBLE MEAN

#01 Ensemble Mean Density
plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of ensemble mean density \n across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
img=ax.pcolor(q.YC, q.Z, q.mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time'), cmap='jet')
plt.colorbar(img)
#ax.set_xlim([-71.8, -71.1])
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/RHO_W'+str(lon)+'_ensemble_mean_1979-2013.png')
plt.show()

#02 Ensemble Mean Zonal Velocity
vmax=np.max([abs(dataU['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time').min()),abs(dataU['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time').max())])
plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of ensemble mean zonal velocity across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
img=ax.pcolor(dataU.YC, dataU.Z, dataU['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar(img)
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/UVEL_W'+str(lon)+'_ensemble_mean_1979-2013.png')
plt.show()

#Difference with IPO
window=1

#ONI=pd.read_pickle('./data/extendedONI_1905-2005.pkl')
#IPO=pd.read_csv('./data/IPO_NOAA_filtered.txt', sep=";")
cindfile='IPO_NOAA_unfiltered.txt' 
path='/data/hpcdata/users/grejan/mitgcm/data/'
cind=pd.read_csv(os.path.join(path, cindfile), names=range(1,13), skiprows=1, skipfooter=11, delim_whitespace=True)
cind=cind.replace(-99.000, np.nan)
cind=cind.replace(-99.99, np.nan)

cind=cind.stack()
cind.index=[pd.to_datetime(str(j)+'-'+str(i)) for i, j in cind.index]

ci=xr.DataArray(data=cind.values, dims=['time'], coords={'time':cind.index})
ci=ci.rolling(time=window).mean(dim='time')

try:
    dataU['time']=dataU.indexes['time'].to_datetimeindex()
except:
    print('Already datetime')

vmax=np.nanmax([abs(dataU['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').where(ci<0).mean(dim='time')-dataU['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').mean(dim='time')).min(),
             abs(dataU['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').where(ci<0).mean(dim='time')-dataU['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').mean(dim='time')).max()])
#difference IPO
plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of difference in zonal velocity during negative IPO vs mean state for ensemble mean  \n  across shelf break at '+str(lon)+' degrees W (mean over 1979-2013; running mean window: '+str(window)+ ' months)')
img=ax.pcolor(dataU.YC, dataU.Z, dataU['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').where(ci<0).mean(dim='time')-dataU['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar(img, label='[m/s]')
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/UVEL_W'+str(lon)+'_negIPOdiff_ensemblemean_1979-2013_rm'+str(window)+'.png')
plt.show()
    
# plt.figure(figsize=(20,10))
# ax=plt.subplot(1,1,1)
# ax.set_title('Longitudinal cross section of ensemble mean potential temperature across shelf break at 115 degrees W (mean over 1979-2013)')
# img=ax.pcolor(dataT.YC, dataT.Z, dataT['THETA'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time'), cmap='seismic', vmin=-2, vmax=2)
# plt.colorbar(img)
# ax.set_xlim([-73, -72])
# ax.set_ylim([-2500, 0])
# ax.set_ylabel('Depth [m]')
# ax.set_xlabel('Latitude [degree N]')
# ax.set_facecolor('grey')
# plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/THETA_W'+str(lon)+'_ensemble_mean_1979-2013.png')
# plt.show()


#----------------------------------------------------------------------------------------------------ENSEMBLE MEMBER
ens=19


vmax=np.nanmax([abs(dataU['UVEL'].sel(ens=ens-1).sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time').min()),abs(dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time').max())])

plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of member '+str(ens)+' zonal velocity across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
img=ax.pcolor(dataU.YC, dataU.Z, dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar(img)
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/UVEL_W'+str(lon)+'_ens'+str(ens)+'_1979-2013.png')
plt.show()


vmax=np.nanmax([abs(dataU['UVEL'].sel(YC=slice(-71.8, -71.1)).sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time')-dataU['UVEL'].mean(dim='ens').sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time')).min(),
             abs(dataU['UVEL'].sel(YC=slice(-71.8, -71.1)).sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time')-dataU['UVEL'].mean(dim='ens').sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time')).max()])

print(vmax)

#difference
plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of difference between mean state of member '+str(ens)+' vs ensemble mean \n of zonal velocity across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
img=ax.pcolor(dataU.YC, dataU.Z, dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time')-dataU['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar(img)
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/UVEL_W'+str(lon)+'_diff_ens'+str(ens)+'_1979-2013.png')
plt.show()

plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of member '+str(ens)+' density across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
img=ax.pcolor(q.YC, q.Z, q.sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time'))
plt.colorbar(img)
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/RHO_W'+str(lon)+'_ens'+str(ens)+'_1979-2013.png')
plt.show()

#difference

vmax=np.nanmax([abs(q.sel(ens=ens-1).sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time')-q.mean(dim='ens').sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time')).min(),
             abs(q.sel(ens=ens-1).sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time')-q.mean(dim='ens').sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time')).max()])

plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of difference of density \n between mean state of member '+str(ens)+' vs ensemble mean across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
img=ax.pcolor(q.YC, q.Z, q.sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time')-q.mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar(img)
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/RHO_W'+str(lon)+'_diff_ens'+str(ens)+'_1979-2013.png')
plt.show()


# #Difference with IPO
# window=60

# #ONI=pd.read_pickle('./data/extendedONI_1905-2005.pkl')
# #IPO=pd.read_csv('./data/IPO_NOAA_filtered.txt', sep=";")
# cindfile='IPO_NOAA_unfiltered.txt' 
# path='/data/hpcdata/users/grejan/mitgcm/data/'
# cind=pd.read_csv(os.path.join(path, cindfile), names=range(1,13), skiprows=1, skipfooter=11, delim_whitespace=True)
# cind=cind.replace(-99.000, np.nan)
# cind=cind.replace(-99.99, np.nan)

# cind=cind.stack()
# cind.index=[pd.to_datetime(str(j)+'-'+str(i)) for i, j in cind.index]

# ci=xr.DataArray(data=cind.values, dims=['time'], coords={'time':cind.index})
# ci=ci.rolling(time=window).mean(dim='time')

try:
    dataU['time']=dataU.indexes['time'].to_datetimeindex()
except:
    print('Already datetime')

vmax=np.nanmax([abs(dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').where(ci<0).mean(dim='time')-dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').mean(dim='time')).min(),
             abs(dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').where(ci<0).mean(dim='time')-dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').mean(dim='time')).max()])


plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of difference in zonal velocity during negative IPO \n vs mean state for member '+str(ens)+' across shelf break at '+str(lon)+' degrees W (mean over 1979-2013; running mean window: '+str(window)+ ' months)')
img=ax.pcolor(dataU.YC, dataU.Z, dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').where(ci<0).mean(dim='time')-dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar(img)
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/UVEL_W'+str(lon)+'_negIPOdiff_ens'+str(ens)+'_1979-2013_rm'+str(window)+'.png')
plt.show()


#----------------------------------------------------------------------------------------------------ENSEMBLE MEMBER
ens=1


vmax=np.nanmax([abs(dataU['UVEL'].sel(ens=ens-1).sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time').min()),abs(dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time').max())])

plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of member '+str(ens)+' zonal velocity across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
img=ax.pcolor(dataU.YC, dataU.Z, dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar(img)
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/UVEL_W'+str(lon)+'_ens'+str(ens)+'_1979-2013.png')
plt.show()


vmax=np.nanmax([abs(dataU['UVEL'].sel(YC=slice(-71.8, -71.1)).sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time')-dataU['UVEL'].mean(dim='ens').sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time')).min(),
             abs(dataU['UVEL'].sel(YC=slice(-71.8, -71.1)).sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time')-dataU['UVEL'].mean(dim='ens').sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time')).max()])

print(vmax)

#difference
plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of difference between mean state of member '+str(ens)+'\n  vs ensemble mean of zonal velocity across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
img=ax.pcolor(dataU.YC, dataU.Z, dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time')-dataU['UVEL'].mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar(img)
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/UVEL_W'+str(lon)+'_diff_ens'+str(ens)+'_1979-2013.png')
plt.show()

plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of member '+str(ens)+' density across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
img=ax.pcolor(q.YC, q.Z, q.sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time'))
plt.colorbar(img)
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/RHO_W'+str(lon)+'_ens'+str(ens)+'_1979-2013.png')
plt.show()

#difference

vmax=np.nanmax([abs(q.sel(ens=ens-1).sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time')-q.mean(dim='ens').sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time')).min(),
             abs(q.sel(ens=ens-1).sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time')-q.mean(dim='ens').sel(YC=slice(-71.8, -71.1)).sel(time=slice('1979', '2013')).mean(dim='time')).max()])

plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of difference of density \n between mean state of member '+str(ens)+' vs ensemble mean across shelf break at '+str(lon)+' degrees W (mean over 1979-2013)')
img=ax.pcolor(q.YC, q.Z, q.sel(ens=ens-1).sel(time=slice('1979', '2013')).mean(dim='time')-q.mean(dim='ens').sel(time=slice('1979', '2013')).mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar(img)
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/RHO_W'+str(lon)+'_diff_ens'+str(ens)+'_1979-2013.png')
plt.show()


# #Difference with IPO
# window=60

# #ONI=pd.read_pickle('./data/extendedONI_1905-2005.pkl')
# #IPO=pd.read_csv('./data/IPO_NOAA_filtered.txt', sep=";")
# cindfile='IPO_NOAA_unfiltered.txt' 
# path='/data/hpcdata/users/grejan/mitgcm/data/'
# cind=pd.read_csv(os.path.join(path, cindfile), names=range(1,13), skiprows=1, skipfooter=11, delim_whitespace=True)
# cind=cind.replace(-99.000, np.nan)
# cind=cind.replace(-99.99, np.nan)

# cind=cind.stack()
# cind.index=[pd.to_datetime(str(j)+'-'+str(i)) for i, j in cind.index]

# ci=xr.DataArray(data=cind.values, dims=['time'], coords={'time':cind.index})
# ci=ci.rolling(time=window).mean(dim='time')

try:
    dataU['time']=dataU.indexes['time'].to_datetimeindex()
except:
    print('Already datetime')

vmax=np.nanmax([abs(dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').where(ci<0).mean(dim='time')-dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').mean(dim='time')).min(),
             abs(dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').where(ci<0).mean(dim='time')-dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').mean(dim='time')).max()])


plt.figure(figsize=(20,10))
ax=plt.subplot(1,1,1)
ax.set_title('Longitudinal cross section of difference in zonal velocity during negative IPO \n vs mean state for member '+str(ens)+' across shelf break at 115 degrees W (mean over 1979-2013; running mean window: '+str(window)+ ' months)')
img=ax.pcolor(dataU.YC, dataU.Z, dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').where(ci<0).mean(dim='time')-dataU['UVEL'].sel(ens=ens-1).sel(time=slice('1979', '2013')).rolling(time=window).mean(dim='time').mean(dim='time'), cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar(img)
ax.set_xlim([-73, -72])
ax.set_ylim([-2500, 0])
ax.set_ylabel('Depth [m]')
ax.set_xlabel('Latitude [degree N]')
ax.set_facecolor('grey')
plt.savefig('/data/hpcdata/users/grejan/mitgcm/03_output/slices/UVEL_W'+str(lon)+'_negIPOdiff_ens'+str(ens)+'_1979-2013_rm'+str(window)+'.png')
plt.show()


