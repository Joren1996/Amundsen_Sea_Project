# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 15:55:00 2022

@author: joren
"""
import os
import sys
sys.path.append('/users/kaight/MITgcm/utils/python/MITgcmutils/')
import mitgcm_python_master.projects.pace as pc

bd='/data/oceans_output/shelf/kaight/mitgcm'
pd=[os.path.join(bd, f) for f in os.listdir(bd) if ('PAS_PACE' in f) & ('hb' not in f)][1:2]
print(pd)
#['PAS_PACE01', 'PAS_PACE02', 'PAS_PACE03', 'PAS_PACE04', 'PAS_PACE05', 'PAS_PACE06', 'PAS_PACE07', 'PAS_PACE08', 'PAS_PACE09', 'PAS_PACE10', 'PAS_PACE11', 'PAS_PACE12', 'PAS_PACE13', 'PAS_PACE14', 'PAS_PACE15', 'PAS_PACE16', 'PAS_PACE17', 'PAS_PACE18', 'PAS_PACE19', 'PAS_PACE20']


#pd=['/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE01', '/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE06', '/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE08', '/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE15']
gp='/data/oceans_output/shelf/kaight/mitgcm/PAS_grid/'
tf='timeseries_final.nc'

#pc.correlate_ismr_forcing (pace_dir, grid_path, timeseries_file='timeseries.nc', fig_dir='./')
#pc.correlate_ismr_forcing (pace_dir=pd, grid_path=gp, timeseries_file=tf, fig_dir='./')
#pc.correlate_melt_with_wind_map(pace_dir=pd, grid_path=gp, timeseries_file=tf, fig_dir='./')
#pc.correlate_timeseries_with_map (pace_dir=pd, grid_path=gp, timeseries_file=tf, fig_dir='./', region='amundsen_sea')
#pc.correlate_ensemble_mean_map (pace_dir=pd, grid_path=gp, timeseries_file=tf, fig_dir='./', region='amundsen_bay')
#pc.correlate_timeseries_with_map_full_ensemble (pace_dir=pd, grid_path=gp, timeseries_file=tf, fig_dir='./', region='amundsen_shelf_break', detrend=True, timeseries=['amundsen_shelf_temp_btw_200_700m'], mapvar='EXFuwind', title='Correlation map of local zonal surface winds with Amundsen Shelf Temperatures (200-700m) \n (both 24 running mean) (All members)')



import matplotlib.pyplot as plt
from mitgcm_python_master.grid import Grid

grid=Grid(gp)
plt.figure()
a=plt.contour(grid.lon_2d, grid.lat_2d, grid.bathy, levels=[-4000])
plt.clabel(a, '-1000m')
plt.savefig('test_bathy.png') 
