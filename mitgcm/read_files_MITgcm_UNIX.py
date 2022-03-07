# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 13:30:03 2022

@author: joren
"""

import xarray as xr
import os
import shutil

real=True

goal='hovmoller_adv.nc'
#goal='timeseries_final.nc'

if real==True:
    mainpath=r'/data/oceans_output/shelf/kaight/mitgcm'
else:
    mainpath=r'C:\Users\joren\OneDrive\Documenten\Master Thesis\Project\MITgcm\data_MITgcm\timeseries'

for folder in os.listdir(mainpath):
    if 'PAS_PACE' in folder:
        if 'hb' in folder:
            print('Skip!')
        else:
            if real==True:
                file=os.path.join(mainpath, folder, 'output', goal)
            else:
                file=os.path.join(mainpath, folder, goal)
            print('Start with {}!'.format(file))
            shutil.copyfile(file, goal[:-3]+folder[-7:]+'.nc')
            print('Finished with {}!'.format(file))
