# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 13:30:03 2022

@author: joren
"""


from mitgcm_python_master.postprocess import precompute_timeseries_coupled

for i in range(1,21):
    print('Start with member: '+str(i))
    precompute_timeseries_coupled (output_dir='/data/oceans_output/shelf/kaight/mitgcm/PAS_PACE{}/output'.format(str(i).zfill(2)), timeseries_file='./data/timeseries_ohc_PACE{}.nc'.format(str(i).zfill(2)), 
                                   hovmoller_file='hovmoller.nc', file_name='output.nc', segment_dir=None, timeseries_types=['amundsen_shelf_ohc_below_0m', 'pine_island_bay_ohc_below_0m'], 
                                   hovmoller_loc=[], key='PAS', time_average=False)

    print('Finished with member: '+str(i))
print('Finished with *EVERYTHING*!')


# import xarray as xr
# import os
# import shutil


# real=True

# goal='hovmoller_adv.nc'
# #goal='timeseries_final.nc'

# if real==True:
#     mainpath=r'/data/oceans_output/shelf/kaight/mitgcm'
# else:
#     mainpath=r'C:\Users\joren\OneDrive\Documenten\Master Thesis\Project\MITgcm\data_MITgcm\timeseries'

# for folder in os.listdir(mainpath):
#     if 'PAS_PACE' in folder:
#         if 'hb' in folder:
#             print('Skip!')
#         else:
#             if real==True:
#                 file=os.path.join(mainpath, folder, 'output', goal)
#             else:
#                 file=os.path.join(mainpath, folder, goal)
#             print('Start with {}!'.format(file))
#             shutil.copyfile(file, goal[:-3]+folder[-7:]+'.nc')
#             print('Finished with {}!'.format(file))
