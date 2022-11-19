# READ ME for  'Amundsen Sea Project'
**Author**: Joren Janzing
**Contact**: gjanzing@mailfence.com

This folder contains the code which I used for my Master Thesis titled: **'Oceanic and Atmospheric Controls on Ice Shelf Basal Melt Variability around the Amundsen Sea'**.

This code is used to create figures based on the output of runs with MITgcm.
Note that some functions rely on access to https://github.com/knaughten/mitgcm_python

The structure is as follows: the main folder is `01_scripts`. This contains all the scripts used for the figures in this project.

The main folder contains several normal python files, which are used to read and process the data in order to create the figures I used in the report:
- `01_main_preprocessing.py`: run the code here to read data from the MITgcm output and save in the 02_data folder
- `02_main_create_indices.py`: run the code here to compute indices used in the report.
- `05_main_correlation_maps.py`: run code here to create the correlation maps
- `09_new_model_run.py`: run code here to show which members were selected for the runs with a new geometry.
- `10_main_convection.py`: the code here is used to check if correlation maps change when convection is not taken into account.
- `11_figures_for_report.py`: most figures in the report can be created here (after running the other files).

Furthermore, it contains several Jupyter Notebook files:
- `xx_baroclinic_feedback.ipynb`: shows that the undercurrent has co-varies in a baroclinic way with the melt.
- `xx_case1940s.ipynb`: shows a case study for the 1940s event.
- `xx_composite_high_and_low_melt.ipynb`: shows the difference in bottom, baroclinic and barotropic flow during high and low melt.
- `xx_connection_wind_and_flow.ipynb`: shows the spread between individual ensemble members in the connection between wind and melt.
- `xx_farFieldForcing.ipynb`: shows the influence of remote forcing on Amundsen Sea
-  `xx_newGeometry.ipynb`: shows the sensitivity of the system to ice shelf geometry
-  `xx_terrain_map.ipynb`: illustrates the domains used in the model

The subfolder `functions` contains python files with functions that are used in the main code.
- `preprocessing.py`: functions related to reading and saving MITgcm output
- `loading_and_processing_data.py`: functions to load the data when it has already been saved during preprocessing.
- `correlation_maps.py`: functions relating to correlation maps
- `convection.py`: functions relating to correlations maps without convection

The folders `02_data` and `03_output` now only contain empty folders, but are used to save data or figures when running the code.