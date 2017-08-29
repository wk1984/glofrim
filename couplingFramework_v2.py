# -*- coding: utf-8 -*-
"""
Introduction:
-------------
This is the key script to couple PCR-GLOBWB with either Delft3D Flexible Mesh ("DFM") by Deltares or
LISFLOOD-FP ("LFP") by University of Bristol.
The coupling is achieved by making use of the Basic Model Interface (BMI) which allows for initializing,
updating, data manipulation, and finalization of models from a shell-like environment.
For couple Delft3D FM or LISFLOOD-FP, the python-module "bmi.wrapper" has to be loaded.
For Delft3D FM, any compiled version (>1.1.201.48898) has already implemented a BMI-compatible structure, and the
required variables accessible to the user.
For LISFLOOD-FP, however, a specifically designed version needs to be compiled which is currently only available for
version 5.9 as the model is not originally BMI-compatible.
Also for PCR-GLOBWB, a BMI-compatible version needs to be used.

Literature and sources:
-----------------------
	BMI         -> https://csdms.colorado.edu/wiki/BMI_Description
				-> http://www.sciencedirect.com/science/article/pii/S0098300412001252
	bmi.wrapper -> https://github.com/openearth/bmi-python

Running the script:
-------------------
To run the script, an ini-file containing the required specifications and paths is necessary.
Using python, run this file along with the ini-file as follows:
	python couplingFramework_v1.py default.set

Disclaimer:
-----------
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

In addition to this disclaimer, the disclaimers of each component involved in this coupling (i.e. PCR-GLOBWB, LIFLOOD-FP, Delft3D Flexible Mesh)
remain valid.
No warranty/responsibility for any outcome of using this coupling script.
Please ensure to cite the models involved when using this coupling script.

Copyright (C) 2017 Jannis Hoch

@author: Jannis Hoch, Department of Physical Geography, Faculty of Geosciences, Utrecht University (j.m.hoch@uu.nl)
@date: 22-05-2017


Progress status 2way-coupling:

09 Aug 2017:	added functionality to activate floodplain infiltration factor
				made changes to logging and printing to console
				successfully ran current version of 2way coupling for NIG, AMA, and LFP_NIG

30 Jun 2017: 	implemented functionality to obtain two separate sets of coupled PCR-FM cells if RFS is on: one for PCR->FM1d, the other for FM2d->PCR;
				if RFS is not on, it's only one array for the 2d-cells;
				running GLOFRIM works
				
TO DO: 	- make script such that not only the set/ini-file for coupling, but also the mdu/par-file for hydrodynamics and ini-file for hydrology can be put in command line -> easier batch-runs
		- to speed up 2way coupling, maybe it's worth investigating whether filtering out all hydrodynamic cells which do not contain any water to reduce number of cells to be coupled?
"""

# -------------------------------------------------------------------------------------------------
# LOAD REQUIRED LIBRARIES
# -------------------------------------------------------------------------------------------------

import netCDF4
from distutils.util import strtobool
import pylab
import matplotlib
import matplotlib.pyplot as plt
import sys
import os
import platform
import numpy as np
import pyproj as pyproj
import datetime
import pdb
import time
import bmi.wrapper
import pcrglobwb_203_30min_2way as pcrglobwb_bmi_v203
from pcrglobwb_203_30min_2way import pcrglobwb_bmi
from pcrglobwb_203_30min_2way import disclaimer
from coupling_PCR_FM_2way import coupling_functions
from coupling_PCR_FM_2way import model_functions
from coupling_PCR_FM_2way import configuration

# -------------------------------------------------------------------------------------------------
# IMPORT MODEL SETTINGS FROM INI-FILE/SET-FILE
# -------------------------------------------------------------------------------------------------

config = configuration.Configuration()
config.parse_configuration_file(sys.argv[1])

# -------------------------------------------------------------------------------------------------
# SPECIFY MODEL SETTINGS
# -------------------------------------------------------------------------------------------------

model_type = config.model_type['model_type']                                                   

use_2way = strtobool(config.general_settings['2way coupling'])
latlon = strtobool(config.general_settings['latlon'])
if latlon == False:
	inProj  = pyproj.Proj(init=config.model_settings['model_projection'])
use_Fluxes = strtobool(config.general_settings['use_Fluxes'])
use_RFS = strtobool(config.general_settings['use_RFS'])
verbose = strtobool(config.general_settings['verbose'])
use_floodplain_infiltration_factor = strtobool(config.general_settings['use_floodplain_infiltration_factor'])
adjust_initial_groundwater = strtobool(config.general_settings['adjust_initial_groundwater'])

# -------------------------------------------------------------------------------------------------
# SPECIFY NUMERICAL SETTINGS
# -------------------------------------------------------------------------------------------------

nr_pcr_timesteps                      = int(config.numerical_settings['number_of_timesteps'])                      
update_step                           = int(config.numerical_settings['update_step'])  
                      
secPerDay                             = 86400.
end_time 							  = nr_pcr_timesteps * secPerDay
fraction_timestep 					  = secPerDay / update_step

threshold_inundated_depth             = float(config.numerical_settings['threshold_inundated_depth'])                                                   

# other
missing_value_landmask                = 255
missing_value_pcr                     = -999

# -------------------------------------------------------------------------------------------------
# PLOT AND PRINT OPTIONS
# -------------------------------------------------------------------------------------------------

# set default figure size
pylab.rcParams['figure.figsize']      = (14.0, 7.0)
# set plot stuff for coloured plots of FM water depths
my_cmap = matplotlib.cm.get_cmap('Blues_r')
my_cmap.set_under('seagreen')
my_cmap.set_bad('seagreen')

# -------------------------------------------------------------------------------------------------
# SET PATHS TO MODELS
# -------------------------------------------------------------------------------------------------

model_dir       	= config.model_settings['model_dir'] 
model_file      	= config.model_settings['model_file']
model_proj			= config.model_settings['model_projection']                                    

config_pcr       	=  config.PCR_settings['config_pcr']
landmask_pcr     	=  config.PCR_settings['landmask_pcr']
clone_pcr        	=  config.PCR_settings['clone_pcr']

adjust_initial_groundwater_file = config.general_settings['adjust_initial_groundwater_file']

# -------------------------------------------------------------------------------------------------
# SET PATHS TO .SO / .DLL FILES
# -------------------------------------------------------------------------------------------------

# these may be changed according to personal file and folder structure
if model_type == 'DFM':
    model_path = config.model_paths['path2DFM']

elif model_type == 'LFP':
    model_path = config.model_paths['path2LFP']

else:
    sys.exit('\nno adequate model defined in configuration file - define either DFM or LFP!\n')

# -------------------------------------------------------------------------------------------------
# INITIALIZE AND SPIN-UP PCR-GLOBWB
# -------------------------------------------------------------------------------------------------
                                  
# initiate logging and define folder for verbose-output
verbose_folder = model_functions.write2log(model_dir, model_file, latlon, use_2way, use_Fluxes, use_RFS, use_floodplain_infiltration_factor, adjust_initial_groundwater, verbose, moment='start')
print 'Model Start-Time: ', datetime.datetime.now()
print ''

# print PCR-GLOBWB disclaimer
print '\n>>> Please consider reading the PCR-GLOBWB disclaimer <<<'
disclaimer.print_disclaimer()
time.sleep(5)

# initiate PCR-GLOBWB
model_pcr = pcrglobwb_bmi_v203.pcrglobwb_bmi.pcrglobwbBMI()
model_pcr.initialize(config_pcr)
print '\n>>> PCR Initialized <<<\n' 

# spin-up PCR-GLOBWB
model_pcr.spinup()

# -------------------------------------------------------------------------------------------------
# INITIALIZING HYDRODYNAMIC MODEL
# -------------------------------------------------------------------------------------------------

# initiate hydraulic model
model_hydr = bmi.wrapper.BMIWrapper(engine = model_path, configfile = (os.path.join(model_dir, model_file)))
model_hydr.initialize()
print '\n>>> ',model_type,' Initialized <<<\n' 

# -------------------------------------------------------------------------------------------------
# EXCTRACTING RELEVANT DATA FROM MODELS
# -------------------------------------------------------------------------------------------------

if model_type == 'DFM':

    #- retrieving data from Delft3D FM    
    x_coords, y_coords, z_coords, bottom_lvl, cell_points_fm, separator_1D, cellAreaSpherical_1way, cellAreaSpherical_2way, xz_coords, yz_coords, modelCoords, modelCoords_2way,\
                cellarea_data_pcr, landmask_data_pcr, clone_data_pcr = model_functions.extractModelData_FM(model_hydr, model_pcr, landmask_pcr, clone_pcr, use_RFS, use_2way)
    coupledFPindices = 0.
    print '\n>>> DFM data retrieved <<<\n'
         
elif model_type == 'LFP':
    
    #- retrieving data from LISFLOOD-FP
    dx, dy, DEM, bottom_lvl, H, waterDepth, rows, cols, \
                list_x_coords, list_x_coords_2way, list_y_coords, list_y_coords_2way, coupledFPindices, coupledFPindices_2way, grid_dA, cellAreaSpherical_1way, SGCQin, separator_1D,\
                cellarea_data_pcr, landmask_data_pcr, clone_data_pcr = model_functions.extractModelData_FP(model_hydr, model_dir, model_pcr, landmask_pcr, clone_pcr, verbose_folder, use_RFS, use_2way, verbose)
                
    cellAreaSpherical_2way = np.copy(cellAreaSpherical_1way)
                
    #- computing FP-coordinates    
    modelCoords = coupling_functions.getVerticesFromMidPoints(list_x_coords, list_y_coords, dx, dy, verbose)
    modelCoords_2way = coupling_functions.getVerticesFromMidPoints(list_x_coords_2way, list_y_coords_2way, dx, dy, verbose)
    print '\n>>> LFP data retrieved <<<\n'

#- computing PCR-coordinates
PCRcoords = coupling_functions.getPCRcoords(landmask_data_pcr)
print '\n>>> PCR data retrieved <<<\n'
		
# -------------------------------------------------------------------------------------------------
# COUPLING THE GRIDS
# -------------------------------------------------------------------------------------------------

# this is only required for plotting later, not for actual coupling process
CoupledCellsInfoAll = coupling_functions.coupleAllCells(modelCoords,PCRcoords)
CoupledCellsInfoAll_2way = coupling_functions.coupleAllCells(modelCoords_2way,PCRcoords)

# converting single indices of coupled PCR cells to double (array,column) indices
if use_2way == False:
    CoupleModel2PCR, CouplePCR2model, CoupledPCRcellIndices = coupling_functions.assignPCR2cells(landmask_pcr, modelCoords, verbose)
elif use_2way == True:
    CoupleModel2PCR, CouplePCR2model, CoupledPCRcellIndices = coupling_functions.assignPCR2cells(landmask_pcr, modelCoords, verbose)
    CoupleModel2PCR_2way, CouplePCR2model_2way, CoupledPCRcellIndices_2way = coupling_functions.assignPCR2cells(landmask_pcr, modelCoords_2way, verbose)

# saving plots of coupled cells to verbose-folder
# currently doesn't work with FM and use_RFS on, due to data structure required (? check this ?)
if (model_type == 'DFM') and (use_RFS == True):
    pass
else:
	if (verbose == True) and (use_2way == False):
		coupling_functions.plotGridfromCoords(PCRcoords, modelCoords)
		plt.savefig(os.path.join(verbose_folder , 'AllCells_1way.png'))
		coupling_functions.plotGridfromCoords(CoupledCellsInfoAll[1],CoupledCellsInfoAll[0])
		plt.savefig(os.path.join(verbose_folder , 'CoupledCells_1way.png'))   
		plt.close('all')
	elif (verbose == True) and (use_2way == True):
		coupling_functions.plotGridfromCoords(PCRcoords, modelCoords)
		plt.savefig(os.path.join(verbose_folder , 'AllCells_1way.png'))
		coupling_functions.plotGridfromCoords(CoupledCellsInfoAll[1],CoupledCellsInfoAll[0])
		plt.savefig(os.path.join(verbose_folder , 'CoupledCells_1way.png'))   
		plt.close('all')    
		coupling_functions.plotGridfromCoords(PCRcoords, modelCoords_2way)
		plt.savefig(os.path.join(verbose_folder , 'AllCells_2way.png'))
		coupling_functions.plotGridfromCoords(CoupledCellsInfoAll_2way[1],CoupledCellsInfoAll_2way[0])
		plt.savefig(os.path.join(verbose_folder , 'CoupledCells_2way.png'))   
		plt.close('all')
    
# -------------------------------------------------------------------------------------------------
# TURNING OFF CHANNELSTORAGE, WATERBODYSTORAGE, WATERBODIES AND RUNOFF TO CHANNELS
# ------------------------------------------------------------------------------------------------- 

model_functions.noStorage(model_pcr, missing_value_pcr, CoupledPCRcellIndices, CouplePCR2model)

# -------------------------------------------------------------------------------------------------
# TURNING OFF ROUTING BY PCR IN COUPLED AREA
# -------------------------------------------------------------------------------------------------

model_functions.noLDD(model_pcr, CoupledPCRcellIndices, verbose_folder, verbose)

# -------------------------------------------------------------------------------------------------
# ACTIVATING A RANGE OF VARIABLES SPECIFICALLY REQUIRED FOR 2WAY-COUPLING
# -------------------------------------------------------------------------------------------------

#TODO: applying the function here leads to tremendous increase of water volume coupled from PCR to DFM
#model_functions.activate2wayVariables(model_pcr, CoupledPCRcellIndices)

model_functions.adjust_iniGR(model_pcr, adjust_initial_groundwater_file, CoupledPCRcellIndices_2way, adjust_initial_groundwater=False) # optional, if adjust_initial_groundwater is set to True

model_functions.activate_floodplain_infiltration_factor(model_pcr, CoupledPCRcellIndices_2way, use_floodplain_infiltration_factor=False) # optional, if use_floodplain_infiltration_factor is set to True

# -------------------------------------------------------------------------------------------------
# CALCULATE DELTA VOLUMES (DAY 1)
# first day outside loop, to make sure PCR time is at start time (timestep 1) and not at end of spin-up (timestep 365)
# ------------------------------------------------------------------------------------------------- 

# using a dummy to avoid feeding back volumes to PCR before 1st update actually took place
water_volume_FM_2_PCR = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)

# setting PCR's topWaterLayer to 0 for all cells coupled to hydrodynamics (1way)
model_functions.set_zeroTopWaterlayer(model_pcr, CoupledPCRcellIndices)

# retrieving PCR-GLOBWB and converting it to m3/d
delta_volume_PCR, delta_volume_PCR_coupled = model_functions.calculateDeltaVolumes(model_pcr, missing_value_pcr, secPerDay, CoupledPCRcellIndices, cellarea_data_pcr, water_volume_FM_2_PCR)

# removing negative delta volumes where necessary
delta_volume_PCR_positiveOnly, delta_volume_PCR_coupled_positiveOnly = model_functions.account4negativeDeltaVolumes(model_hydr, model_type, CoupledPCRcellIndices, CoupledPCRcellIndices_2way, CouplePCR2model_2way, delta_volume_PCR, cellAreaSpherical_2way)

# dividing delta volume from PCR-GLOBWB over hydraulic cells, depending on model specifications
delta_water_fm, verbose_volume = model_functions.calculateDeltaWater(model_hydr, CoupleModel2PCR, CouplePCR2model, delta_volume_PCR_coupled_positiveOnly, cellAreaSpherical_1way, fraction_timestep, model_type, use_Fluxes)

#pdb.set_trace()

# saving PCR-GLOBWB output volumes and volumes used as input to hydraulic models to verbose-folder
if verbose == True:
	# initial file objects
    fo_PCR_V_tot = open(os.path.join(verbose_folder, 'delta_volume_PCR_coupled.txt'), 'w')
    fo_verbose_volume = open(os.path.join(verbose_folder, 'verbose_volume.txt'), 'w')
    # aggregate daily totals
    delta_volume_PCR_total_aggr = np.sum(delta_volume_PCR_coupled)
    verbose_volume_aggr = np.sum(verbose_volume)
    # write daily totals to file
    fo_PCR_V_tot.write(str(delta_volume_PCR_total_aggr) + os.linesep)
    fo_verbose_volume.write(str(verbose_volume_aggr) + os.linesep)

# check to ensure that volumes are equal, i.e. no errors in water balance
if np.round(np.sum(np.asfarray(delta_volume_PCR_coupled)),4) != np.round(np.sum(np.asfarray(verbose_volume)), 4):
    print 'PCR volume out: ', np.sum(np.asfarray(delta_volume_PCR_coupled))
    print 'FM volume in: ', np.sum(np.asfarray(verbose_volume))
    sys.exit('\nERROR: input volumes to hydrodynamics do not match PCR output volumes!\n')

# reshaping data for LISFLOOD-FP from list to arrays
if model_type == 'LFP':
    delta_water_fm = model_functions.fillLFPgrid(model_hydr, coupledFPindices, delta_water_fm, DEM, verbose_folder, verbose)
  
# -------------------------------------------------------------------------------------------------
# FIRST UPDATE (DAY 1)
# -------------------------------------------------------------------------------------------------

print '\n>>> update 1 started <<<\n'

# updating arrays with computed additional volumes; array used depends on model specifications
if (model_type == 'LFP'):
    model_functions.updateModel(model_hydr, delta_water_fm, update_step, separator_1D, use_Fluxes, use_RFS, model_type, verbose)

while model_hydr.get_current_time() < (model_pcr.get_time_step() * secPerDay):
    
    # updating FM on user-specified time step
    if (model_type == 'DFM'):
        model_functions.updateModel(model_hydr, delta_water_fm, update_step, separator_1D, use_Fluxes, use_RFS, model_type, verbose)
    
    # updating FM or FP on daily time step
    if (model_type == 'LFP'):
        model_hydr.update()

# ----------------------------------------------------------------------------------------------------
# UPDATE FM FOR THE REST OF THE MODEL PERIOD
# ----------------------------------------------------------------------------------------------------

while model_pcr.get_time_step() < nr_pcr_timesteps:
	
    # COMPUTING INUNDATION AREA AND AS FRACTION OF PCR-CELL
    inundated_area_FM_2_PCR_coupled, inundated_fraction_FM_2_PCR =  model_functions.determine_InundationArea_Hydrodynamics(model_type, model_hydr, CouplePCR2model_2way, CoupledPCRcellIndices_2way, threshold_inundated_depth, cellAreaSpherical_2way, cellarea_data_pcr, landmask_pcr, missing_value_landmask)

    # COMPUTING WATER DEPTH AND VOLUME TO BE COUPLED BACK TO PCR
    water_depths_FM_2_PCR, water_volume_FM_2_PCR = model_functions.determine_InundationDepth_Hydrodynamics(model_type, model_hydr, landmask_pcr, missing_value_landmask, inundated_area_FM_2_PCR_coupled, CouplePCR2model_2way, CoupledPCRcellIndices_2way, cellAreaSpherical_1way)

    # DETERMINING NEW STORAGE IN PCR-CHANNELS
    new_storage_pcr = model_functions.determine_new_channelStoragePCR(model_pcr, landmask_pcr, missing_value_landmask, water_volume_FM_2_PCR)
    
    # UPDATE VARIABLES IN PCR
    model_functions.updateHydrologicVariables(model_pcr, water_depths_FM_2_PCR, inundated_fraction_FM_2_PCR, new_storage_pcr)
    
    # retrieving PCR-GLOBWB and converting it to m3/d
    delta_volume_PCR, delta_volume_PCR_coupled = model_functions.calculateDeltaVolumes(model_pcr, missing_value_pcr, secPerDay, CoupledPCRcellIndices, cellarea_data_pcr, water_volume_FM_2_PCR)                                                                                                  
    
    # removing negative delta volumes where necessary
    delta_volume_PCR_positiveOnly, delta_volume_PCR_coupled_positiveOnly = model_functions.account4negativeDeltaVolumes(model_hydr, model_type, CoupledPCRcellIndices, CoupledPCRcellIndices_2way, CouplePCR2model_2way, delta_volume_PCR, cellAreaSpherical_2way)
				                 
    # dividing delta volume from PCR-GLOBWB over hydraulic cells, depending on model specifications
    delta_water_fm, verbose_volume = model_functions.calculateDeltaWater(model_hydr, CoupleModel2PCR, CouplePCR2model, delta_volume_PCR_coupled_positiveOnly, cellAreaSpherical_1way, fraction_timestep, model_type, use_Fluxes)
    
    # saving PCR-GLOBWB output volumes and volumes used as input to hydraulic models to verbose-folder
    if verbose == True:
        # aggregate daily totals
        delta_volume_PCR_total_aggr = np.sum(delta_volume_PCR_coupled)
        verbose_volume_aggr = np.sum(verbose_volume)
        # write daily totals to file
        fo_PCR_V_tot.write(str(delta_volume_PCR_total_aggr) + os.linesep)
        fo_verbose_volume.write(str(verbose_volume_aggr) + os.linesep)
    
    # reshaping data for LISFLOOD-FP from list to arrays
    if model_type == 'LFP':
        delta_water_fm = model_functions.fillLFPgrid(model_hydr, coupledFPindices, delta_water_fm, DEM, verbose_folder, verbose)  
    
    # updating arrays with computed additional volumes; array used depends on model specifications
    if (model_type == 'LFP'):
        model_functions.updateModel(model_hydr, delta_water_fm, update_step, separator_1D, use_Fluxes, use_RFS, model_type, verbose)      

    # update FM unless it has has reached the same time as PCR
    while model_hydr.get_current_time() < (model_pcr.get_time_step() * secPerDay):
        
        # updating FM on user-specified time step
        if (model_type == 'DFM'):
            model_functions.updateModel(model_hydr, delta_water_fm, update_step, separator_1D, use_Fluxes, use_RFS, model_type, verbose)      
        
        # updating FM or FP on daily time step
        if (model_type == 'LFP'):
            model_hydr.update()   
        
# ----------------------------------------------------------------------------------------------------
# END OF MODEL PERIOD REACHED
# ----------------------------------------------------------------------------------------------------
    
# update and finalize logging
model_functions.write2log(model_dir, model_file, latlon, use_2way, use_Fluxes, use_RFS, use_floodplain_infiltration_factor, adjust_initial_groundwater, verbose, moment='end') 
# close files
if verbose == True:
    fo_PCR_V_tot.close()
    fo_verbose_volume.close()

#- finalizing hydrodynamic model to properly end execution
model_hydr.finalize()
print '\nModel End-Time: ', datetime.datetime.now()
