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
				
TO DO: 	- make script such that not only the set/ini-file for coupling, but also the mdu/par-file for hydrodynamicModel and ini-file for hydrologicModel can be put in command line -> easier batch-runs
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
if use_RFS == False:
	sys.exit('running GLOFRIM with RFS deactivated not possible anymore... Bummer!')
verbose = strtobool(config.general_settings['verbose'])
use_floodplain_infiltration_factor = strtobool(config.general_settings['use_floodplain_infiltration_factor'])
adjust_initial_groundwater = strtobool(config.general_settings['adjust_initial_groundwater'])

# -------------------------------------------------------------------------------------------------
# SPECIFY NUMERICAL SETTINGS
# -------------------------------------------------------------------------------------------------

nr_model_timesteps                    = int(config.numerical_settings['number_of_timesteps'])                      
update_step                           = int(config.numerical_settings['update_step'])  
                      
secPerDay                             = 86400.
end_time 							  = nr_model_timesteps * secPerDay
fraction_timestep 					  = secPerDay / update_step

threshold_inundated_depth             = float(config.numerical_settings['threshold_inundated_depth'])                                                   

# other
missing_value_landmask                = 255
missing_value_pcr                     = -999

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
hydrologicModel = pcrglobwb_bmi_v203.pcrglobwb_bmi.pcrglobwbBMI()
hydrologicModel.initialize(config_pcr)
print '\n>>> PCR Initialized <<<\n' 

# spin-up PCR-GLOBWB
hydrologicModel.spinup()

# -------------------------------------------------------------------------------------------------
# INITIALIZING HYDRODYNAMIC MODEL
# -------------------------------------------------------------------------------------------------

# initiate hydraulic model
hydrodynamicModel = bmi.wrapper.BMIWrapper(engine = model_path, configfile = (os.path.join(model_dir, model_file)))
hydrodynamicModel.initialize()
print '\n>>> ',model_type,' Initialized <<<\n' 

# -------------------------------------------------------------------------------------------------
# EXCTRACTING RELEVANT DATA FROM MODELS
# -------------------------------------------------------------------------------------------------

if model_type == 'DFM':

    #- retrieving data from Delft3D FM    
    x_coords, y_coords, z_coords, bottom_lvl, cell_points_fm, separator_1D, cellAreaSpherical_1D, cellAreaSpherical_2D, xz_coords, yz_coords, hydrodynamic_coords_1D, hydrodynamic_coords_2D,\
                cellarea_data_pcr, landmask_data_pcr, clone_data_pcr = model_functions.extractModelData_FM(hydrodynamicModel, hydrologicModel, landmask_pcr, clone_pcr, use_RFS, use_2way)
    coupledFPindices = 0.
    print '\n>>> DFM data retrieved <<<\n'
         
elif model_type == 'LFP':
    
    #- retrieving data from LISFLOOD-FP
    dx, dy, DEM, bottom_lvl, H, waterDepth, rows, cols, \
                list_x_coords, list_x_coords_2way, list_y_coords, list_y_coords_2way, coupledFPindices, coupledFPindices_2way, grid_dA, cellAreaSpherical_1D, SGCQin, separator_1D,\
                cellarea_data_pcr, landmask_data_pcr, clone_data_pcr = model_functions.extractModelData_FP(hydrodynamicModel, model_dir, hydrologicModel, landmask_pcr, clone_pcr, verbose_folder, use_RFS, use_2way, verbose)
                
    cellAreaSpherical_2D = np.copy(cellAreaSpherical_1D)
                
    #- computing FP-coordinates    
    hydrodynamic_coords_1D = coupling_functions.getVerticesFromMidPoints(list_x_coords, list_y_coords, dx, dy, verbose)
    hydrodynamic_coords_2D = coupling_functions.getVerticesFromMidPoints(list_x_coords_2way, list_y_coords_2way, dx, dy, verbose)
    print '\n>>> LFP data retrieved <<<\n'

#- computing PCR-coordinates
PCRcoords = coupling_functions.getPCRcoords(landmask_data_pcr)
print '\n>>> PCR data retrieved <<<\n'
		
# -------------------------------------------------------------------------------------------------
# COUPLING THE GRIDS
# -------------------------------------------------------------------------------------------------

# this is only required for plotting later, not for actual coupling process
# remove once code works! not needed for actual applicability

CoupledCellsInfoAll = coupling_functions.coupleAllCells(hydrodynamic_coords_1D,PCRcoords)
CoupledCellsInfoAll_2way = coupling_functions.coupleAllCells(hydrodynamic_coords_2D,PCRcoords)

# converting single indices of coupled PCR cells to double (array,column) indices

if use_2way == False:
	# coupled_hydrodynamicModel_2_hydrologicModel: list containing indices of all coupled hydrodynamic cells with their corresponding index of coupled hydrologic cell (one-to-one)
	# coupled_hydrologicModel_2_hydrodynamicModel: list containing indices of all coupled hydrologic cells with their corresponding indices of coupled hydrodynamic cells (one-to-many_
	# coupled_hydrologicModel_indices: list containing indices of all coupled PCR cells
	# linking PCR-cells and 1d-channels (DFM) or river cells (LFP)
    couple_hydrodynamicModel_2_hydrologicModel, couple_hydrologicModel_2_hydrodynamicModel, coupled_hydrologicModel_indices = coupling_functions.assignPCR2cells(landmask_pcr, hydrodynamic_coords_1D, verbose)
    
elif use_2way == True:
	# idem
	# all variables "*_2way" are obtained by linking PCR-cells with 2d-cells (DFM) or all cells not defined as river cells (LFP)
    couple_hydrodynamicModel_2_hydrologicModel, couple_hydrologicModel_2_hydrodynamicModel, coupled_hydrologicModel_indices = coupling_functions.assignPCR2cells(landmask_pcr, hydrodynamic_coords_1D, verbose)
    couple_hydrodynamicModel_2_hydrologicModel_2way, couple_hydrologicModel_2_hydrodynamicModel_2way, coupled_hydrologicModel_indices_2way = coupling_functions.assignPCR2cells(landmask_pcr, hydrodynamic_coords_2D, verbose)

# saving plots of coupled cells to verbose-folder
# currently doesn't work with FM and use_RFS on, due to data structure required (? check this ?)
# remove once code works! not needed for actual applicability

if (model_type == 'DFM') and (use_RFS == True):
    pass
    
else:
	if (verbose == True) and (use_2way == False):
		coupling_functions.plotGridfromCoords(PCRcoords, hydrodynamic_coords_1D)
		plt.savefig(os.path.join(verbose_folder , 'AllCells_1way.png'))
		coupling_functions.plotGridfromCoords(CoupledCellsInfoAll[1],CoupledCellsInfoAll[0])
		plt.savefig(os.path.join(verbose_folder , 'CoupledCells_1way.png'))   
		plt.close('all')
	elif (verbose == True) and (use_2way == True):
		coupling_functions.plotGridfromCoords(PCRcoords, hydrodynamic_coords_1D)
		plt.savefig(os.path.join(verbose_folder , 'AllCells_1way.png'))
		coupling_functions.plotGridfromCoords(CoupledCellsInfoAll[1],CoupledCellsInfoAll[0])
		plt.savefig(os.path.join(verbose_folder , 'CoupledCells_1way.png'))   
		plt.close('all')    
		coupling_functions.plotGridfromCoords(PCRcoords, hydrodynamic_coords_2D)
		plt.savefig(os.path.join(verbose_folder , 'AllCells_2way.png'))
		coupling_functions.plotGridfromCoords(CoupledCellsInfoAll_2way[1],CoupledCellsInfoAll_2way[0])
		plt.savefig(os.path.join(verbose_folder , 'CoupledCells_2way.png'))   
		plt.close('all')
    
# -------------------------------------------------------------------------------------------------
# TURNING OFF CHANNELSTORAGE AND WATERBODYSTORAGE
# 
# for PCR <-> 1D 
# this is required to avoid the entire water volume in channels and waterbodies to empty at once into DFM/LFP
# besides, coupling is activated for relevant sections (has to be de-activated initially, otherwise spin-up does not work)
# ------------------------------------------------------------------------------------------------- 

model_functions.noStorage(hydrologicModel, missing_value_pcr, coupled_hydrologicModel_indices, couple_hydrologicModel_2_hydrodynamicModel)

# -------------------------------------------------------------------------------------------------
# TURNING OFF ROUTING BY PCR IN COUPLED AREA
# 
# for PCR <-> 1D 
# this is required to avoid double-routing by DFM/LFP and PCR (as we want to perform the routing by hydrodynamicModel only!)
# -------------------------------------------------------------------------------------------------

model_functions.noLDD(hydrologicModel, coupled_hydrologicModel_indices, verbose_folder, verbose)

# -------------------------------------------------------------------------------------------------
# ACTIVATING A RANGE OF VARIABLES SPECIFICALLY REQUIRED FOR 2WAY-COUPLING
#
# for PCR <-> 2D
# optional
# keep it turned off until remaining coupling scheme works!
# -------------------------------------------------------------------------------------------------

model_functions.adjust_iniGR(hydrologicModel, adjust_initial_groundwater_file, coupled_hydrologicModel_indices_2way, adjust_initial_groundwater=False) # optional

model_functions.activate_floodplain_infiltration_factor(hydrologicModel, coupled_hydrologicModel_indices_2way, use_floodplain_infiltration_factor=False) # optional

# -------------------------------------------------------------------------------------------------
# CALCULATE DELTA VOLUMES (DAY 1)
#
# first day outside loop, to make sure PCR time is at start time (timestep 1) and not at end of spin-up (timestep 365)
# ------------------------------------------------------------------------------------------------- 

# using a dummy here (creates 0's for all entries in array) to avoid feeding back volumes to PCR before 1st update actually took place
water_volume_FM_2_PCR = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)

# setting PCR's topWaterLayer to 0 for all PCR cells coupled to hydrodynamicModel (1way)
model_functions.set_zeroTopWaterlayer(hydrologicModel, coupled_hydrologicModel_indices)

# retrieving PCR-GLOBWB and converting it to m3/d
# for PCR <-> 1D
# delta volume = volume PCR - water_volume_FM_2_PCR
delta_volume_PCR, delta_volume_PCR_1way = model_functions.calculateDeltaVolumes(hydrologicModel, missing_value_pcr, secPerDay, coupled_hydrologicModel_indices, cellarea_data_pcr, water_volume_FM_2_PCR)

# dividing delta volume from PCR-GLOBWB over hydraulic cells
# for PCR <-> 1D
# if applied to DFM/LFP as flux, then delta_volume_PCR_1way is used
# if applied to DFM/LFP as state, then delta_volume_PCR_1way needs to be converted to delta_water_DFM_1way
delta_water_DFM_1way, verbose_volume_DFM_1way = model_functions.calculateDeltaWater(hydrodynamicModel, couple_hydrodynamicModel_2_hydrologicModel, couple_hydrologicModel_2_hydrodynamicModel, delta_volume_PCR_1way, cellAreaSpherical_1D, fraction_timestep, model_type, use_Fluxes)

# saving PCR-GLOBWB output volumes and volumes used as input to hydraulic models to verbose-folder
if verbose == True:
	# initial file objects
    fo_PCR_V_tot = open(os.path.join(verbose_folder, 'delta_volume_PCR_1way.txt'), 'w')
    fo_verbose_volume_DFM_1way = open(os.path.join(verbose_folder, 'verbose_volume_DFM_1way.txt'), 'w')
    # aggregate daily totals
    delta_volume_PCR_total_aggr = np.sum(delta_volume_PCR_1way)
    verbose_volume_DFM_1way_aggr = np.sum(verbose_volume_DFM_1way)
    # write daily totals to file
    fo_PCR_V_tot.write(str(delta_volume_PCR_total_aggr) + os.linesep)
    fo_verbose_volume_DFM_1way.write(str(verbose_volume_DFM_1way_aggr) + os.linesep)

# check to ensure that volumes are equal, i.e. no errors in water balance
if np.round(np.sum(np.asfarray(delta_volume_PCR_1way)),4) != np.round(np.sum(np.asfarray(verbose_volume_DFM_1way)), 4):
    print 'PCR volume out: ', np.sum(np.asfarray(delta_volume_PCR_1way))
    print 'FM volume in: ', np.sum(np.asfarray(verbose_volume_DFM_1way))
    sys.exit('\nERROR: input volumes to hydrodynamicModel do not match PCR output volumes!\n')

# reshaping data for LISFLOOD-FP from list to arrays
# from flattened array to actual 2D-array
if model_type == 'LFP':
    delta_water_DFM_1way = model_functions.fillLFPgrid(hydrodynamicModel, coupledFPindices, delta_water_DFM_1way, DEM, verbose_folder, verbose)
  
# -------------------------------------------------------------------------------------------------
# FIRST UPDATE (DAY 1)
# -------------------------------------------------------------------------------------------------

# overwriting variables in LFP
# separated from DFM function because LFP can only be updated on daily time step
if (model_type == 'LFP'):
    model_functions.updateModel(hydrodynamicModel, delta_water_DFM_1way, update_step, separator_1D, use_Fluxes, use_RFS, model_type, verbose)

while hydrodynamicModel.get_current_time() < (hydrologicModel.get_time_step() * secPerDay):
    
    # overwriting variables and updating DFM on user-specified time step
    if (model_type == 'DFM'):
        model_functions.updateModel(hydrodynamicModel, delta_water_DFM_1way, update_step, separator_1D, use_Fluxes, use_RFS, model_type, verbose)
    
    # updating LFP on daily time step
    if (model_type == 'LFP'):
        hydrodynamicModel.update()

# ----------------------------------------------------------------------------------------------------
# UPDATING REMAINING MODELLING TIME
# ----------------------------------------------------------------------------------------------------

while hydrologicModel.get_time_step() < nr_model_timesteps:
	
    # COMPUTING INUNDATION AS AREA AND AS FRACTION OF PCR-CELL
    # for PCR <-> 1D
    # for PCR <-> 2D
    # returns both as list for coupled PCR cells only ('*_coupled') and as 2d-arrays to be set back into PCR via BMI
    inundated_area_FM_2_PCR_coupled_1D, inundated_area_FM_2_PCR_coupled, inundated_area_FM_2_PCR, inundated_fraction_FM_2_PCR_coupled, inundated_fraction_FM_2_PCR = model_functions.determine_InundationArea_Hydrodynamics(model_type, 
                                                                                                                                                                                        hydrodynamicModel, 
                                                                                                                                                                                        couple_hydrologicModel_2_hydrodynamicModel,
                                                                                                                                                                                        coupled_hydrologicModel_indices,
																																														couple_hydrologicModel_2_hydrodynamicModel_2way, 
																																														coupled_hydrologicModel_indices_2way, 
																																														threshold_inundated_depth, 
																																														cellAreaSpherical_1D,
																																														cellAreaSpherical_2D, 
																																														cellarea_data_pcr, 
																																														landmask_pcr, 
																																														missing_value_landmask)
    
    # testing output, but only for one time step
    if hydrologicModel.get_time_step() == 2.:
		plt.figure()
		plt.imshow(inundated_area_FM_2_PCR)
		plt.savefig(os.path.join(verbose_folder , 'inundated_area_FM_2_PCR.png'))
		plt.figure()
		plt.imshow(inundated_fraction_FM_2_PCR)
		plt.savefig(os.path.join(verbose_folder , 'inundated_fraction_FM_2_PCR.png'))
		plt.close('all')
    

    # COMPUTING WATER VOLUME AND DEPTH TO BE COUPLED BACK TO PCR
    # for PCR <-> 1D
    # for PCR <-> 2D
    # returns only 2d-arrays to be set back into PCR via BMI
    water_volume_FM_2_PCR_1D, water_depths_FM_2_PCR_1D, water_volume_FM_2_PCR, water_depths_FM_2_PCR  = model_functions.determine_InundationDepth_Hydrodynamics(model_type, 
																											hydrodynamicModel, 
																											landmask_pcr, 
																											missing_value_landmask,
																											inundated_area_FM_2_PCR_coupled_1D, 
																											inundated_area_FM_2_PCR_coupled,
																											couple_hydrologicModel_2_hydrodynamicModel,
																											coupled_hydrologicModel_indices, 
																											couple_hydrologicModel_2_hydrodynamicModel_2way, 
																											coupled_hydrologicModel_indices_2way)

    # DETERMINING NEW STORAGE IN PCR-CHANNELS
    # for PCR <-> 1D
    # returns 2d-array to be set back into PCR via BMI
    new_storage_pcr = model_functions.determine_new_channelStoragePCR(hydrologicModel, 
                                                                      landmask_pcr, 
                                                                      missing_value_landmask, 
                                                                      water_volume_FM_2_PCR_1D)
    
    # UPDATE VARIABLES IN PCR
    # for PCR <-> 1D
    # for PCR <-> 2D
    # doesn't return anything; updates several variables in PCR based on previously determined 2d-arrays
    # REMARK: is it actually necessary to update channel storage in PCR?
    model_functions.updateHydrologicVariables(hydrologicModel, 
                                              water_depths_FM_2_PCR, 
                                              inundated_fraction_FM_2_PCR, 
                                              new_storage_pcr)
    
    # retrieving PCR-GLOBWB and converting it to m3/d
    # for PCR <-> 1D
    delta_volume_PCR, delta_volume_PCR_coupled_1D = model_functions.calculateDeltaVolumes(hydrologicModel, 
																						  missing_value_pcr, 
                                                                                          secPerDay, 
                                                                                          coupled_hydrologicModel_indices, 
                                                                                          cellarea_data_pcr, 
                                                                                          water_volume_FM_2_PCR) 
	
    # removing negative delta volumes where necessary
#    delta_volume_PCR_positiveOnly, delta_volume_PCR_1way_positiveOnly = model_functions.account4negativeDeltaVolumes(hydrodynamicModel, 
#                                                                                                                     model_type, 
#                                                                                                                     coupled_hydrologicModel_indices, 
#                                                                                                                     coupled_hydrologicModel_indices_2way, 
#                                                                                                                     couple_hydrologicModel_2_hydrodynamicModel_2way, 
#                                                                                                                     delta_volume_PCR, 
#                                                                                                                     cellAreaSpherical_2D)
				                 
    # dividing delta volume from PCR-GLOBWB over hydraulic cells, depending on model specifications
    delta_water_DFM_1way, verbose_volume_DFM_1way = model_functions.calculateDeltaWater(hydrodynamicModel, 
                                                                                        couple_hydrodynamicModel_2_hydrologicModel, 
                                                                                        couple_hydrologicModel_2_hydrodynamicModel, 
                                                                                        delta_volume_PCR_1way_positiveOnly, 
                                                                                        cellAreaSpherical_1D, 
                                                                                        fraction_timestep, 
                                                                                        model_type, 
                                                                                        use_Fluxes)
    
    # saving PCR-GLOBWB output volumes and volumes used as input to hydraulic models to verbose-folder
    if verbose == True:
        # aggregate daily totals
        delta_volume_PCR_total_aggr = np.sum(delta_volume_PCR_1way)
        verbose_volume_DFM_1way_aggr = np.sum(verbose_volume_DFM_1way)
        # write daily totals to file
        fo_PCR_V_tot.write(str(delta_volume_PCR_total_aggr) + os.linesep)
        fo_verbose_volume_DFM_1way.write(str(verbose_volume_DFM_1way_aggr) + os.linesep)
    
    # reshaping data for LISFLOOD-FP from list to arrays
    if model_type == 'LFP':
        delta_water_DFM_1way = model_functions.fillLFPgrid(hydrodynamicModel, coupledFPindices, delta_water_DFM_1way, DEM, verbose_folder, verbose)  
    
    # updating arrays with computed additional volumes; array used depends on model specifications
    if (model_type == 'LFP'):
        model_functions.updateModel(hydrodynamicModel, delta_water_DFM_1way, update_step, separator_1D, use_Fluxes, use_RFS, model_type, verbose)      

    # update FM unless it has has reached the same time as PCR
    while hydrodynamicModel.get_current_time() < (hydrologicModel.get_time_step() * secPerDay):
        
        # updating FM on user-specified time step
        if (model_type == 'DFM'):
            model_functions.updateModel(hydrodynamicModel, delta_water_DFM_1way, update_step, separator_1D, use_Fluxes, use_RFS, model_type, verbose)      
        
        # updating FM or FP on daily time step
        if (model_type == 'LFP'):
            hydrodynamicModel.update()   
        
# ----------------------------------------------------------------------------------------------------
# END OF MODEL PERIOD REACHED
# ----------------------------------------------------------------------------------------------------
    
# update and finalize logging
model_functions.write2log(model_dir, model_file, latlon, use_2way, use_Fluxes, use_RFS, use_floodplain_infiltration_factor, adjust_initial_groundwater, verbose, moment='end') 

# close files
if verbose == True:
    fo_PCR_V_tot.close()
    fo_verbose_volume_DFM_1way.close()

#- finalizing hydrodynamic model
hydrodynamicModel.finalize()
print '\nModel End-Time: ', datetime.datetime.now()
