# -*- coding: utf-8 -*-
"""
New approach for two-way coupling: first update hydrodynamics, then align hydrology.
This reduces the risk of numerical instabilites.
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
if use_2way == True:
	use_RFS = True

verbose = strtobool(config.general_settings['verbose'])

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

# hydrodynamics
model_dir       	= config.model_settings['model_dir']
model_file      	= config.model_settings['model_file']
model_proj			= config.model_settings['model_projection']

# hydrology
config_pcr       	=  config.PCR_settings['config_pcr']
configPCR           = configuration.Configuration()
configPCR.parse_configuration_file(config_pcr )
inputDIR 			= configPCR.globalOptions['inputDir']
clone_pcr 			= os.path.join(inputDIR, configPCR.globalOptions['cloneMap'])
landmask_pcr 		= os.path.join(inputDIR, configPCR.globalOptions['landmask'])

# -------------------------------------------------------------------------------------------------
# SET PATHS TO .SO / .DLL FILES
# -------------------------------------------------------------------------------------------------

# these may be changed according to personal file and folder structure
if model_type == 'DFM':
    model_path = '/home/jannis/Programmes/DFLOWFM/lib/libdflowfm.so'

elif model_type == 'LFP':
    model_path = '/home/jannis/Programmes/LISFLOODFP/lisflood-bmi-v5.9/liblisflood.so'

else:
    sys.exit('\nno adequate model defined in configuration file - define either DFM or LFP!\n')

# -------------------------------------------------------------------------------------------------
# INITIALIZE AND SPIN-UP PCR-GLOBWB
# -------------------------------------------------------------------------------------------------

# initiate logging and define folder for verbose-output
verbose_folder = model_functions.write2log(model_dir, model_file, latlon, use_2way, use_Fluxes, use_RFS, use_floodplain_infiltration_factor=False, adjust_initial_groundwater=False, verbose=False, moment='start')
print 'Model Start-Time: ', datetime.datetime.now()
print ''

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
#

if model_type == 'DFM':

    #- retrieving data from Delft3D FM
    bottom_lvl, cell_points_fm, separator_1D, cellAreaSpherical_1D, cellAreaSpherical_2D, xz_coords, yz_coords, \
					hydrodynamic_coords_1D, hydrodynamic_coords_2D = model_functions.extractModelData_DFM(hydrodynamicModel, use_RFS, use_2way)

    coupledFPindices = 0.
    print '\n>>> DFM data retrieved <<<\n'

elif model_type == 'LFP':

    #- retrieving data from LISFLOOD-FP
    dx, dy, DEM, bottom_lvl, H, waterDepth, rows, cols, \
                list_x_coords, list_x_coords_2way, list_y_coords, list_y_coords_2way, coupledFPindices, coupledFPindices_2way, grid_dA, cellAreaSpherical_1D, \
				SGCQin, separator_1D, = model_functions.extractModelData_LFP(hydrodynamicModel, model_dir, verbose_folder, use_RFS, use_2way, verbose)

    cellAreaSpherical_2D = np.copy(cellAreaSpherical_1D)

    #- computing FP-coordinates
    hydrodynamic_coords_1D = coupling_functions.getVerticesFromMidPoints(list_x_coords, list_y_coords, dx, dy, verbose)
    hydrodynamic_coords_2D = coupling_functions.getVerticesFromMidPoints(list_x_coords_2way, list_y_coords_2way, dx, dy, verbose)
    print '\n>>> LFP data retrieved <<<\n'

#- computing PCR-coordinates
cellarea_data_pcr, landmask_data_pcr, clone_data_pcr = model_functions.extractModelData_PCR(hydrologicModel, landmask_pcr, clone_pcr)
PCRcoords = coupling_functions.getPCRcoords(landmask_data_pcr)
print '\n>>> PCR data retrieved <<<\n'

# -------------------------------------------------------------------------------------------------
# ACTIVATING COUPLING FOR RELEVANT PCR SECTIONS
# -------------------------------------------------------------------------------------------------

hydrologicModel.set_var(('grassland','ActivateCoupling'), 'True') #2way
hydrologicModel.set_var(('forest','ActivateCoupling'), 'True') #2way
hydrologicModel.set_var(('routing','ActivateCoupling'), 'True')
hydrologicModel.set_var(('WaterBodies', 'ActivateCoupling'), 'True')

# -------------------------------------------------------------------------------------------------
# COUPLING THE GRIDS
# -------------------------------------------------------------------------------------------------

# converting single indices of coupled PCR cells to double (array,column) indices
# couple_HDYN_2_HLOG: list containing indices of all coupled hydrodynamic cells with their corresponding index of coupled hydrologic cell (one-to-one)
# couple_HLOG_2_HDYN: list containing indices of all coupled hydrologic cells with their corresponding indices of coupled hydrodynamic cells (one-to-many_
# coupled_HLOG_indices: list containing indices of all coupled PCR cells

if use_2way == False:
	# linking PCR-cells with 1D hydrodynamic cells (if use_RFS=True) or 2D hydrodynamic cells (if use_RFS=False)
    couple_HDYN_2_HLOG, couple_HLOG_2_HDYN, coupled_HLOG_indices = coupling_functions.assignPCR2cells(landmask_pcr, hydrodynamic_coords_1D, verbose)

elif use_2way == True:
	# linking PCR-cells with 1D hydrodynamic cells (since use_RFS=True)
    couple_HDYN_2_HLOG, couple_HLOG_2_HDYN, coupled_HLOG_indices = coupling_functions.assignPCR2cells(landmask_pcr, hydrodynamic_coords_1D, verbose)
	# linking PCR-cells with 2D hydrodynamic cells
    couple_HDYN_2_HLOG_2way, couple_HLOG_2_HDYN_2way, coupled_HLOG_indices_2way = coupling_functions.assignPCR2cells(landmask_pcr, hydrodynamic_coords_2D, verbose)

# -------------------------------------------------------------------------------------------------
# TURNING OFF ROUTING BY PCR IN COUPLED AREA
#
# for 1D <-> PCR
# this is required to avoid double-routing by HDYN and HLOG (as we want to perform the routing by HDYN only!)
# -------------------------------------------------------------------------------------------------

model_functions.noLDD(hydrologicModel, coupled_HLOG_indices, verbose_folder, verbose)

# -------------------------------------------------------------------------------------------------
# COMPUTING INUNDATION AS AREA AND AS FRACTION OF PCR-CELL
#
# inundatedArea_HDYN1D_2_HLOG			= list with upscaled inundated area of all 1D hydrodynamic cells within a hydrologic cell
# inundatedArea_HDYN2D_2_HLOG			= list with upscaled inundated area of all 2D hydrodynamic cells within a hydrologic cell
# inundatedArea_HDYN2D_2_HLOG_BMI		= array of upscaled inundated area of all 2D hydrodynamic cells within a hydrologic cell; to be used with BMI-command
# inundatedFraction_HDYN2D_2_HLOG		= list with fractional inundated area of all 2D hydrodynamic cells per hydrologic cell
# inundatedFraction_HDYN2D_2_HLOG_BMI	= array of fractional inundated area of all 2D hydrodynamic cells per hydrologic cell; to be used with BMI-command
# -------------------------------------------------------------------------------------------------

inundatedArea_HDYN1D_2_HLOG, inundatedArea_HDYN2D_2_HLOG, inundatedArea_HDYN2D_2_HLOG_BMI, \
			inundatedFraction_HDYN2D_2_HLOG, inundatedFraction_HDYN2D_2_HLOG_BMI = model_functions.determine_InundationArea_Hydrodynamics(model_type,
                                                                                                                                        hydrodynamicModel,
                                                                                                                                        couple_HLOG_2_HDYN,
                                                                                                                                        coupled_HLOG_indices,
																																		couple_HLOG_2_HDYN_2way,
																																		coupled_HLOG_indices_2way,
																																		threshold_inundated_depth,
																																		cellAreaSpherical_1D,
																																		cellAreaSpherical_2D,
																																		cellarea_data_pcr,
																																		landmask_pcr,
																																		missing_value_landmask)
plt.figure()
plt.imshow(inundatedArea_HDYN2D_2_HLOG_BMI)
plt.colorbar()
plt.savefig(os.path.join(verbose_folder,'inundatedArea_HDYN2D_2_HLOG_BMI.png'), dpi=300)
plt.figure()
plt.imshow(inundatedFraction_HDYN2D_2_HLOG_BMI)
plt.colorbar()
plt.savefig(os.path.join(verbose_folder,'inundatedFraction_HDYN2D_2_HLOG_BMI.png'), dpi=300)

# -------------------------------------------------------------------------------------------------
# COMPUTING WATER VOLUME AND DEPTH TO BE COUPLED BACK TO PCR
#
# waterVolume_HDYN1D_2_HLOG_BMI		= array of accumulated water volume of all 1D hydrodynamic cells within a hydrologic cell
# waterVolume_HDYN2D_2_HLOG_BMI		= array of accumulated water volume of all 2D hydrodynamic cells within a hydrologic cell
# waterDepths_HDYN2D_2_HLOG_BMI		= array of accumulated water depth of all 2D hydrodynamic cells within a hydrologic cell (waterdepth = Vhdyn/Ainund.hdyn)
# -------------------------------------------------------------------------------------------------

waterVolume_HDYN1D_2_HLOG_BMI, waterVolume_HDYN2D_2_HLOG_BMI, waterDepths_HDYN2D_2_HLOG_BMI  = model_functions.determine_InundationDepth_Hydrodynamics(model_type,
																										hydrodynamicModel,
																										landmask_pcr,
																										missing_value_landmask,
																										inundatedArea_HDYN1D_2_HLOG,
																										inundatedArea_HDYN2D_2_HLOG,
																										couple_HLOG_2_HDYN,
																										coupled_HLOG_indices,
																										couple_HLOG_2_HDYN_2way,
																										coupled_HLOG_indices_2way)
plt.figure()
plt.imshow(waterVolume_HDYN1D_2_HLOG_BMI)
plt.colorbar()
plt.savefig(os.path.join(verbose_folder,'waterVolume_HDYN1D_2_HLOG_BMI.png'), dpi=300)
plt.figure()
plt.imshow(hydrologicModel.get_var('channelStorage'))
plt.colorbar()
plt.savefig(os.path.join(verbose_folder,'initialChannelStoragePCR.png'), dpi=300)
plt.figure()
plt.imshow(inundatedFraction_HDYN2D_2_HLOG_BMI)
plt.colorbar()
plt.savefig(os.path.join(verbose_folder,'waterVolume_HDYN2D_2_HLOG_BMI.png'), dpi=300)
plt.figure()
plt.imshow(inundatedFraction_HDYN2D_2_HLOG_BMI)
plt.colorbar()
plt.savefig(os.path.join(verbose_folder,'waterDepths_HDYN2D_2_HLOG_BMI.png'), dpi=300)

# -------------------------------------------------------------------------------------------------
# UPDATING CHANNELSTORAGE AND TURNING OFF WATERBODYSTORAGE
#
# for 1D HDYN <-> HLOG
# this is required to avoid the entire water volume in channels and waterbodies to empty at once into DFM/LFP
# also, channel storage in HLOG will be updated based on HDYN after first update
# -------------------------------------------------------------------------------------------------

model_functions.updateStorage(hydrologicModel, landmask_pcr, missing_value_pcr, missing_value_landmask, coupled_HLOG_indices, couple_HLOG_2_HDYN, waterVolume_HDYN1D_2_HLOG_BMI)

plt.figure()
plt.imshow(hydrologicModel.get_var('channelStorage'))
plt.colorbar()
plt.savefig(os.path.join(verbose_folder,'updatedChannelStoragePCR.png'), dpi=300)

### - model works until here - ###
pdb.set_trace()

# -------------------------------------------------------------------------------------------------
# CALCULATE DELTA VOLUMES (DAY 1)
#
# first day outside loop, to make sure PCR time is at start time (timestep 1) and not at end of spin-up (timestep 365)
# -------------------------------------------------------------------------------------------------

# using a dummy here (creates 0's for all entries in array) to avoid feeding back volumes to PCR before 1st update actually took place
waterVolume_HDYN2D_2_HLOG_BMI = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)

# setting PCR's topWaterLayer to 0 for all PCR cells coupled to hydrodynamicModel (1way)
model_functions.set_zeroTopWaterlayer(hydrologicModel, coupled_HLOG_indices)

# retrieving PCR-GLOBWB and converting it to m3/d
# for PCR <-> 1D
# delta volume = volume PCR - waterVolume_HDYN2D_2_HLOG_BMI
delta_volume_PCR, delta_volume_PCR_1way = model_functions.calculateDeltaVolumes(hydrologicModel, missing_value_pcr, secPerDay, coupled_HLOG_indices, cellarea_data_pcr, waterVolume_HDYN2D_2_HLOG_BMI)

# dividing delta volume from PCR-GLOBWB over hydraulic cells
# for PCR <-> 1D
# if applied to DFM/LFP as flux, then delta_volume_PCR_1way is used
# if applied to DFM/LFP as state, then delta_volume_PCR_1way needs to be converted to delta_water_DFM_1way
delta_water_DFM_1way, verbose_volume_DFM_1way = model_functions.calculateDeltaWater(hydrodynamicModel, couple_HDYN_2_HLOG, couple_HLOG_2_HDYN, delta_volume_PCR_1way, cellAreaSpherical_1D, fraction_timestep, model_type, use_Fluxes)

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






    # UPDATE VARIABLES IN PCR
    # for PCR <-> 1D
    # for PCR <-> 2D
    # doesn't return anything; updates several variables in PCR based on previously determined 2d-arrays
    # REMARK: is it actually necessary to update channel storage in PCR?
    # QUESTION: why doing this before actually accounting for negative water volumes? wouldn't it make more sense to first adapt water levels in DFM and then determine water depth and inundation fraction?
	#           and these water levels/fraction would then be set back to PCR to force MODFLOW....
    model_functions.updateHydrologicVariables(hydrologicModel,
                                              waterDepths_HDYN2D_2_HLOG_BMI,
                                              inundatedFraction_HDYN2D_2_HLOG_BMI,
                                              new_storage_pcr,
                                              couple_channel_storage)

    # retrieving PCR-GLOBWB and converting it to m3/d
    # for PCR <-> 1D
    delta_volume_PCR, delta_volume_PCR_coupled_1D = model_functions.calculateDeltaVolumes(hydrologicModel,
																						  missing_value_pcr,
                                                                                          secPerDay,
                                                                                          coupled_HLOG_indices,
                                                                                          cellarea_data_pcr,
                                                                                          waterVolume_HDYN2D_2_HLOG_BMI)

    # removing negative delta volumes where necessary
    # QUESTION: would it be possible to just adapt volumes in DFM based on delta volume? what would happen then with water levels if volume is added/removed? are they automatically increased/lowered?
    delta_volume_PCR_positiveOnly, delta_volume_PCR_1way_positiveOnly = model_functions.account4negativeDeltaVolumes(hydrodynamicModel,
                                                                                                                     model_type,
                                                                                                                     coupled_HLOG_indices,
                                                                                                                     coupled_HLOG_indices_2way,
                                                                                                                     couple_HLOG_2_HDYN_2way,
                                                                                                                     delta_volume_PCR,
                                                                                                                     cellAreaSpherical_2D)

    # dividing delta volume from PCR-GLOBWB over hydraulic cells, depending on model specifications
    delta_water_DFM_1way, verbose_volume_DFM_1way = model_functions.calculateDeltaWater(hydrodynamicModel,
                                                                                        couple_HDYN_2_HLOG,
                                                                                        couple_HLOG_2_HDYN,
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
