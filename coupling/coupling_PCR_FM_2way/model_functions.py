# -*- coding: utf-8 -*-
"""
Set of functions related to the coupling of PCR-GLOBWB with Delft 3D Flexible Mesh and LISFLOOD-FP.
Unlike the "coupling_functions" where primarily the grids are spatially coupled,
this set of functions executes the actual model functions to manipulate model
arrays and to compute the required volumes and depths, as well as updating
the models. 

Please note: currently only to be used for 1way-coupling!!!
Development for 2way-coupling is ongoing.

@author: Jannis Hoch M.Sc., Department of Physical Geography, Faculty of Geosciences, Utrecht University (j.m.hoch@uu.nl)
@date: 08-05-2017
"""

import numpy as np
import datetime
import os
import pdb
import sys
import pcraster as pcr
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from coupling_PCR_FM import coupling_functions

# =============================================================================

def set_values_in_array(vals, idx, update_vals):
    """
    Sets new values in an existing array at specified locations. Locations are provided by a list of (y, x) indices
    :param vals: array (numpy) of values that will be updated
    :param idx: list of (y, x) indices
    :param update_vals: single value or list of values (equal size as idx) for updating
    :return: vals: updated list
    """
    vals[zip(*idx)] = update_vals
    return vals

# =============================================================================

def write2log(model_dir, model_file, latlon, use_2way, useFluxes, use_RFS, use_floodplain_infiltration_factor, verbose, moment):
    """
    Writing model settings/paths/etc to a txt-file in directory.
    Note that PCR-GLOBWB is additionally writing its own log-file and Delft3D DFM adds model
    diagnosis and settings to the dia-file.
    Also verbose-output (if verbose == True) is stored in the created directory.
    
    TO DO: avoid printing of model settings at end of run if not total number of timesteps is specified in set-file
    """
    
    #- creating folder to store verbose and log output
    folder_name = model_dir + model_file + '_verboseOut'
        
    #-write end time of simulation finished when finished    
    if moment == 'start':
        #- print info to console at start of simulation
        print '\n##############################'
        print '### MODEL COUPLING STARTED ###'
        print '##############################'
        print '\nmodel file chosen: ', model_file
        print '2way coupling activated: ', bool(use_2way)
        print 'lat-lon on: ', bool(latlon)
        print 'fluxes on: ', bool(useFluxes)
        print 'RFS on: ', bool(use_RFS)
        print 'use floodplain infiltration factor: ', bool(use_floodplain_infiltration_factor)
        print 'verbose mode on: ', bool(verbose)
        print '\nVerbose Output and Log-File saved in: ', folder_name + os.linesep
		
    elif moment == 'end':
        #- check if folder has already been created
        if os.path.exists(folder_name):
            pass
        else:
            os.mkdir(folder_name)
			
		#- crate log-file
        fo_name = model_file + '_couplingLog.txt'
        fo = open((os.path.join(folder_name, fo_name)),'w')
		
        #- write info to log-file
        fo.write('\nmodel_file chosen: ')
        fo.write(model_file + os.linesep)
        fo.write('model start-time: ')
        fo.write(str(datetime.datetime.now()) + os.linesep)
        fo.write('2way activated: ')
        fo.write(str(bool(use_2way)) + os.linesep)
        fo.write('lat-lon activated: ')
        fo.write(str(bool(latlon)) + os.linesep)
        fo.write('forcing by fluxes activated: ')
        fo. write(str(bool(useFluxes)) + os.linesep)
        fo.write('river-floodplain-scheme activated: ')
        fo.write(str(bool(use_RFS)) + os.linesep)
        fo.write('model end-time:')
        fo.write(str(datetime.datetime.now()) + os.linesep)
		#- close log-file    
        fo.close()
    
    return folder_name

# =============================================================================

def extractModelData_FM(model, model_pcr, landmask_pcr, clone_pcr, useRFS, use_2way):
    """
    Extracting data by using BMI-command "get_var", and preparing it depending 
    on model specification from initialized Delft3D FM model to be used later in the model run.
    There is a wider range of flexible mesh variables exposed which can be made use of
    if required. The number of variables exposed depends on DFM version in use.
    """
	
    cell_points_fm           = model.get_var('flowelemnode')	# list with number of 2D cell points
    
    # define separator between 2D and 1D parts of arrays
    if useRFS == True:
        separator = len(cell_points_fm)
    else:
        separator = 0.   
		 
    x_coords                 = model.get_var('xk') 				# x-coords of each cell corner point
    y_coords                 = model.get_var('yk') 				# y-coords of each cell corner point
    z_coords                 = model.get_var('zk') 				# elevation value of each cell corner point
    bottom_lvl               = model.get_var('bl') 				# surface elevation of each cell point
    bottom_lvl 				 = bottom_lvl[separator:]
    cellAreaSpherical        = model.get_var('ba') 				# cell area
    cellAreaSpherical 		 = cellAreaSpherical[separator:]
    xz_coords	             = model.get_var('xz')				# x-coords of each cell centre point
    xz_coords 				 = xz_coords[separator:]
    yz_coords          		 = model.get_var('yz')				# y-coords of each cell centre point
    yz_coords 				 = yz_coords[separator:]
    
    # preparing the DFM-coords as tuples containing x,y coords
    # for 1D part, x- and y-coords of centre points just have to be combined
    if useRFS == True:
        modelCoords = []
        for i in xrange(len(xz_coords)):
            x_i = xz_coords[i]
            y_i = yz_coords[i]
            xy_coords = (x_i, y_i)
            modelCoords.append([xy_coords])
        # TODO: perhaps use shorter function: alternative way to pair
        # modelCoords = zip(xz_coords, yz_coords)
    # for 2D part, x- and y-coords of corner points have to be allocated to cells
    elif useRFS == False:
		modelCoords = coupling_functions.getFMcoords(cell_points_fm, x_coords, y_coords)
    if use_2way == True:
        modelCoords_2way = coupling_functions.getFMcoords(cell_points_fm, x_coords, y_coords)
    elif use_2way == False:
        modelCoords_2way = []

    # retrieve PCR-data
    cellarea_data_pcr    = model_pcr.get_var('cellArea')
    landmask_data_pcr    = pcr.readmap(landmask_pcr)
    clone_data_pcr       = pcr.readmap(clone_pcr)
    
    return x_coords, y_coords, z_coords, bottom_lvl, cell_points_fm, separator, cellAreaSpherical, xz_coords, yz_coords, modelCoords, modelCoords_2way, \
                cellarea_data_pcr, landmask_data_pcr, clone_data_pcr
                
# =============================================================================
    
def extractModelData_FP(model, model_dir, model_pcr, landmask_pcr, clone_pcr, verbose_folder, use_RFS, use_2way, verbose):
    """
    Extracting data by using BMI-command "get_var", and preparing it depending 
    on model specification from initialized LISFLOOD-FP model to be used later in the model run.
    List of variables can be extended, but needs to be specified in 'lib_bmi.cpp' first with
    correct data type declaration.
    
    TO DO: get value specified as missing value in DEM automatically!
    """
	
    BLx                      = model.get_var('blx') 			# x-coord of bottom left corner of grid
    BLy                      = model.get_var('bly') 			# y-coord of bottom left corner of grid
    dx                       = model.get_var('dx') 				# incremental distance in x-direction
    dy                       = model.get_var('dy') 				# incremental distance in y-direction
    grid_dA                  = model.get_var('dA') 				# if lat/lon: array with cell area values; else: one uniform value
    DEM                      = model.get_var('DEM') 			# array with surface elevation values
    H                        = model.get_var('H') 				# array with water depth values
    SGCwidth                 = model.get_var('SGCwidth') 		# array with water depth values
    SGCQin                   = model.get_var('SGCQin') 			# array with discharge to be added to model (better than changing H)
    
    # getting dimensions of model grid
    rows, cols = DEM.shape 
    
    # converting arrays to flattened array to be processed later
    # required since Delft3D DFM stores data in flattened arrays too
    list_dA = grid_dA.ravel() 
    bottom_lvl = DEM.ravel()
    waterDepth = H.ravel()
    
    # initialize empty grids for model area to be filled later
    grid_x_coords = np.zeros([rows, cols])
    grid_y_coords = np.zeros([rows, cols]) 
    
    # computing centre points x- and y-coordinates of each grid cell
    for i in xrange(rows):
        grid_x_coords[i, :] = BLx + (np.arange(cols) + 0.5) * dx
    for i in xrange(cols):
        grid_y_coords[:, i] = BLy + (np.arange(rows) + 0.5) * dy
    # flipping array required since we start from bottom left
    grid_y_coords = np.flipud(grid_y_coords) 
    
    # compute list with centre point coords of each LFP-cell to be coupled to PCR-GLOBWB 
    # if RFS active, mask LFP-cells only to those with channel data and without missing values
    if use_RFS == True:
        i, j = np.where(np.logical_and(SGCwidth > 0., DEM != -9999))
        list_x_coords = grid_x_coords[i, j]
        list_y_coords = grid_y_coords[i, j]
        coupledFPindices = zip(i, j)
    elif use_RFS == False:
        i, j = np.where(DEM != -9999)
        list_x_coords = grid_x_coords[i, j]
        list_y_coords = grid_y_coords[i, j]
        coupledFPindices = zip(i, j)    
    if use_2way == True:
        i, j = np.where(DEM != -9999)
        list_x_coords_2way = grid_x_coords[i, j]
        list_y_coords_2way = grid_y_coords[i, j]
        coupledFPindices_2way = zip(i, j)
    elif use_2way == False:
		list_x_coords_2way = []
		list_y_coords_2way = []
		coupledFPindices_2way = []
		        
    # print and save verbose output
    if verbose == True:
        print 'x, y of corners of overall grid'
        for i in (0,-1):
            for j in (0,-1):
                print grid_x_coords[i, j], grid_y_coords[i, j]
        print 'rows', rows
        print 'cols', cols
        print 'number coupled FP-cells 1-D: ', len(coupledFPindices)
        print 'number coupled FP-cells 2-D: ', len(coupledFPindices_2way)
        fig = plt.figure()
        plt.imshow(np.ma.masked_outside(DEM, -1000,1000), cmap='terrain')
        plt.colorbar()
        plt.savefig(os.path.join(verbose_folder, 'DEM.png')) 
        plt.close(fig)
                
    # retrieve PCR-data
    cellarea_data_pcr    = model_pcr.get_var('cellArea')
    landmask_data_pcr    = pcr.readmap(landmask_pcr)
    clone_data_pcr       = pcr.readmap(clone_pcr)
    
    separator_1D = 0. # setting separator between 1-D and 2-D to 0 as only used for DFM
    
    return dx, dy, DEM, bottom_lvl, H, waterDepth, rows, cols, \
                list_x_coords, list_x_coords_2way, list_y_coords, list_y_coords_2way, coupledFPindices, coupledFPindices_2way, grid_dA, list_dA, SGCQin, separator_1D,\
                cellarea_data_pcr, landmask_data_pcr, clone_data_pcr 
    
# =============================================================================
    
def fillLFPgrid(model, indices_list, value_list, array_in, verbose_folder, verbose):
    """
    Reshaping lists used for volumes coupling into array as used by LISFLOOD-LFP. 
    
    Input: 
		- indices list: list with index-pointers to coupled LFP-grids in array
        - value_list: list with values per coupled LFP-grid specified in indices_list
        - array_in: map with extent of model extend (e.g. DEM array)
    Output: 
		- array with dimension of LFP-domain, filled with values for coupled LFP-cells
    """              
    
    #- creating files based on model extent as specified by DEM array
    filled_map = np.zeros([len(array_in), len(array_in[0])])
    test_map = np.copy(filled_map)

    #- loop through list of coupled LFP-cells and assign values
    for i in range(len(indices_list)):
        test_map[indices_list[i]] = 1. #- creating boolean map to quickly assess output
        filled_map[indices_list[i]] = value_list[i]  #- filling in actual values from list
        
    #- saving boolean map to verbose folder, if specified
    if (verbose == True) and (model.get_current_time() == 0.):
        print 'number of cells in filled LFP-grid that were filled', len((np.where(filled_map > 0.0))[1])
        fig = plt.figure()
        plt.imshow(test_map)
        plt.savefig(os.path.join(verbose_folder, 'filledFPgrid_testGrid.png'))
        plt.close(fig) 

    return filled_map
    
# =============================================================================

def activate2wayVariables(model_pcr, CoupledPCRcellIndices):
	
    # get variable/map from PCR
    new_preventRunoffToDischarge = model_pcr.get_var('preventRunoffToDischarge')
    new_controlDynamicFracWat = model_pcr.get_var('controlDynamicFracWat')
    new_waterBodyIdsAdjust = model_pcr.get_var(('WaterBodies', 'waterBodyIdsAdjust'))

    # adjust variable/map so it contains zeros at coupled cell locations
    new_preventRunoffToDischarge = set_values_in_array(new_preventRunoffToDischarge, CoupledPCRcellIndices, 0.)
    new_controlDynamicFracWat = set_values_in_array(new_controlDynamicFracWat, CoupledPCRcellIndices, 0.)
    new_waterBodyIdsAdjust = set_values_in_array(new_waterBodyIdsAdjust, CoupledPCRcellIndices, 0.)
    
    return new_preventRunoffToDischarge, new_controlDynamicFracWat, new_waterBodyIdsAdjust
    
# =============================================================================

def determine_InundationArea_Hydrodynamics(CouplePCR2model_2way, CoupledPCRcellIndices_2way, threshold_inundated_depth_floodplains, FMcellAreaSpherical, cellarea_data_pcr):
	"""
	get total area of all flooded FM cells coupled to current PCR cell
	"""
            
    # temporary variable for storing total area of flooded FM cells coupled to current PCR cell
    temp_inundated_area_FM = 0.
    
    # loop over all coupled FM cells
    for j in range(len(CouplePCR2model_2way[i][1])):
		# get current FM cell index
        current_FM_cell_index = CouplePCR2model_2way[i][1][j]
        # check if water depth of current FM cell is above chosen threshold
        if current_water_depth_fm[current_FM_cell_index] > threshold_inundated_depth_floodplains:
			# if so, add cell area to temporary variable
            temp_inundated_area_FM += FMcellAreaSpherical[current_FM_cell_index]
    # at end of loop, assign temporary variable to array storing inundated areas
    inundated_area_FM_2_PCR_coupled[i] = temp_inundated_area_FM
    # divide this value by PCR cell area to obtain flooded area fraction
    inundated_fraction_FM_2_PCR_coupled[i] = inundated_area_FM_2_PCR_coupled[i] / \
                                                        cellarea_data_pcr[CoupledPCRcellIndices_2way[i]]                                                        
    # assign value of current cell to zero array to be used with BMI
    inundated_fraction_FM_2_PCR[CoupledPCRcellIndices[i]] = inundated_fraction__FM_2_PCR_coupled[i]
	
	return inundated_area_FM_2_PCR_coupled, inundated_fraction_FM_2_PCR

# =============================================================================

def determine_InundationDepth_Hydrodynamics(model, landmask, missing_value_landmask, inundated_area_FM_2_PCR_coupled, CouplePCR2model_2way,CoupledPCRcellIndices_2way):
    
    # get water volume of all FM cells
    current_volume_fm = model.get_var('vol1')
    
    # create zero array for filling in total volumes and water depths for each coupled PCR cell
    water_volume_FM_2_PCR_coupled = np.zeros(len(CouplePCR2model_2way))
    water_depths_FM_2_PCR_coupled = np.zeros(len(CouplePCR2model_2way))
    # create zero array of appropriate size for using BMI function 'set_var'
    water_depths_FM_2_PCR = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
    
    # loop over all coupled PCR cells and fill in and/or calculate values
    for i in range(len(CouplePCR2model_2way)):
		# get total volume of all FM cells coupled to current PCR cell
        water_volume_FM_2_PCR_coupled[i] = np.sum(current_volume_fm[CouplePCR2model_2way[i][1]])
        # divide total volume by inundated FM area to obtain water depth
        if inundated_area_FM_2_PCR_coupled[i] > 0. :
            water_depths_FM_2_PCR_coupled[i] = water_volume_FM_2_PCR_coupled[i] / inundated_area_FM_2_PCR_coupled[i]
        else:
            water_depths_floodplains_FM_2_PCR_coupled[i] = 0.
        # assign value of current cell to zero array to be used with BMI
        water_depths_FM_2_PCR[CoupledPCRcellIndices_2way[i]] = water_depths_FM_2_PCR_coupled[i]

    return water_depths_FM_2_PCR 
      
# =============================================================================
    
def noStorage(model_pcr, missing_value_pcr, CoupledPCRcellIndices, CouplePCR2model):
    """
    This is done to prevent the lakes/reservoirs in PCR from draining all at once,
    introducing a huge flood wave at certain areas, e.g. reservoirs.
    By means of this function, waterbody and channel storage in PCR-GLOBWB are removed
    for all coupled PCR-cells.
    """
    
    # get required variables from PCR-GLOBWB
    current_channel_storage_pcr     = model_pcr.get_var('channelStorage')
    current_waterbody_storage_pcr   = model_pcr.get_var(('routing', 'waterBodyStorage'))
    
    # no channel storage    
    new_channel_storage_pcr = set_values_in_array(current_channel_storage_pcr, CoupledPCRcellIndices, 0.)
    # # no waterbody storage
    new_waterbody_storage_pcr = set_values_in_array(current_waterbody_storage_pcr, CoupledPCRcellIndices, 0.)    

    # activating coupling for relevant sections
    model_pcr.set_var(('grassland','ActivateCoupling'), 'True') #2way
    model_pcr.set_var(('forest','ActivateCoupling'), 'True') #2way
    model_pcr.set_var(('routing','ActivateCoupling'), 'True')
    model_pcr.set_var(('WaterBodies', 'ActivateCoupling'), 'True')
    
    # overwriting variables with new values
    model_pcr.set_var('channelStorage', new_channel_storage_pcr, missing_value_pcr)
    model_pcr.set_var(('routing','waterBodyStorage'), new_waterbody_storage_pcr)
    
    return
    
# =============================================================================
    
def updateHydrologicVariables(model_pcr, new_preventRunoffToDischarge, new_controlDynamicFracWat, new_waterBodyIdsAdjust, water_depths_floodplains_FM_2_PCR, \
            inundated_fraction_floodplains_FM_2_PCR, new_channelStorage_pcr, inundated_fraction_rivers_FM_2_PCR, new_controlFloodplainFactor, use_floodplain_infiltration_factor):
    """
    This functions update variables in the hydrologic models based on hydrodynamics
    """
    # adjusting maps controlling the locations where certain PCR variables should be updated
    # NOTE: THIS ONLY HAS TO BE DONE ONCE!
    if model_pcr.get_current_timestep() == 1.:
        model_pcr.set_var(('routing','preventRunoffToDischarge'), new_preventRunoffToDischarge)
        model_pcr.set_var(('routing','controlDynamicFracWat'), new_controlDynamicFracWat)
        model_pcr.set_var(('WaterBodies', 'waterBodyIdsAdjust'), new_waterBodyIdsAdjust)

    elif model_pcr.get_current_timestep() > 1.:
        # add water from FM floodplains back to PCR
        model_pcr.set_var(('grassland','floodplainWaterLayer'), water_depths_floodplains_FM_2_PCR)
        model_pcr.set_var(('forest','floodplainWaterLayer'), water_depths_floodplains_FM_2_PCR)
    
        # set the variable for dealing with floodplain inundated area fraction in PCR
        model_pcr.set_var(('grassland','inundatedFraction'), inundated_fraction_floodplains_FM_2_PCR)
        model_pcr.set_var(('forest','inundatedFraction'), inundated_fraction_floodplains_FM_2_PCR)
    
        # add water from FM rivers back to PCR
        model_pcr.set_var(('routing','channelStorage'), new_channelStorage_pcr)
    
        # set the variable for dealing with river (water bodies) area fraction in PCR
        model_pcr.set_var(('routing','waterBodyFractionFM'), inundated_fraction_rivers_FM_2_PCR)
    
        # floodplain scaling should only be changed using BMI if specified at model initialization
        if use_floodplain_infiltration_factor == True:
        
            # activate the scaling factor functions
            model_pcr.set_var(('forest','ActivateFactorInfiltrationFloodplain'), 'True')
            model_pcr.set_var(('grassland','ActivateFactorInfiltrationFloodplain'), 'True')
        
            # set the variable controlling the floodplain infiltration scaling factor
            model_pcr.set_var(('forest','controlFloodplainFactor'), new_controlFloodplainFactor)
            model_pcr.set_var(('grassland','controlFloodplainFactor'), new_controlFloodplainFactor)
    
    return
    
# =============================================================================
    
def noLDD(model_pcr, CoupledPCRcellIndices, verbose_folder, verbose):
    """
    The LDD of PCR needs to have 'pits' at all coupled cells to make sure that the routing is handled by DFM only.
    Adjusting the LDD to have a pit at every coupled cell (i.e. value 5), to make sure discharge values are not fed into DFM multiple times.
    
    Input:
		- List with  indices pointing to coupled PCR-cells
    """    
    # retrieving current LDD map
    LDD_PCR_new = np.copy(model_pcr.get_var(('routing', 'lddMap')))
    # replace LDD values within the hydrodynamic model domain to 5 (pit cell)
    LDD_PCR_new = set_values_in_array(LDD_PCR_new, CoupledPCRcellIndices, 5)

    # saving boolean map with locations of PCR-cells with pits to verbose folder, if specified
    if verbose == True:
        locCoupledPCRCells = checkLocationCoupledPCRcells(model_pcr, CoupledPCRcellIndices)
        fig = plt.figure()
        plt.imshow(locCoupledPCRCells)
        plt.savefig(os.path.join(verbose_folder, 'locCoupledPCRCells.png'))
        plt.close(fig)
        
    # overwriting current with new LDD information
    model_pcr.set_var(('routing', 'lddMap'), LDD_PCR_new, 255)
    return

# =============================================================================

def calculateDeltaVolumes(model_pcr, missing_value_pcr, secPerDay, CoupledPCRcellIndices, cellarea_data_pcr):
    """
    Calculating the delta volumes [m3/d] for all coupled PCR-cells.
    Delta volumes are based on discharge, surfaceRunoff, and topWaterLayer (only 2way-coupling) in PCR-GLOBWB.
    
    Input:
        - list with indexes pointing to coupled PCR-cells
        - PCR landmask data
        
    Output:
        - if RFS active, two arrays with delta volumes for river and floodplain cells, respectively
        - if RFS not active, one array with aggregated delta volumes
        - all outputs are in m3/day
    """
    #- update PCR for one day to get values
    model_pcr.update(1)

    #- retrieve data from PCR-GLOBWB
    current_discharge_pcr  = model_pcr.get_var('discharge')
    current_runoff_pcr     = model_pcr.get_var('landSurfaceRunoff')
    current_waterlayer_pcr = model_pcr.get_var('topWaterLayer')

    # # TODO: remove below when everything works
    # clone = model_pcr.get_var(('routing', 'lddMap'))
    # len_y = clone.shape[0]
    # len_x = clone.shape[1]
    # current_discharge_pcr = np.random.rand(len_y, len_x)
    # current_runoff_pcr = np.random.rand(len_y, len_x)
    # current_waterlayer_pcr  = np.random.rand(len_y, len_x)

    # 1a. Discharge
    
    # prepare empty array for all PCR-cells
    water_volume_PCR_rivers = np.zeros([len(current_discharge_pcr),len(current_discharge_pcr[0])])
    
    # loop over current discharge and convert to m3/d; missing values are replaced with zero
    water_volume_PCR_rivers = current_discharge_pcr * secPerDay
    water_volume_PCR_rivers[current_discharge_pcr==missing_value_pcr] = 0.

	# prepare empty array for coupled PCR-cells
    # get daily discharge volumes for all coupled PCR-cells [m3/day]
    water_volume_PCR_rivers_coupled = water_volume_PCR_rivers[zip(*CoupledPCRcellIndices)]
    # 1b. Runoff and Waterlayer
    water_volume_PCR_runoff = current_runoff_pcr * cellarea_data_pcr
    water_volume_PCR_runoff[current_runoff_pcr==missing_value_pcr] = 0.
    water_volume_PCR_waterlayer = current_waterlayer_pcr * cellarea_data_pcr
    water_volume_PCR_waterlayer[current_waterlayer_pcr==missing_value_pcr] = 0.
    # # loop over current runoff and waterlayer and convert to m3/d; missing values are replaced with zero

    water_volume_PCR_runoff_coupled = water_volume_PCR_runoff[zip(*CoupledPCRcellIndices)]
    water_volume_PCR_waterlayer_coupled = water_volume_PCR_waterlayer[zip(*CoupledPCRcellIndices)]
    water_volume_PCR_floodplains_coupled = water_volume_PCR_runoff_coupled + water_volume_PCR_waterlayer_coupled
    water_volume_PCR_coupled = water_volume_PCR_floodplains_coupled + water_volume_PCR_rivers_coupled

    delta_volume_PCR_coupled = water_volume_PCR_coupled - 0.

    return delta_volume_PCR_coupled
    
# =============================================================================

def calculateDeltaWater(CouplePCR2model, CoupleModel2PCR, delta_volume_PCR_coupled, CellAreaSpherical, fraction_timestep, model_type, useFluxes):
    """
    In this function the calculated daily delta volumes [m3/d] is translated to suitable units later to be used in the updating step.
    The input volumes of PCR-GLOBWB are here divided over the number of hydrodynamic cells within each PCR-cell.
    When fluxes are used, daily volumes are simply translated to averaged flux of m3/s per day.
    When states are used, daily volumes are first divided by the cell area of each coupled hydrodynamic cell, and then divided by the
    user-specified timestep, if provided.
    
    Input:
        - delta volumes
        - cell areas
        - model settings
    Output:
        - delta state [m/d] or [m/timestep]
        - delta flux [m3/s]
    """
    
    #- creating zero arrays to be populated with data    
    additional_water_level = np.zeros(len(CoupleModel2PCR))
    additional_water_volume = np.zeros(len(CoupleModel2PCR))
    verbose_volume = np.zeros(len(CoupleModel2PCR))

    # loop over all coupled PCR cells
    for i in range(len(CouplePCR2model)):
        
        # check if delta volume is positive
        if delta_volume_PCR_coupled[i] >= 0.:
				
            # divided volume by number of hydrodynamic cells in each PCR-cell
            temp_water_volume = delta_volume_PCR_coupled[i] / len(CouplePCR2model[i][1])
            
            # at each PCR-cell...
            for j in range(len(CouplePCR2model[i][1])):
                
                # ...get current hydrodynamic-cell index
                current_cell_index = CouplePCR2model[i][1][j]
                   
                # ...calculate additional water levels [m/day] for current index
                additional_water_level[current_cell_index] = temp_water_volume / CellAreaSpherical[current_cell_index]
                		
		# ...calculate additional water volume [m3/day] for current index
                additional_water_volume[current_cell_index] = temp_water_volume  / 1.
						
                verbose_volume[current_cell_index] = temp_water_volume  / 1.
                    
        else:
            os.sys.exit('delta volume PCR coupled is negative, should not be the case!')
                   
    # calculate additional water levels or fluxes based on chosen settings
    if (useFluxes == False) and (model_type == 'DFM'):
        delta_water = additional_water_level / fraction_timestep        # [m/timestep]
        
    elif (useFluxes == False) and (model_type == 'LFP'):
        delta_water = np.copy(additional_water_level)                   # [m/day]
   
    elif (useFluxes == True) and (model_type == 'DFM'):              
        delta_water = additional_water_level * 1000.                    # [mm/day]
        
    elif (useFluxes == True) and (model_type == 'LFP'):              
        delta_water = additional_water_volume / 86400.                  # [m3/s]
        
    return delta_water, verbose_volume
    
# =============================================================================
    
def updateModel(model, delta_water, update_step, separator, useFluxes, use_RFS, model_type, verbose):
    """
    Calculating the new water depth based on current depth and previously computed delta water depth (see calculateDeltaWater).
    Returning the new water depth to the hydrodynamic model and updating it according with user-specified time step.
    Alternatively, additional water input can be provided using fluxes.
    
    NOTE: For coupling with LFP, updating the model is done outside this function.
    NOTE: For coupling wiht LFP, it is generally recommended to use fluxes to avoid instabilites and increase computational speed. 
    
    Input:
		- delta_water: list or array with states or fluxes for all coupled DFM- or LFP-cells
		- update_step: user-specified update step (only for DFM applicable)
    Output:
		- check_array_in: can be used to assess whether delta_volumes and delta_water are correctly set into hydrodynamic model
    """
    
    if useFluxes == False:
		
        if model_type == 'DFM':
			
            # get current water level in DFM
            current_water_levels = np.copy(model.get_var('s1'))
            
            if use_RFS == True:
				current_water_levels_2d = current_water_levels[:separator]
				current_water_levels = current_water_levels[separator:]            
			      
            # create empty array for filling in values
            new_water_levels = np.zeros(len(current_water_levels))
    
            # loop over all DFM-cells
            for i in range(len(new_water_levels)):
				
                # if at "known" cells:
                if i < len(delta_water):
                    new_water_levels[i] = current_water_levels[i] + delta_water[i]
                # for cells outside the range of "known" cells (e.g. 1D channels):
                else:
                    new_water_levels[i] = current_water_levels[i]
            
            if use_RFS == True:
				new_water_levels = np.append(current_water_levels_2d , new_water_levels)
            
            # overwriting current with new water levels
            model.set_var('s1', new_water_levels)
            
            # compute output array
            check_array_in = model.get_var('s1')[separator:] - current_water_levels
            
            # update DFM with user-specified timestep
            model.update(update_step)
                    
        elif model_type == 'LFP':
			
			# get current water level in LFP
            current_water_levels = np.copy(model.get_var('H'))
            
            # create empty array for filling in values
            new_water_levels = np.zeros([len(current_water_levels), len(current_water_levels[0])])
            
            # compute new water levels
            new_water_levels = current_water_levels + delta_water
			
            # overwriting current with new water levels
            model.get_var('H')[:] = new_water_levels 
            
    elif useFluxes == True:

        if model_type == 'DFM': 
                    
            # creating zero-array for 2D part
            zeroValuesFor2D  = np.zeros(len(model.get_var('flowelemnode')))
            
            # resetting rain-array to zero before every update (otherwise it's added!)
            model.get_var('rain')[:] = np.zeros(len(model.get_var('rain')))
            
            # creating array combined of 2D and 1D part
            if use_RFS == True:
				delta_water = np.append(zeroValuesFor2D, delta_water)

            # adding new information
            if use_RFS == True:
                model.get_var('rain')[:] = delta_water
            elif use_RFS == False:
                model.get_var('rain')[:len(model.get_var('flowelemnode'))] = delta_water
            
            # updating model for one day
            model.update(86400.)
            
        elif model_type == 'LFP':
                        
            # overwriting current with new subgrid inflow discharge
            model.get_var('SGCQin')[:] = delta_water

    return
    
# =============================================================================  

def checkLocationCoupledPCRcells(model_pcr, CoupledPCRcellIndices):
    """
    This is mainly a debug-function to check whether PCR-cells are corretly coupled.
    
    Output:
		- boolean map that allows to assess the locations of the coupled PCR-GLOBWB cells.    
    """              
    
    # retrieving array representing PCR-GLOBWB data extent
    shapemap = np.copy(model_pcr.get_var(('routing', 'lddMap')))
    
    # creating empty array with PCR-GLOBWB extent 
    check_map = np.zeros([len(shapemap), len(shapemap[0])])

    # assigning True value for all coupled PCR-GLOBWB cells
    for i in range(len(CoupledPCRcellIndices)):
        check_map[CoupledPCRcellIndices[i]] = 1.    

    return check_map
