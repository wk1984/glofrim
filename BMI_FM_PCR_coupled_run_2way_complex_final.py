
# specify which PCR-GLOBWB model version to use
# 0  = 2.0.1, 30min
# 1  = 2.0.1, 05min
# 10 = 2.0.3, 30min, accuTravel
# 11 = 2.0.3, 05min, accuTravel
# 20 = 2.0.3, 30min, kinematicWave
# 21 = 2.0.3, 05min, kinematicWave
use_pcr_version = 10

# specify whether to use upstream GRDC discharge data instead of PCR-GLOBWB output
# 0 = forcing by PCR-GLOBWB only (no GRDC data)
# 1 = GRDC discharge data at specified upstream points
use_GRDC_discharge_data = 0

print ''
if use_pcr_version < 10:
    print 'PCR-GLOBWB version:  2.0.1'
else:
    print 'PCR-GLOBWB version:  2.0.3'
if (use_pcr_version == 0) or (use_pcr_version == 10) or (use_pcr_version == 20):
    print 'model resolution:    30min'
else:
    print 'model resolution:    05min'
if use_pcr_version < 20:
    print 'routing method:      accuTravel'
else:
    print 'routing method:      kinematicWave'
print ''
if use_GRDC_discharge_data == 0:
    print 'upstream discharge:  PCR-GLOBWB'
else:
    print 'upstream discharge:  GRDC data'
print ''

# -------------------------------------------------------------------------------------------------
# Required libraries
# -------------------------------------------------------------------------------------------------

import os
import datetime
import bmi.wrapper
import pcraster as pcr
if use_pcr_version == 0:
    import pcrglobwb_2way_simple
    from pcrglobwb_2way_simple import pcrglobwb_bmi
elif use_pcr_version == 1:
    import pcrglobwb_2way_05min
    from pcrglobwb_2way_05min import pcrglobwb_bmi
elif use_pcr_version == 10:
    import pcrglobwb_203_30min_2way
    from pcrglobwb_203_30min_2way import pcrglobwb_bmi
elif use_pcr_version == 11:
    import pcrglobwb_203_05min_2way
    from pcrglobwb_203_05min_2way import pcrglobwb_bmi
from coupling_PCR_FM import coupling_functions

# -------------------------------------------------------------------------------------------------
# Model settings
# -------------------------------------------------------------------------------------------------

# specify model timesteps
end_time_pcr                          = 365    # [days] (1995 = 365, 1995-2001 = 2557)
update_step_fm                        = 1800   # [seconds]

# specify whether to adjust initial groundwater stores based on previously calculated recharge
# 0 = no, 1 = yes
adjust_initial_groundwater            = 1

# specify whether to use the river/floodplain cells scheme
# 0 = no, 1 = yes
use_river_floodplain_scheme           = 1

# specify whether to use scaling factor for floodplain infiltration once 2-way coupling begins
# 0 = no, 1 = yes
use_floodplain_infiltration_factor    = 0

# set threshold above which a cell is considered 'inundated'
threshold_inundated_depth             = 0.0    # [m] (for when river/floodplain scheme is not used)
threshold_inundated_depth_rivers      = 0.0    # [m] (for river cells only)
threshold_inundated_depth_floodplains = 0.0    # [m] (for floodplain cells only)

# specify whether (and how much) burning of river cells should be done
# (set burn depth to zero for no burning, set s1 negative to create initial water levels)
burn_depth_bl                         = 0.     # [m] (variable 'bl' = bottom level)
burn_depth_s1                         = -2.    # [m] (variable 's1' = water level)
limit_bl_to_zero                      = 1      # prevent negative values yes (1) or no (0)
limit_s1_to_zero                      = 1      # prevent negative values yes (1) or no (0)

# other
if (use_pcr_version == 0) or (use_pcr_version == 10) or (use_pcr_version == 20):
    missing_value_landmask            = 255    # 30min resolution landmask
else:
    missing_value_landmask            = 255    # 05min resolution landmask
missing_value_pcr                     = -999
secPerDay                             = 86400. # seconds per day

# specify whether to save arrays with daily values of all FM water levels and depths to lists
# 0 = no, 1 = yes
save_fm_lists = 0

# -------------------------------------------------------------------------------------------------
# Settings for using GRDC discharge data as input
# -------------------------------------------------------------------------------------------------

# coordinates of inflow points
if (use_pcr_version == 0) or (use_pcr_version == 10):
    inflow_point_Niger                = [(12.25, -8.25)]
    inflow_point_Bani                 = [(12.75, -6.75)]
elif (use_pcr_version == 1) or (use_pcr_version == 11):
    inflow_point_Niger                = [(12.04166698, -8.29166698)]
    inflow_point_Bani                 = [(12.54166698, -6.70833349), (12.54166698, -6.54166651)]

# coordinates of upstream cells of inflow points
if (use_pcr_version == 0) or (use_pcr_version == 10):
    upstream_points_Niger             = [(11.75,-8.25), (11.75,-8.75)]
    upstream_points_Bani              = [(12.25,-6.75), (12.25,-7.25)]
elif (use_pcr_version == 1) or (use_pcr_version == 11):
    upstream_points_Niger             = [(11.95833302, -8.375), (11.95833302, -8.29166698), (11.95833302, -8.20833302)]
    upstream_points_Bani              = [(12.54166698, -6.79166651), (12.54166698, -6.45833349), (12.45833302, -6.45833349), \
                                         (12.45833302, -6.54166651)]

# path GRDC data
path_GRDC_data = r"C:\Users\haag\Files\Subversion\PCR-GLOBWB-FM\data\GRDC" # work laptop
#r"H:\MSc thesis\Data_and_files\Discharge_Niger\GRDC" # external disk
#r"C:\Users\Arjen.Arjen-PC\Documents\MSc thesis\Data_and_files\Discharge_Niger\GRDC" # home laptop
#r"C:\Users\haag\Files\Subversion\PCR-GLOBWB-FM\data\GRDC" # node

# file names
file_Niger = 'Koulikoro_1995_2001_adjusted.txt'    #'1134100.day' = original data (with possible missing values)
file_Bani  = 'Douna_1995_2001_adjusted.txt'        #'1134300.day' = original data (with possible missing values)

# missing values option (-999 in raw data)
# 0 = do nothing, 1 = set missing values to zero
mv_option = 1

# make a plot of GRDC discharge input
make_plot_GRDC_discharge = 1

# data reading options
# number of rows to skip in header, whether to use names of columns, delimiter, 
# column names for time, discharge and modified discharge and used time string
header_skip              = 40
column_names             = True
file_delimiter           = ';'
GRDC_data_name_time      = 'YYYYMMDD'
GRDC_data_name           = 'Original'    # other option: 'Calculated'
file_time_str            = '%Y-%m-%d'

# -------------------------------------------------------------------------------------------------
# Region (e.g. Inner Delta) to save total inundated area and volume
# -------------------------------------------------------------------------------------------------

# area of interest
area_pcr_lon_min = -5.0
area_pcr_lon_max = -3.0
area_pcr_lat_min = 13.5
area_pcr_lat_max = 17.0

# threshold to include a PCR cell in area of interest [fraction]
# (if fraction of total area of FM cells exceeds this value, the PCR cell is included)
area_pcr_threshold_fraction_FM_in_PCR_cell = 0.3

# threshold below which an area is not stored in the threshold list [m]
area_pcr_threshold_save_to_list = 0.1

# lists to store data in
area_pcr_total_inundated_area_floodplains = []
area_pcr_total_inundated_area_threshold   = []
area_pcr_total_flooded_volume_floodplains = []

# path, including file names, to save data to files
area_pcr_path_inundated_area      = r"C:\Users\haag\Files\Subversion\PCR-GLOBWB-FM\data\temp_output\inundated_area_daily.txt"
area_pcr_path_inundated_threshold = r"C:\Users\haag\Files\Subversion\PCR-GLOBWB-FM\data\temp_output\inundated_threshold_daily.txt"
area_pcr_path_flooded_volume      = r"C:\Users\haag\Files\Subversion\PCR-GLOBWB-FM\data\temp_output\flooded_volume_daily.txt"

# -------------------------------------------------------------------------------------------------
# Plot and print options
# -------------------------------------------------------------------------------------------------

# set default figure size
pylab.rcParams['figure.figsize']      = (20.0, 10.0)

# set plot stuff for coloured plots of FM water depths
my_cmap = matplotlib.cm.get_cmap('Blues_r')
my_cmap.set_under('seagreen')
my_cmap.set_bad('seagreen')

# print timesteps to screen (to check speed of calculations)
# 0 = no, 1 = yes
print_timesteps                       = 0

# control whether to print remaining delta volume info during model updates
# 0 = no, 1 = yes
print_delta_volume_info               = 0

# -------------------------------------------------------------------------------------------------
# Set paths (based on system architecture)
# -------------------------------------------------------------------------------------------------

config_fm    = os.path.abspath(r'..\data\fm\model_test_epsg4326.mdu')

if (use_pcr_version == 0) or (use_pcr_version == 10) or (use_pcr_version == 20):
    landmask_pcr    = os.path.abspath(r'..\data\pcrglobwb\input30min\Niger_landmask_30min.map')
    clone_pcr       = os.path.abspath(r'..\data\pcrglobwb\input30min\Niger_clone_30min.map')
    if adjust_initial_groundwater == 1:
        GW_recharge_pcr = os.path.abspath(r'..\data\pcrglobwb\input30min\GWR_daily_mean_1995_2001.txt')
    if use_pcr_version == 0:
        config_pcr  = os.path.abspath(r'..\data\pcrglobwb\niger_initialstates_2way.ini')
    elif use_pcr_version == 10:
        config_pcr  = os.path.abspath(r'..\data\pcrglobwb\niger_initialstates_2way_203_30min_accutravel.ini')
    elif use_pcr_version == 20:
        config_pcr  = os.path.abspath(r'..\data\pcrglobwb\niger_initialstates_2way_203_30min_kinematic.ini')
else:
    landmask_pcr    = os.path.abspath(r'..\data\pcrglobwb\input05min\Niger_landmask_05min_2.map') # clone_niger_05min_2_filled_2
    clone_pcr       = os.path.abspath(r'..\data\pcrglobwb\input05min\domain_niger_05min_2.map')
    if adjust_initial_groundwater == 1:
        GW_recharge_pcr = os.path.abspath(r'..\data\pcrglobwb\input05min\GWR_daily_mean_1979_1980.txt')
    if use_pcr_version == 1:
        config_pcr  = os.path.abspath(r'..\data\pcrglobwb\niger_initialstates_2way_05min.ini')
    elif use_pcr_version == 11:
        config_pcr  = os.path.abspath(r'..\data\pcrglobwb\niger_initialstates_2way_203_05min_accutravel.ini')
    elif use_pcr_version == 21:
        config_pcr  = os.path.abspath(r'..\data\pcrglobwb\niger_initialstates_2way_203_05min_kinematic.ini')

# System architecture (32bit or 64 bit), important for dynamic-link library (dll) of FM 
is64 = sys.maxsize > 2**32
if is64:
    dll_path = os.path.abspath(r'..\bin\fm-x64\unstruc.dll')
else:
    dll_path = os.path.abspath(r'..\bin\fm-w32\unstruc.dll')

# -------------------------------------------------------------------------------------------------
# Calculate relevant values from chosen settings
# -------------------------------------------------------------------------------------------------

end_time_fm            = end_time_pcr * secPerDay    # FM uses seconds for counting its current model time

# calculate fraction of update value per day (to divide water levels per update)
fraction_timestep_fm   = secPerDay / update_step_fm

# create zero arrays for filling in the inundated area and volume within the specified region (above)
saved_inundated_area   = np.zeros(end_time_pcr)
saved_inundated_volume = np.zeros(end_time_pcr)

# -------------------------------------------------------------------------------------------------
# Initialize models
# -------------------------------------------------------------------------------------------------

# initialize D-FLow FM
model_fm = bmi.wrapper.BMIWrapper(engine=dll_path, configfile=config_fm)
model_fm.initialize()

# initialize PCR-GLOBWB
model_pcr = pcrglobwb_bmi.pcrglobwbBMI()
model_pcr.initialize(config_pcr)

# spinning-up PCR
model_pcr.spinup()

# Extracting relevant data from models:

# D-Flow FM
x_fm           = model_fm.get_var('xk')
y_fm           = model_fm.get_var('yk')
z_fm           = model_fm.get_var('zk')
bottom_lvl_fm  = model_fm.get_var('bl')
cell_points_fm = model_fm.get_var('flowelemnode')
# NOTE: the coupling from data (earlier notebooks) used 'netelemnode' for 'cell_points_fm', just like Gena's script for plotting data. 
# While this 'netelemnode' is also available as a variable now, it gives different values and results in errors in the code below. 
# It turns out that 'flowelemnode' actually gives 100% the same values as the previously used 'netelemnode', 
# which is why it is used to get the 'cell_points_fm' in the code above.

# PCR-GLOBWB
cellarea_data_pcr    = model_pcr.get_var('cellArea')
landmask_data_pcr    = pcr.readmap(landmask_pcr)
clone_data_pcr       = pcr.readmap(clone_pcr)

# Coupling the grids

# getting coordinates of all PCR cells
PCRcoords = coupling_functions.getPCRcoords(landmask_data_pcr)

# getting coordinates of all FM cells
FMcoords = coupling_functions.getFMcoords(cell_points_fm, x_fm, y_fm)

# coupling the grids, getting:
# - coordinates of all coupled FM cells
# - coordinates of all coupled PCR cells
# - list with indices of all FM cells and their coupled PCR cells (if any)
# - list with indices of all unique coupled PCR cells including array of all their coupled FM cells
# - list with indices of all unique coupled PCR cells only
CoupledCellsInfoAll = coupling_functions.coupleAllCells(FMcoords,PCRcoords)

# extract relevant lists (easier to read code of coupled run further down this notebook)
CoupleFM2PCR = CoupledCellsInfoAll[2]    # list of all FM cells with coupled PCR cells (if any)
CouplePCR2FM = CoupledCellsInfoAll[3]    # list of all coupled PCR cells and their coupled FM cells

# converting single indices of coupled PCR cells to double (array,column) indices
CoupledPCRcellIndices = coupling_functions.convertPCRindices(CoupledCellsInfoAll[4], landmask_pcr, missingValues=missing_value_landmask)
#CoupledPCRcellIndices = coupling_functions.convertPCRindicesOLD(CoupledCellsInfoAll[4], landmask_pcr)    # old version    

# calculating cell area for all FM cells
FMcellAreaSpherical = model_fm.get_var('ba')
#FMcellAreaSpherical = coupling_functions.calculateAllCellAreaSpherical(FMcoords)    # old version, before I could get 'ba' with BMI

# plotting the grids to quickly check results

# complete grids
coupling_functions.plotGridfromCoords(PCRcoords, FMcoords)

# coupled grids
coupling_functions.plotGridfromCoords(CoupledCellsInfoAll[1],CoupledCellsInfoAll[0])

if use_river_floodplain_scheme == 1:
    # creating seperate lists for FM river cells and floodplain cells, getting:
    # - list with FM cell indices representing river cells
    # - list with FM cell indices representing floodplain cells
    # - list with coordinates of FM cells representing river cells
    # - list with coordinates of FM cells representing floodplain cells
    # - list with a list of FM river cells per coupled PCR cell
    # - boolean array indicating whether coupled PCR cells have FM river cells (True) or not (False)
    # - list with a list of FM floodplain cells per coupled PCR cell
    FM_rivers_floodplains = coupling_functions.getFMcellsRiverFloodplain(CouplePCR2FM, FMcoords, FMcellAreaSpherical, \
                                                                         area_threshold='all', \
                                                                         threshold_factor_3vertices=1.5, \
                                                                         threshold_factor_4vertices=3.0)
    
    FM_river_cells                           = FM_rivers_floodplains[0]
    FM_floodplain_cells                      = FM_rivers_floodplains[1]
    FM_river_coords                          = FM_rivers_floodplains[2]
    FM_floodplain_coords                     = FM_rivers_floodplains[3]
    FM_river_cells_per_coupled_PCR_cell      = FM_rivers_floodplains[4]
    boolean_river_cell_in_coupled_PCR_cell   = FM_rivers_floodplains[5]
    FM_floodplain_cells_per_coupled_PCR_cell = FM_rivers_floodplains[6]
    
    # printing and plotting output to quickly check results
    print 'number of river cells found:', len(FM_river_cells)
    coupling_functions.plotGridfromCoords(FM_river_coords)
    
else:
    print 'original cell coupling chosen, no river/floodplain stuff loaded'

# burning of 'bl' (and 's1' to prevent initial water levels) using BMI

# check if burning can be done (requires identified river cells)
if use_river_floodplain_scheme == 1:
    # burn 'bl' if specified depth > 0
    if (burn_depth_bl > 0.):
        print "BMI burning 'bl' of", burn_depth_bl, 'meter(s)'
        burned_bottom_lvl_fm = np.copy(bottom_lvl_fm)
        burned_bottom_lvl_fm[np.array(FM_river_cells)] -= burn_depth_bl
        # set negative values to zero if specified
        if limit_bl_to_zero == 1:
            burned_bottom_lvl_fm[np.where(burned_bottom_lvl_fm<0.)] = 0.
        # pass to FM using BMI's 'set_var'
        model_fm.set_var('bl', burned_bottom_lvl_fm)
        # and update variable used in this notebook
        bottom_lvl_fm = model_fm.get_var('bl')
    else:
        print "No burning of 'bl'"
    # burn 's1' if specified depth > 0, or add initial water levels if depth < 0 ("reversed burning")
    if (burn_depth_s1 != 0.):
        print "BMI burning 's1' of", burn_depth_s1, 'meter(s)'
        updated_s1_fm = np.copy(model_fm.get_var('s1'))
        updated_s1_fm[np.array(FM_river_cells)] -= burn_depth_s1
        # set negative values to zero if specified
        if limit_s1_to_zero == 1:
            updated_s1_fm[np.where(updated_s1_fm<0.)] = 0.
        # pass to FM using BMI's 'set_var'
        model_fm.set_var('s1', updated_s1_fm)
    else:
        print "No burning of 's1'"

# get all latitudes (xcoords) and longitudes (ycoords) from clone map
pcr_clone_xcoords = pcr.pcr2numpy(pcr.xcoordinate(clone_data_pcr),-999)
pcr_clone_ycoords = pcr.pcr2numpy(pcr.ycoordinate(clone_data_pcr),-999)

pcr_clone_lons = pcr_clone_xcoords[0]
pcr_clone_lats = pcr_clone_ycoords[:,0]

# calculate (approximate) dimension of data
area_pcr_dimension = abs(pcr_clone_lats[0]-pcr_clone_lats[1])

# obtain indices of study area
area_pcr_index_lon_1 = abs(pcr_clone_lons-(area_pcr_lon_min+(area_pcr_dimension/2))).argmin()
area_pcr_index_lon_2 = abs(pcr_clone_lons-(area_pcr_lon_max-(area_pcr_dimension/2))).argmin()
area_pcr_index_lat_1 = abs(pcr_clone_lats-(area_pcr_lat_min+(area_pcr_dimension/2))).argmin()
area_pcr_index_lat_2 = abs(pcr_clone_lats-(area_pcr_lat_max-(area_pcr_dimension/2))).argmin()

# make sure index 1 comes before 2 (depends on how data is stored)
if area_pcr_index_lon_1 > area_pcr_index_lon_2:
    area_pcr_index_lon_1, area_pcr_index_lon_2 = area_pcr_index_lon_2, area_pcr_index_lon_1
if area_pcr_index_lat_1 > area_pcr_index_lat_2:
    area_pcr_index_lat_1, area_pcr_index_lat_2 = area_pcr_index_lat_2, area_pcr_index_lat_1

# convert to arrays
area_pcr_lon_array = np.arange(area_pcr_index_lon_1, area_pcr_index_lon_2 + 1)
area_pcr_lat_array = np.arange(area_pcr_index_lat_1, area_pcr_index_lat_2 + 1)

# construct array with all PCR (lat,lon) indices of chosen area
area_pcr_lat_lon_array = []

for i in range(len(area_pcr_lat_array)):
    for j in range(len(area_pcr_lon_array)):
        area_pcr_lat_lon_array.append((area_pcr_lat_array[i], area_pcr_lon_array[j]))
        
# find all indices that are coupled to FM cells
area_pcr_lat_lon_array_coupled = []
for i in range(len(area_pcr_lat_lon_array)):
    if area_pcr_lat_lon_array[i] in CoupledPCRcellIndices:
        area_pcr_lat_lon_array_coupled.append(area_pcr_lat_lon_array[i])
        
# find indices of coupled lat/lon array within 'CoupledPCRcellIndices'
# which will correspond to indices within 'CouplePCR2FM' (so that I can calculate area of FM cells per coupled PCR cell)
area_pcr_indices_lat_lon_array_coupled_in_CoupledPCRcellIndices = []
area_pcr_indices_to_get_coords_1 = []
for i in range(len(area_pcr_lat_lon_array_coupled)):
    for j in range(len(CoupledPCRcellIndices)):
        if area_pcr_lat_lon_array_coupled[i] == CoupledPCRcellIndices[j]:
            area_pcr_indices_lat_lon_array_coupled_in_CoupledPCRcellIndices.append(j)
            area_pcr_indices_to_get_coords_1.append(CouplePCR2FM[j][0])

# get from the coupled lat/lon array only the indices which satify the area fraction
area_pcr_indices_in_CoupledPCRcellIndices_final = []
area_pcr_lat_lon_array_coupled_AreaFraction = []
area_pcr_indices_to_get_coords_2 = []
for i in range(len(area_pcr_lat_lon_array_coupled)):
    area_pcr_temp_area_FM  = np.sum(FMcellAreaSpherical[CouplePCR2FM[area_pcr_indices_lat_lon_array_coupled_in_CoupledPCRcellIndices[i]][1]])
    area_pcr_temp_area_PCR = cellarea_data_pcr[area_pcr_lat_lon_array_coupled[i]]
    area_pcr_temp_fraction_area_FM_in_PCR_cell = area_pcr_temp_area_FM / area_pcr_temp_area_PCR
    if area_pcr_temp_fraction_area_FM_in_PCR_cell >= area_pcr_threshold_fraction_FM_in_PCR_cell:
        area_pcr_indices_in_CoupledPCRcellIndices_final.append(area_pcr_indices_lat_lon_array_coupled_in_CoupledPCRcellIndices[i])
        area_pcr_lat_lon_array_coupled_AreaFraction.append(area_pcr_lat_lon_array_coupled[i])
        area_pcr_indices_to_get_coords_2.append(area_pcr_indices_to_get_coords_1[i])
        
# get corresponding FM floodplain cell indices (to store inundated area and volume)
area_pcr_indices_FM_floodplain_cells = []
for i in range(len(CouplePCR2FM)):
    if CouplePCR2FM[i][0] in area_pcr_indices_to_get_coords_2:
        for j in range(len(FM_floodplain_cells_per_coupled_PCR_cell[i])):
            area_pcr_indices_FM_floodplain_cells.append(FM_floodplain_cells_per_coupled_PCR_cell[i][j])
            
# remove possible 'None' values from FM list
list_indices_to_remove_because_of_None = []
for i in range(len(area_pcr_indices_FM_floodplain_cells)):
    if type(area_pcr_indices_FM_floodplain_cells[i]) != numpy.int64:
        list_indices_to_remove_because_of_None.append(i)
list_indices_to_remove_because_of_None.reverse()
for i in range(len(list_indices_to_remove_because_of_None)):
    del area_pcr_indices_FM_floodplain_cells[list_indices_to_remove_because_of_None[i]]
        
# get PCR coords of calculated areas (could be useful for plotting/checking)
PCRcoords_InnerDelta_coupled_area = []
for i in range(len(area_pcr_indices_to_get_coords_2)):
    PCRcoords_InnerDelta_coupled_area.append(PCRcoords[area_pcr_indices_to_get_coords_2[i]])
    
# get FM coords of calculated areas (could be useful for plotting/checking)
FMcoords_InnerDelta_coupled_area = []
for i in range(len(area_pcr_indices_FM_floodplain_cells)):
    FMcoords_InnerDelta_coupled_area.append(FMcoords[area_pcr_indices_FM_floodplain_cells[i]])

# plot grids of specified area

# calculate 'buffer' to properly scale plot
extra_buffer_for_plot = max(0., (area_pcr_lat_max - area_pcr_lat_min) / (area_pcr_lon_max - area_pcr_lon_min))

# full FM grid with PCR cells that represent the specified area
coupling_functions.plotGridfromCoords(PCRcoords_InnerDelta_coupled_area, FMcoords, linewidthGrid1=3, linewidthGrid2=1)

# scale axis
plt.axis([area_pcr_lon_min - extra_buffer_for_plot - 0.5, area_pcr_lon_max + extra_buffer_for_plot + 0.5, \
          area_pcr_lat_min - 0.1, area_pcr_lat_max + 0.1])
plt.show()

# show number of PCR cells
print 'number of PCR cells:', len(area_pcr_indices_to_get_coords_2)

# FM (floodplain) cells and PCR cells of the specified area only
coupling_functions.plotGridfromCoords(PCRcoords_InnerDelta_coupled_area, FMcoords_InnerDelta_coupled_area,\
                                      linewidthGrid1=2, linewidthGrid2=1)

# scale axis
plt.axis([area_pcr_lon_min - extra_buffer_for_plot - 0.5, area_pcr_lon_max + extra_buffer_for_plot + 0.5, \
          area_pcr_lat_min - 0.1, area_pcr_lat_max + 0.1])
plt.show()

# show number of FM cells
print 'number of FM cells:', len(area_pcr_indices_FM_floodplain_cells)

# save water levels and depths to lists for time 0 (after initialization, before any updates)
if save_fm_lists == 1:
    output_water_levels_fm = []
    output_water_depths_fm = []
    output_water_levels_fm.append(model_fm.get_var('s1'))
    output_water_depths_fm.append(model_fm.get_var('s1') - bottom_lvl_fm)

# create a plot to see water depths of FM after initialization (before any updates)
f, ax = plt.subplots(1, 1)
coupling_functions.plotGridValues(ax, FMcoords, (model_fm.get_var('s1') - bottom_lvl_fm), model_name='FM', limits = None, \
                                  cell_cmap=my_cmap, extend_colorbar='both')

# and a plot of channel storage and water body storage values of PCR after spin-up

# channelStorage
current_storage_pcr = model_pcr.get_var('channelStorage')
f, ax = plt.subplots(1, 1)
coupling_functions.plotGridValues(ax, PCRcoords, current_storage_pcr, model_name='PCR', cell_cmap=my_cmap, \
                                  extend_colorbar='both')

# waterBodyStorage
# it does not have missing values, instead it is just zero at all locations without storage
# but if I want to plot it using my coupling functions, it needs missing values (as the code assumes it has), so:
if use_pcr_version < 10:
    waterBodyStorage_with_mv = model_pcr.get_var(('WaterBodies', 'waterBodyStorage'))
else:
    waterBodyStorage_with_mv = model_pcr.get_var(('routing', 'waterBodyStorage'))
for i in range(len(current_storage_pcr)):
    for j in range(len(current_storage_pcr[0])):
        if current_storage_pcr[i][j] == -999:
            waterBodyStorage_with_mv[i][j] = -999

# plot it
f, ax = plt.subplots(1, 1)
coupling_functions.plotGridValues(ax, PCRcoords, waterBodyStorage_with_mv, model_name='PCR', cell_cmap=my_cmap, \
                                  extend_colorbar='both')

# setting channelStorage of coupled cells to zero before coupled run
current_storage_pcr   = model_pcr.get_var('channelStorage')

new_storage_pcr = np.zeros([len(current_storage_pcr),len(current_storage_pcr[0])])
for i in range(len(new_storage_pcr)):
    for j in range(len(new_storage_pcr[0])):
        new_storage_pcr[i][j] = current_storage_pcr[i][j]

for i in range(len(CoupledPCRcellIndices)):
    new_storage_pcr[CoupledPCRcellIndices[i]] = 0.

model_pcr.set_var('channelStorage', new_storage_pcr, missing_value_pcr)

# setting waterBodyStorage of coupled cells to zero before coupled run
if use_pcr_version < 10:
    current_waterbody_storage_pcr = model_pcr.get_var(('WaterBodies', 'waterBodyStorage'))
else:
    current_waterbody_storage_pcr = model_pcr.get_var(('routing', 'waterBodyStorage'))

new_waterbody_storage_pcr = np.zeros([len(current_waterbody_storage_pcr),len(current_waterbody_storage_pcr[0])])
for i in range(len(new_waterbody_storage_pcr)):
    for j in range(len(new_waterbody_storage_pcr[0])):
        new_waterbody_storage_pcr[i][j] = current_waterbody_storage_pcr[i][j]

for i in range(len(CoupledPCRcellIndices)):
    new_waterbody_storage_pcr[CoupledPCRcellIndices[i]] = 0.

if use_pcr_version < 10:
    model_pcr.set_var(('WaterBodies','waterBodyStorage'), new_waterbody_storage_pcr)
else:
    model_pcr.set_var(('routing','waterBodyStorage'), new_waterbody_storage_pcr)

# plot to inspect resulting channelStorage and waterBodyStorage

# channelStorage
new_storage_pcr = model_pcr.get_var('channelStorage')
f, ax = plt.subplots(1, 1)
coupling_functions.plotGridValues(ax, PCRcoords, new_storage_pcr, model_name='PCR', cell_cmap=my_cmap, \
                                  extend_colorbar='both')

# waterBodyStorage
# it does not have missing values, instead it is just zero at all locations without storage
# but if I want to plot it using my coupling functions, it needs missing values (as the code assumes it has), so:
if use_pcr_version < 10:
    new_waterBodyStorage_with_mv = model_pcr.get_var(('WaterBodies', 'waterBodyStorage'))
else:
    new_waterBodyStorage_with_mv = model_pcr.get_var(('routing', 'waterBodyStorage'))
for i in range(len(new_storage_pcr)):
    for j in range(len(new_storage_pcr[0])):
        if new_storage_pcr[i][j] == -999:
            new_waterBodyStorage_with_mv[i][j] = -999

# plot it
f, ax = plt.subplots(1, 1)
coupling_functions.plotGridValues(ax, PCRcoords, new_waterBodyStorage_with_mv, model_name='PCR', cell_cmap=my_cmap, \
                                  extend_colorbar='both')

###--- added to new function "adjust_iniGR" ---###
#
#if adjust_initial_groundwater == 1:
#    # load groundwater recharge data from text file
#    GW_average_recharge = np.loadtxt(GW_recharge_pcr)
#    # get recession coefficient using BMI
#    GW_coefficient      = model_pcr.get_var('recessionCoeff')
#    # calculate initial groundwater stores
#    GW_initial_stores   = GW_average_recharge / GW_coefficient
#    # get current groundwater storage using BMI
#    GW_current_storage  = model_pcr.get_var('storGroundwater')
#    # create an array that will be used to set the new groundwater storage
#    GW_new_storage      = np.copy(GW_current_storage)
#    # adjust this array, using the maximum value of the current and new initial stores
#    for i in range(len(CoupledPCRcellIndices)):
#        if GW_current_storage[CoupledPCRcellIndices[i]] != -999:
#            if GW_initial_stores[CoupledPCRcellIndices[i]] > GW_current_storage[CoupledPCRcellIndices[i]]:
#                GW_new_storage[CoupledPCRcellIndices[i]] = GW_initial_stores[CoupledPCRcellIndices[i]]
#    # set this into the model using BMI
#    model_pcr.set_var('storGroundwater', GW_new_storage)
#    
###--- added to new function "adjust_iniGR" ---###

# Adjusting the LDD to have a pit at every coupled cell, to make sure discharge values are not fed into FM multiple times

# make a copy of the current LDD
LDD_PCR_new = np.copy(model_pcr.get_var(('routing', 'lddMap')))

# set all values at coupled PCR cells to 5 (pit/sink)
for i in range(len(CoupledPCRcellIndices)):
    LDD_PCR_new[CoupledPCRcellIndices[i]] = 5.

# set the adjusted LDD into the model using BMI's 'set_var'
model_pcr.set_var(('routing', 'lddMap'), LDD_PCR_new, 255)

if use_GRDC_discharge_data == 1:
    
    # Finding the indices corresponding to the chosen inflow and upstream points
    # -------------------------------------------------------------------------------------
    # get coordinates of PCR centroids
    x_coord_centroids_PCR = pcr.xcoordinate(landmask_data_pcr)
    x_coord_centroids_np  = pcr.pcr2numpy(x_coord_centroids_PCR, -999)
    
    y_coord_centroids_PCR = pcr.ycoordinate(landmask_data_pcr)
    y_coord_centroids_np  = pcr.pcr2numpy(y_coord_centroids_PCR, -999)
    
    # find where these match the specified points
    indices_inflow_points  = []
    indices_upstream_points = []
    
    # Niger inflow point:
    for i in range(len(inflow_point_Niger)):
        temp_indices_array = np.where((x_coord_centroids_np == inflow_point_Niger[i][1]) & \
                                      (y_coord_centroids_np == inflow_point_Niger[i][0]))
        indices_inflow_points.append((temp_indices_array[0][0], temp_indices_array[1][0]))
    
    # Bani inflow point:
    for i in range(len(inflow_point_Bani)):
        temp_indices_array = np.where((x_coord_centroids_np == inflow_point_Bani[i][1]) & \
                                      (y_coord_centroids_np == inflow_point_Bani[i][0]))
        indices_inflow_points.append((temp_indices_array[0][0], temp_indices_array[1][0]))
        
    # adjust inflow points for 05min version (Bani has 2, the first one will be removed)
    if len(indices_inflow_points) > 2:
        del indices_inflow_points[1]
    
    # Niger upstream points:
    for i in range(len(upstream_points_Niger)):
        temp_indices_array = np.where((x_coord_centroids_np == upstream_points_Niger[i][1]) & \
                                      (y_coord_centroids_np == upstream_points_Niger[i][0]))
        indices_upstream_points.append((temp_indices_array[0][0], temp_indices_array[1][0]))
        
    # Bani upstream points:
    for i in range(len(upstream_points_Bani)):
        temp_indices_array = np.where((x_coord_centroids_np == upstream_points_Bani[i][1]) & \
                                      (y_coord_centroids_np == upstream_points_Bani[i][0]))
        indices_upstream_points.append((temp_indices_array[0][0], temp_indices_array[1][0]))
    # -------------------------------------------------------------------------------------
    
    # Loading data
    # -------------------------------------------------------------------------------------
    # get total paths
    path_GRDC_data_Niger      = os.path.join(path_GRDC_data, file_Niger)
    path_GRDC_data_Bani       = os.path.join(path_GRDC_data, file_Bani)
    
    # loading data for Niger
    if 'adjusted' in file_Niger:
        GRDC_data_Niger_discharge = np.genfromtxt(path_GRDC_data_Niger)
    else:
        GRDC_data_Niger_total     = np.genfromtxt(path_GRDC_data_Niger, skip_header=header_skip, names=column_names, \
                                                  delimiter=file_delimiter, dtype=None)
        GRDC_data_Niger_time      = GRDC_data_Niger_total[GRDC_data_name_time]
        GRDC_data_Niger_discharge = GRDC_data_Niger_total[GRDC_data_name]
    
    # loading data for Bani
    if 'adjusted' in file_Bani:
        GRDC_data_Bani_discharge  = np.genfromtxt(path_GRDC_data_Bani)
    else:
        GRDC_data_Bani_total      = np.genfromtxt(path_GRDC_data_Bani, skip_header=header_skip, names=column_names, \
                                                  delimiter=file_delimiter, dtype=None)
        GRDC_data_Bani_time       = GRDC_data_Bani_total[GRDC_data_name_time]
        GRDC_data_Bani_discharge  = GRDC_data_Bani_total[GRDC_data_name]
    # -------------------------------------------------------------------------------------
    
    # Creating arrays of data for correct time period
    # -------------------------------------------------------------------------------------
    # find start time
    start_date_as_datetime = model_pcr.get_start_time()
    start_date_as_string   = start_date_as_datetime.strftime("%Y-%m-%d")
    
    # find indices of start date
    if 'adjusted' in file_Niger:
        start_index_Niger = 0
    else:
        start_index_Niger = np.where(GRDC_data_Niger_time==start_date_as_string)[0]
    if 'adjusted' in file_Bani:
        start_index_Bani = 0
    else:
        start_index_Bani  = np.where(GRDC_data_Bani_time==start_date_as_string)[0]
    
    # create separate arrays
    GRDC_discharge_Niger = GRDC_data_Niger_discharge[start_index_Niger:(start_index_Niger+end_time_pcr)]
    GRDC_discharge_Bani  = GRDC_data_Bani_discharge[start_index_Bani:(start_index_Bani+end_time_pcr)]
    
    # handling of missing values
    if mv_option == 1:
        GRDC_discharge_Niger = GRDC_discharge_Niger.clip(min=0)
        GRDC_discharge_Bani  = GRDC_discharge_Bani.clip(min=0)
    else:
        pass
    
    # create single list with arrays
    GRDC_discharge_input = zip(GRDC_discharge_Niger,GRDC_discharge_Bani)
    # -------------------------------------------------------------------------------------
    
    # Adjusting the LDD to have a pit at every identified upstream cell, to allow input of GRDC data without also using PCR-GLOBWB data
    # -------------------------------------------------------------------------------------
    # make a copy of the current LDD
    LDD_PCR_new_2 = np.copy(model_pcr.get_var('lddMap'))
    
    # set all values at identified upstream locations to 5 (pit/sink)
    for i in range(len(indices_upstream_points)):
        LDD_PCR_new_2[indices_upstream_points[i]] = 5.
    
    # set the adjusted LDD into the model using BMI's 'set_var'
    model_pcr.set_var('lddMap', LDD_PCR_new_2, 255)
    # -------------------------------------------------------------------------------------
    
    # check if plots are to be made
    if make_plot_GRDC_discharge == 1:
        print 'Niger plot:'
        plt.plot(GRDC_discharge_Niger)
        plt.show()
        print 'Bani plot:'
        plt.plot(GRDC_discharge_Bani)
        plt.show()
    
    print 'WARNING: using GRDC discharge data at specified inflow points!'
    
else:
    print "'normal' model run, not using GRDC data"

# 1. Create map that turns PCR water bodies off for coupled cells

# get variable/map from PCR
new_waterBodyIdsAdjust = model_pcr.get_var(('WaterBodies', 'waterBodyIdsAdjust'))

# adjust variable/map so it contains zeros at coupled cell locations
for i in range(len(CouplePCR2FM)):
    new_waterBodyIdsAdjust[CoupledPCRcellIndices[i]] = 0.0

# 2. Create map to prevent runoff from entering channels at coupled cells

# OLD METHOD (before adjustment, commented out since new method is better)
# ------------------------------------------------------------------------
#new_preventRunoffToDischarge = np.ones([len(current_storage_pcr),len(current_storage_pcr[0])])
#
# assign missing values to correct locations
#for i in range(len(new_preventRunoffToDischarge)):
#    for j in range(len(new_preventRunoffToDischarge[0])):
#        if current_storage_pcr[i][j] == missing_value_pcr:
#            new_preventRunoffToDischarge[i][j] = missing_value_pcr
# ------------------------------------------------------------------------

# get variable/map from PCR
new_preventRunoffToDischarge = model_pcr.get_var('preventRunoffToDischarge')

# adjust variable/map so it contains zeros at coupled cell locations
for i in range(len(CouplePCR2FM)):
    new_preventRunoffToDischarge[CoupledPCRcellIndices[i]] = 0.0

# 3. Create map that keeps PCR's 'dynamicFracWat' in non-coupled cells and allows overwriting in coupled cells

# get variable/map from PCR
new_controlDynamicFracWat = model_pcr.get_var('controlDynamicFracWat')

# adjust variable/map so it contains zeros at coupled cell locations
# OLD METHOD (commented out, new method is better)
# ------------------------------------------------------------------------
#for i in range(len(new_controlDynamicFracWat)):
#    for j in range(len(new_controlDynamicFracWat[0])):
#        if (i,j) in CoupledPCRcellIndices:
#            new_controlDynamicFracWat[i][j] = 0.
# ------------------------------------------------------------------------
for i in range(len(CouplePCR2FM)):
    new_controlDynamicFracWat[CoupledPCRcellIndices[i]] = 0.0

# 1. calculating fraction and area of (inundated) water per coupled PCR cell

# temporary variable to store total inundated area to-be-stored to list
temp_inundated_area_threshold_list = 0.

# get water depth of all FM cells
current_water_depth_fm = model_fm.get_var('s1') - bottom_lvl_fm

# if the original method (not the river/floodplain scheme) is to be used:
if use_river_floodplain_scheme == 0:
    
    # create zero array for filling in total area and fraction of flooded FM cells for each coupled PCR cell
    inundated_area_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
    inundated_fraction_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
    
    # create zero array of appropriate size for using BMI function 'set_var'
    inundated_fraction_FM_2_PCR = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
    
    # loop over all coupled PCR cells and fill in and/or calculate values
    for i in range(len(CouplePCR2FM)):
        
        # get total area of all flooded FM cells coupled to current PCR cell
        
        # temporary variable for storing total area of flooded FM cells coupled to current PCR cell
        temp_inundated_area_FM = 0.
        
        # loop over all coupled FM cells
        for j in range(len(CouplePCR2FM[i][1])):
            
            # get current FM cell index
            current_FM_cell_index = CouplePCR2FM[i][1][j]
            
            # check if water depth of current FM cell is above chosen threshold
            if current_water_depth_fm[current_FM_cell_index] > threshold_inundated_depth:
                
                # if so, add cell area to temporary variable
                temp_inundated_area_FM += FMcellAreaSpherical[current_FM_cell_index]
                
                # check if current FM cell falls within specified region and if water depth is above threshold
                if (current_FM_cell_index in area_pcr_indices_FM_floodplain_cells) and \
                   (current_water_depth_fm[current_FM_cell_index] >= area_pcr_threshold_save_to_list):
                
                    temp_inundated_area_threshold_list += FMcellAreaSpherical[current_FM_cell_index]
                
        # at end of loop, assign temporary variable to array storing inundated areas
        inundated_area_FM_2_PCR_coupled[i] = temp_inundated_area_FM
        
        # limit inundated area to area of coupled PCR cell
        #inundated_area_FM_2_PCR_coupled[i] = min([inundated_area_FM_2_PCR_coupled[i], cellarea_data_pcr[CoupledPCRcellIndices[i]]])
        
        # divide this value by PCR cell area to obtain flooded area fraction
        inundated_fraction_FM_2_PCR_coupled[i] = inundated_area_FM_2_PCR_coupled[i] / cellarea_data_pcr[CoupledPCRcellIndices[i]]
        
        # limit inundated fraction between 0 and 1
        #inundated_fraction_FM_2_PCR_coupled[i] = max([0., min([inundated_fraction_FM_2_PCR_coupled[i], 1.])])
        
        # assign value of current cell to zero array to be used with BMI
        inundated_fraction_FM_2_PCR[CoupledPCRcellIndices[i]] = inundated_fraction_FM_2_PCR_coupled[i]

# or else, if the river/floodplain scheme is to be used:
elif use_river_floodplain_scheme == 1:
    
    # create zero array for filling in total area and fraction of flooded FM cells for each coupled PCR cell
    inundated_area_rivers_FM_2_PCR_coupled          = np.zeros(len(CouplePCR2FM))
    inundated_fraction_rivers_FM_2_PCR_coupled      = np.zeros(len(CouplePCR2FM))
    
    inundated_area_floodplains_FM_2_PCR_coupled     = np.zeros(len(CouplePCR2FM))
    inundated_fraction_floodplains_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
    
    # create zero array of appropriate size for using BMI function 'set_var'
    inundated_fraction_rivers_FM_2_PCR      = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
    inundated_fraction_floodplains_FM_2_PCR = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
    
    # loop over all coupled PCR cells and fill in and/or calculate values
    for i in range(len(CouplePCR2FM)):
        
        # check if there are any FM river cells for this PCR cell, in which case a distinction between river/floodplain cells is made
        if boolean_river_cell_in_coupled_PCR_cell[i]:
            
            # 1a. calculating fraction and area of wet FM river cells per coupled PCR cell
            
            # temporary variable for storing total area of wet FM river cells coupled to current PCR cell
            temp_inundated_area_FM = 0.
            
            # loop over all river cells coupled to current PCR cell
            for j in range(len(FM_river_cells_per_coupled_PCR_cell[i])):
            
                # get current FM cell index
                current_FM_cell_index = FM_river_cells_per_coupled_PCR_cell[i][j]
                
                # check if water depth of current FM cell is above chosen threshold
                if current_water_depth_fm[current_FM_cell_index] > threshold_inundated_depth_rivers:
                    
                    # if so, add cell area to temporary variable
                    temp_inundated_area_FM += FMcellAreaSpherical[current_FM_cell_index]
            
            # at end of loop, assign temporary variable to array storing inundated areas
            inundated_area_rivers_FM_2_PCR_coupled[i] = temp_inundated_area_FM
            
            # limit inundated area to area of coupled PCR cell
            #inundated_area_rivers_FM_2_PCR_coupled[i] = min([inundated_area_rivers_FM_2_PCR_coupled[i], \
            #                                                      cellarea_data_pcr[CoupledPCRcellIndices[i]]])
            
            # divide this value by PCR cell area to obtain flooded area fraction
            inundated_fraction_rivers_FM_2_PCR_coupled[i] = inundated_area_rivers_FM_2_PCR_coupled[i] / \
                                                                    cellarea_data_pcr[CoupledPCRcellIndices[i]]
            
            # limit inundated fraction between 0 and 1
            #inundated_fraction_rivers_FM_2_PCR_coupled[i] = max([0., min([inundated_fraction_rivers_FM_2_PCR_coupled[i], 1.])])
            
            # assign value of current cell to zero array to be used with BMI
            inundated_fraction_rivers_FM_2_PCR[CoupledPCRcellIndices[i]] = inundated_fraction_rivers_FM_2_PCR_coupled[i]
            
            # 1b. calculating fraction and area of inundated FM floodplain cells per coupled PCR cell
            
            # temporary variable for storing total area of inundated FM floodplain cells coupled to current PCR cell
            temp_inundated_area_FM = 0.
            
            # loop over all floodplain cells coupled to current PCR cell
            for j in range(len(FM_floodplain_cells_per_coupled_PCR_cell[i])):
            
                # get current FM cell index
                current_FM_cell_index = FM_floodplain_cells_per_coupled_PCR_cell[i][j]
                
                # check if water depth of current FM cell is above chosen threshold
                if current_water_depth_fm[current_FM_cell_index] > threshold_inundated_depth_floodplains:
                    
                    # if so, add cell area to temporary variable
                    temp_inundated_area_FM += FMcellAreaSpherical[current_FM_cell_index]
                    
                    # check if current FM cell falls within specified region and if water depth is above threshold
                    if (current_FM_cell_index in area_pcr_indices_FM_floodplain_cells) and \
                       (current_water_depth_fm[current_FM_cell_index] >= area_pcr_threshold_save_to_list):
                    
                        temp_inundated_area_threshold_list += FMcellAreaSpherical[current_FM_cell_index]
            
            # at end of loop, assign temporary variable to array storing inundated areas
            inundated_area_floodplains_FM_2_PCR_coupled[i] = temp_inundated_area_FM
            
            # limit inundated area to area of coupled PCR cell
            #inundated_area_floodplains_FM_2_PCR_coupled[i] = min([inundated_area_floodplains_FM_2_PCR_coupled[i], \
            #                                                      cellarea_data_pcr[CoupledPCRcellIndices[i]]])
            
            # divide this value by PCR cell area to obtain flooded area fraction
            inundated_fraction_floodplains_FM_2_PCR_coupled[i] = inundated_area_floodplains_FM_2_PCR_coupled[i] / \
                                                                    cellarea_data_pcr[CoupledPCRcellIndices[i]]
            
            # limit inundated fraction between 0 and 1
            #inundated_fraction_floodplains_FM_2_PCR_coupled[i] = max([0., min([inundated_fraction_floodplains_FM_2_PCR_coupled[i], 1.])])
            
            # assign value of current cell to zero array to be used with BMI
            inundated_fraction_floodplains_FM_2_PCR[CoupledPCRcellIndices[i]] = inundated_fraction_floodplains_FM_2_PCR_coupled[i]
            
        # or if there are no FM river cells, use the same same procedure as in the original method
        else:
                
            # get total area of all flooded FM cells coupled to current PCR cell
            
            # temporary variable for storing total area of flooded FM cells coupled to current PCR cell
            temp_inundated_area_FM = 0.
            
            # loop over all coupled FM cells
            for j in range(len(CouplePCR2FM[i][1])):
                
                # get current FM cell index
                current_FM_cell_index = CouplePCR2FM[i][1][j]
                
                # check if water depth of current FM cell is above chosen threshold
                if current_water_depth_fm[current_FM_cell_index] > threshold_inundated_depth_floodplains:
                    
                    # if so, add cell area to temporary variable
                    temp_inundated_area_FM += FMcellAreaSpherical[current_FM_cell_index]
                    
                    # check if current FM cell falls within specified region and if water depth is above threshold
                    if (current_FM_cell_index in area_pcr_indices_FM_floodplain_cells) and \
                       (current_water_depth_fm[current_FM_cell_index] >= area_pcr_threshold_save_to_list):
                    
                        temp_inundated_area_threshold_list += FMcellAreaSpherical[current_FM_cell_index]
            
            # at end of loop, assign temporary variable to array storing inundated areas
            inundated_area_floodplains_FM_2_PCR_coupled[i] = temp_inundated_area_FM
            
            # limit inundated area to area of coupled PCR cell
            #inundated_area_floodplains_FM_2_PCR_coupled[i] = min([inundated_area_floodplains_FM_2_PCR_coupled[i], \
            #                                                      cellarea_data_pcr[CoupledPCRcellIndices[i]]])
            
            # divide this value by PCR cell area to obtain flooded area fraction
            inundated_fraction_floodplains_FM_2_PCR_coupled[i] = inundated_area_floodplains_FM_2_PCR_coupled[i] / \
                                                        cellarea_data_pcr[CoupledPCRcellIndices[i]]
            
            # limit inundated fraction between 0 and 1
            #inundated_fraction_floodplains_FM_2_PCR_coupled[i] = max([0., min([inundated_fraction_floodplains_FM_2_PCR_coupled[i], 1.])])
            
            # assign value of current cell to zero array to be used with BMI
            inundated_fraction_floodplains_FM_2_PCR[CoupledPCRcellIndices[i]] = inundated_fraction_floodplains_FM_2_PCR_coupled[i]

# 2. extracting FM water volume for each coupled PCR cell and converting these to water depths to add back to PCR

# get water volume of all FM cells
current_volume_fm = model_fm.get_var('vol1')

# if the original method (not the river/floodplain scheme) is to be used:
if use_river_floodplain_scheme == 0:
    
    # create zero array for filling in total volumes and water depths for each coupled PCR cell
    water_volume_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
    water_depths_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
    
    # create zero array of appropriate size for using BMI function 'set_var'
    water_depths_FM_2_PCR = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
    
    # loop over all coupled PCR cells and fill in and/or calculate values
    for i in range(len(CouplePCR2FM)):
        
        # get total volume of all FM cells coupled to current PCR cell
        water_volume_FM_2_PCR_coupled[i] = np.sum(current_volume_fm[CouplePCR2FM[i][1]])
        
        # OLD: divide total volume by PCR cell area to obtain water depth
        #water_depths_FM_2_PCR_coupled[i] = water_volume_FM_2_PCR_coupled[i] / cellarea_data_pcr[CoupledPCRcellIndices[i]]
        
        # NEW: divide total volume by inundated FM area to obtain water depth
        if inundated_area_FM_2_PCR_coupled[i] > 0. :
            water_depths_FM_2_PCR_coupled[i] = water_volume_FM_2_PCR_coupled[i] / inundated_area_FM_2_PCR_coupled[i]
        else:
            water_depths_FM_2_PCR_coupled[i] = 0.
        
        # assign value of current cell to zero array to be used with BMI
        water_depths_FM_2_PCR[CoupledPCRcellIndices[i]] = water_depths_FM_2_PCR_coupled[i]

# or else, if the river/floodplain scheme is to be used:
elif use_river_floodplain_scheme == 1:
    
    # create zero array for filling in total volumes and water depths for each coupled PCR cell
    water_volume_rivers_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
    
    water_volume_floodplains_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
    water_depths_floodplains_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
    
    # create zero array of appropriate size for using BMI function 'set_var'
    water_volume_rivers_FM_2_PCR      = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
    water_depths_floodplains_FM_2_PCR = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
    
    # loop over all coupled PCR cells and fill in and/or calculate values
    for i in range(len(CouplePCR2FM)):
        
        # check if there are any FM river cells for this PCR cell, in which case a distinction between river/floodplain cells is made
        if boolean_river_cell_in_coupled_PCR_cell[i]:
            
            # 2a. extracting water volume from FM river cells
            
            # get total volume of all FM cells coupled to current PCR cell
            water_volume_rivers_FM_2_PCR_coupled[i] = np.sum(current_volume_fm[FM_river_cells_per_coupled_PCR_cell[i]])
            
            # assign value of current cell to zero array to be used with BMI
            water_volume_rivers_FM_2_PCR[CoupledPCRcellIndices[i]] = water_volume_rivers_FM_2_PCR_coupled[i]
            
            # 2b. extracting water volume from FM floodplain cells
            
            # get total volume of all FM cells coupled to current PCR cell
            water_volume_floodplains_FM_2_PCR_coupled[i] = np.sum(current_volume_fm[FM_floodplain_cells_per_coupled_PCR_cell[i]])
            
            # divide total volume by inundated FM area to obtain water depth
            if inundated_area_floodplains_FM_2_PCR_coupled[i] > 0. :
                water_depths_floodplains_FM_2_PCR_coupled[i] = water_volume_floodplains_FM_2_PCR_coupled[i] / \
                                                                inundated_area_floodplains_FM_2_PCR_coupled[i]
            else:
                water_depths_floodplains_FM_2_PCR_coupled[i] = 0.
            
            # assign value of current cell to zero array to be used with BMI
            water_depths_floodplains_FM_2_PCR[CoupledPCRcellIndices[i]] = water_depths_floodplains_FM_2_PCR_coupled[i]            
            
        # or if there are no FM river cells, use the same same procedure as in the original method
        else:
            
            # get total volume of all FM cells coupled to current PCR cell
            water_volume_floodplains_FM_2_PCR_coupled[i] = np.sum(current_volume_fm[CouplePCR2FM[i][1]])
            
            # OLD: divide total volume by PCR cell area to obtain water depth
            #water_depths_floodplains_FM_2_PCR_coupled[i] = water_volume_floodplains_FM_2_PCR_coupled[i] / \
            #                                                cellarea_data_pcr[CoupledPCRcellIndices[i]]
            
            # NEW: divide total volume by inundated FM area to obtain water depth
            if inundated_area_floodplains_FM_2_PCR_coupled[i] > 0. :
                water_depths_floodplains_FM_2_PCR_coupled[i] = water_volume_floodplains_FM_2_PCR_coupled[i] / \
                                                                inundated_area_floodplains_FM_2_PCR_coupled[i]
            else:
                water_depths_floodplains_FM_2_PCR_coupled[i] = 0.
            
            # assign value of current cell to zero array to be used with BMI
            water_depths_floodplains_FM_2_PCR[CoupledPCRcellIndices[i]] = water_depths_floodplains_FM_2_PCR_coupled[i]

# save the total inundated area and volume within the specified area
area_pcr_total_inundated_area_floodplains.append(np.sum(inundated_area_floodplains_FM_2_PCR_coupled[area_pcr_indices_in_CoupledPCRcellIndices_final]))
area_pcr_total_inundated_area_threshold.append(temp_inundated_area_threshold_list)
area_pcr_total_flooded_volume_floodplains.append(np.sum(water_depths_floodplains_FM_2_PCR_coupled[area_pcr_indices_in_CoupledPCRcellIndices_final]))

# check if floodplain scaling factor is to be activated
if use_floodplain_infiltration_factor == 1:
    
    # obtain default scaling factor
    new_controlFloodplainFactor = np.copy(model_pcr.get_var(('forest','controlFloodplainFactor')))
    
    # assign zero value to coupled PCR cells with FM river cells
    for i in range(len(CoupledPCRcellIndices)):
        if boolean_river_cell_in_coupled_PCR_cell[i]:
            new_controlFloodplainFactor[CoupledPCRcellIndices[i]] = 0.
            
    print 'Using floodplain scaling factor!'
else:
    print 'NOT using floodplain scaling factor!'

current_channelStorage_pcr  = model_pcr.get_var(('routing','channelStorage'))

new_channelStorage_pcr = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)

for i in range(len(new_channelStorage_pcr)):
    for j in range(len(new_channelStorage_pcr[0])):
        if current_channelStorage_pcr[i][j] != -999:
            new_channelStorage_pcr[i][j] = current_channelStorage_pcr[i][j] + water_volume_rivers_FM_2_PCR[i][j]

# use BMI functions for coupling FM stuff back to PCR

# activate 2-way coupling functions
# NOTE: THIS ONLY HAS TO BE DONE ONCE!
model_pcr.set_var(('grassland','ActivateCoupling'), 'True') # added to function "noStorage"
model_pcr.set_var(('forest','ActivateCoupling'), 'True') # added to function "noStorage"
model_pcr.set_var(('routing','ActivateCoupling'), 'True')
model_pcr.set_var(('WaterBodies', 'ActivateCoupling'), 'True')

###--- added to new function "updateHydrologicVariables" ---###

## adjusting maps controlling the locations where certain PCR variables should be updated
## NOTE: THIS ONLY HAS TO BE DONE ONCE!
#model_pcr.set_var(('routing','preventRunoffToDischarge'), new_preventRunoffToDischarge)
#model_pcr.set_var(('routing','controlDynamicFracWat'), new_controlDynamicFracWat)
#model_pcr.set_var(('WaterBodies', 'waterBodyIdsAdjust'), new_waterBodyIdsAdjust)

## add water from FM floodplains back to PCR
#model_pcr.set_var(('grassland','floodplainWaterLayer'), water_depths_floodplains_FM_2_PCR)
#model_pcr.set_var(('forest','floodplainWaterLayer'), water_depths_floodplains_FM_2_PCR)

## set the variable for dealing with floodplain inundated area fraction in PCR
#model_pcr.set_var(('grassland','inundatedFraction'), inundated_fraction_floodplains_FM_2_PCR)
#model_pcr.set_var(('forest','inundatedFraction'), inundated_fraction_floodplains_FM_2_PCR)

## add water from FM rivers back to PCR
#model_pcr.set_var(('routing','channelStorage'), new_channelStorage_pcr)

## set the variable for dealing with river (water bodies) area fraction in PCR
#model_pcr.set_var(('routing','waterBodyFractionFM'), inundated_fraction_rivers_FM_2_PCR)

## floodplain scaling should only be changed using BMI if specified at model initialization
#if use_floodplain_infiltration_factor == 1:
    
    ## activate the scaling factor functions
    #model_pcr.set_var(('forest','ActivateFactorInfiltrationFloodplain'), 'True')
    #model_pcr.set_var(('grassland','ActivateFactorInfiltrationFloodplain'), 'True')
    
    ## set the variable controlling the floodplain infiltration scaling factor
    #model_pcr.set_var(('forest','controlFloodplainFactor'), new_controlFloodplainFactor)
    #model_pcr.set_var(('grassland','controlFloodplainFactor'), new_controlFloodplainFactor)
    
###--- added to new function "updateHydrologicVariables" ---###

# COUPLED RUN (day 1)
# first day outside loop, to make sure PCR time is at start time (timestep 1) and not at end of spin-up (timestep 365)

# Update PCR and obtain variables
# ----------------------------------------------------------------------------------------------------
# update PCR
model_pcr.update(1)

# get relevant variables at current timestep
current_discharge_pcr  = model_pcr.get_var('discharge')
current_runoff_pcr     = model_pcr.get_var('landSurfaceRunoff')
current_waterlayer_pcr = model_pcr.get_var('topWaterLayer')

# check settings of using GRDC discharge data
if use_GRDC_discharge_data == 1:
    for i in range(len(indices_inflow_points)):
        current_discharge_pcr[indices_inflow_points[i]] += GRDC_discharge_input[0][i]
# ----------------------------------------------------------------------------------------------------

# Calculate delta volumes for all coupled PCR cells (output PCR - input FM)
# ----------------------------------------------------------------------------------------------------

# 1. Rivers

# convert PCR discharge [m3/s] to water volume for river cells [m3/day]
water_volume_PCR_rivers = np.zeros([len(current_discharge_pcr),len(current_discharge_pcr[0])])
for i in range(len(current_discharge_pcr)):
    for j in range(len(current_discharge_pcr[0])):
        if current_discharge_pcr[i][j] != missing_value_pcr:
            water_volume_PCR_rivers[i][j] = current_discharge_pcr[i][j] * secPerDay
        else:
            water_volume_PCR_rivers[i][j] = 0.
        
# get water volume for rivers for all coupled cells [m3/day]
water_volume_PCR_rivers_coupled = np.zeros(len(CoupledPCRcellIndices))
for i in range(len(CoupledPCRcellIndices)):
    water_volume_PCR_rivers_coupled[i] = water_volume_PCR_rivers[CoupledPCRcellIndices[i]]

# 2. Floodplains

# convert PCR runoff and top water layer [m/day] to water volumes [m3/day]
water_volume_PCR_runoff      = np.zeros([len(current_runoff_pcr),len(current_runoff_pcr[0])])
water_volume_PCR_waterlayer  = np.zeros([len(current_waterlayer_pcr),len(current_waterlayer_pcr[0])])
for i in range(len(current_runoff_pcr)):
    for j in range(len(current_runoff_pcr[0])):
        if current_runoff_pcr[i][j] != missing_value_pcr:
            water_volume_PCR_runoff[i][j]     = current_runoff_pcr[i][j] * cellarea_data_pcr[i][j]
            water_volume_PCR_waterlayer[i][j] = current_waterlayer_pcr[i][j] * cellarea_data_pcr[i][j]
        else:
            water_volume_PCR_runoff[i][j]     = 0.
            water_volume_PCR_waterlayer[i][j] = 0.
            
# get water volumes for all coupled cells [m3/day]
water_volume_PCR_runoff_coupled      = np.zeros(len(CoupledPCRcellIndices))
water_volume_PCR_waterlayer_coupled  = np.zeros(len(CoupledPCRcellIndices))
water_volume_PCR_floodplains_coupled = np.zeros(len(CoupledPCRcellIndices))
for i in range(len(CoupledPCRcellIndices)):
    water_volume_PCR_runoff_coupled[i]      = water_volume_PCR_runoff[CoupledPCRcellIndices[i]]
    water_volume_PCR_waterlayer_coupled[i]  = water_volume_PCR_waterlayer[CoupledPCRcellIndices[i]]
    water_volume_PCR_floodplains_coupled[i] = water_volume_PCR_runoff_coupled[i] + water_volume_PCR_waterlayer_coupled[i]

# calculate delta volume for all coupled cells [m3/day]
delta_volume_PCR_rivers_coupled      = water_volume_PCR_rivers_coupled - water_volume_rivers_FM_2_PCR_coupled
delta_volume_PCR_floodplains_coupled = water_volume_PCR_floodplains_coupled - water_volume_floodplains_FM_2_PCR_coupled

# 3. set PCR's 'topWaterLayer' to zero

# get variable/map from PCR
new_topWaterLayer_forest = model_pcr.get_var(('forest', 'topWaterLayer'))
new_topWaterLayer_grass  = model_pcr.get_var(('grassland', 'topWaterLayer'))

# adjust variable/map so it contains zeros at coupled cell locations
for i in range(len(CouplePCR2FM)):
    new_topWaterLayer_forest[CoupledPCRcellIndices[i]] = 0.0
    new_topWaterLayer_grass[CoupledPCRcellIndices[i]]  = 0.0

# set in PCR model for each land cover type using BMI function 'set_var'
model_pcr.set_var(('forest','topWaterLayer'), new_topWaterLayer_forest)
model_pcr.set_var(('grassland','topWaterLayer'), new_topWaterLayer_grass)
# ----------------------------------------------------------------------------------------------------

# Calculate variables to update FM water levels based on delta volumes
# NOTE: CURRENTLY USING ALGORITHM 2 TO HANDLE NEGATIVE DELTA VOLUMES
# ----------------------------------------------------------------------------------------------------
# get current FM water levels
new_water_levels_fm_1 = np.copy(model_fm.get_var('s1'))

# initiate additional water levels array for each FM cell (same indices as FM cell index to easily add up later) [m/day]
additional_water_level_FM = np.zeros(len(CoupleFM2PCR))

# loop over all coupled PCR cells
for i in range(len(CouplePCR2FM)):
    
    # if the original method (not the river/floodplain scheme) is to be used:
    if use_river_floodplain_scheme == 0:
        
        # check if delta volume is positive, and if so, perform the usual procedure that is also used in 1-way coupling
        if delta_volume_PCR_coupled[i] > 0:
            
            temp_water_volume_FM = delta_volume_PCR_coupled[i]/len(CouplePCR2FM[i][1])
            
            for j in range(len(CouplePCR2FM[i][1])):
                
                # get current FM cell index
                current_fm_cell_index = CouplePCR2FM[i][1][j]
                # calculate additional water level for current FM cell
                additional_water_level_FM[current_fm_cell_index] = temp_water_volume_FM / FMcellAreaSpherical[current_fm_cell_index]
            
        # or, if delta volume is negative, perform the new procedure which will remove water from FM instead of adding to it
        # (and if delta volume exactly equals zero, do nothing)
        elif delta_volume_PCR_coupled[i] < 0:
            
            # get all coupled FM cells of the current coupled PCR cell
            current_fm_cell_indices = CouplePCR2FM[i][1]
            
            # get current water depths of these cells
            current_water_depths_FM = model_fm.get_var('s1')[current_fm_cell_indices] - bottom_lvl_fm[current_fm_cell_indices]
            
            # set variable for keeping track of water depths
            new_water_depths_FM = np.copy(current_water_depths_FM)
            
            # remaining 'delta volume'
            remaining_delta_water_volume_FM = delta_volume_PCR_coupled[i]
            
            # while there is 'delta volume' remaining, remove water from wet cells
            while remaining_delta_water_volume_FM < 0:
                
                # find remaining wet cells
                remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                
                # if no remaining wet cells are found, no more water can be removed and the while-loop should stop
                # this probably means that the remaining delta volume is very low (but still just above zero)
                if len(remaining_wet_cells_FM_temp_indices) != 0:
                    pass
                else:
                    if print_delta_volume_info == 1:
                        print 'No more wet cells found! Info:'
                        print 'PCR cell single index: ', CouplePCR2FM[i][0]
                        print 'Remaining delta volume:', remaining_delta_water_volume_FM
                    break
                
                # find cell with least amount of water
                min_water_depth_FM_value = np.min(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                min_water_depth_FM_index = np.argmin(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                
                # calculate the total volume that is about to be removed from FM
                total_volume_about_to_be_removed = np.sum(min_water_depth_FM_value * \
                                                          FMcellAreaSpherical[remaining_wet_cells_FM_total_indices])
                    
                # check if removing this from all wet cells will not exceed the remaining 'delta volume'
                if (-1 * remaining_delta_water_volume_FM) > total_volume_about_to_be_removed:
                    
                    # remove water from all wet cells
                    new_water_depths_FM[remaining_wet_cells_FM_temp_indices] -= min_water_depth_FM_value
                    
                    # remove corresponding volume from remaining 'delta volume'
                    remaining_delta_water_volume_FM += total_volume_about_to_be_removed
                    
                # or, if volume to be removed would exceed 'delta volume', remove less than this instead
                # for this I identified 2 possible solutions (see notebook 'BMI_FM_PCR_coupled_run_2way_test_delta'), using 1 of those here
                else:
                    
                    # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                    remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                    
                    # remove water from all wet cells
                    new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                    
                    # check if there are negative water depths, and if so, repeat process
                    # -------------------------------------------------------------------
                    while any(new_water_depths_FM < 0):
                        
                        # reset negative water depths to zero and count the 'missing' volume
                        temp_vol_was_negative = 0
                        for k in range(len(new_water_depths_FM)):
                            if new_water_depths_FM[k] < 0:
                                temp_vol_was_negative -= (new_water_depths_FM[k] * -1) * FMcellAreaSpherical[current_fm_cell_indices[k]]
                                new_water_depths_FM[k] = 0.
                                
                        # this represents the volume that is being 'added back' to the delta volume to then be removed again from wet cells
                        remaining_delta_water_volume_FM = temp_vol_was_negative
                        
                        # find remaining wet cells
                        remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                        remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                        
                        # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                        remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                    FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                            
                        # remove water from all wet cells
                        new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                    # -------------------------------------------------------------------
                    
                    # set 'delta volume' to zero (so while-loop will end)
                    remaining_delta_water_volume_FM = 0.
        
            # adjusting array to be used with BMI 'set_var' to remove water from FM for negative delta volumes
            new_water_levels_fm_1[current_fm_cell_indices] = bottom_lvl_fm[current_fm_cell_indices] + new_water_depths_FM
    
    # or else, if the river/floodplain scheme is to be used:
    elif use_river_floodplain_scheme == 1:
        
        # check if there are any FM river cells for this PCR cell, in which case a distinction between river/floodplain cells is made
        if boolean_river_cell_in_coupled_PCR_cell[i]:
            
            # first up, river cells:
            
            # check if the river delta volume is positive, and if so, add water to FM's river cells
            if delta_volume_PCR_rivers_coupled[i] > 0:
                
                temp_water_volume_river_FM = delta_volume_PCR_rivers_coupled[i]/len(FM_river_cells_per_coupled_PCR_cell[i])
                
                for j in range(len(FM_river_cells_per_coupled_PCR_cell[i])):
                    
                    # get current FM cell index
                    current_fm_cell_index = FM_river_cells_per_coupled_PCR_cell[i][j]
                    # calculate additional water level for current FM cell
                    additional_water_level_FM[current_fm_cell_index] = temp_water_volume_river_FM / FMcellAreaSpherical[current_fm_cell_index]
            
            # or, if the river delta volume is negative, remove water from FM's river cells instead of adding to it
            # (and if delta volume exactly equals zero, do nothing)
            elif delta_volume_PCR_rivers_coupled[i] < 0:
                
                # get all coupled FM cells of the current coupled PCR cell
                current_fm_cell_indices = np.array(FM_river_cells_per_coupled_PCR_cell[i])
                
                # get current water depths of these cells
                current_water_depths_FM = model_fm.get_var('s1')[current_fm_cell_indices] - bottom_lvl_fm[current_fm_cell_indices]
                
                # set variable for keeping track of water depths
                new_water_depths_FM = np.copy(current_water_depths_FM)
                
                # remaining 'delta volume'
                remaining_delta_water_volume_FM = delta_volume_PCR_rivers_coupled[i]
                
                # while there is 'delta volume' remaining, remove water from wet cells
                while remaining_delta_water_volume_FM < 0:
                    
                    # find remaining wet cells
                    remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                    remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                    
                    # if no remaining wet cells are found, no more water can be removed and the while-loop should stop
                    # this probably means that the remaining delta volume is very low (but still just above zero)
                    if len(remaining_wet_cells_FM_temp_indices) != 0:
                        pass
                    else:
                        if print_delta_volume_info == 1:
                            print 'No more wet cells found! Info:'
                            print 'PCR cell single index: ', CouplePCR2FM[i][0]
                            print 'Remaining delta volume:', remaining_delta_water_volume_FM
                        break
                    
                    # find cell with least amount of water
                    min_water_depth_FM_value = np.min(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                    min_water_depth_FM_index = np.argmin(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                    
                    # calculate the total volume that is about to be removed from FM
                    total_volume_about_to_be_removed = np.sum(min_water_depth_FM_value * \
                                                              FMcellAreaSpherical[remaining_wet_cells_FM_total_indices])
                        
                    # check if removing this from all wet cells will not exceed the remaining 'delta volume'
                    if (-1 * remaining_delta_water_volume_FM) > total_volume_about_to_be_removed:
                        
                        # remove water from all wet cells
                        new_water_depths_FM[remaining_wet_cells_FM_temp_indices] -= min_water_depth_FM_value
                        
                        # remove corresponding volume from remaining 'delta volume'
                        remaining_delta_water_volume_FM += total_volume_about_to_be_removed
                        
                    # or, if volume to be removed would exceed 'delta volume', remove less than this instead
                    # for this I identified 2 possible solutions (see notebook 'BMI_FM_PCR_coupled_run_2way_test_delta'), using 1 of those here
                    else:
                        
                        # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                        remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                    FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                        
                        # remove water from all wet cells
                        new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                        
                        # check if there are negative water depths, and if so, repeat process
                        # -------------------------------------------------------------------
                        while any(new_water_depths_FM < 0):
                            
                            # reset negative water depths to zero and count the 'missing' volume
                            temp_vol_was_negative = 0
                            for k in range(len(new_water_depths_FM)):
                                if new_water_depths_FM[k] < 0:
                                    temp_vol_was_negative -= (new_water_depths_FM[k] * -1) * FMcellAreaSpherical[current_fm_cell_indices[k]]
                                    new_water_depths_FM[k] = 0.
                                    
                            # this represents the volume that is being 'added back' to the delta volume to then be removed again from wet cells
                            remaining_delta_water_volume_FM = temp_vol_was_negative
                            
                            # find remaining wet cells
                            remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                            remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                            
                            # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                            remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                        FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                                
                            # remove water from all wet cells
                            new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                        # -------------------------------------------------------------------
                        
                        # set 'delta volume' to zero (so while-loop will end)
                        remaining_delta_water_volume_FM = 0.
            
                # adjusting array to be used with BMI 'set_var' to remove water from FM for negative delta volumes
                new_water_levels_fm_1[current_fm_cell_indices] = bottom_lvl_fm[current_fm_cell_indices] + new_water_depths_FM
            
            # and then, floodplain cells:
            
            # check if the floodplain delta volume is positive, and if so, add water to FM's floodplain cells
            if delta_volume_PCR_floodplains_coupled[i] > 0:
                
                temp_water_volume_river_FM = delta_volume_PCR_floodplains_coupled[i]/len(FM_floodplain_cells_per_coupled_PCR_cell[i])
                
                for j in range(len(FM_floodplain_cells_per_coupled_PCR_cell[i])):
                    
                    # get current FM cell index
                    current_fm_cell_index = FM_floodplain_cells_per_coupled_PCR_cell[i][j]
                    # calculate additional water level for current FM cell
                    additional_water_level_FM[current_fm_cell_index] = temp_water_volume_river_FM / FMcellAreaSpherical[current_fm_cell_index]
            
            # or, if the floodplain delta volume is negative, remove water from FM's floodplain cells instead of adding to it
            # (and if delta volume exactly equals zero, do nothing)
            elif delta_volume_PCR_floodplains_coupled[i] < 0:
                
                # get all coupled FM cells of the current coupled PCR cell
                current_fm_cell_indices = np.array(FM_floodplain_cells_per_coupled_PCR_cell[i])
                
                # get current water depths of these cells
                current_water_depths_FM = model_fm.get_var('s1')[current_fm_cell_indices] - bottom_lvl_fm[current_fm_cell_indices]
                
                # set variable for keeping track of water depths
                new_water_depths_FM = np.copy(current_water_depths_FM)
                
                # remaining 'delta volume'
                remaining_delta_water_volume_FM = delta_volume_PCR_floodplains_coupled[i]
                
                # while there is 'delta volume' remaining, remove water from wet cells
                while remaining_delta_water_volume_FM < 0:
                    
                    # find remaining wet cells
                    remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                    remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                    
                    # if no remaining wet cells are found, no more water can be removed and the while-loop should stop
                    # this probably means that the remaining delta volume is very low (but still just above zero)
                    if len(remaining_wet_cells_FM_temp_indices) != 0:
                        pass
                    else:
                        if print_delta_volume_info == 1:
                            print 'No more wet cells found! Info:'
                            print 'PCR cell single index: ', CouplePCR2FM[i][0]
                            print 'Remaining delta volume:', remaining_delta_water_volume_FM
                        break
                    
                    # find cell with least amount of water
                    min_water_depth_FM_value = np.min(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                    min_water_depth_FM_index = np.argmin(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                    
                    # calculate the total volume that is about to be removed from FM
                    total_volume_about_to_be_removed = np.sum(min_water_depth_FM_value * \
                                                              FMcellAreaSpherical[remaining_wet_cells_FM_total_indices])
                        
                    # check if removing this from all wet cells will not exceed the remaining 'delta volume'
                    if (-1 * remaining_delta_water_volume_FM) > total_volume_about_to_be_removed:
                        
                        # remove water from all wet cells
                        new_water_depths_FM[remaining_wet_cells_FM_temp_indices] -= min_water_depth_FM_value
                        
                        # remove corresponding volume from remaining 'delta volume'
                        remaining_delta_water_volume_FM += total_volume_about_to_be_removed
                        
                    # or, if volume to be removed would exceed 'delta volume', remove less than this instead
                    # for this I identified 2 possible solutions (see notebook 'BMI_FM_PCR_coupled_run_2way_test_delta'), using 1 of those here
                    else:
                        
                        # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                        remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                    FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                        
                        # remove water from all wet cells
                        new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                        
                        # check if there are negative water depths, and if so, repeat process
                        # -------------------------------------------------------------------
                        while any(new_water_depths_FM < 0):
                            
                            # reset negative water depths to zero and count the 'missing' volume
                            temp_vol_was_negative = 0
                            for k in range(len(new_water_depths_FM)):
                                if new_water_depths_FM[k] < 0:
                                    temp_vol_was_negative -= (new_water_depths_FM[i] * -1) * FMcellAreaSpherical[current_fm_cell_indices[k]]
                                    new_water_depths_FM[k] = 0.
                                    
                            # this represents the volume that is being 'added back' to the delta volume to then be removed again from wet cells
                            remaining_delta_water_volume_FM = temp_vol_was_negative
                            
                            # find remaining wet cells
                            remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                            remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                            
                            # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                            remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                        FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                                
                            # remove water from all wet cells
                            new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                        # -------------------------------------------------------------------
                        
                        # set 'delta volume' to zero (so while-loop will end)
                        remaining_delta_water_volume_FM = 0.
            
                # adjusting array to be used with BMI 'set_var' to remove water from FM for negative delta volumes
                new_water_levels_fm_1[current_fm_cell_indices] = bottom_lvl_fm[current_fm_cell_indices] + new_water_depths_FM
        
        # or if there are no FM river cells, use the same same procedure as in the original method
        else:
            
            # add the river (from PCR 'discharge') and floodplain (from PCR 'runoff') volumes to obain the total delta volume
            delta_volume_PCR_total_coupled = delta_volume_PCR_rivers_coupled[i] + delta_volume_PCR_floodplains_coupled[i]
            
            # check if the total delta volume is positive, and if so, add water to the relevant FM cells
            if delta_volume_PCR_total_coupled > 0:
                
                temp_water_volume_FM = delta_volume_PCR_total_coupled/len(CouplePCR2FM[i][1])
                
                for j in range(len(CouplePCR2FM[i][1])):
                    
                    # get current FM cell index
                    current_fm_cell_index = CouplePCR2FM[i][1][j]
                    # calculate additional water level for current FM cell
                    additional_water_level_FM[current_fm_cell_index] = temp_water_volume_FM / FMcellAreaSpherical[current_fm_cell_index]
            
            # or, if delta volume is negative, perform the new procedure which will remove water from FM instead of adding to it
            # (and if delta volume exactly equals zero, do nothing)
            elif delta_volume_PCR_total_coupled < 0:
                
                # get all coupled FM cells of the current coupled PCR cell
                current_fm_cell_indices = CouplePCR2FM[i][1]
                
                # get current water depths of these cells
                current_water_depths_FM = model_fm.get_var('s1')[current_fm_cell_indices] - bottom_lvl_fm[current_fm_cell_indices]
                
                # set variable for keeping track of water depths
                new_water_depths_FM = np.copy(current_water_depths_FM)
                
                # remaining 'delta volume'
                remaining_delta_water_volume_FM = delta_volume_PCR_total_coupled
                
                # while there is 'delta volume' remaining, remove water from wet cells
                while remaining_delta_water_volume_FM < 0:
                    
                    # find remaining wet cells
                    remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                    remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                    
                    # if no remaining wet cells are found, no more water can be removed and the while-loop should stop
                    # this probably means that the remaining delta volume is very low (but still just above zero)
                    if len(remaining_wet_cells_FM_temp_indices) != 0:
                        pass
                    else:
                        if print_delta_volume_info == 1:
                            print 'No more wet cells found! Info:'
                            print 'PCR cell single index: ', CouplePCR2FM[i][0]
                            print 'Remaining delta volume:', remaining_delta_water_volume_FM
                        break
                    
                    # find cell with least amount of water
                    min_water_depth_FM_value = np.min(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                    min_water_depth_FM_index = np.argmin(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                    
                    # calculate the total volume that is about to be removed from FM
                    total_volume_about_to_be_removed = np.sum(min_water_depth_FM_value * \
                                                              FMcellAreaSpherical[remaining_wet_cells_FM_total_indices])
                        
                    # check if removing this from all wet cells will not exceed the remaining 'delta volume'
                    if (-1 * remaining_delta_water_volume_FM) > total_volume_about_to_be_removed:
                        
                        # remove water from all wet cells
                        new_water_depths_FM[remaining_wet_cells_FM_temp_indices] -= min_water_depth_FM_value
                        
                        # remove corresponding volume from remaining 'delta volume'
                        remaining_delta_water_volume_FM += total_volume_about_to_be_removed
                        
                    # or, if volume to be removed would exceed 'delta volume', remove less than this instead
                    # for this I identified 2 possible solutions (see notebook 'BMI_FM_PCR_coupled_run_2way_test_delta'), using 1 of those here
                    else:
                        
                        # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                        remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                    FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                        
                        # remove water from all wet cells
                        new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                        
                        # check if there are negative water depths, and if so, repeat process
                        # -------------------------------------------------------------------
                        while any(new_water_depths_FM < 0):
                            
                            # reset negative water depths to zero and count the 'missing' volume
                            temp_vol_was_negative = 0
                            for k in range(len(new_water_depths_FM)):
                                if new_water_depths_FM[k] < 0:
                                    temp_vol_was_negative -= (new_water_depths_FM[k] * -1) * FMcellAreaSpherical[current_fm_cell_indices[k]]
                                    new_water_depths_FM[k] = 0.
                                    
                            # this represents the volume that is being 'added back' to the delta volume to then be removed again from wet cells
                            remaining_delta_water_volume_FM = temp_vol_was_negative
                            
                            # find remaining wet cells
                            remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                            remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                            
                            # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                            remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                        FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                                
                            # remove water from all wet cells
                            new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                        # -------------------------------------------------------------------
                        
                        # set 'delta volume' to zero (so while-loop will end)
                        remaining_delta_water_volume_FM = 0.
            
                # adjusting array to be used with BMI 'set_var' to remove water from FM for negative delta volumes
                new_water_levels_fm_1[current_fm_cell_indices] = bottom_lvl_fm[current_fm_cell_indices] + new_water_depths_FM

# adjusting 's1' for negative delta volumes
model_fm.set_var('s1', new_water_levels_fm_1)

# calculate additional water levels for every FM update step (based on positive delta volumes)
delta_water_levels_fm = additional_water_level_FM / fraction_timestep_fm
# ----------------------------------------------------------------------------------------------------

# Update FM
# ----------------------------------------------------------------------------------------------------
# update FM unless it has has reached the same time as PCR
while model_fm.get_current_time() < (model_pcr.get_time_step() * secPerDay):
    
    # print FM model time to screen, to check progress and speed of model
    if print_timesteps == 1:
        print 'current FM model time: ', model_fm.get_current_time()
    
    # get current FM water levels
    current_water_levels_fm = model_fm.get_var('s1')
    
    # calculate new water levels for fm
    
    # empty array for filling in values
    new_water_levels_fm_2 = np.zeros(len(current_water_levels_fm))
    
    # loop over cells
    for i in range(len(new_water_levels_fm_2)):
        
        # if at "known" cells:
        if i < len(delta_water_levels_fm):
            
            new_water_levels_fm_2[i] = current_water_levels_fm[i] + delta_water_levels_fm[i]
        
        # for cells outside the range of "known" cells:
        else:
            
            new_water_levels_fm_2[i] = current_water_levels_fm[i]
    
    # set new water levels in fm model
    model_fm.set_var('s1', new_water_levels_fm_2)
    
    # update FM
    model_fm.update(update_step_fm)
    
    # save water levels and water depths as output
    # only every full day
    if save_fm_lists == 1:
        if model_fm.get_current_time()%secPerDay == 0.0:
            # calculate water depths
            water_depths_fm = model_fm.get_var('s1') - bottom_lvl_fm
            # append to lists
            output_water_levels_fm.append(model_fm.get_var('s1'))
            output_water_depths_fm.append(water_depths_fm)
            print 'saving fm output at time:', model_fm.get_current_time()
# ----------------------------------------------------------------------------------------------------

# COUPLED RUN (full run)

# coupled update
while model_pcr.get_time_step() < end_time_pcr:
    
    # CALCULATIONS FOR 2-WAY COUPLING
    # ---------------------------------------------------------------------------------------------------------------
    
    # 1. calculating fraction and area of (inundated) water per coupled PCR cell
    
    # temporary variable to store total inundated area to-be-stored to list
    temp_inundated_area_threshold_list = 0.
    
    # get water depth of all FM cells
    current_water_depth_fm = model_fm.get_var('s1') - bottom_lvl_fm
    
    # if the original method (not the river/floodplain scheme) is to be used:
    if use_river_floodplain_scheme == 0:
        
        # create zero array for filling in total area and fraction of flooded FM cells for each coupled PCR cell
        inundated_area_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
        inundated_fraction_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
        
        # create zero array of appropriate size for using BMI function 'set_var'
        inundated_fraction_FM_2_PCR = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
        
        # loop over all coupled PCR cells and fill in and/or calculate values
        for i in range(len(CouplePCR2FM)):
            
            # get total area of all flooded FM cells coupled to current PCR cell
            
            # temporary variable for storing total area of flooded FM cells coupled to current PCR cell
            temp_inundated_area_FM = 0.
            
            # loop over all coupled FM cells
            for j in range(len(CouplePCR2FM[i][1])):
                
                # get current FM cell index
                current_FM_cell_index = CouplePCR2FM[i][1][j]
                
                # check if water depth of current FM cell is above chosen threshold
                if current_water_depth_fm[current_FM_cell_index] > threshold_inundated_depth:
                    
                    # if so, add cell area to temporary variable
                    temp_inundated_area_FM += FMcellAreaSpherical[current_FM_cell_index]
                    
                    # check if current FM cell falls within specified region and if water depth is above threshold
                    if (current_FM_cell_index in area_pcr_indices_FM_floodplain_cells) and \
                       (current_water_depth_fm[current_FM_cell_index] >= area_pcr_threshold_save_to_list):
                    
                        temp_inundated_area_threshold_list += FMcellAreaSpherical[current_FM_cell_index]
            
            # at end of loop, assign temporary variable to array storing inundated areas
            inundated_area_FM_2_PCR_coupled[i] = temp_inundated_area_FM
            
            # limit inundated area to area of coupled PCR cell
            #inundated_area_FM_2_PCR_coupled[i] = min([inundated_area_FM_2_PCR_coupled[i], \
            #                                          cellarea_data_pcr[CoupledPCRcellIndices[i]]])
            
            # divide this value by PCR cell area to obtain flooded area fraction
            inundated_fraction_FM_2_PCR_coupled[i] = inundated_area_FM_2_PCR_coupled[i] / \
                                                        cellarea_data_pcr[CoupledPCRcellIndices[i]]
            
            # limit inundated fraction between 0 and 1
            #inundated_fraction_FM_2_PCR_coupled[i] = max([0., min([inundated_fraction_FM_2_PCR_coupled[i], 1.])])
            
            # assign value of current cell to zero array to be used with BMI
            inundated_fraction_FM_2_PCR[CoupledPCRcellIndices[i]] = inundated_fraction_FM_2_PCR_coupled[i]
    
    # or else, if the river/floodplain scheme is to be used:
    elif use_river_floodplain_scheme == 1:
        
        # create zero array for filling in total area and fraction of flooded FM cells for each coupled PCR cell
        inundated_area_rivers_FM_2_PCR_coupled          = np.zeros(len(CouplePCR2FM))
        inundated_fraction_rivers_FM_2_PCR_coupled      = np.zeros(len(CouplePCR2FM))
        
        inundated_area_floodplains_FM_2_PCR_coupled     = np.zeros(len(CouplePCR2FM))
        inundated_fraction_floodplains_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
        
        # create zero array of appropriate size for using BMI function 'set_var'
        inundated_fraction_rivers_FM_2_PCR      = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
        inundated_fraction_floodplains_FM_2_PCR = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
        
        # loop over all coupled PCR cells and fill in and/or calculate values
        for i in range(len(CouplePCR2FM)):
            
            # check if there are any FM river cells for this PCR cell
            # (in which case a distinction between river/floodplain cells is made)
            if boolean_river_cell_in_coupled_PCR_cell[i]:
                
                # 1a. calculating fraction and area of wet FM river cells per coupled PCR cell
                
                # temporary variable for storing total area of wet FM river cells coupled to current PCR cell
                temp_inundated_area_FM = 0.
                
                # loop over all river cells coupled to current PCR cell
                for j in range(len(FM_river_cells_per_coupled_PCR_cell[i])):
                
                    # get current FM cell index
                    current_FM_cell_index = FM_river_cells_per_coupled_PCR_cell[i][j]
                    
                    # check if water depth of current FM cell is above chosen threshold
                    if current_water_depth_fm[current_FM_cell_index] > threshold_inundated_depth_rivers:
                        
                        # if so, add cell area to temporary variable
                        temp_inundated_area_FM += FMcellAreaSpherical[current_FM_cell_index]
                
                # at end of loop, assign temporary variable to array storing inundated areas
                inundated_area_rivers_FM_2_PCR_coupled[i] = temp_inundated_area_FM
                
                # limit inundated area to area of coupled PCR cell
                #inundated_area_rivers_FM_2_PCR_coupled[i] = min([inundated_area_rivers_FM_2_PCR_coupled[i], \
                #                                                      cellarea_data_pcr[CoupledPCRcellIndices[i]]])
                
                # divide this value by PCR cell area to obtain flooded area fraction
                inundated_fraction_rivers_FM_2_PCR_coupled[i] = inundated_area_rivers_FM_2_PCR_coupled[i] / \
                                                                        cellarea_data_pcr[CoupledPCRcellIndices[i]]
                
                # limit inundated fraction between 0 and 1
                #inundated_fraction_rivers_FM_2_PCR_coupled[i] = max([0., min([inundated_fraction_rivers_FM_2_PCR_coupled[i], 1.])])
                
                # assign value of current cell to zero array to be used with BMI
                inundated_fraction_rivers_FM_2_PCR[CoupledPCRcellIndices[i]] = inundated_fraction_rivers_FM_2_PCR_coupled[i]
                
                # 1b. calculating fraction and area of inundated FM floodplain cells per coupled PCR cell
                
                # temporary variable for storing total area of inundated FM floodplain cells coupled to current PCR cell
                temp_inundated_area_FM = 0.
                
                # loop over all floodplain cells coupled to current PCR cell
                for j in range(len(FM_floodplain_cells_per_coupled_PCR_cell[i])):
                
                    # get current FM cell index
                    current_FM_cell_index = FM_floodplain_cells_per_coupled_PCR_cell[i][j]
                    
                    # check if water depth of current FM cell is above chosen threshold
                    if current_water_depth_fm[current_FM_cell_index] > threshold_inundated_depth_floodplains:
                        
                        # if so, add cell area to temporary variable
                        temp_inundated_area_FM += FMcellAreaSpherical[current_FM_cell_index]
                        
                        # check if current FM cell falls within specified region and if water depth is above threshold
                        if (current_FM_cell_index in area_pcr_indices_FM_floodplain_cells) and \
                           (current_water_depth_fm[current_FM_cell_index] >= area_pcr_threshold_save_to_list):
                        
                            temp_inundated_area_threshold_list += FMcellAreaSpherical[current_FM_cell_index]
                
                # at end of loop, assign temporary variable to array storing inundated areas
                inundated_area_floodplains_FM_2_PCR_coupled[i] = temp_inundated_area_FM
                
                # limit inundated area to area of coupled PCR cell
                #inundated_area_floodplains_FM_2_PCR_coupled[i] = min([inundated_area_floodplains_FM_2_PCR_coupled[i], \
                #                                                      cellarea_data_pcr[CoupledPCRcellIndices[i]]])
                
                # divide this value by PCR cell area to obtain flooded area fraction
                inundated_fraction_floodplains_FM_2_PCR_coupled[i] = inundated_area_floodplains_FM_2_PCR_coupled[i] / \
                                                                        cellarea_data_pcr[CoupledPCRcellIndices[i]]
                
                # limit inundated fraction between 0 and 1
                #inundated_fraction_floodplains_FM_2_PCR_coupled[i] = max([0., min([inundated_fraction_floodplains_FM_2_PCR_coupled[i], \
                #                                                                   1.])])
                
                # assign value of current cell to zero array to be used with BMI
                inundated_fraction_floodplains_FM_2_PCR[CoupledPCRcellIndices[i]] = inundated_fraction_floodplains_FM_2_PCR_coupled[i]
                
            # or if there are no FM river cells, use the same same procedure as in the original method
            else:
                    
                # get total area of all flooded FM cells coupled to current PCR cell
                
                # temporary variable for storing total area of flooded FM cells coupled to current PCR cell
                temp_inundated_area_FM = 0.
                
                # loop over all coupled FM cells
                for j in range(len(CouplePCR2FM[i][1])):
                    
                    # get current FM cell index
                    current_FM_cell_index = CouplePCR2FM[i][1][j]
                    
                    # check if water depth of current FM cell is above chosen threshold
                    if current_water_depth_fm[current_FM_cell_index] > threshold_inundated_depth_floodplains:
                        
                        # if so, add cell area to temporary variable
                        temp_inundated_area_FM += FMcellAreaSpherical[current_FM_cell_index]
                        
                        # check if current FM cell falls within specified region and if water depth is above threshold
                        if (current_FM_cell_index in area_pcr_indices_FM_floodplain_cells) and \
                           (current_water_depth_fm[current_FM_cell_index] >= area_pcr_threshold_save_to_list):
                        
                            temp_inundated_area_threshold_list += FMcellAreaSpherical[current_FM_cell_index]
                
                # at end of loop, assign temporary variable to array storing inundated areas
                inundated_area_floodplains_FM_2_PCR_coupled[i] = temp_inundated_area_FM
                
                # limit inundated area to area of coupled PCR cell
                #inundated_area_floodplains_FM_2_PCR_coupled[i] = min([inundated_area_floodplains_FM_2_PCR_coupled[i], \
                #                                                      cellarea_data_pcr[CoupledPCRcellIndices[i]]])
                
                # divide this value by PCR cell area to obtain flooded area fraction
                inundated_fraction_floodplains_FM_2_PCR_coupled[i] = inundated_area_floodplains_FM_2_PCR_coupled[i] / \
                                                            cellarea_data_pcr[CoupledPCRcellIndices[i]]
                
                # limit inundated fraction between 0 and 1
                #inundated_fraction_floodplains_FM_2_PCR_coupled[i] = max([0., min([inundated_fraction_floodplains_FM_2_PCR_coupled[i], \
                #                                                                   1.])])
                
                # assign value of current cell to zero array to be used with BMI
                inundated_fraction_floodplains_FM_2_PCR[CoupledPCRcellIndices[i]] = inundated_fraction_floodplains_FM_2_PCR_coupled[i]
    
    # 2. extracting FM water volume for each coupled PCR cell and converting these to water depths to add back to PCR
    
    # get water volume of all FM cells
    current_volume_fm = model_fm.get_var('vol1')
    
    # if the original method (not the river/floodplain scheme) is to be used:
    if use_river_floodplain_scheme == 0:
        
        # create zero array for filling in total volumes and water depths for each coupled PCR cell
        water_volume_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
        water_depths_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
        
        # create zero array of appropriate size for using BMI function 'set_var'
        water_depths_FM_2_PCR = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
        
        # loop over all coupled PCR cells and fill in and/or calculate values
        for i in range(len(CouplePCR2FM)):
            
            # get total volume of all FM cells coupled to current PCR cell
            water_volume_FM_2_PCR_coupled[i] = np.sum(current_volume_fm[CouplePCR2FM[i][1]])
            
            # OLD: divide total volume by PCR cell area to obtain water depth
            #water_depths_FM_2_PCR_coupled[i] = water_volume_FM_2_PCR_coupled[i] / cellarea_data_pcr[CoupledPCRcellIndices[i]]
            
            # NEW: divide total volume by inundated FM area to obtain water depth
            if inundated_area_FM_2_PCR_coupled[i] > 0. :
                water_depths_FM_2_PCR_coupled[i] = water_volume_FM_2_PCR_coupled[i] / inundated_area_FM_2_PCR_coupled[i]
            else:
                water_depths_FM_2_PCR_coupled[i] = 0.
            
            # assign value of current cell to zero array to be used with BMI
            water_depths_FM_2_PCR[CoupledPCRcellIndices[i]] = water_depths_FM_2_PCR_coupled[i]
    
    # or else, if the river/floodplain scheme is to be used:
    elif use_river_floodplain_scheme == 1:
        
        # create zero array for filling in total volumes and water depths for each coupled PCR cell
        water_volume_rivers_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
        
        water_volume_floodplains_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
        water_depths_floodplains_FM_2_PCR_coupled = np.zeros(len(CouplePCR2FM))
        
        # create zero array of appropriate size for using BMI function 'set_var'
        water_volume_rivers_FM_2_PCR      = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
        water_depths_floodplains_FM_2_PCR = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
        
        # loop over all coupled PCR cells and fill in and/or calculate values
        for i in range(len(CouplePCR2FM)):
            
            # check if there are any FM river cells for this PCR cell, in which case a distinction between river/floodplain cells is made
            if boolean_river_cell_in_coupled_PCR_cell[i]:
                
                # 2a. extracting water volume from FM river cells
                
                # get total volume of all FM cells coupled to current PCR cell
                water_volume_rivers_FM_2_PCR_coupled[i] = np.sum(current_volume_fm[FM_river_cells_per_coupled_PCR_cell[i]])
                
                # assign value of current cell to zero array to be used with BMI
                water_volume_rivers_FM_2_PCR[CoupledPCRcellIndices[i]] = water_volume_rivers_FM_2_PCR_coupled[i]
                
                # 2b. extracting water volume from FM floodplain cells
                
                # get total volume of all FM cells coupled to current PCR cell
                water_volume_floodplains_FM_2_PCR_coupled[i] = np.sum(current_volume_fm[FM_floodplain_cells_per_coupled_PCR_cell[i]])
                
                # divide total volume by inundated FM area to obtain water depth
                if inundated_area_floodplains_FM_2_PCR_coupled[i] > 0. :
                    water_depths_floodplains_FM_2_PCR_coupled[i] = water_volume_floodplains_FM_2_PCR_coupled[i] / \
                                                                    inundated_area_floodplains_FM_2_PCR_coupled[i]
                else:
                    water_depths_floodplains_FM_2_PCR_coupled[i] = 0.
                
                # assign value of current cell to zero array to be used with BMI
                water_depths_floodplains_FM_2_PCR[CoupledPCRcellIndices[i]] = water_depths_floodplains_FM_2_PCR_coupled[i]            
                
            # or if there are no FM river cells, use the same same procedure as in the original method
            else:
                
                # get total volume of all FM cells coupled to current PCR cell
                water_volume_floodplains_FM_2_PCR_coupled[i] = np.sum(current_volume_fm[CouplePCR2FM[i][1]])
                
                # divide total volume by inundated FM area to obtain water depth
                if inundated_area_floodplains_FM_2_PCR_coupled[i] > 0. :
                    water_depths_floodplains_FM_2_PCR_coupled[i] = water_volume_floodplains_FM_2_PCR_coupled[i] / \
                                                                    inundated_area_floodplains_FM_2_PCR_coupled[i]
                else:
                    water_depths_floodplains_FM_2_PCR_coupled[i] = 0.
                
                # assign value of current cell to zero array to be used with BMI
                water_depths_floodplains_FM_2_PCR[CoupledPCRcellIndices[i]] = water_depths_floodplains_FM_2_PCR_coupled[i]
                
    # save the total inundated area and volume within the specified area
    area_pcr_total_inundated_area_floodplains.append(np.sum(\
                                        inundated_area_floodplains_FM_2_PCR_coupled[area_pcr_indices_in_CoupledPCRcellIndices_final]))
    area_pcr_total_inundated_area_threshold.append(temp_inundated_area_threshold_list)
    area_pcr_total_flooded_volume_floodplains.append(np.sum(\
                                        water_depths_floodplains_FM_2_PCR_coupled[area_pcr_indices_in_CoupledPCRcellIndices_final]))
    
    # 3. Add the calculated volume of river water in FM to the current channel storage in PCR
    
    # get current channelStorage
    current_channelStorage_pcr  = model_pcr.get_var(('routing','channelStorage'))
    
    # initialize new map
    new_channelStorage_pcr = coupling_functions.zeroMapArray(landmask_pcr, missingValuesMap=missing_value_landmask)
    
    # loop over all indices and add values
    for i in range(len(new_channelStorage_pcr)):
        for j in range(len(new_channelStorage_pcr[0])):
            if current_channelStorage_pcr[i][j] != -999:
                new_channelStorage_pcr[i][j] = current_channelStorage_pcr[i][j] + water_volume_rivers_FM_2_PCR[i][j]
    # ---------------------------------------------------------------------------------------------------------------
    
    # USE BMI FUNCTIONS TO ADD RELEVANT VARIABLES FROM FM TO PCR
    # ---------------------------------------------------------------------------------------------------------------
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
    # ---------------------------------------------------------------------------------------------------------------
    
    # UPDATE PCR AND OBTAIN DELTA VOLUME
    # ---------------------------------------------------------------------------------------------------------------
    # update PCR
    model_pcr.update(1)
    
    # get relevant variables at current timestep
    current_discharge_pcr  = model_pcr.get_var('discharge')
    current_runoff_pcr     = model_pcr.get_var('landSurfaceRunoff')
    current_waterlayer_pcr = model_pcr.get_var('topWaterLayer')
    
    # check settings of using GRDC discharge data
    if use_GRDC_discharge_data == 1:
        current_timestep_pcr = model_pcr.get_time_step() - 1
        for i in range(len(indices_inflow_points)):
            current_discharge_pcr[indices_inflow_points[i]] += GRDC_discharge_input[current_timestep_pcr][i]
           
    # 1. Water volume for river cells
    
    # convert PCR discharge [m3/s] to water volume for river cells [m3/day]
    water_volume_PCR_rivers = np.zeros([len(current_discharge_pcr),len(current_discharge_pcr[0])])
    for i in range(len(current_discharge_pcr)):
        for j in range(len(current_discharge_pcr[0])):
            if current_discharge_pcr[i][j] != missing_value_pcr:
                water_volume_PCR_rivers[i][j] = current_discharge_pcr[i][j] * secPerDay
            else:
                water_volume_PCR_rivers[i][j] = 0.
            
    # get water volume for rivers for all coupled cells [m3/day]
    water_volume_PCR_rivers_coupled = np.zeros(len(CoupledPCRcellIndices))
    for i in range(len(CoupledPCRcellIndices)):
        water_volume_PCR_rivers_coupled[i] = water_volume_PCR_rivers[CoupledPCRcellIndices[i]]
    
    # 2. Water volume for floodplain cells
    
    # convert PCR runoff and top water layer [m/day] to water volumes [m3/day]
    water_volume_PCR_runoff      = np.zeros([len(current_runoff_pcr),len(current_runoff_pcr[0])])
    water_volume_PCR_waterlayer  = np.zeros([len(current_waterlayer_pcr),len(current_waterlayer_pcr[0])])
    for i in range(len(current_runoff_pcr)):
        for j in range(len(current_runoff_pcr[0])):
            if current_runoff_pcr[i][j] != missing_value_pcr:
                water_volume_PCR_runoff[i][j]     = current_runoff_pcr[i][j] * cellarea_data_pcr[i][j]
                water_volume_PCR_waterlayer[i][j] = current_waterlayer_pcr[i][j] * cellarea_data_pcr[i][j]
            else:
                water_volume_PCR_runoff[i][j]     = 0.
                water_volume_PCR_waterlayer[i][j] = 0.
                
    # get water volumes for all coupled cells [m3/day]
    water_volume_PCR_runoff_coupled      = np.zeros(len(CoupledPCRcellIndices))
    water_volume_PCR_waterlayer_coupled  = np.zeros(len(CoupledPCRcellIndices))
    water_volume_PCR_floodplains_coupled = np.zeros(len(CoupledPCRcellIndices))
    for i in range(len(CoupledPCRcellIndices)):
        water_volume_PCR_runoff_coupled[i]      = water_volume_PCR_runoff[CoupledPCRcellIndices[i]]
        water_volume_PCR_waterlayer_coupled[i]  = water_volume_PCR_waterlayer[CoupledPCRcellIndices[i]]
        water_volume_PCR_floodplains_coupled[i] = water_volume_PCR_runoff_coupled[i] + water_volume_PCR_waterlayer_coupled[i]
    
    # 3. Calculate delta volume for all coupled cells [m3/day]
    delta_volume_PCR_rivers_coupled      = water_volume_PCR_rivers_coupled - water_volume_rivers_FM_2_PCR_coupled
    delta_volume_PCR_floodplains_coupled = water_volume_PCR_floodplains_coupled - water_volume_floodplains_FM_2_PCR_coupled
    
    # 4. set PCR's 'topWaterLayer' to zero
    
    # get variable/map from PCR
    new_topWaterLayer_forest = model_pcr.get_var(('forest', 'topWaterLayer'))
    new_topWaterLayer_grass  = model_pcr.get_var(('grassland', 'topWaterLayer'))
    
    # adjust variable/map so it contains zeros at coupled cell locations
    for i in range(len(CouplePCR2FM)):
        new_topWaterLayer_forest[CoupledPCRcellIndices[i]] = 0.0
        new_topWaterLayer_grass[CoupledPCRcellIndices[i]]  = 0.0
    
    # set in PCR model for each land cover type using BMI function 'set_var'
    model_pcr.set_var(('forest','topWaterLayer'), new_topWaterLayer_forest)
    model_pcr.set_var(('grassland','topWaterLayer'), new_topWaterLayer_grass)
    # ---------------------------------------------------------------------------------------------------------------
    
    # CONVERT DELTA VOLUME TO NEW WATER LEVELS FOR FM
    # NOTE: CURRENTLY USING ALGORITHM 2 TO HANDLE NEGATIVE DELTA VOLUMES, MAYBE I ALSO HAVE TO CONSIDER ALGORITHM 1?!
    # ---------------------------------------------------------------------------------------------------------------
    # get current FM water levels
    new_water_levels_fm_1 = np.copy(model_fm.get_var('s1'))
    
    # initiate additional water levels array for each FM cell (same indices as FM cell index to easily add up later) [m/day]
    additional_water_level_FM = np.zeros(len(CoupleFM2PCR))
    
    # loop over all coupled PCR cells
    for i in range(len(CouplePCR2FM)):
        
        # if the original method (not the river/floodplain scheme) is to be used:
        if use_river_floodplain_scheme == 0:
            
            # check if delta volume is positive, and if so, perform the usual procedure that is also used in 1-way coupling
            if delta_volume_PCR_coupled[i] > 0:
                
                temp_water_volume_FM = delta_volume_PCR_coupled[i]/len(CouplePCR2FM[i][1])
                
                for j in range(len(CouplePCR2FM[i][1])):
                    
                    # get current FM cell index
                    current_fm_cell_index = CouplePCR2FM[i][1][j]
                    # calculate additional water level for current FM cell
                    additional_water_level_FM[current_fm_cell_index] = temp_water_volume_FM / FMcellAreaSpherical[current_fm_cell_index]
                
            # or, if delta volume is negative, perform the new procedure which will remove water from FM instead of adding to it
            # (and if delta volume exactly equals zero, do nothing)
            elif delta_volume_PCR_coupled[i] < 0:
                
                # get all coupled FM cells of the current coupled PCR cell
                current_fm_cell_indices = CouplePCR2FM[i][1]
                
                # get current water depths of these cells
                current_water_depths_FM = model_fm.get_var('s1')[current_fm_cell_indices] - bottom_lvl_fm[current_fm_cell_indices]
                
                # set variable for keeping track of water depths
                new_water_depths_FM = np.copy(current_water_depths_FM)
                
                # remaining 'delta volume'
                remaining_delta_water_volume_FM = delta_volume_PCR_coupled[i]
                
                # while there is 'delta volume' remaining, remove water from wet cells
                while remaining_delta_water_volume_FM < 0:
                    
                    # find remaining wet cells
                    remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                    remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                    
                    # if no remaining wet cells are found, no more water can be removed and the while-loop should stop
                    # this probably means that the remaining delta volume is very low (but still just above zero)
                    if len(remaining_wet_cells_FM_temp_indices) != 0:
                        pass
                    else:
                        if print_delta_volume_info == 1:
                            print 'No more wet cells found! Info:'
                            print 'PCR cell single index: ', CouplePCR2FM[i][0]
                            print 'Remaining delta volume:', remaining_delta_water_volume_FM
                        break
                    
                    # find cell with least amount of water
                    min_water_depth_FM_value = np.min(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                    min_water_depth_FM_index = np.argmin(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                    
                    # calculate the total volume that is about to be removed from FM
                    total_volume_about_to_be_removed = np.sum(min_water_depth_FM_value * \
                                                              FMcellAreaSpherical[remaining_wet_cells_FM_total_indices])
                        
                    # check if removing this from all wet cells will not exceed the remaining 'delta volume'
                    if (-1 * remaining_delta_water_volume_FM) > total_volume_about_to_be_removed:
                        
                        # remove water from all wet cells
                        new_water_depths_FM[remaining_wet_cells_FM_temp_indices] -= min_water_depth_FM_value
                        
                        # remove corresponding volume from remaining 'delta volume'
                        remaining_delta_water_volume_FM += total_volume_about_to_be_removed
                        
                    # or, if volume to be removed would exceed 'delta volume', remove less than this instead
                    # (for this I identified 2 possible solutions [see notebook 'BMI_FM_PCR_coupled_run_2way_test_delta']
                    #  using 1 of those here)
                    else:
                        
                        # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                        remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                    FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                        
                        # remove water from all wet cells
                        new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                        
                        # check if there are negative water depths, and if so, repeat process
                        # -------------------------------------------------------------------
                        while any(new_water_depths_FM < 0):
                            
                            # reset negative water depths to zero and count the 'missing' volume
                            temp_vol_was_negative = 0
                            for k in range(len(new_water_depths_FM)):
                                if new_water_depths_FM[k] < 0:
                                    temp_vol_was_negative -= (new_water_depths_FM[k] * -1) * \
                                                                FMcellAreaSpherical[current_fm_cell_indices[k]]
                                    new_water_depths_FM[k] = 0.
                                    
                            # this represents the volume that is being 'added back' to the delta volume
                            remaining_delta_water_volume_FM = temp_vol_was_negative
                            
                            # find remaining wet cells
                            remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                            remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                            
                            # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                            remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                        FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                                
                            # remove water from all wet cells
                            new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                        # -------------------------------------------------------------------
                        
                        # set 'delta volume' to zero (so while-loop will end)
                        remaining_delta_water_volume_FM = 0.
            
                # adjusting array to be used with BMI 'set_var' to remove water from FM for negative delta volumes
                new_water_levels_fm_1[current_fm_cell_indices] = bottom_lvl_fm[current_fm_cell_indices] + new_water_depths_FM
        
        # or else, if the river/floodplain scheme is to be used:
        elif use_river_floodplain_scheme == 1:
            
            # check if there are any FM river cells for this PCR cell
            # (in which case a distinction between river/floodplain cells is made)
            if boolean_river_cell_in_coupled_PCR_cell[i]:
                
                # first up, river cells:
                
                # check if the river delta volume is positive, and if so, add water to FM's river cells
                if delta_volume_PCR_rivers_coupled[i] > 0:
                    
                    temp_water_volume_river_FM = delta_volume_PCR_rivers_coupled[i]/len(FM_river_cells_per_coupled_PCR_cell[i])
                    
                    for j in range(len(FM_river_cells_per_coupled_PCR_cell[i])):
                        
                        # get current FM cell index
                        current_fm_cell_index = FM_river_cells_per_coupled_PCR_cell[i][j]
                        # calculate additional water level for current FM cell
                        additional_water_level_FM[current_fm_cell_index] = temp_water_volume_river_FM / \
                                                                            FMcellAreaSpherical[current_fm_cell_index]
                
                # or, if the river delta volume is negative, remove water from FM's river cells instead of adding to it
                # (and if delta volume exactly equals zero, do nothing)
                elif delta_volume_PCR_rivers_coupled[i] < 0:
                    
                    # get all coupled FM cells of the current coupled PCR cell
                    current_fm_cell_indices = np.array(FM_river_cells_per_coupled_PCR_cell[i])
                    
                    # get current water depths of these cells
                    current_water_depths_FM = model_fm.get_var('s1')[current_fm_cell_indices] - bottom_lvl_fm[current_fm_cell_indices]
                    
                    # set variable for keeping track of water depths
                    new_water_depths_FM = np.copy(current_water_depths_FM)
                    
                    # remaining 'delta volume'
                    remaining_delta_water_volume_FM = delta_volume_PCR_rivers_coupled[i]
                    
                    # while there is 'delta volume' remaining, remove water from wet cells
                    while remaining_delta_water_volume_FM < 0:
                        
                        # find remaining wet cells
                        remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                        remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                        
                        # if no remaining wet cells are found, no more water can be removed and the while-loop should stop
                        # this probably means that the remaining delta volume is very low (but still just above zero)
                        if len(remaining_wet_cells_FM_temp_indices) != 0:
                            pass
                        else:
                            if print_delta_volume_info == 1:
                                print 'No more wet cells found! Info:'
                                print 'PCR cell single index: ', CouplePCR2FM[i][0]
                                print 'Remaining delta volume:', remaining_delta_water_volume_FM
                            break
                        
                        # find cell with least amount of water
                        min_water_depth_FM_value = np.min(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                        min_water_depth_FM_index = np.argmin(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                        
                        # calculate the total volume that is about to be removed from FM
                        total_volume_about_to_be_removed = np.sum(min_water_depth_FM_value * \
                                                                  FMcellAreaSpherical[remaining_wet_cells_FM_total_indices])
                            
                        # check if removing this from all wet cells will not exceed the remaining 'delta volume'
                        if (-1 * remaining_delta_water_volume_FM) > total_volume_about_to_be_removed:
                            
                            # remove water from all wet cells
                            new_water_depths_FM[remaining_wet_cells_FM_temp_indices] -= min_water_depth_FM_value
                            
                            # remove corresponding volume from remaining 'delta volume'
                            remaining_delta_water_volume_FM += total_volume_about_to_be_removed
                            
                        # or, if volume to be removed would exceed 'delta volume', remove less than this instead
                        # (for this I identified 2 possible solutions [see notebook 'BMI_FM_PCR_coupled_run_2way_test_delta'] 
                        #  using 1 of those here)
                        else:
                            
                            # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                            remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                        FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                            
                            # remove water from all wet cells
                            new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                            
                            # check if there are negative water depths, and if so, repeat process
                            # -------------------------------------------------------------------
                            while any(new_water_depths_FM < 0):
                                
                                # reset negative water depths to zero and count the 'missing' volume
                                temp_vol_was_negative = 0
                                for k in range(len(new_water_depths_FM)):
                                    if new_water_depths_FM[k] < 0:
                                        temp_vol_was_negative -= (new_water_depths_FM[k] * -1) * \
                                                                    FMcellAreaSpherical[current_fm_cell_indices[k]]
                                        new_water_depths_FM[k] = 0.
                                        
                                # this represents the volume that is being 'added back' to the delta volume
                                remaining_delta_water_volume_FM = temp_vol_was_negative
                                
                                # find remaining wet cells
                                remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                                remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                                
                                # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                                remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                            FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                                    
                                # remove water from all wet cells
                                new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                            # -------------------------------------------------------------------
                            
                            # set 'delta volume' to zero (so while-loop will end)
                            remaining_delta_water_volume_FM = 0.
                
                    # adjusting array to be used with BMI 'set_var' to remove water from FM for negative delta volumes
                    new_water_levels_fm_1[current_fm_cell_indices] = bottom_lvl_fm[current_fm_cell_indices] + new_water_depths_FM
                
                # and then, floodplain cells:
                
                # check if the floodplain delta volume is positive, and if so, add water to FM's floodplain cells
                if delta_volume_PCR_floodplains_coupled[i] > 0:
                    
                    temp_water_volume_river_FM = delta_volume_PCR_floodplains_coupled[i]/len(FM_floodplain_cells_per_coupled_PCR_cell[i])
                    
                    for j in range(len(FM_floodplain_cells_per_coupled_PCR_cell[i])):
                        
                        # get current FM cell index
                        current_fm_cell_index = FM_floodplain_cells_per_coupled_PCR_cell[i][j]
                        # calculate additional water level for current FM cell
                        additional_water_level_FM[current_fm_cell_index] = temp_water_volume_river_FM / \
                                                                            FMcellAreaSpherical[current_fm_cell_index]
                
                # or, if the floodplain delta volume is negative, remove water from FM's floodplain cells instead of adding to it
                # (and if delta volume exactly equals zero, do nothing)
                elif delta_volume_PCR_floodplains_coupled[i] < 0:
                    
                    # get all coupled FM cells of the current coupled PCR cell
                    current_fm_cell_indices = np.array(FM_floodplain_cells_per_coupled_PCR_cell[i])
                    
                    # get current water depths of these cells
                    current_water_depths_FM = model_fm.get_var('s1')[current_fm_cell_indices] - bottom_lvl_fm[current_fm_cell_indices]
                    
                    # set variable for keeping track of water depths
                    new_water_depths_FM = np.copy(current_water_depths_FM)
                    
                    # remaining 'delta volume'
                    remaining_delta_water_volume_FM = delta_volume_PCR_floodplains_coupled[i]
                    
                    # while there is 'delta volume' remaining, remove water from wet cells
                    while remaining_delta_water_volume_FM < 0:
                        
                        # find remaining wet cells
                        remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                        remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                        
                        # if no remaining wet cells are found, no more water can be removed and the while-loop should stop
                        # this probably means that the remaining delta volume is very low (but still just above zero)
                        if len(remaining_wet_cells_FM_temp_indices) != 0:
                            pass
                        else:
                            if print_delta_volume_info == 1:
                                print 'No more wet cells found! Info:'
                                print 'PCR cell single index: ', CouplePCR2FM[i][0]
                                print 'Remaining delta volume:', remaining_delta_water_volume_FM
                            break
                        
                        # find cell with least amount of water
                        min_water_depth_FM_value = np.min(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                        min_water_depth_FM_index = np.argmin(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                        
                        # calculate the total volume that is about to be removed from FM
                        total_volume_about_to_be_removed = np.sum(min_water_depth_FM_value * \
                                                                  FMcellAreaSpherical[remaining_wet_cells_FM_total_indices])
                            
                        # check if removing this from all wet cells will not exceed the remaining 'delta volume'
                        if (-1 * remaining_delta_water_volume_FM) > total_volume_about_to_be_removed:
                            
                            # remove water from all wet cells
                            new_water_depths_FM[remaining_wet_cells_FM_temp_indices] -= min_water_depth_FM_value
                            
                            # remove corresponding volume from remaining 'delta volume'
                            remaining_delta_water_volume_FM += total_volume_about_to_be_removed
                            
                        # or, if volume to be removed would exceed 'delta volume', remove less than this instead
                        # (for this I identified 2 possible solutions [see notebook 'BMI_FM_PCR_coupled_run_2way_test_delta']
                        #  using 1 of those here)
                        else:
                            
                            # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                            remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                        FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                            
                            # remove water from all wet cells
                            new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                            
                            # check if there are negative water depths, and if so, repeat process
                            # -------------------------------------------------------------------
                            while any(new_water_depths_FM < 0):
                                
                                # reset negative water depths to zero and count the 'missing' volume
                                temp_vol_was_negative = 0
                                for k in range(len(new_water_depths_FM)):
                                    if new_water_depths_FM[k] < 0:
                                        temp_vol_was_negative -= (new_water_depths_FM[k] * -1) * \
                                                                    FMcellAreaSpherical[current_fm_cell_indices[k]]
                                        new_water_depths_FM[k] = 0.
                                        
                                # this represents the volume that is being 'added back' to the delta volume
                                remaining_delta_water_volume_FM = temp_vol_was_negative
                                
                                # find remaining wet cells
                                remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                                remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                                
                                # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                                remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                            FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                                    
                                # remove water from all wet cells
                                new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                            # -------------------------------------------------------------------
                            
                            # set 'delta volume' to zero (so while-loop will end)
                            remaining_delta_water_volume_FM = 0.
                
                    # adjusting array to be used with BMI 'set_var' to remove water from FM for negative delta volumes
                    new_water_levels_fm_1[current_fm_cell_indices] = bottom_lvl_fm[current_fm_cell_indices] + new_water_depths_FM
            
            # or if there are no FM river cells, use the same same procedure as in the original method
            else:
                
                # add the river (from PCR 'discharge') and floodplain (from PCR 'runoff') volumes to obain the total delta volume
                delta_volume_PCR_total_coupled = delta_volume_PCR_rivers_coupled[i] + delta_volume_PCR_floodplains_coupled[i]
                
                # check if the total delta volume is positive, and if so, add water to the relevant FM cells
                if delta_volume_PCR_total_coupled > 0:
                    
                    temp_water_volume_FM = delta_volume_PCR_total_coupled/len(CouplePCR2FM[i][1])
                    
                    for j in range(len(CouplePCR2FM[i][1])):
                        
                        # get current FM cell index
                        current_fm_cell_index = CouplePCR2FM[i][1][j]
                        # calculate additional water level for current FM cell
                        additional_water_level_FM[current_fm_cell_index] = temp_water_volume_FM / FMcellAreaSpherical[current_fm_cell_index]
                
                # or, if delta volume is negative, perform the new procedure which will remove water from FM instead of adding to it
                # (and if delta volume exactly equals zero, do nothing)
                elif delta_volume_PCR_total_coupled < 0:
                    
                    # get all coupled FM cells of the current coupled PCR cell
                    current_fm_cell_indices = CouplePCR2FM[i][1]
                    
                    # get current water depths of these cells
                    current_water_depths_FM = model_fm.get_var('s1')[current_fm_cell_indices] - bottom_lvl_fm[current_fm_cell_indices]
                    
                    # set variable for keeping track of water depths
                    new_water_depths_FM = np.copy(current_water_depths_FM)
                    
                    # remaining 'delta volume'
                    remaining_delta_water_volume_FM = delta_volume_PCR_total_coupled
                    
                    # while there is 'delta volume' remaining, remove water from wet cells
                    while remaining_delta_water_volume_FM < 0:
                        
                        # find remaining wet cells
                        remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                        remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                        
                        # if no remaining wet cells are found, no more water can be removed and the while-loop should stop
                        # this probably means that the remaining delta volume is very low (but still just above zero)
                        if len(remaining_wet_cells_FM_temp_indices) != 0:
                            pass
                        else:
                            if print_delta_volume_info == 1:
                                print 'No more wet cells found! Info:'
                                print 'PCR cell single index: ', CouplePCR2FM[i][0]
                                print 'Remaining delta volume:', remaining_delta_water_volume_FM
                            break
                        
                        # find cell with least amount of water
                        min_water_depth_FM_value = np.min(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                        min_water_depth_FM_index = np.argmin(new_water_depths_FM[remaining_wet_cells_FM_temp_indices])
                        
                        # calculate the total volume that is about to be removed from FM
                        total_volume_about_to_be_removed = np.sum(min_water_depth_FM_value * \
                                                                  FMcellAreaSpherical[remaining_wet_cells_FM_total_indices])
                            
                        # check if removing this from all wet cells will not exceed the remaining 'delta volume'
                        if (-1 * remaining_delta_water_volume_FM) > total_volume_about_to_be_removed:
                            
                            # remove water from all wet cells
                            new_water_depths_FM[remaining_wet_cells_FM_temp_indices] -= min_water_depth_FM_value
                            
                            # remove corresponding volume from remaining 'delta volume'
                            remaining_delta_water_volume_FM += total_volume_about_to_be_removed
                            
                        # or, if volume to be removed would exceed 'delta volume', remove less than this instead
                        # (for this I identified 2 possible solutions [see notebook 'BMI_FM_PCR_coupled_run_2way_test_delta']
                        #  using 1 of those here)
                        else:
                            
                            # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                            remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                        FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                            
                            # remove water from all wet cells
                            new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                            
                            # check if there are negative water depths, and if so, repeat process
                            # -------------------------------------------------------------------
                            while any(new_water_depths_FM < 0):
                                
                                # reset negative water depths to zero and count the 'missing' volume
                                temp_vol_was_negative = 0
                                for k in range(len(new_water_depths_FM)):
                                    if new_water_depths_FM[k] < 0:
                                        temp_vol_was_negative -= (new_water_depths_FM[k] * -1) * \
                                                                    FMcellAreaSpherical[current_fm_cell_indices[k]]
                                        new_water_depths_FM[k] = 0.
                                        
                                # this represents the volume that is being 'added back' to the delta volume
                                remaining_delta_water_volume_FM = temp_vol_was_negative
                                
                                # find remaining wet cells
                                remaining_wet_cells_FM_temp_indices  = np.where(new_water_depths_FM > 0.)[0]
                                remaining_wet_cells_FM_total_indices = current_fm_cell_indices[remaining_wet_cells_FM_temp_indices]
                                
                                # calculate additional water depths to be removed from all wet cells to reach remaining delta volume
                                remove_water_depths_extra = (remaining_delta_water_volume_FM/len(remaining_wet_cells_FM_temp_indices)) / \
                                                            FMcellAreaSpherical[remaining_wet_cells_FM_total_indices]
                                    
                                # remove water from all wet cells
                                new_water_depths_FM[remaining_wet_cells_FM_temp_indices] += remove_water_depths_extra
                            # -------------------------------------------------------------------
                            
                            # set 'delta volume' to zero (so while-loop will end)
                            remaining_delta_water_volume_FM = 0.
                
                    # adjusting array to be used with BMI 'set_var' to remove water from FM for negative delta volumes
                    new_water_levels_fm_1[current_fm_cell_indices] = bottom_lvl_fm[current_fm_cell_indices] + new_water_depths_FM
    
    # adjusting 's1' for negative delta volumes
    model_fm.set_var('s1', new_water_levels_fm_1)
    
    # calculate additional water levels for every FM update step (based on positive delta volumes)
    delta_water_levels_fm = additional_water_level_FM / fraction_timestep_fm
    # ---------------------------------------------------------------------------------------------------------------
    
    # UPDATE FM, ADDING WATER AT SPECIFIED UPDATE STEPS
    # ---------------------------------------------------------------------------------------------------------------
    # update FM unless it has has reached the same time as PCR
    while model_fm.get_current_time() < (model_pcr.get_time_step() * secPerDay):
        
        # print FM model time to screen, to check progress and speed of model
        if print_timesteps == 1:
            print 'current FM model time: ', model_fm.get_current_time()
        
        # get current FM water levels
        current_water_levels_fm = model_fm.get_var('s1')
        
        # calculate new water levels for fm
        
        # empty array for filling in values
        new_water_levels_fm_2 = np.zeros(len(current_water_levels_fm))
        
        # loop over cells
        for i in range(len(new_water_levels_fm_2)):
            
            # if at "known" cells:
            if i < len(delta_water_levels_fm):
                
                new_water_levels_fm_2[i] = current_water_levels_fm[i] + delta_water_levels_fm[i]
            
            # for cells outside the range of "known" cells:
            else:
                
                new_water_levels_fm_2[i] = current_water_levels_fm[i]
        
        # set new water levels in fm model
        model_fm.set_var('s1', new_water_levels_fm_2)
        
        # update FM
        model_fm.update(update_step_fm)
        
        # save water levels and water depths as output
        # only every full day
        if save_fm_lists == 1:
            if model_fm.get_current_time()%secPerDay == 0.0:
                # calculate water depths
                water_depths_fm = model_fm.get_var('s1') - bottom_lvl_fm
                # append to lists
                output_water_levels_fm.append(model_fm.get_var('s1'))
                output_water_depths_fm.append(water_depths_fm)
                print 'saving fm output at time:', model_fm.get_current_time()
    # ---------------------------------------------------------------------------------------------------------------
    
# at the end of the model run, save daily values of inundated area and volume within a specified region to individual files
with open(area_pcr_path_inundated_area, 'w') as f:
    for daily_value in area_pcr_total_inundated_area_floodplains:
        f.write(str(daily_value) + '\n')
        
with open(area_pcr_path_inundated_threshold, 'w') as f:
    for daily_value in area_pcr_total_inundated_area_threshold:
        f.write(str(daily_value) + '\n')
        
with open(area_pcr_path_flooded_volume, 'w') as f:
    for daily_value in area_pcr_total_flooded_volume_floodplains:
        f.write(str(daily_value) + '\n')

# end
print ''
print ''
print '>>>   MODEL RUN COMPLETED   <<<'
print ''
print ''

# create plots to see water levels
show_day_of_model_run = 2

if save_fm_lists == 1:
    f, ax = plt.subplots(1, 1)
    coupling_functions.plotGridValues(ax, FMcoords, output_water_depths_fm[show_day_of_model_run], model_name='FM', \
                                      limits=([0.001,10]), cell_cmap=my_cmap, extend_colorbar='both')
else:
    print 'no water levels and/or depths saved, cannot plot values!'
    print 'check output files instead!'

# check inundated area and volume
plt.plot(area_pcr_total_inundated_area_floodplains)
plt.show()
plt.plot(area_pcr_total_inundated_area_threshold)
plt.show()
plt.plot(area_pcr_total_flooded_volume_floodplains)
plt.show()


