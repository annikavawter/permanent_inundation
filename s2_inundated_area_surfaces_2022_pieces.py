from cgi import test
import arcpy
from arcpy import env
from arcpy.sa import *
import os
import glob
import numpy
from numpy import genfromtxt
arcpy.CheckOutExtension("Spatial")

arcpy.env.overwriteOutput = True

print ('STATUS 1: imports successful')

# ----------------------------------------------------------------------------------------------------
# SET REGIONS, PROJECTIONS, YEARS : toggle on/off
#     Regions: west_coast, east_coast, pr_usvi_coast, hi_coast, as_coast, gu_coast
#     Years: 2022, 2030, 2050, 2100 ...more? 2045
#     Scenarios: high, int_low ...more?

# Set vars to run
#RUN west_coast high_2022, high_2100, int_low_2022, int_low_2100
#RUN east_coast -
region = 'west_coast'
year = '2065'          #NOT RUN ['2035','2060','2080','2100']
projection = 'int_low'
flood_frequency = '26' #keep 26

print ('STATUS 2: vars set')
print ('Running region ' + str(region) + ' & year ' + str(year) + ' & projection ' + str(projection))

# ----------------------------------------------------------------------------------------------------
# 
gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
arcpy.env.workspace = gdb

# get the specific WLS from the gdb. name: water_level_surface_wrt_navd_26x_2022_high_west_coast
water_level_surface = arcpy.ListRasters('water_level_surface_wrt_navd_26x_{0}_{1}_{2}'.format(year,projection,region))
print (water_level_surface)

# ----------------------------------------------------------------------------------------------------
# west_dems  = ['CA_ChannelIslands_GCS_5m_NAVDm','CA_EKA1_GCS_5m_NAVD88m','CA_EKA2_GCS_5m_NAVD88m','CA_LOX1_GCS_5m_NAVD88m','CA_LOX2_GCS_5m_NAVD88m','CA_MTR1_GCS_5m_NAVD88m','CA_MTR2_GCS_5m_NAVD88m','CA_MTR3_GCS_5m_NAVD88m','CA_SGX_GCS_5m_NAVD88m','OR_MFR_GCS_5m_NAVD88m','OR_PQR1_GCS_5m_NAVD88m_TillamookFix','OR_PQR2_GCS_5m_NAVD88m','WA_PQR_GCS_5m_NAVD88m','WA_PugetSoundNW_GCS_3m_NAVDm','WA_PugetSoundSW_GCS_3m_NAVDm','WA_SEW1_GCS_5m_NAVD88m','WA_SEW1b_GCS_5m_NAVD88m','WA_SEW2_GCS_5m_NAVDm']
# dem = west_dems[17] #run 0 - 17 (18 total) #3-4min per run

# east_dems = ['AL_GCS_3m_NAVDm','CT_GCS_3m_NAVDm','DC_GCS_3m_NAVDm','DE_GCS_3m_NAVDm','FL_A3_GCS_5m_NAVDm','FL_JAX_1_GCS_5m_NAVD88m','FL_JAX_2_GCS_5m_NAVD88m','FL_MFL_1_GCS_5m_NAVD88m','FL_MFL_2_GCS_5m_NAVD88m','FL_MLB_1_GCS_5m_NAVD88m','FL_MLB_2_GCS_5m_NAVD88m','FL_Pan_East_GCS_3m_NAVDm','FL_Pan_West_GCS_3m_NAVDm','FL_TBW_1_GCS_5m_NAVD88m','FL_TBW_2_GCS_5m_NAVD88m','GA_CHS_GCS_5m_NAVDm','GA_JAX_GCS_5m_NAVDm','LA_Central_GCS_3m_NAVDm','LA_CentralEast_GCS_3m_NAVDm','LA_CentralNorth_GCS_3m_NAVDm','LA_Delta_GCS_3m_NAVDm','LA_LakePontchartrain_GCS_3m_NAVDm','LA_West_GCS_3m_NAVDm','MA_GCS_3m_NAVDm','MD_East_GCS_3m_NAVDm','MD_North_GCS_3m_NAVDm','MD_Southeast_GCS_3m_NAVDm','MD_Southwest_GCS_3m_NAVDm','MD_West_GCS_3m_NAVDm','ME_East_GCS_3m_NAVDm','ME_West_GCS_3m_NAVDm','MS_GCS_3m_NAVDm','NC_Middle1_GCS_3m_NAVDm','NC_Middle2_GCS_3m_NAVDm','NC_Northern_GCS_3m_NAVDm','NC_Southern1_GCS_3m_NAVDm','NC_Southern2_GCS_3m_NAVDm','NH_GCS_3m_NAVDm','NJ_Middle_GCS_3m_NAVDm_FY21Update','NJ_Northern_GCS_3m_NAVDm','NJ_Southern_GCS_3m_NAVDm_FY21Update','NY_Hudson_GCS_3m_NAVDm','NY_Metro_GCS_3m_NAVDm','NY_Suffolk_GCS_3m_NAVDm','PA_GCS_3m_NAVDm','RI_GCS_3m_NAVDm','SC_Central_GCS_3m_NAVDm','SC_North_GCS_3m_NAVDm','SC_South_GCS_3m_NAVDm','TX_Central_GCS_3m_NAVDm','TX_North1_GCS_3m_NAVDm','TX_North2_GCS_3m_NAVDm','TX_South1_GCS_3m_NAVDm','TX_South2_GCS_3m_NAVDm','VA_EasternShore_GCS_3m_NAVDm','VA_Middle_GCS_3m_NAVDm','VA_Northern_GCS_3m_NAVDm','VA_Southern_GCS_3m_NAVDm']
#east_dems_notevil = ['AL_GCS_3m_NAVDm','CT_GCS_3m_NAVDm','DC_GCS_3m_NAVDm','DE_GCS_3m_NAVDm','FL_A3_GCS_5m_NAVDm','FL_JAX_2_GCS_5m_NAVD88m','FL_MFL_2_GCS_5m_NAVD88m','FL_MLB_2_GCS_5m_NAVD88m','FL_Pan_East_GCS_3m_NAVDm','FL_Pan_West_GCS_3m_NAVDm','LA_Central_GCS_3m_NAVDm','LA_CentralEast_GCS_3m_NAVDm','LA_CentralNorth_GCS_3m_NAVDm','LA_Delta_GCS_3m_NAVDm','LA_LakePontchartrain_GCS_3m_NAVDm','LA_West_GCS_3m_NAVDm','MA_GCS_3m_NAVDm','MD_East_GCS_3m_NAVDm','MD_North_GCS_3m_NAVDm','MD_Southeast_GCS_3m_NAVDm','MD_Southwest_GCS_3m_NAVDm','MD_West_GCS_3m_NAVDm','ME_East_GCS_3m_NAVDm','ME_West_GCS_3m_NAVDm','MS_GCS_3m_NAVDm','NC_Middle1_GCS_3m_NAVDm','NC_Middle2_GCS_3m_NAVDm','NC_Northern_GCS_3m_NAVDm','NC_Southern1_GCS_3m_NAVDm','NC_Southern2_GCS_3m_NAVDm','NH_GCS_3m_NAVDm','NJ_Middle_GCS_3m_NAVDm_FY21Update','NJ_Northern_GCS_3m_NAVDm','NJ_Southern_GCS_3m_NAVDm_FY21Update','NY_Hudson_GCS_3m_NAVDm','NY_Metro_GCS_3m_NAVDm','NY_Suffolk_GCS_3m_NAVDm','PA_GCS_3m_NAVDm','RI_GCS_3m_NAVDm','SC_Central_GCS_3m_NAVDm','SC_North_GCS_3m_NAVDm','SC_South_GCS_3m_NAVDm','TX_Central_GCS_3m_NAVDm','TX_North1_GCS_3m_NAVDm','TX_North2_GCS_3m_NAVDm','TX_South1_GCS_3m_NAVDm','TX_South2_GCS_3m_NAVDm','VA_EasternShore_GCS_3m_NAVDm','VA_Middle_GCS_3m_NAVDm','VA_Northern_GCS_3m_NAVDm','VA_Southern_GCS_3m_NAVDm',]
# east_dems_muchevil = ['FL_TBW_1_GCS_5m_NAVD88m','FL_TBW_2_GCS_5m_NAVD88m','FL_MFL_1_GCS_5m_NAVD88m','FL_MLB_1_GCS_5m_NAVD88m','FL_JAX_1_GCS_5m_NAVD88m','GA_JAX_GCS_5m_NAVDm','GA_CHS_GCS_5m_NAVDm']

# dem = east_dems[11] #run 0 - 57. (58 total)  #runtime
#     #need 12 - end for east coast high 2100 except the special ones in ladems and njdems

# dem = east_dems_muchevil[0] #run 0 - ?. ( total)  #runtime

# -------------
ladems = ['LA_Central_GCS_3m_NAVDm','LA_CentralEast_GCS_3m_NAVDm','LA_CentralNorth_GCS_3m_NAVDm','LA_Delta_GCS_3m_NAVDm','LA_LakePontchartrain_GCS_3m_NAVDm','LA_West_GCS_3m_NAVDm']
njdems = ['NJ_Middle_GCS_3m_NAVDm_FY21Update','NJ_Northern_GCS_3m_NAVDm','NJ_Southern_GCS_3m_NAVDm_FY21Update']
dem =  'CA_SGX_GCS_5m_NAVD88m'  #njdems[0]   #ladems[3]

# ----------------------------------------------------------------------------------------------------
# set mask to dem extent
arcpy.env.mask = dem
arcpy.env.cellSize = "MINOF"

# Compare DEM and water level surface. Where WLS >= DEM, set value of 1 ; i.e. water above land is TRUE
#print ('DEM is: ' + str(dem))
#print ('Creating inundated area surface')
# water_level_surface_setnulls = Con(Raster(water_level_surface) > (-100), Raster(water_level_surface)) # sets -9999 vals as nodata
# inundated_area_surface = Con(Raster(water_level_surface_setnulls) >= Raster(dem), 1) # creates flat inundation area raster
inundated_area_surface = Con(Raster(water_level_surface) >= Raster(dem), 1) # creates flat inundation area raster

print ('Created inundated area surface for ' + str(dem))
outname_inundated_area_surface = 'inundated_area_surface_%sx_%s_%s_%s' %(flood_frequency, year, projection, region) + '_' + dem[:20]  #previously dem[:-10]
print('New layer is: ' + str(outname_inundated_area_surface))
inundated_area_surface.save(outname_inundated_area_surface)

#TURN OFF DEPTH
#inundated_area_surface_depth = Con(Raster(water_level_surface)* .3048 >= Raster(dem), Raster(water_level_surface) - Raster(dem)) # creates depth raster
#outname_inundated_area_surface_depth = 'depth_inundated_area_surface_%sx_%s_%s_' %(flood_frequency, year, projection) + '_' + dem[:-10]
#inundated_area_surface_depth.save(outname_inundated_area_surface_depth)

# inundated_area_surfaces_raw.append(outname_inundated_area_surface)

# return inundated_area_surfaces_raw
# print ('STATUS 2: functions set up successfully')

print ('STATUS 3: inundation surface run successfully')
print('end')