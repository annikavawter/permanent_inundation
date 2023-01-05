print('FIXING DEMS')

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

# Set vars to run
regions = ['east_coast','west_coast','pr_coast','vi_coast','gu_coast','hi_coast','as_coast']
regions = ['pr_coast','vi_coast','gu_coast','hi_coast','as_coast']
regions = ['east_coast','west_coast']
# regions = ['east_coast'] #test

print ('STATUS 2: vars set')

# ----------------------------------------------------------------------------------------------------
# SET GDB

gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
arcpy.env.workspace = gdb

# ----------------------------------------------------------------------------------------------------
# SET DEMS

west_dems  = ['CA_ChannelIslands_GCS_5m_NAVDm','CA_EKA1_GCS_5m_NAVD88m','CA_EKA2_GCS_5m_NAVD88m','CA_LOX1_GCS_5m_NAVD88m','CA_LOX2_GCS_5m_NAVD88m','CA_MTR1_GCS_5m_NAVD88m','CA_MTR2_GCS_5m_NAVD88m','CA_MTR3_GCS_5m_NAVD88m','CA_SGX_GCS_5m_NAVD88m','OR_MFR_GCS_5m_NAVD88m','OR_PQR1_GCS_5m_NAVD88m_TillamookFix','OR_PQR2_GCS_5m_NAVD88m','WA_PQR_GCS_5m_NAVD88m','WA_PugetSoundNW_GCS_3m_NAVDm','WA_PugetSoundSW_GCS_3m_NAVDm','WA_SEW1_GCS_5m_NAVD88m','WA_SEW1b_GCS_5m_NAVD88m','WA_SEW2_GCS_5m_NAVDm']
east_dems = ['AL_GCS_3m_NAVDm','CT_GCS_3m_NAVDm','DC_GCS_3m_NAVDm','DE_GCS_3m_NAVDm','FL_A3_GCS_5m_NAVDm','FL_JAX_1_GCS_5m_NAVD88m','FL_JAX_2_GCS_5m_NAVD88m','FL_MFL_1_GCS_5m_NAVD88m','FL_MFL_2_GCS_5m_NAVD88m','FL_MLB_1_GCS_5m_NAVD88m','FL_MLB_2_GCS_5m_NAVD88m','FL_Pan_East_GCS_3m_NAVDm','FL_Pan_West_GCS_3m_NAVDm','FL_TBW_1_GCS_5m_NAVD88m','FL_TBW_2_GCS_5m_NAVD88m','GA_CHS_GCS_5m_NAVDm','GA_JAX_GCS_5m_NAVDm','LA_Central_GCS_3m_NAVDm','LA_CentralEast_GCS_3m_NAVDm','LA_CentralNorth_GCS_3m_NAVDm','LA_Delta_GCS_3m_NAVDm','LA_LakePontchartrain_GCS_3m_NAVDm','LA_West_GCS_3m_NAVDm','MA_GCS_3m_NAVDm','MD_East_GCS_3m_NAVDm','MD_North_GCS_3m_NAVDm','MD_Southeast_GCS_3m_NAVDm','MD_Southwest_GCS_3m_NAVDm','MD_West_GCS_3m_NAVDm','ME_East_GCS_3m_NAVDm','ME_West_GCS_3m_NAVDm','MS_GCS_3m_NAVDm','NC_Middle1_GCS_3m_NAVDm','NC_Middle2_GCS_3m_NAVDm','NC_Northern_GCS_3m_NAVDm','NC_Southern1_GCS_3m_NAVDm','NC_Southern2_GCS_3m_NAVDm','NH_GCS_3m_NAVDm','NJ_Middle_GCS_3m_NAVDm_FY21Update','NJ_Northern_GCS_3m_NAVDm','NJ_Southern_GCS_3m_NAVDm_FY21Update','NY_Hudson_GCS_3m_NAVDm','NY_Metro_GCS_3m_NAVDm','NY_Suffolk_GCS_3m_NAVDm','PA_GCS_3m_NAVDm','RI_GCS_3m_NAVDm','SC_Central_GCS_3m_NAVDm','SC_North_GCS_3m_NAVDm','SC_South_GCS_3m_NAVDm','TX_Central_GCS_3m_NAVDm','TX_North1_GCS_3m_NAVDm','TX_North2_GCS_3m_NAVDm','TX_South1_GCS_3m_NAVDm','TX_South2_GCS_3m_NAVDm','VA_EasternShore_GCS_3m_NAVDm','VA_Middle_GCS_3m_NAVDm','VA_Northern_GCS_3m_NAVDm','VA_Southern_GCS_3m_NAVDm']
# east_dems_notevil = ['AL_GCS_3m_NAVDm','CT_GCS_3m_NAVDm','DC_GCS_3m_NAVDm','DE_GCS_3m_NAVDm','FL_A3_GCS_5m_NAVDm','FL_JAX_2_GCS_5m_NAVD88m','FL_MFL_2_GCS_5m_NAVD88m','FL_MLB_2_GCS_5m_NAVD88m','FL_Pan_East_GCS_3m_NAVDm','FL_Pan_West_GCS_3m_NAVDm','LA_Central_GCS_3m_NAVDm','LA_CentralEast_GCS_3m_NAVDm','LA_CentralNorth_GCS_3m_NAVDm','LA_Delta_GCS_3m_NAVDm','LA_LakePontchartrain_GCS_3m_NAVDm','LA_West_GCS_3m_NAVDm','MA_GCS_3m_NAVDm','MD_East_GCS_3m_NAVDm','MD_North_GCS_3m_NAVDm','MD_Southeast_GCS_3m_NAVDm','MD_Southwest_GCS_3m_NAVDm','MD_West_GCS_3m_NAVDm','ME_East_GCS_3m_NAVDm','ME_West_GCS_3m_NAVDm','MS_GCS_3m_NAVDm','NC_Middle1_GCS_3m_NAVDm','NC_Middle2_GCS_3m_NAVDm','NC_Northern_GCS_3m_NAVDm','NC_Southern1_GCS_3m_NAVDm','NC_Southern2_GCS_3m_NAVDm','NH_GCS_3m_NAVDm','NJ_Middle_GCS_3m_NAVDm_FY21Update','NJ_Northern_GCS_3m_NAVDm','NJ_Southern_GCS_3m_NAVDm_FY21Update','NY_Hudson_GCS_3m_NAVDm','NY_Metro_GCS_3m_NAVDm','NY_Suffolk_GCS_3m_NAVDm','PA_GCS_3m_NAVDm','RI_GCS_3m_NAVDm','SC_Central_GCS_3m_NAVDm','SC_North_GCS_3m_NAVDm','SC_South_GCS_3m_NAVDm','TX_Central_GCS_3m_NAVDm','TX_North1_GCS_3m_NAVDm','TX_North2_GCS_3m_NAVDm','TX_South1_GCS_3m_NAVDm','TX_South2_GCS_3m_NAVDm','VA_EasternShore_GCS_3m_NAVDm','VA_Middle_GCS_3m_NAVDm','VA_Northern_GCS_3m_NAVDm','VA_Southern_GCS_3m_NAVDm',]
# east_dems_muchevil = ['FL_TBW_1_GCS_5m_NAVD88m','FL_TBW_2_GCS_5m_NAVD88m','FL_MFL_1_GCS_5m_NAVD88m','FL_MLB_1_GCS_5m_NAVD88m','FL_JAX_1_GCS_5m_NAVD88m','GA_JAX_GCS_5m_NAVDm','GA_CHS_GCS_5m_NAVDm']
hi_dems = ['HI_HFO_Hawaii_UTM_3m_LMSLm','HI_HFO_Kauai_UTM_3m_LMSLm','HI_HFO_Lanai_UTM_3m_LMSLm','HI_HFO_Maui_UTM_3m_LMSLm','HI_HFO_Molokai_UTM_3m_LMSLm','HI_Oahu_GCS_3m_LMSLm']
gu_dems = ['Guam_GCS_3m_GUVD04m']
as_dems = ['AS_OfuOlosega_DEM_v2_UTM','AS_Tau_DEM_v2_UTM','AS_Tutuila_DEM_v2_UTM']
pr_dems = ['PR_GCS_3m_PRVD02m']
vi_dems = ['USVI_GCS_3m_NAVDm']

# ----------------------------------------------------------------------------------------------------
# FIX ERROR IN DEMS

for region in regions:
    #set save location
    outfolder_str_dem = 'C:/Users/annik/Documents/UCS/Data/2022/DEMs_fix/{0}/{0}_dem_fix.gdb' .format(region)

    print('Region is: ' + str(region))

    if region == 'west_coast':
        dems = west_dems #18 total
    elif region == 'east_coast':
        dems = east_dems #58 total
    elif region == 'hi_coast':
        dems = hi_dems 
    elif region == 'gu_coast':
        dems = gu_dems
    elif region == 'as_coast':
        dems = as_dems
    elif region == 'pr_coast':
        dems = pr_dems
    elif region == 'vi_coast':
        dems = vi_dems

    for dem in dems:
        # set mask to dem extent
        arcpy.env.mask = dem
        arcpy.env.cellSize = "MINOF"

        print ('DEM is: ' + str(dem))

        # demObj = arcpy.Raster(dem)
        # print('properties:')
        # nodata = arcpy.GetRasterProperties_management(demObj,"ANYNODATA").getOutput(0) #error

        #fix error in dems (-9999 or -99 vs NoData) and save in separate regional gdb
        print('...calculating...')
        # dem_fixnulls = Con(Raster(dem) > (-50), (dem)) # sets -9999 and -99 vals as nodata
        dem_fixnulls = Con(Raster(dem) > (-1000), (dem)) # sets -9999 as nodata

        # outname_dem_fixnulls = str(outfolder_str_dem) + '/' + str(dem) + '_CONNECT'
        outname_dem_fixnulls = str(outfolder_str_dem) + '/' + str(dem)
        #print(outname_dem_fixnulls)

        print('...saving...')
        dem_fixnulls.save(outname_dem_fixnulls)
        print('New DEM created: ' + str(outname_dem_fixnulls))
        

print ('STATUS 3: DEMs fixed successfully')
print('end')