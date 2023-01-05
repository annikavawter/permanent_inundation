print('RUN INUNDATED SURFACES')

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
#     Years: 2022, 2030, 2035, 2045, 2050, 2060, 2065, 2080, 2100 ...more?
#     Scenarios: high, int_low ...more?

# Set vars to run
# region = 'east_coast'
# regions = ['west_coast','east_coast']
# regions = ['east_coast']
regions = ['west_coast']
years = ['2045','2060','2065','2080']     #['2022','2100','2060','2045','2035','2080','2050','2030','2065', ...]
projection = 'high'  #['high','int_low']
flood_frequency = '26' #keep 26

print ('STATUS 2: vars set')
# print ('Running region ' + str(region) + ' & year ' + str(year) + ' & projection ' + str(projection))

# ----------------------------------------------------------------------------------------------------
# SET WORKSPACE

gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
arcpy.env.workspace = gdb

# ----------------------------------------------------------------------------------------------------
# SET UP REGIONAL DEMS

#dem gdb: C:\Users\annik\Documents\UCS\Data\2022\DEMs_regional\west_coast\west_coast_dem.gdb
# dem_gdb = 'C:/Users/annik/Documents/UCS/Data/2022/DEMs_regional/west_coast/west_coast_dem.gdb'
# dem_gdb = 'C:/Users/annik/Documents/UCS/Data/2022/DEMs_fix/west_coast/west_coast_dem_fix.gdb'

west_dems = ['CA_ChannelIslands_GCS_5m_NAVDm','CA_EKA1_GCS_5m_NAVD88m','CA_EKA2_GCS_5m_NAVD88m','CA_LOX1_GCS_5m_NAVD88m','CA_LOX2_GCS_5m_NAVD88m','CA_MTR1_GCS_5m_NAVD88m','CA_MTR2_GCS_5m_NAVD88m','CA_MTR3_GCS_5m_NAVD88m','CA_SGX_GCS_5m_NAVD88m','OR_MFR_GCS_5m_NAVD88m','OR_PQR1_GCS_5m_NAVD88m_TillamookFix','OR_PQR2_GCS_5m_NAVD88m','WA_PQR_GCS_5m_NAVD88m','WA_PugetSoundNW_GCS_3m_NAVDm','WA_PugetSoundSW_GCS_3m_NAVDm','WA_SEW1_GCS_5m_NAVD88m','WA_SEW1b_GCS_5m_NAVD88m','WA_SEW2_GCS_5m_NAVDm']
east_dems = ['AL_GCS_3m_NAVDm','CT_GCS_3m_NAVDm','DC_GCS_3m_NAVDm','DE_GCS_3m_NAVDm','FL_A3_GCS_5m_NAVDm','FL_JAX_1_GCS_5m_NAVD88m','FL_JAX_2_GCS_5m_NAVD88m','FL_MFL_1_GCS_5m_NAVD88m','FL_MFL_2_GCS_5m_NAVD88m','FL_MLB_1_GCS_5m_NAVD88m','FL_MLB_2_GCS_5m_NAVD88m','FL_Pan_East_GCS_3m_NAVDm','FL_Pan_West_GCS_3m_NAVDm','FL_TBW_1_GCS_5m_NAVD88m','FL_TBW_2_GCS_5m_NAVD88m','GA_CHS_GCS_5m_NAVDm','GA_JAX_GCS_5m_NAVDm','LA_Central_GCS_3m_NAVDm','LA_CentralEast_GCS_3m_NAVDm','LA_CentralNorth_GCS_3m_NAVDm','LA_Delta_GCS_3m_NAVDm','LA_LakePontchartrain_GCS_3m_NAVDm','LA_West_GCS_3m_NAVDm','MA_GCS_3m_NAVDm','MD_East_GCS_3m_NAVDm','MD_North_GCS_3m_NAVDm','MD_Southeast_GCS_3m_NAVDm','MD_Southwest_GCS_3m_NAVDm','MD_West_GCS_3m_NAVDm','ME_East_GCS_3m_NAVDm','ME_West_GCS_3m_NAVDm','MS_GCS_3m_NAVDm','NC_Middle1_GCS_3m_NAVDm','NC_Middle2_GCS_3m_NAVDm','NC_Northern_GCS_3m_NAVDm','NC_Southern1_GCS_3m_NAVDm','NC_Southern2_GCS_3m_NAVDm','NH_GCS_3m_NAVDm','NJ_Middle_GCS_3m_NAVDm_FY21Update','NJ_Northern_GCS_3m_NAVDm','NJ_Southern_GCS_3m_NAVDm_FY21Update','NY_Hudson_GCS_3m_NAVDm','NY_Metro_GCS_3m_NAVDm','NY_Suffolk_GCS_3m_NAVDm','PA_GCS_3m_NAVDm','RI_GCS_3m_NAVDm','SC_Central_GCS_3m_NAVDm','SC_North_GCS_3m_NAVDm','SC_South_GCS_3m_NAVDm','TX_Central_GCS_3m_NAVDm','TX_North1_GCS_3m_NAVDm','TX_North2_GCS_3m_NAVDm','TX_South1_GCS_3m_NAVDm','TX_South2_GCS_3m_NAVDm','VA_EasternShore_GCS_3m_NAVDm','VA_Middle_GCS_3m_NAVDm','VA_Northern_GCS_3m_NAVDm','VA_Southern_GCS_3m_NAVDm']
hi_dems = ['HI_HFO_Hawaii_UTM_3m_LMSLm','HI_HFO_Kauai_UTM_3m_LMSLm','HI_HFO_Lanai_UTM_3m_LMSLm','HI_HFO_Maui_UTM_3m_LMSLm','HI_HFO_Molokai_UTM_3m_LMSLm','HI_Oahu_GCS_3m_LMSLm']
gu_dems = ['Guam_GCS_3m_GUVD04m']
as_dems = ['AS_OfuOlosega_DEM_v2_UTM','AS_Tau_DEM_v2_UTM','AS_Tutuila_DEM_v2_UTM']
pr_dems = ['PR_GCS_3m_PRVD02m']
vi_dems = ['USVI_GCS_3m_NAVDm']

# ----------------------------------------------------------------------------------------------------
# TESTING
# west_dems  = ['CA_ChannelIslands_GCS_5m_NAVDm'] #testing
# west_dems  = ['CA_SGX_GCS_5m_NAVD88m'] #testing
# east_dems  = ['GA_CHS_GCS_5m_NAVDm'] #testing, corrected dem
# west_dems = ['WA_PugetSoundSW_GCS_3m_NAVDm','WA_SEW1_GCS_5m_NAVD88m','WA_SEW1b_GCS_5m_NAVD88m','WA_SEW2_GCS_5m_NAVDm']

# ----------------------------------------------------------------------------------------------------
# COMPARE FIXED DEMS TO WLS & OUTPUT RASTER  

for region in regions:
    for year in years:
        print ('Running region ' + str(region) + ' & year ' + str(year) + ' & projection ' + str(projection))
        # print('Region is: ' + str(region))

        # ----------------------------------------------------------------
        # get the specific WLS from the gdb
        #name: water_level_surface_wrt_navd_26x_2022_high_west_coast
        water_level_surface = arcpy.ListRasters('water_level_surface_wrt_navd_26x_{0}_{1}_{2}'.format(year,projection,region))
        print (water_level_surface)

        # ----------------------------------------------------------------
        # set save location
        # outfolder_str = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
        # print (outfolder_str)

        # ----------------------------------------------------------------
        # get dems from corrected regional dems

        # gdb_dem = 'C:/Users/annik/Documents/UCS/Data/2022/DEMs_regional/{0}/{0}_dem.gdb' #previous
        gdb_dem = 'C:/Users/annik/Documents/UCS/Data/2022/DEMs_fix/{0}/{0}_dem_fix.gdb' .format(region)

        if region == 'west_coast':
            dems = west_dems #18 total
        elif region == 'east_coast':
            dems = east_dems  #58 total
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

        # for dem in dems:
        #     # print (dems)
        #     dem_full = str(gdb_dem) + '/' + str(dem)
        #     print(dem_full)

        # ----------------------------------------------------------------
        # RUN DEMS LOOP

        for dem in dems:

            #set path for pulling in corrected dems
            dem_path = str(gdb_dem) + '/' + str(dem)

            #set path for saving inundated raster
            gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

            # set mask to dem extent
            arcpy.env.mask = dem_path
            arcpy.env.cellSize = "MINOF"

            # Compare DEM and water level surface. Where WLS >= DEM, set value of 1 ; i.e. water above land is TRUE
            print ('DEM is: ' + str(dem))
            #print ('Creating inundated area surface')
            inundated_area_surface = Con(Raster(water_level_surface) >= Raster(dem_path), 1) # creates flat inundation area raster
            
            print ('Created inundated area surface for ' + str(dem))
            outname_inundated_area_surface = 'inundated_area_surface_%sx_%s_%s_%s' %(flood_frequency, year, projection, region) + '_' + dem[:20]  #previously dem[:-10]
            outname_inundated_area_surface_full = str(gdb_out) + '/' + str(outname_inundated_area_surface)

            print('New layer is: ' + str(outname_inundated_area_surface))
            #print('Full save path is: ' + str(outname_inundated_area_surface_full))

            # save raster
            inundated_area_surface.save(outname_inundated_area_surface_full)

        print ('STATUS 3: inundation surface run successfully for ' + str(projection) + ' ' + str(year))

print('end')