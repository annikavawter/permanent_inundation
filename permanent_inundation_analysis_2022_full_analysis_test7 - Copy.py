print('SCRIPT ALL INUNDATION THRU POLYGONS')

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

print ('STATUS: imports successful')

import datetime
t1 = datetime.datetime.now() ; #print('Time start: ' + str(t1))

# ----------------------------------------------------------------------------------------------------
# # SET INPUTS
flood_frequency = '26'
## regions = ['west_coast']     #['east_coast','west_coast','pr_coast','vi_coast','gu_coast','hi_coast','as_coast']
# years = ['2060']     #['2022','2100','2060','2045','2035','2080','2050','2030','2065', ...]
# projections = ['int_low']  #['high','int_low']

# ----------------------------------------------------------------------------------------------------
# SET WORKSPACE
gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
arcpy.env.workspace = gdb

# Set off 0 | on 1
runscratch=1
runscratchout=1
showtimes=0

# ----------------------------------------------------------------------------------------------------
# DEM ERROR FIX
def dem_fix(regions):
    print('RUNNING FUNCTION FIXING DEM ERRORS')
    print('Setting all -9999 values in DEMS to NoData and resaving in new regional geodatabases')

    # set project gdb
    gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'

    # set gdb to pull dems from, and set as workspace
    # gdb_dem_orig = 'C:/Users/annik/Documents/UCS/Data/NOAA_DEM.gdb'
    # arcpy.env.workspace = gdb_dem_orig

    # SET DEMS: assigns original DEM file names to regions to be pulled from project gdb (done manually due to scripting issues)
    west_dems = ['CA_ChannelIslands_GCS_5m_NAVDm','CA_EKA1_GCS_5m_NAVD88m','CA_EKA2_GCS_5m_NAVD88m','CA_LOX1_GCS_5m_NAVD88m','CA_LOX2_GCS_5m_NAVD88m','CA_MTR1_GCS_5m_NAVD88m','CA_MTR2_GCS_5m_NAVD88m','CA_MTR3_GCS_5m_NAVD88m','CA_SGX_GCS_5m_NAVD88m','OR_MFR_GCS_5m_NAVD88m','OR_PQR1_GCS_5m_NAVD88m_TillamookFix','OR_PQR2_GCS_5m_NAVD88m','WA_PQR_GCS_5m_NAVD88m','WA_PugetSoundNW_GCS_3m_NAVDm','WA_PugetSoundSW_GCS_3m_NAVDm','WA_SEW1_GCS_5m_NAVD88m','WA_SEW1b_GCS_5m_NAVD88m','WA_SEW2_GCS_5m_NAVDm']
    east_dems = ['AL_GCS_3m_NAVDm','CT_GCS_3m_NAVDm','DC_GCS_3m_NAVDm','DE_GCS_3m_NAVDm','FL_A3_GCS_5m_NAVDm','FL_JAX_1_GCS_5m_NAVD88m','FL_JAX_2_GCS_5m_NAVD88m','FL_MFL_1_GCS_5m_NAVD88m','FL_MFL_2_GCS_5m_NAVD88m','FL_MLB_1_GCS_5m_NAVD88m','FL_MLB_2_GCS_5m_NAVD88m','FL_Pan_East_GCS_3m_NAVDm','FL_Pan_West_GCS_3m_NAVDm','FL_TBW_1_GCS_5m_NAVD88m','FL_TBW_2_GCS_5m_NAVD88m','GA_CHS_GCS_5m_NAVDm','GA_JAX_GCS_5m_NAVDm','LA_Central_GCS_3m_NAVDm','LA_CentralEast_GCS_3m_NAVDm','LA_CentralNorth_GCS_3m_NAVDm','LA_Delta_GCS_3m_NAVDm','LA_LakePontchartrain_GCS_3m_NAVDm','LA_West_GCS_3m_NAVDm','MA_GCS_3m_NAVDm','MD_East_GCS_3m_NAVDm','MD_North_GCS_3m_NAVDm','MD_Southeast_GCS_3m_NAVDm','MD_Southwest_GCS_3m_NAVDm','MD_West_GCS_3m_NAVDm','ME_East_GCS_3m_NAVDm','ME_West_GCS_3m_NAVDm','MS_GCS_3m_NAVDm','NC_Middle1_GCS_3m_NAVDm','NC_Middle2_GCS_3m_NAVDm','NC_Northern_GCS_3m_NAVDm','NC_Southern1_GCS_3m_NAVDm','NC_Southern2_GCS_3m_NAVDm','NH_GCS_3m_NAVDm','NJ_Middle_GCS_3m_NAVDm_FY21Update','NJ_Northern_GCS_3m_NAVDm','NJ_Southern_GCS_3m_NAVDm_FY21Update','NY_Hudson_GCS_3m_NAVDm','NY_Metro_GCS_3m_NAVDm','NY_Suffolk_GCS_3m_NAVDm','PA_GCS_3m_NAVDm','RI_GCS_3m_NAVDm','SC_Central_GCS_3m_NAVDm','SC_North_GCS_3m_NAVDm','SC_South_GCS_3m_NAVDm','TX_Central_GCS_3m_NAVDm','TX_North1_GCS_3m_NAVDm','TX_North2_GCS_3m_NAVDm','TX_South1_GCS_3m_NAVDm','TX_South2_GCS_3m_NAVDm','VA_EasternShore_GCS_3m_NAVDm','VA_Middle_GCS_3m_NAVDm','VA_Northern_GCS_3m_NAVDm','VA_Southern_GCS_3m_NAVDm']
    hi_dems = ['HI_HFO_Hawaii_UTM_3m_LMSLm','HI_HFO_Kauai_UTM_3m_LMSLm','HI_HFO_Lanai_UTM_3m_LMSLm','HI_HFO_Maui_UTM_3m_LMSLm','HI_HFO_Molokai_UTM_3m_LMSLm','HI_Oahu_GCS_3m_LMSLm']
    gu_dems = ['Guam_GCS_3m_GUVD04m']
    as_dems = ['AS_OfuOlosega_DEM_v2_UTM','AS_Tau_DEM_v2_UTM','AS_Tutuila_DEM_v2_UTM']
    pr_dems = ['PR_GCS_3m_PRVD02m']
    vi_dems = ['USVI_GCS_3m_NAVDm']

    for region in regions:
        print('Region is: ' + str(region))

        # set save location
        outfolder_str_dem = 'C:/Users/annik/Documents/UCS/Data/2022/DEMs_fix/{0}/{0}_dem_fix.gdb' .format(region)

        #assign dem set based on region
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

            # check if dems already fixed and saved; create if not
            arcpy.env.workspace = outfolder_str_dem
            check_demfix = arcpy.ListRasters(dem) ; #print(check_demfix) ; print(len(check_demfix))
            arcpy.env.workspace = gdb

            if len(check_demfix)!=0:
                print('DEM ' + str(dem) + ' already fixed')

            elif len(check_demfix)==0:
                arcpy.env.mask = dem  #set mask to dem extent
                arcpy.env.cellSize = "MINOF"

                print ('DEM is: ' + str(dem))

                # fix error in dems (-9999 or -99 vs NoData) and save in separate regional gdb
                print('...calculating...')
                dem_fixnulls = Con(Raster(dem) > (-1000), (dem)) # sets -9999 as nodata; leaves -99 water values in rasters

                outname_dem_fixnulls = str(outfolder_str_dem) + '/' + str(dem)

                print('...saving raster...')
                dem_fixnulls.save(outname_dem_fixnulls)
                print('New DEM created: ' + str(outname_dem_fixnulls))
            
        print ('All DEMs fixed for ' + region)
    print('STATUS: DEMs fixed successfully')

def prep_noaa_mhhw_layers (regions):  ## THE MAKING PART HAS NOT BEEN TESTED< MAKE NEW GDB OUT AND RUN
    for region in regions:
        print('RUNNING PREP MHHW LAYERS FUNCTION')
        print ('...running ' + str(region) + '...')

        # ----------------------------------------------------------------------------------------
        # set gdb in and out
        gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
        # set path for loading mhhw data
        gdb_mhhw = 'C:/Users/annik/Documents/UCS/Data/MHHW/NOAA_OCM_MHHW.gdb'
        # set path for saving
        gdb_out = gdb

        # set workspace
        arcpy.env.workspace = gdb

        if region=='east_coast' or region=='west_coast':
            exit()

            # ----------------------------------------------------------------------------------------
            # mosaic mhhw rasters if not already done
            
            outname_mhhw_mos = 'all_noaa_mhhw_navd_mosaic'

            mhhw_mos_check = arcpy.ListRasters(outname_mhhw_mos)
            # print(mhhw_mos_check)
            # print(str(len(mhhw_mos_check)))
            
            if len(mhhw_mos_check)>0:
                print('NOAA MHHW data already prepped')

            else:                 
                # grab all noaa mhhw rasters from downloaded gdb
                arcpy.env.workspace = gdb_mhhw
                mhhw_rasters = arcpy.ListRasters('*')

                # grab regional rasters, assigned by state, for masking later  ##EDIT THIS
                mhhw_rasters_regional = []
                for raster in mhhw_rasters:               
                    #set states in regions from mhhw file names
                    if region=='east_coast':
                        #state_set = ['TX','LA','MS','AL','FL','GA','SC','NC','VA','DC','MD','DE','PA','NJ','NY','CT','RI','MA','NH','ME'] #east
                        state_set = ['TX','LA','MSAL','FL','GA','SC','NC','MDVA','NJDEPA','NYCT','NewEngland'] #east
                    elif region=='west_coast':
                        state_set = ['WA','OR','CA'] #west
                    elif region=='pr_coast' or region=='vi_coast':
                        state_set = ['PR_USVI'] #pr

                    #name format: CA_MHHW_UTM_100m_NAVDm
                    mhhw_rasters_regional = arcpy.ListRasters('{0}_MHHW_UTM_*_NAVDm*' .format(state_set))  #added * at end due to manual name adjustment


                arcpy.env.workspace = gdb

                mhhw_paths = []
                print('Rasters to mosaic:')
                for raster in mhhw_rasters:
                    print(str(raster))
                    mhhw_paths.append(str(gdb_mhhw) + '/' + str(raster))

                print('...creating mosaic...')
                arcpy.MosaicToNewRaster_management(mhhw_paths, gdb_out, outname_mhhw_mos,"","32_BIT_FLOAT","",1,"MAXIMUM","FIRST")   #auto cellsize; too large and tool won't run. NOT cellsz_pt000045
                print('Created NOAA MHHW mosaic')

                #from Pro tool:
                #arcpy.management.MosaicToNewRaster(r"C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\CA_MHHW_UTM_100m_NAVDm_west_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\FL_MHHW_UTM_50m_NAVDm_east_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\GA_MHHW_UTM_100m_NAVDm_east_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\LA_MHHW_UTM_100m_NAVDm_east_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\MDVA_MHHW_UTM_100m_NAVDm_east_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\MSAL_MHHW_UTM_100m_NAVDm_east_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\NC_MHHW_UTM_50m_NAVDm_east_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\NewEngland_MHHW_UTM_100m_NAVDm_east_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\NJDEPA_MHHW_UTM_100m_NAVDm_east_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\NYCT_MHHW_UTM_100m_NAVDm_east_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\OR_MHHW_UTM_100m_NAVDm_west_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\PR_USVI_MHHW_UTM_50m_NAVDm_prvi_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\SC_MHHW_UTM_100m_NAVDm_east_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\TX_MHHW_UTM_100m_NAVDm_east_coast;C:\Users\annik\Documents\UCS\Data\MHHW\NOAA_OCM_MHHW.gdb\WA_MHHW_UTM_100m_NAVDm_west_coast", r"C:\Users\annik\Documents\ArcGIS\Projects\UCS_Flood_2022\UCS_Flood_2022.gdb", "all_noaa_mhhw_navd_mosaic", 'GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]', "32_BIT_FLOAT", None, 1, "MAXIMUM", "FIRST")

            # ----------------------------------------------------------------------------------------
            #  create mhhw regional polygons if not already done

            arcpy.env.workspace = gdb

            outname_mhhw_poly = 'all_noaa_mhhw_navd_mosaic_polygon_multi'
            outname_mhhw_poly_reg = 'all_noaa_mhhw_navd_mosaic_polygon_dissolve_{0}' .format(region)

            mhhw_poly_checkMain = arcpy.ListFeatureClasses(outname_mhhw_poly)
            mhhw_poly_checkReg = arcpy.ListFeatureClasses(outname_mhhw_poly_reg)
            if len(mhhw_poly_checkMain)>0 and len(mhhw_poly_checkReg)>0:
                print('NOAA MHHW polygon already created for ' + region)
            
            else:
                if len(mhhw_poly_checkMain)>0:  #is this right?  CHECK

                    # make raster values compatible with feature class conversion (integer); actual values not used henceforth
                    mhhw_raster_mos_calc = arcpy.sa.RasterCalculator(' "all_noaa_mhhw_navd_mosaic" *100000000 ')
                    mhhw_raster_mos_calc.save(str(outname_mhhw_mos) + '_cellX10e8')

                    mhhw_raster_mos_calc_int = Int(mhhw_raster_mos_calc)
                    mhhw_raster_mos_calc_int.save(str(mhhw_raster_mos_calc) + '_int')

                    # make polygon from mosaiced raster
                    # raster_in = arcpy.ListRasters(mhhw_raster_mos_calc_int)[0]
                    with arcpy.EnvManager(outputZFlag="Disabled", outputMFlag="Disabled"):
                        # arcpy.conversion.RasterToPolygon(r"C:\Users\annik\Documents\ArcGIS\Projects\UCS_Flood_2022\UCS_Flood_2022.gdb\all_noaa_mhhw_navd_mosaic_cellX10e8_int", r"C:\Users\annik\Documents\ArcGIS\Projects\UCS_Flood_2022\UCS_Flood_2022.gdb\all_noaa_mhhw_navd_mosaic_polygon", "SIMPLIFY", "Value", "SINGLE_OUTER_PART", None)
                        mhhw_poly = arcpy.conversion.RasterToPolygon(mhhw_raster_mos_calc_int, outname_mhhw_poly, "SIMPLIFY", "Value", "SINGLE_OUTER_PART", None)
                
                    #reset gridcode values
                    mhhw_poly_calc1 = arcpy.management.CalculateField(mhhw_poly, "gridcode_float", "!gridcode! / 100000000", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")

                    #add fields for regions
                    with arcpy.EnvManager(extent="-129.639770225289 31.6125128491048 -111.698078318046 49.8281237167632"): #west coast
                        mhhw_poly_calc2 = arcpy.management.CalculateField(mhhw_poly_calc1, "region1", "'west_coast'", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")
                    with arcpy.EnvManager(extent="-103.188871709047 23.6185868603371 -64.1317607620386 48.7390194518053"): #east coast
                        mhhw_poly_calc3 = arcpy.management.CalculateField(mhhw_poly_calc2, "region1", "'east_coast'", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")
                    with arcpy.EnvManager(extent="-70.408042225505 14.9827777633202 -60.1390883470893 21.5874793767505"): #prvi
                        mhhw_poly_calc4 = arcpy.management.CalculateField(mhhw_poly_calc3, "region1", "'prvi_coast'", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")
                    
                    mhhw_poly = mhhw_poly_calc4

                if len(mhhw_poly_checkReg)>0:

                    #dissolve polygon by regions
                    mhhw_diss = arcpy.management.Dissolve(mhhw_poly, "all_noaa_mhhw_navd_mosaic_polygon_dissolve", "region", None, "MULTI_PART", "DISSOLVE_LINES")
                    
                    #create regional polygon
                    arcpy.conversion.FeatureClassToFeatureClass(mhhw_diss, gdb_out, "all_noaa_mhhw_navd_mosaic_polygon_dissolve_{0}" .format(region), "region LIKE '%east%'", 'region "region1" true true false 255 Text 0 0,First,#,all_noaa_mhhw_navd_mosaic_polygon_dissolve,region,0,255;Shape_Length "Shape_Length" false true true 8 Double 0 0,First,#,all_noaa_mhhw_navd_mosaic_polygon_dissolve,Shape_Length,-1,-1;Shape_Area "Shape_Area" false true true 8 Double 0 0,First,#,all_noaa_mhhw_navd_mosaic_polygon_dissolve,Shape_Area,-1,-1', '')

                print('Created polygons for MHHW layer in ' + region)

def prep_aev_mhhw_layers (region): #GOOD
    print('RUNNING PREP MHHW LAYERS FUNCTION')
    print ('...running ' + str(region) + '...')

    # ----------------------------------------------------------------------------------------
    # set gdb in and out
    gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
    # set gdb for tidal surfaces not made by NOAA matching local DEM vertical datums, made by AEV
    gdb_tidal_aev = 'C:/Users/annik/Documents/UCS/Data/MHHW_Ext/tidal_surf_aev/tidal_surf_aev.gdb'
    # set path for saving
    gdb_out = gdb

    if region=='pr_coast' or region=='vi_coast' or region=='gu_coast' or region=='hi_coast':
        
        # set workspace
        arcpy.env.workspace = gdb_tidal_aev

        mhhw_raster_strings = arcpy.ListRasters('*{0}*' .format(region)) ; print(mhhw_raster_strings)
        mhhw_poly_strings = arcpy.ListFeatureClasses('*{0}*' .format(region)) ; print(mhhw_poly_strings)

        # check if tidal surfaces have been imported to project gdb
        arcpy.env.workspace = gdb

        for tidal_raster in mhhw_raster_strings:
            if arcpy.Exists(tidal_raster)==True:
                print('Tidal surface raster {0} already in project geodatabase' .format(str(tidal_raster)))
            elif arcpy.Exists(tidal_raster)==False:
                print('Copying raster...')
                arcpy.env.workspace = gdb_tidal_aev
                arcpy.conversion.RasterToGeodatabase(tidal_raster,gdb_out)
                print('Tidal surface raster copied to project geodatabase: {0}' .format(str(tidal_raster)))
                                
        for tidal_poly in mhhw_poly_strings:
            if arcpy.Exists(tidal_poly)==True:
                print('Tidal surface polygon {0} already in project geodatabase' .format(str(tidal_poly)))
            elif arcpy.Exists(tidal_poly)==False:
                print('Copying polygon...')
                arcpy.env.workspace = gdb_tidal_aev
                arcpy.conversion.FeatureClassToGeodatabase(tidal_poly,gdb_out)
                print('Tidal surface polygon copied to project geodatabase: {0}' .format(str(tidal_poly)))

    else:
        print('ERROR: No tidal surface layer for ' + region + ' region') ; exit()

# BEFORE CONTINUING: 
# create transect points layers 
#   (make polylines feature class, draw lines along coast through TG stations, and generate points along lines), 
# generate near table using TG station list so Near_FID field == the nearby TG Station ID,
# join TG station list, add Station_ID field (type long), and calculate the field from the TG station list Station ID field,
# then remove the join.
# This will create a transect points layer with near values and matching station IDs 
# which will allow water level projections to be joined by station ID.

# This requires the near table to be generated and Station_ID field to be added before running; should be in projected coordinate system
def prep_gauge_data_for_interpolation (region,projection):
            
    print('RUNNING PREP GAUGE DATA FUNCTION')
    print ('...running ' + str(region) + ' ' + str(projection) + '...')
    if showtimes==1: tf1 = datetime.datetime.now() ; #print('Time start: ' + str(tf1))

    # ----------------------------------------------------------------------------------------
    # Set gdb, env
    gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
    arcpy.env.workspace = gdb

    # ----------------------------------------------------------------------------------------
    # Make layers from empty transect points file and file with gauges and water levels for each year/projection

    if region=='pr_coast' or region=='vi_coast': region_str = 'pr_usvi_coast'
    else: region_str = region

    gauges_name = 'all_{0}_stations_26x_{1}_proj' .format(region_str, projection) #from csv from tg gauges analysis
    if arcpy.Exists(gauges_name)==True:
        gauges_layer = arcpy.ListFeatureClasses(gauges_name)[0]
    elif arcpy.Exists(gauges_name)==False:
        print('Import tide gauge analysis data') ; exit()
    print('Gauge layer is: ' + str(gauges_layer))

    transect_points_name = '{0}_empty_points_proj_2022_join_{1}' .format(region_str,projection) #2022 here refers to year of project, not analysis/projection
    check_transect_layer = arcpy.ListFeatureClasses(transect_points_name)
    if len(check_transect_layer)!=0:
        transect_points_layer = arcpy.ListFeatureClasses(transect_points_name)[0]
    elif len(check_transect_layer)==0:
        transect_points_layer = arcpy.MakeFeatureLayer_management('{0}_empty_points_proj_2022'.format(region_str),transect_points_name) 

        # Join empty points and gauges file
        arcpy.AddJoin_management(transect_points_layer,'Station_ID', gauges_layer, 'Station_ID','KEEP_COMMON')

        print('Joined data to points for ' + str(projection) + ' projection')
    print('Transect points layer is: ' + str(transect_points_layer))
    
    print ('STATUS: gauge data prep run successfully for ' + str(region) + ' ' + str(projection))
    if showtimes==1: tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2)) 

    print('END PREP GAUGE DATA FUNCTION; RETURNING TO INTERPOLATE AND CREATE WLS FUNCTION')
    return transect_points_layer          

def interpolate_and_create_water_level_surfaces (regions,projections,years, flood_frequency=flood_frequency):  #WORKING, double check in full run
    print('RUNNING CREATE WATER LEVEL SURFACES FUNCTION')
    for region in regions:
        for projection in projections:

            print ('...running ' + region + ' ' + projection + '...')
            if showtimes==1: tf1 = datetime.datetime.now() ; print('Time start: ' + str(tf1))

            # ----------------------------------------------------------------------------------------
            # Set gdb in and out
            gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
            # set path for saving
            gdb_out = gdb

            # Set workspace
            arcpy.env.workspace = gdb

            # ----------------------------------------------------------------------------------------
            # Run prep functions

            print('...running prep functions...')
            #prep_noaa_mhhw_layers(region)  #function not tested
            prep_aev_mhhw_layers(region)
            transect_points_layer = prep_gauge_data_for_interpolation(region, projection)

            # Comment out following lines:
            # number_of_points = arcpy.GetCount_management('transect_points_layer')
            # count = int(number_of_points.getOutput(0))
            # print(count)

            if region=='pr_coast' or region=='vi_coast': region_str = 'pr_usvi_coast'
            else: region_str = region

            # print('List of fields in transect points layer:')
            # fields = arcpy.ListFields(transect_points_layer)
            # for field in fields:
            #     print(field.name)

            fields = arcpy.ListFields(transect_points_layer)
            check_join = []
            for field in fields:
                if str(field.name)=='all_{0}_stations_26x_{1}_proj.{1}_{2}' .format(region_str,projection,'2022'):
                    check_join.append(str(field.name))
            if len(check_join)!=0:
                print('Tide gauge data successfully joined to points layer')

            print('Prep functions completed')

            # ----------------------------------------------------------------------------------------
            # Load regional mhhw data
            print('...loading regional tidal surface water data...')

            if region=='west_coast' or region=='east_coast':
                mhhw_raster_string = 'all_noaa_mhhw_navd_mosaic'
                mhhw_poly_string = 'all_noaa_mhhw_navd_mosaic_polygon_dissolve_{0}' .format(region)
            elif region=='pr_coast':
                mhhw_raster_string = 'aev_tidal_surface_pr_coast_mhhw_prvd'
                mhhw_poly_string = 'aev_tidal_surface_pr_coast_mhhw_prvd_polygon'
            elif region == 'vi_coast':
                mhhw_raster_string = 'aev_tidal_surface_vi_coast_mhhw_vivd'
                mhhw_poly_string = 'aev_tidal_surface_vi_coast_mhhw_vivd_polygon'
            elif region == 'gu_coast':
                mhhw_raster_string = 'aev_tidal_surface_gu_coast_mhhw_guvd'
                mhhw_poly_string = 'aev_tidal_surface_gu_coast_mhhw_guvd_polygon'
            elif region == 'hi_coast':
                mhhw_raster_string = 'aev_tidal_surface_hi_coast_mhhw_msl'
                mhhw_poly_string = 'aev_tidal_surface_hi_coast_mhhw_msl_polygon'                
            else:
                print('ERROR: No tidal surface layer for ' + region + ' region') ; exit()

            mhhw_data = arcpy.ListRasters(mhhw_raster_string)[0]
            mhhw_surface = Raster(mhhw_data)
                
            print('Regional tidal surface is: ' + str(mhhw_surface))
            print('Tidal surface data loaded')
            
            # Set interpolation limit as polygon created from regional coastal counties
            mhhw_mask = arcpy.ListFeatureClasses(mhhw_poly_string)[0]
            arcpy.env.mask = mhhw_mask
            print('Mask is: ' + str(mhhw_mask))

            # Get cell size properties of mhhw_layer to use in interpolation between transect points
            cellsizes = []
            size_tmp = arcpy.GetRasterProperties_management(mhhw_surface, "CELLSIZEX")
            size = size_tmp.getOutput(0)
            cellsizes.append(size)
            cellsizes = [float(item) for item in cellsizes]
            cellsize = cellsizes[0]
            if region=='pr_coast' or region=='vi_coast':
                cellsize = 0.00122148208269197  #to match other regions; interpolation covers both pr and vi areas
            print('Cell size of mhhw layer is: ' + str(cellsize)) ; #cellsize = cellsize*4
            
            # ----------------------------------------------------------------------------------------
            # For each year/proj horizon, interpolate empty points file based on field from joined table
            for year in years:

                print('Region is: ' + region + ' and projection is: ' + projection + ' and year is: ' + year)

                # Interpolate between points
                print('...creating interpolated surface...')
                outname_interpol = 'interpol_points_{0}_proj_{1}_{2}' .format(region_str,projection,year)

                if arcpy.Exists(outname_interpol)==True:
                    print('Interpolated surface already created')
                    interpolated_surface_wrt_mhhw = outname_interpol ###arcpy.ListRasters(outname_interpol)[0]

                elif arcpy.Exists(outname_interpol)==False:
                    interpolated_surface_wrt_mhhw = NaturalNeighbor(transect_points_layer,str('all_{0}_stations_{1}x_{2}_proj.{3}_{4}' .format(region_str,flood_frequency,projection,  projection,year) ), cellsize)
                    print('Created interpolated surface')

                    print('...saving interpolated surface...')
                    interpolated_surface_wrt_mhhw.save(outname_interpol)
                    print('Saved interpolated surface layer')

                # Add mhhw water surface to interpolated projection to get total water level wrt navd or other vertical orthometric datum, in meters
                if region=='east_coast' or region=='west_coast': vdatum_name = 'navd'
                elif region=='pr_coast': vdatum_name = 'prvd'
                elif region=='vi_coast': vdatum_name = 'vivd'
                elif region=='gu_coast' : vdatum_name = 'guvd'
                elif region=='hi_coast' : vdatum_name = 'msl'

                outname_wls = 'water_level_surface_wrt_{0}_{1}x_{2}_{3}_{4}' .format(vdatum_name,flood_frequency,year,projection,region)  
                
                if arcpy.Exists(outname_wls)==True:
                    print('Water level surface already created')

                elif arcpy.Exists(outname_wls)==False:
                    print('...creating water level surface...')
                    
                    water_level_surface =  Raster(mhhw_surface) + Raster(interpolated_surface_wrt_mhhw)  #GOOD
                    print('Created water level surface')

                    print('Outname is: ' + outname_wls)
            
                    print('...saving water level surface...')
                    water_level_surface.save(outname_wls)
                    print('Saved water level surface layer')

                # outfind = arcpy.ListRasters(outname_wls) ; print(len(outfind))
                
                print ('STATUS: interpolation and create water level surfaces run successfully for ' + region + ' ' + projection + ' ' + year)
                if showtimes==1: tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2)) 

    print('END CREATE WATER LEVEL SURFACES FUNCTION')

# ----------------------------------------------------------------------------------------------------
# INUNDATED AREAS               # 8 hours east
def inundated_areas (regions,projections,years, flood_frequency=flood_frequency):
    print('RUNNING INUNDATED SURFACES FUNCTION')
    for region in regions:
        for projection in projections:
            for year in years:
                
                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                if showtimes==1: tf1 = datetime.datetime.now() ; print('Time start: ' + str(tf1))

                # ----------------------------------------------------------------------------------------
                # Set gdb in and out
                gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                # Set path for loading dems
                gdb_dem = 'C:/Users/annik/Documents/UCS/Data/2022/DEMs_fix/{0}/{0}_dem_fix.gdb' .format(region)  #
                # Set path for saving inundated raster
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                if runscratch==1:
                    gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/scratch.gdb'

                # Set workspace
                arcpy.env.workspace = gdb

                # ----------------------------------------------------------------------------------------
                # Set vdatum by region
                if region=='east_coast' or region=='west_coast': vdatum_name = 'navd' ; print('Contiguous US vertical datum is NAVD88')
                elif region=='pr_coast': vdatum_name = 'prvd' ; print('Puerto Rico vertical datum is PRVD02')
                elif region=='vi_coast': 
                    vdatum_name = 'vivd' ; print('US Virgin Isles vertical datum is VIVD09')
                    print('  NOTE: USVI DEM name is misprinted; actual vertical datum is VIVD, not NAVD')
                elif region=='gu_coast': vdatum_name = 'guvd' ; print('Guam vertical datum is GUVD04')
                elif region=='hi_coast': 
                    vdatum_name = 'msl' ; print('Hawaii vertical datum is LMSL ([Local] Mean Sea Level; tidal)')
                    print('  NOTE: Hawaii does not yet have a standardized orthometric vertical datum')

                # ----------------------------------------------------------------------------------------
                # Get the specific WLS from the gdb
                water_level = arcpy.ListRasters('water_level_surface_wrt_{0}_26x_{1}_{2}_{3}'.format(vdatum_name,year,projection,region))
                water_level_surface = Raster(water_level)
                print('Water level surface is: ' + str(water_level_surface))

                # ----------------------------------------------------------------------------------------
                # Get dems
                rundemmanual=0 #manual dem assignment turned off; auto assigned via regional dem gdbs
                if rundemmanual==1:
                    west_dems = ['CA_ChannelIslands_GCS_5m_NAVDm','CA_EKA1_GCS_5m_NAVD88m','CA_EKA2_GCS_5m_NAVD88m','CA_LOX1_GCS_5m_NAVD88m','CA_LOX2_GCS_5m_NAVD88m','CA_MTR1_GCS_5m_NAVD88m','CA_MTR2_GCS_5m_NAVD88m','CA_MTR3_GCS_5m_NAVD88m','CA_SGX_GCS_5m_NAVD88m','OR_MFR_GCS_5m_NAVD88m','OR_PQR1_GCS_5m_NAVD88m_TillamookFix','OR_PQR2_GCS_5m_NAVD88m','WA_PQR_GCS_5m_NAVD88m','WA_PugetSoundNW_GCS_3m_NAVDm','WA_PugetSoundSW_GCS_3m_NAVDm','WA_SEW1_GCS_5m_NAVD88m','WA_SEW1b_GCS_5m_NAVD88m','WA_SEW2_GCS_5m_NAVDm']
                    east_dems = ['AL_GCS_3m_NAVDm','CT_GCS_3m_NAVDm','DC_GCS_3m_NAVDm','DE_GCS_3m_NAVDm','FL_A3_GCS_5m_NAVDm','FL_JAX_1_GCS_5m_NAVD88m','FL_JAX_2_GCS_5m_NAVD88m','FL_MFL_1_GCS_5m_NAVD88m','FL_MFL_2_GCS_5m_NAVD88m','FL_MLB_1_GCS_5m_NAVD88m','FL_MLB_2_GCS_5m_NAVD88m','FL_Pan_East_GCS_3m_NAVDm','FL_Pan_West_GCS_3m_NAVDm','FL_TBW_1_GCS_5m_NAVD88m','FL_TBW_2_GCS_5m_NAVD88m','GA_CHS_GCS_5m_NAVDm','GA_JAX_GCS_5m_NAVDm','LA_Central_GCS_3m_NAVDm','LA_CentralEast_GCS_3m_NAVDm','LA_CentralNorth_GCS_3m_NAVDm','LA_Delta_GCS_3m_NAVDm','LA_LakePontchartrain_GCS_3m_NAVDm','LA_West_GCS_3m_NAVDm','MA_GCS_3m_NAVDm','MD_East_GCS_3m_NAVDm','MD_North_GCS_3m_NAVDm','MD_Southeast_GCS_3m_NAVDm','MD_Southwest_GCS_3m_NAVDm','MD_West_GCS_3m_NAVDm','ME_East_GCS_3m_NAVDm','ME_West_GCS_3m_NAVDm','MS_GCS_3m_NAVDm','NC_Middle1_GCS_3m_NAVDm','NC_Middle2_GCS_3m_NAVDm','NC_Northern_GCS_3m_NAVDm','NC_Southern1_GCS_3m_NAVDm','NC_Southern2_GCS_3m_NAVDm','NH_GCS_3m_NAVDm','NJ_Middle_GCS_3m_NAVDm_FY21Update','NJ_Northern_GCS_3m_NAVDm','NJ_Southern_GCS_3m_NAVDm_FY21Update','NY_Hudson_GCS_3m_NAVDm','NY_Metro_GCS_3m_NAVDm','NY_Suffolk_GCS_3m_NAVDm','PA_GCS_3m_NAVDm','RI_GCS_3m_NAVDm','SC_Central_GCS_3m_NAVDm','SC_North_GCS_3m_NAVDm','SC_South_GCS_3m_NAVDm','TX_Central_GCS_3m_NAVDm','TX_North1_GCS_3m_NAVDm','TX_North2_GCS_3m_NAVDm','TX_South1_GCS_3m_NAVDm','TX_South2_GCS_3m_NAVDm','VA_EasternShore_GCS_3m_NAVDm','VA_Middle_GCS_3m_NAVDm','VA_Northern_GCS_3m_NAVDm','VA_Southern_GCS_3m_NAVDm']
                    hi_dems = ['HI_HFO_Hawaii_UTM_3m_LMSLm','HI_HFO_Kauai_UTM_3m_LMSLm','HI_HFO_Lanai_UTM_3m_LMSLm','HI_HFO_Maui_UTM_3m_LMSLm','HI_HFO_Molokai_UTM_3m_LMSLm','HI_Oahu_GCS_3m_LMSLm']
                    gu_dems = ['Guam_GCS_3m_GUVD04m']
                    as_dems = ['AS_OfuOlosega_DEM_v2_UTM','AS_Tau_DEM_v2_UTM','AS_Tutuila_DEM_v2_UTM']
                    pr_dems = ['PR_GCS_3m_PRVD02m']
                    vi_dems = ['USVI_GCS_3m_NAVDm']

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

                # Reset workspace
                arcpy.env.workspace = gdb_dem

                dems = arcpy.ListRasters('*')
                print(dems)

                num_elem = len(dems)
                print (str(num_elem) + ' inundation layers to create')

                # Reset workspace
                arcpy.env.workspace = gdb

                # ----------------------------------------------------------------------------------------
                # Compare dems to water level surface and output inundated raster
                for dem in dems:

                    print ('DEM is: ' + str(dem))

                    # Set path for pulling in corrected dems
                    dem_path = str(gdb_dem) + '/' + str(dem) ; #print(dem_path)

                    # Set path out
                    outname_inundated_area_surface = 'inundated_area_surface_{0}x_{1}_{2}_{3}'.format(flood_frequency, year, projection, region) + '_' + dem[:20]
                    outname_inundated_area_surface_full = str(gdb_out) + '/' + str(outname_inundated_area_surface)

                    # Check if inundation file already exists, create if not
                    arcpy.env.workspace = gdb_out
                    inundated_check = arcpy.ListRasters(outname_inundated_area_surface) ; #print(len(inundated_check))
                    arcpy.env.workspace = gdb
                    if len(inundated_check)==0:  #if not already created

                        # Set mask to dem extent
                        arcpy.env.mask = dem_path
                        arcpy.env.cellSize = "MINOF"

                        # Compare DEM and water level surface. Where WLS >= DEM, set value of 1 ; i.e. water above land is TRUE
                        print ('...creating inundated area surface...')
                        inundated_area_surface = Con(Raster(water_level_surface) >= Raster(dem_path), 1) # creates flat inundation area raster
                        # print ('Created inundated area surface for ' + str(dem))

                        # Save raster
                        print('...saving layer...')
                        inundated_area_surface.save(outname_inundated_area_surface_full)
                        # print('Saved layer')

                    print ('Created inundated area surface: ' + str(outname_inundated_area_surface))

                print ('STATUS: inundated areas run successfully for ' + str(projection) + ' ' + str(year))
                if showtimes==1: tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))       
    
    print('END INUNDATED SURFACES FUNCTION')

def inundated_areas_select_dems (regions,projections,years,dems, flood_frequency=flood_frequency):  #for individual dems
    print('RUNNING INUNDATED SURFACES (SELECT DEMS) FUNCTION')
    for region in regions:
        for projection in projections:
            for year in years:
                    
                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                if showtimes==1: tf1 = datetime.datetime.now() ; print('Time start: ' + str(tf1))

                # ----------------------------------------------------------------------------------------
                # Set gdb in and out
                gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                # Set path for loading dems
                gdb_dem = 'C:/Users/annik/Documents/UCS/Data/2022/DEMs_fix/{0}/{0}_dem_fix.gdb' .format(region)  #
                # Set path for saving inundated raster
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

                # Set workspace
                arcpy.env.workspace = gdb

                # ----------------------------------------------------------------------------------------
                # Get the specific WLS from the gdb
                water_level = arcpy.ListRasters('water_level_surface_wrt_navd_26x_{0}_{1}_{2}'.format(year,projection,region))
                water_level_surface = Raster(water_level)
                print('Water level surface is: ' + str(water_level_surface))

                # ----------------------------------------------------------------------------------------
                # Get dems
                # rundemmanual=0
                # if rundemmanual==1:
                #     west_dems = ['CA_ChannelIslands_GCS_5m_NAVDm','CA_EKA1_GCS_5m_NAVD88m','CA_EKA2_GCS_5m_NAVD88m','CA_LOX1_GCS_5m_NAVD88m','CA_LOX2_GCS_5m_NAVD88m','CA_MTR1_GCS_5m_NAVD88m','CA_MTR2_GCS_5m_NAVD88m','CA_MTR3_GCS_5m_NAVD88m','CA_SGX_GCS_5m_NAVD88m','OR_MFR_GCS_5m_NAVD88m','OR_PQR1_GCS_5m_NAVD88m_TillamookFix','OR_PQR2_GCS_5m_NAVD88m','WA_PQR_GCS_5m_NAVD88m','WA_PugetSoundNW_GCS_3m_NAVDm','WA_PugetSoundSW_GCS_3m_NAVDm','WA_SEW1_GCS_5m_NAVD88m','WA_SEW1b_GCS_5m_NAVD88m','WA_SEW2_GCS_5m_NAVDm']
                #     east_dems = ['AL_GCS_3m_NAVDm','CT_GCS_3m_NAVDm','DC_GCS_3m_NAVDm','DE_GCS_3m_NAVDm','FL_A3_GCS_5m_NAVDm','FL_JAX_1_GCS_5m_NAVD88m','FL_JAX_2_GCS_5m_NAVD88m','FL_MFL_1_GCS_5m_NAVD88m','FL_MFL_2_GCS_5m_NAVD88m','FL_MLB_1_GCS_5m_NAVD88m','FL_MLB_2_GCS_5m_NAVD88m','FL_Pan_East_GCS_3m_NAVDm','FL_Pan_West_GCS_3m_NAVDm','FL_TBW_1_GCS_5m_NAVD88m','FL_TBW_2_GCS_5m_NAVD88m','GA_CHS_GCS_5m_NAVDm','GA_JAX_GCS_5m_NAVDm','LA_Central_GCS_3m_NAVDm','LA_CentralEast_GCS_3m_NAVDm','LA_CentralNorth_GCS_3m_NAVDm','LA_Delta_GCS_3m_NAVDm','LA_LakePontchartrain_GCS_3m_NAVDm','LA_West_GCS_3m_NAVDm','MA_GCS_3m_NAVDm','MD_East_GCS_3m_NAVDm','MD_North_GCS_3m_NAVDm','MD_Southeast_GCS_3m_NAVDm','MD_Southwest_GCS_3m_NAVDm','MD_West_GCS_3m_NAVDm','ME_East_GCS_3m_NAVDm','ME_West_GCS_3m_NAVDm','MS_GCS_3m_NAVDm','NC_Middle1_GCS_3m_NAVDm','NC_Middle2_GCS_3m_NAVDm','NC_Northern_GCS_3m_NAVDm','NC_Southern1_GCS_3m_NAVDm','NC_Southern2_GCS_3m_NAVDm','NH_GCS_3m_NAVDm','NJ_Middle_GCS_3m_NAVDm_FY21Update','NJ_Northern_GCS_3m_NAVDm','NJ_Southern_GCS_3m_NAVDm_FY21Update','NY_Hudson_GCS_3m_NAVDm','NY_Metro_GCS_3m_NAVDm','NY_Suffolk_GCS_3m_NAVDm','PA_GCS_3m_NAVDm','RI_GCS_3m_NAVDm','SC_Central_GCS_3m_NAVDm','SC_North_GCS_3m_NAVDm','SC_South_GCS_3m_NAVDm','TX_Central_GCS_3m_NAVDm','TX_North1_GCS_3m_NAVDm','TX_North2_GCS_3m_NAVDm','TX_South1_GCS_3m_NAVDm','TX_South2_GCS_3m_NAVDm','VA_EasternShore_GCS_3m_NAVDm','VA_Middle_GCS_3m_NAVDm','VA_Northern_GCS_3m_NAVDm','VA_Southern_GCS_3m_NAVDm']
                #     hi_dems = ['HI_HFO_Hawaii_UTM_3m_LMSLm','HI_HFO_Kauai_UTM_3m_LMSLm','HI_HFO_Lanai_UTM_3m_LMSLm','HI_HFO_Maui_UTM_3m_LMSLm','HI_HFO_Molokai_UTM_3m_LMSLm','HI_Oahu_GCS_3m_LMSLm']
                #     gu_dems = ['Guam_GCS_3m_GUVD04m']
                #     as_dems = ['AS_OfuOlosega_DEM_v2_UTM','AS_Tau_DEM_v2_UTM','AS_Tutuila_DEM_v2_UTM']
                #     pr_dems = ['PR_GCS_3m_PRVD02m']
                #     vi_dems = ['USVI_GCS_3m_NAVDm']

                # Reset workspace
                arcpy.env.workspace = gdb_dem
                dems = arcpy.ListRasters('*')
                ## dems = ['NJ_Middle_GCS_3m_NAVDm_FY21Update'] #running single dem
                print(dems)

                num_elem = len(dems)
                print (str(num_elem) + ' inundation layers to create')

                # ----------------------------------------------------------------------------------------
                # Compare dems to water level surface and output inundated raster
                for dem in dems:

                    print ('DEM is: ' + str(dem))

                    # Set path for pulling in corrected dems
                    dem_path = str(gdb_dem) + '/' + str(dem) ; #print(dem_path)

                    # Set path out
                    outname_inundated_area_surface = 'inundated_area_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region) + '_' + dem[:20]
                    outname_inundated_area_surface_full = str(gdb_out) + '/' + str(outname_inundated_area_surface)
                    print('Outname is: ' + str(outname_inundated_area_surface_full))

                    # Check if inundation file already exists, create if not
                    arcpy.env.workspace = gdb_out #gdb to check
                    inundated_check = arcpy.ListRasters(outname_inundated_area_surface) ; #print(len(inundated_check))
                    arcpy.env.workspace = gdb  #switch back to project gdb
                    if len(inundated_check)==0:
                        # Set mask to dem extent
                        arcpy.env.mask = dem_path
                        arcpy.env.cellSize = "MINOF"

                        # Compare DEM and water level surface. Where WLS >= DEM, set value of 1 ; i.e. 1 = water above land is TRUE
                        print ('...creating inundated area surface...')
                        inundated_area_surface = Con(Raster(water_level_surface) >= Raster(dem_path), 1) # creates flat inundation area raster
                        print ('Created inundated area surface for ' + str(dem))

                        # Save raster
                        print('...saving layer...')
                        inundated_area_surface.save(outname_inundated_area_surface_full)
                        print('Saved layer')

                    print ('Created inundated area surface: ' + str(outname_inundated_area_surface))

                print ('STATUS: inundated areas run successfully for ' + str(projection) + ' ' + str(year))
                if showtimes==1: tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))       
    
    print('END INUNDATED SURFACES (SELECT DEMS) FUNCTION')

# ----------------------------------------------------------------------------------------
# MOSAIC                        #5 hours east -- to 11          ***** ADD HAWAII PIECE *****
def mosaic_chunks (regions,projections,years, flood_frequency=flood_frequency, areas='all'):
    print('RUNNING MOSAIC INUNDATED SURFACES FUNCTION')
    for region in regions:
        for projection in projections:
            for year in years:
                
                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                if showtimes==1: tf1 = datetime.datetime.now() ; #print('Time start: ' + str(tf1))

                # ----------------------------------------------------------------------------------------
                # Set gdb in
                gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                if runscratch==1:
                    gdb_scratch = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/scratch.gdb'
                    gdb_in = gdb_scratch

                # Set workspace
                arcpy.env.workspace = gdb_in

                # ----------------------------------------------------------------------------------------
                # Set path out
                # gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                gdb_out = gdb_in
                if runscratchout==1:
                    gdb_scratch = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/scratch.gdb'
                    gdb_out = gdb_scratch

                # ----------------------------------------------------------------------------------------
                # Standardize area names
                if len(areas)==1 and areas[0]=='all' or (len(areas)==0) or 'all' in areas:
                    if region=='west_coast': areas = ['pnw_area','pacific_area']
                    if region=='east_coast': areas = ['gulf_area','atlantic_area']
                    if region=='hi_coast' :  areas = ['hawaii_isle','kauai_isle','lanai_isle','maui_isle','molokai_isle']
                else:
                    set_areas = []
                    for area in areas:
                        if region=='west_coast':
                            if 'pnw' in area.lower():
                                set_areas.append('pnw_area')
                            elif 'pacific' in area.lower():
                                set_areas.append('pacific_area')
                        if region=='east_coast':
                            if 'gulf' in area.lower():
                                set_areas.append('gulf_area')
                            elif 'atlantic' in area.lower():
                                set_areas.append('atlantic_area')
                        # pr and vi don't use areas; can remain 'all'
                    if len(set_areas)>0:
                        areas = set_areas

                # ----------------------------------------------------------------------------------------
                # Resample in certain regions due to processing speed, then mosaic

                if region == 'east_coast':
                    ##NOTES: resample 3 min per element
                    ##NOTES: full resample and mosaic, east coast ~5 hrs to run (pre-area break)
                    
                    # Set cellsize based on speed and quality testing
                    cellsize_out = "0.000045"
                    print('Cellsize out is ' + str(cellsize_out))
                    
                    # ----------------------------------------------------------------------------------
                    for area in areas:
                        print('...running ' + area + '...')
                        if showtimes==1: tta = datetime.datetime.now()
                        
                        # Set states to areas
                        if area=='gulf_area':
                            area_set = ['TX','LA','MS','AL','FL','GA_JAX'] #gulf
                        elif area=='atlantic_area':
                            area_set = ['FL_JAX_1','GA','SC','NC','VA','DC','MD','DE','PA','NJ','NY','CT','RI','MA','NH','ME'] #atlantic

                        # Get inundation raster surfaces raster list
                        raster_inund_files = []
                        for state in area_set:
                            rasters = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}_{4}*' .format(flood_frequency,year,projection,region,state))
                            for file in rasters:
                                raster_inund_files.append(file)

                        print('Files to resample and mosaic:')        
                        for file in raster_inund_files:
                            print(file)
                        
                        # Check area files
                        if area=='gulf_area':
                            if len(raster_inund_files)!=25:
                                print('Wrong number of resampled elements for ' + area + ' - need 25!'); exit()
                            else:
                                print (str(len(raster_inund_files)) + ' files total to resample in ' + area)  
                        elif area=='atlantic_area':
                            if len(raster_inund_files)!=35:
                                print('Wrong number of resampled elements for ' + area + ' - need 35!'); exit()
                            else:
                                print (str(len(raster_inund_files)) + ' files total to resample in ' + area)           
                        
                        # ----------------------------------------------------------------------------------
                        # Resample inundation rasters and save        #NOTES 3-5 min each raster !!!
                        list_resampled_rast_paths = []
                        for raster_inund_file in raster_inund_files:
                            # Set path for pulling in inundation raster
                            rast_path = str(gdb_in) + '/' + str(raster_inund_file)
                            
                            # Set resampled raster file name
                            outname_resample = 'scratch_' + str(raster_inund_file) + '_resamp_temp'  #cellsz_pt000045
                            outname_resample_path = str(gdb_out) + '/' + str(outname_resample)  

                            # Delete old resampled files
                            delweird=1
                            if delweird==1:
                                resample_check_TODEL = arcpy.ListRasters('scratch_' + str(raster_inund_file) + '_resamp*')
                                if len(resample_check_TODEL)!=0:
                                    for del_file in resample_check_TODEL:
                                        print('...deleting old resampled file...')
                                        arcpy.Delete_management(del_file)
                                        print('Deleted old resampled file: ' + str(del_file))

                            # Resample file, if resampled file doesn't already exist
                            resample_check = arcpy.ListRasters(outname_resample)
                            if len(resample_check)==0:  #if not already resampled
                                print ('File to resample: ' + str(raster_inund_file))

                                print ('Resample outname is: ' + outname_resample)

                                # getcellsize_temp = arcpy.GetRasterProperties_management((rast_test),"CELLSIZEX")
                                # cellsize_temp = getcellsize_temp.getOutput(0)
                                # print("Original cellsize is: " + str(cellsize_temp))

                                # Resample
                                print('...resampling...')
                                rast_resamp = arcpy.Resample_management(rast_path,outname_resample_path,cellsize_out,"NEAREST")  #NEAREST does not change cell values
                            print('Resampled ' + str(outname_resample))
                            
                            # Make new set of resampled raster path names to mosaic
                            list_resampled_rast_paths.append(outname_resample_path)
                        
                        if showtimes==1: ttb = datetime.datetime.now() ; print('Area resampling ' + area + ' Time start: ' + str(tta) + '  &  Time end: ' + str(ttb))
                    
                        # ----------------------------------------------------------------------------------
                        # Combine chunks and create mosaic, separated by area
                        print ('...resampling done, setting up mosaic...')

                        # for file in list_resampled_rast_paths:
                        #     print(file)

                        # Check number of files to be mosaiced
                        num_elem = len(list_resampled_rast_paths)
                        print (str(num_elem) + ' elements total to mosaic in ' + area)
                        if area=='gulf_area':
                            if num_elem!=25:
                                print('Wrong number of resampled elements for ' + area + ' - need 25!'); exit()
                        elif area=='atlantic_area':
                            if num_elem!=35:
                                print('Wrong number of resampled elements for ' + area + ' - need 35!'); exit()

                        outname_mosaic = 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area)
                        print ('Mosaic outname is: ' + outname_mosaic)

                        if showtimes==1: tta = datetime.datetime.now() ; print('Start time: ' + str(tta))

                        print ('...creating ' + area + ' mosaic...')
                        arcpy.MosaicToNewRaster_management(list_resampled_rast_paths, gdb_out, outname_mosaic,"","",cellsize_out,1,"MAXIMUM")   #cellsz_pt000045
                        # with arcpy.EnvManager(extent="-98.0633652753512 24.393838762558 -66.8803114143969 46.3976170001712"): #east coast extent
                        #     arcpy.MosaicToNewRaster_management(list_resampled_rast_paths, gdb_out, outname_mosaic,"","",cellsize_out,1,"MAXIMUM")   #cellsz_pt000045
                        print('Created mosaic for ' + area)

                        moscreated=1
                        if moscreated==1:
                            #delete resampled files
                            for del_file in list_resampled_rast_paths:
                                arcpy.Delete_management(del_file)
                            print('Deleted temporary resampled files')
                        
                        if showtimes==1: ttb = datetime.datetime.now() ; print('Area mosaic ' + area + ' Time start: ' + str(tta) + '  &  Time end: ' + str(ttb))
                        
                        print('Mosaicing ' + str(region) + ' ' + str(projection) + ' ' + str(year) + ' area ' + area + ' complete')
                    
                    print('Mosaicing complete for ' + str(region) + ' ' + str(projection) + ' ' + str(year) + ' ' + str(areas))

                elif region == 'west_coast':
                    runareas=0
                    if runareas==1:

                        for area in areas:
                            print('...running ' + area + '...')
                                                    
                            # Set states to areas
                            if area=='pnw_area':
                                area_set = ['WA_PugetSoundNW','WA_PugetSoundSW','WA_SEW2'] #pacific northwest, pnw
                            elif area=='pacific_area':
                                area_set = ['WA_SEW1','WA_SEW1b','WA_PQR','OR','CA'] #pacific

                            # Get inundation raster surfaces raster list
                            raster_inund_files = []
                            for state in area_set:
                                rasters = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}_{4}*' .format(flood_frequency,year,projection,region,state))
                                for file in rasters:
                                    raster_inund_files.append(file)

                            print('Files to resample and mosaic:')        
                            for file in raster_inund_files:
                                print(file)
                            
                            # Create list with full paths of inundation raster surfaces rasters
                            list_rast_paths = []
                            for raster_inund_file in raster_inund_files:
                                # Set path for pulling in inundation rasters
                                rast_path = str(gdb_in) + '/' + str(raster_inund_file)
                                list_rast_paths.append(rast_path)
                            print(list_rast_paths)

                            # Check number of files to be mosaiced
                            num_elem = len(list_rast_paths)
                            print (str(num_elem) + ' elements total to mosaic in ' + area)
                            if area=='pnw_area':
                                if len(raster_inund_files)!=3:
                                    print('Wrong number of resampled elements for ' + area + ' - need 3!'); exit()
                                else:
                                    print (str(len(raster_inund_files)) + ' files total to resample in ' + area)  
                            elif area=='pacific_area':
                                if len(raster_inund_files)!=15:
                                    print('Wrong number of resampled elements for ' + area + ' - need 15!'); exit()
                                else:
                                    print (str(len(raster_inund_files)) + ' files total to resample in ' + area)           
                        
                            if showtimes==1: tta = datetime.datetime.now() ; print('Start time: ' + str(tta))    

                            outname_merged = 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency, year, projection, region, area)  #change test nums
                            print ('Outname is: ' + outname_merged)

                            # Combine chunks and create mosaic
                            print ('...creating ' + area + ' mosaic...')
                            arcpy.MosaicToNewRaster_management(list_rast_paths, gdb_out, outname_merged,"","","",1,"MAXIMUM") 
                            print('Created mosaic for ' + area)
                                                       
                            if showtimes==1: ttb = datetime.datetime.now() ; print('Area mosaic ' + area + ' Time start: ' + str(tta) + '  &  Time end: ' + str(ttb))
                            print('Mosaicing area ' + area + ' complete')
                        
                        print('Mosaicing complete for ' + str(region) + ' ' + str(projection) + ' ' + str(year) + ' ' + str(areas))
                    
                    elif runareas==0:
                        # Create list with full paths of inundation raster surfaces rasters
                        list_rast_paths = []
                        for raster_inund_file in raster_inund_files:
                            # Set path for pulling in inundation rasters
                            rast_path = str(gdb_in) + '/' + str(raster_inund_file)
                            list_rast_paths.append(rast_path)
                        print(list_rast_paths)

                        outname_merged = 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region)
                        print ('Outname is: ' + outname_merged)

                        # Combine chunks and create mosaic
                        print ('...creating mosaic...')                    
                        arcpy.MosaicToNewRaster_management(list_rast_paths, gdb_out, outname_merged,"","","",1,"MAXIMUM") 
                        print('Created mosaic')     

                elif region=='pr_coast' or region=='vi_coast' or region=='gu_coast':
                    # no mosaics needed in these regions; simply rename inundation file

                    # Get inundation raster surfaces list
                    raster_inund_files = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}*' .format(flood_frequency,year,projection,region))

                    if len(raster_inund_files)==1:
                        print('No mosaic needed; only one inundated surface in {0}' .format(region))

                        print('Files to copy and rename:')        
                        for file in raster_inund_files:
                            print(file)

                        # Create list with full paths of inundation raster surfaces rasters
                        list_rast_paths = []
                        for raster_inund_file in raster_inund_files:
                            # Set path for pulling in inundation rasters
                            rast_path = str(gdb_in) + '/' + str(raster_inund_file)
                            list_rast_paths.append(rast_path)
                        #print(list_rast_paths)

                        outname_merged = 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region)  #change test nums
                        print ('Outname is: ' + outname_merged)

                        # Rename and save
                        if len(arcpy.ListRasters(outname_merged))!=0:
                            print('File already copied and saved')
                        elif len(arcpy.ListRasters(outname_merged))==0:
                            print ('...copying and renaming...')                    
                            arcpy.management.CopyRaster(list_rast_paths[0], outname_merged) 
                            print('Copied and saved inundated surface')   

                else:
                    print('region not built for mosaicing!')
                    exit()

                print ('STATUS: mosaic run successfully for ' + str(projection) + ' ' + str(year))
                if showtimes==1: tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))
    ## NOTES: original script had cellSize=10, bandsNum=1.  This script and ArcGIS Pro tool give errors with cellSize=10.
    #   Cellsize "0.000045" selected based on cell size of existing layers, speed testing, and quality control

    print('END MOSAIC INUNDATED SURFACES FUNCTION')

# CLEAR SCRATCH GDB
# Run if scratch gdb used during inundation or mosaic steps
def shift_rasters_from_scratch_gdb (regions,projections,years, flood_frequency=flood_frequency):  
    print('RUNNING FUNCTION SHIFTING INUNDATION AND MOSAIC RASTERS FROM SCRATCH TO MAIN GDB')
    for region in regions:
        for projection in projections:
            for year in years:
                if runscratchout==1 or runscratch==1:
                    print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                    if showtimes==1: tf1 = datetime.datetime.now() ; #print('Time start: ' + str(tf1))

                    # ----------------------------------------------------------------------------------------
                    # Set gdb in and out
                    gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                    gdb_scratch = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/scratch.gdb'
                    gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)                   
                    print('GDB out is: ' + gdb_out)

                    # Set workspace
                    arcpy.env.workspace = gdb_scratch
                    
                    # ----------------------------------------------------------------------------------------
                    # Grab rasters in scratch gdb and move to out gdb                
                    
                    rasters_to_shift = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}*' .format(flood_frequency,year,projection,region))
                    
                    if len(rasters_to_shift)>0:
                        for raster in rasters_to_shift:
                            print('Raster to shift: ' + str(raster))

                            print('...shifting...')
                            arcpy.conversion.RasterToGeodatabase(raster, gdb_out, '')
                            #print('Complete')
                        print('Inundation files shifted from scratch gdb')

                        for del_file in rasters_to_shift:
                            arcpy.Delete_management(del_file)
                        print('Inundation files deleted from scratch gdb')

                    elif len(rasters_to_shift)==0:
                        print('No inundation rasters found in scratch gdb')

                    rasters_to_shift = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}*' .format(flood_frequency,year,projection,region))
                    
                    if len(rasters_to_shift)>0:
                        for raster in rasters_to_shift:
                            print('Raster to shift: ' + str(raster))

                            print('...shifting...')
                            arcpy.conversion.RasterToGeodatabase(raster, gdb_out, '')
                            #print('Complete')
                        print('Mosaic files shifted from scratch gdb')

                        for del_file in rasters_to_shift:
                            arcpy.Delete_management(del_file)
                        print('Mosaic files deleted from scratch gdb')

                    elif len(rasters_to_shift)==0:
                        print('No mosaic rasters found in scratch gdb')

    print('STATUS: rasters shifted successfully')
    if showtimes==1: tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))
    print('END SHIFT RASTERS FROM SCRATCH GDB FUNCTION')

# ----------------------------------------------------------------------------------------
# REGION GROUP                  #3 - 4:30 hours east            ***** ADD HAWAII PIECE *****
def group_regions (regions,projections,years, flood_frequency=flood_frequency, areas='all'):
    print('RUNNING REGION GROUPING FUNCTION FOR CONNECTIVITY')
    for region in regions:
        for projection in projections:
            for year in years:

                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                if showtimes==1: tf1 = datetime.datetime.now() ; #print('Time start: ' + str(tf1))

                # ----------------------------------------------------------------------------------------
                # Set gdb in and out
                gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

                # Set workspace
                arcpy.env.workspace = gdb_in

                # ----------------------------------------------------------------------------------------
                # Standardize area names
                if len(areas)==1 and areas[0]=='all' or (len(areas)==0) or 'all' in areas:
                    if region=='west_coast': areas = ['pnw_area','pacific_area']
                    if region=='east_coast': areas = ['gulf_area','atlantic_area','RI_block_isle']
                    if region=='hi_coast' :  areas = ['hawaii_isle','kauai_isle','lanai_isle','maui_isle','molokai_isle']
                else:
                    set_areas = []
                    for area in areas:
                        if region=='west_coast':
                            if 'pnw' in area.lower():
                                set_areas.append('pnw_area')
                            elif 'pacific' in area.lower():
                                set_areas.append('pacific_area')
                        if region=='east_coast':
                            if 'gulf' in area.lower():
                                set_areas.append('gulf_area')
                            elif 'atlantic' in area.lower():
                                set_areas.append('atlantic_area')
                            elif 'ri' or 'block' or 'rhode' in area.lower():
                                set_areas.append('RI_block_isle')
                        # pr and vi don't use areas; can remain 'all'
                    if len(set_areas)>0:
                        areas = set_areas
                
                # ----------------------------------------------------------------------------------------
                # Region group
                if region=='west_coast':
                    print('Region grouping for both mainland Pacific coast and PNW area')

                    # Perform region group for each area: pacific and pnw
                    for area in areas:
                        if area=='pacific_area':
                            # get merged raster
                            merged_raw_surface = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region))[0]
                            ## print(str(merged_raw_surface))

                            # ----------------------------------------------------------------------------------------
                            # Perform region group for mainland
                            print('Pacific raw surface name is: ' + str(merged_raw_surface))
                            surface_path = gdb_in + '/' +  str(merged_raw_surface) ; #print(surface_path)
                            ##filename = os.path.basename(fullname)

                            fullname = str(merged_raw_surface)
                            outname_rg = 'rg_' + fullname
                            if len(arcpy.ListRasters(outname_rg))==0:
                                print('...region grouping ' + fullname + '...')
                                outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                                outRegionGrp.save(outname_rg)
                                print('Region grouped ' + merged_raw_surface)  
                        
                        # ----------------------------------------------------------------------------------------
                        elif area=='pnw_area':
                            # get merged raster
                            merged_raw_surface = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW_mosaicAll' .format(flood_frequency, year, projection, region))[0]
                            ## print(str(merged_raw_surface))

                            # ----------------------------------------------------------------------------------------
                            # Perform region group for PNW
                            print('PNW raw surface name is: ' + str(merged_raw_surface))
                            surface_path_pnw = gdb_in + '/' +  str(merged_raw_surface) ; #print(surface_path_pnw)
                            ##filename = os.path.basename(fullname)

                            fullname = str(merged_raw_surface)
                            outname_rg = 'rg_' + 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW' .format(flood_frequency, year, projection, region)
                            if len(arcpy.ListRasters(outname_rg))==0:
                                print('...region grouping ' + fullname + '...')
                                outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                                outRegionGrp.save(outname_rg)
                                print('Region grouped ' + merged_raw_surface)  
                            
                    print ('STATUS: region group run successfully for mainland & PNW ' + str(region) + ' ' + str(projection) + ' ' + str(year))
                        
                elif region=='east_coast':
                    print('Region grouping for three east coast areas: Gulf Coast, Atlantic Coast, and RI Block Island')
                    # Perform region group for each area: mainland Gulf+Florida, Atlantic, and RI island
                    # Gulf and Atlantic have overlapping DEM area to allow for proper connectivity assessment at edge
                    
                    for area in areas:
                        if area=='gulf_area' or area=='atlantic_area':

                            # Get merged raster
                            merged_raw_surface = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area))[0]
                            #print(str(merged_raw_surface))

                            # ----------------------------------------------------------------------------------------
                            # Perform region group
                            print(area + ' raw surface name is: ' + str(merged_raw_surface))
                            surface_path = gdb_in + '/' +  str(merged_raw_surface) ; #print(surface_path)
                            #filename = os.path.basename(fullname)
                            fullname = str(merged_raw_surface)

                            # Check that rg doesn't already exist
                            check_rg = arcpy.ListRasters('rg_' + fullname)
                            if len(check_rg)>0:
                                print('Region group already performed for ' + area)
                            elif len(check_rg)==0:
                                print('...region grouping ' + fullname + '...')
                                outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                                outname_rg = 'rg_' + fullname
                                outRegionGrp.save(outname_rg)
                                print('Region grouped ' + str(merged_raw_surface))  
                        # ----------------------------------------------------------------------------------------
                        elif area=='RI_block_isle':

                            # Get RI inundation raster, NOT merged raster
                            raw_surface = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}_RI*' .format(flood_frequency,year,projection,region))[0]
                            #print(str(raw_surface))

                            # ----------------------------------------------------------------------------------------
                            # Perform region group
                            print(area + ' raw surface name is: ' + str(raw_surface))
                            surface_path_riisle = gdb_in + '/' +  str(raw_surface) ; #print(surface_path_riisle)
                            #filename = os.path.basename(fullname)
                            fullname = str(raw_surface)

                            # Check that rg doesn't already exist
                            check_rg = arcpy.ListRasters('rg_' + fullname)
                            if len(check_rg)>0:
                                print('Region group already performed for ' + area)
                            elif len(check_rg)==0:
                                print('...region grouping ' + fullname + '...')
                                outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                                outname_rg = 'rg_' + 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area)
                                outRegionGrp.save(outname_rg)
                                print('Region grouped ' + str(raw_surface)) 
                    
                    print ('STATUS: region group run successfully for ' + str(projection) + ' ' + str(year))
                
                elif region=='pr_coast' or region=='vi_coast' or region=='gu_coast' or region=='hi_coast':
                    print('Region grouping for {0}' .format(region))
                    
                    # Get merged raster
                    merged_raw_surface = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region))[0]

                    print('Raw surface name is: ' + str(merged_raw_surface))
                    surface_path = gdb_in + '/' +  str(merged_raw_surface) ; #print(surface_path)

                    # Set mask
                    arcpy.env.workspace = gdb
                    mask_name = arcpy.ListFeatureClasses('aev_tidal_surface_{0}_mhhw_*_polygon' .format(region))[0]
                    mask_path = gdb + '/' + mask_name
                    arcpy.env.workspace = gdb_in

                    # Perform region group; islands each have one DEM
                    fullname = str(merged_raw_surface)
                    outname_rg = 'rg_' + fullname
                    print('Outname is: ' + outname_rg)
                    
                    if len(arcpy.ListRasters(outname_rg))==0:
                        print('...region grouping ' + fullname + '...')
                        arcpy.env.mask = mask_path ; #print(mask_path)
                        outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                        outRegionGrp.save(outname_rg)
                        print('Region grouped ' + merged_raw_surface)  

                    print ('STATUS: region group run successfully for ' + str(region) + ' ' + str(projection) + ' ' + str(year))

                # Other regions
                else:
                    print ('Region not built yet!') ; exit()
                    print ('STATUS: region group run successfully for ' + str(projection) + ' ' + str(year))

                if showtimes==1: tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))

    print('END REGION GROUPING FUNCTION FOR CONNECTIVITY')

# ----------------------------------------------------------------------------------------
# EXTRACT CONNECTED REGIONS     #3 hours east
def extract_connected (regions,projections,years, flood_frequency=flood_frequency, areas='all'):
    print('RUNNING EXTRACTING CONNECTED REGION GROUPS FUNCTION')
    for region in regions:
        for projection in projections:
            for year in years:

                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                if showtimes==1: tf1 = datetime.datetime.now() ; #print('Time start: ' + str(tf1))

                # ----------------------------------------------------------------------------------------
                # Set gdb in and out
                gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

                # Set workspace
                arcpy.env.workspace = gdb_in

                # ----------------------------------------------------------------------------------------
                # Standardize area names
                if len(areas)==1 and areas[0]=='all' or (len(areas)==0) or 'all' in areas:
                    if region=='west_coast': areas = ['pnw_area','pacific_area']
                    if region=='east_coast': areas = ['gulf_area','atlantic_area','RI_block_isle']
                else:
                    set_areas = []
                    for area in areas:
                        if region=='west_coast':
                            if 'pnw' in area.lower():
                                set_areas.append('pnw_area')
                            elif 'pacific' in area.lower():
                                set_areas.append('pacific_area')
                        if region=='east_coast':
                            if 'gulf' in area.lower():
                                set_areas.append('gulf_area')
                            elif 'atlantic' in area.lower():
                                set_areas.append('atlantic_area')
                            elif 'ri' or 'block' or 'rhode' in area.lower():
                                set_areas.append('RI_block_isle')
                        # pr and vi don't use areas; can remain 'all'
                    if len(set_areas)>0:
                        areas = set_areas

                # ----------------------------------------------------------------------------------------
                # Extract
                if region=='west_coast':
                    print('Extracting both mainland Pacific west coast and PNW area')

                    # ----------------------------------------------------------------------------------------
                    # Perform region group for each area: pacific and pnw
                    for area in areas:
                        if area=='pacific_area':
                            print('Extracting west coast' + str.lower(area))

                            # Get region group raster
                            rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region))[0]
                            rg_surface_path = gdb_in + '/' + rg_surface
                            print(area + ' file to extract is ' + rg_surface)

                            # Set outname
                            outname_extract = 'extract_' + rg_surface

                            # Check that extract doesn't already exist
                            check_extract = arcpy.ListRasters('outname_extract')
                            if len(check_extract)!=0: print('Extraction already performed for ' + area)
                            
                            elif len(check_extract)==0:
                                # Set areas for extraction - mainland (maximum[0]) and Channel Isles (maximum[2])
                                print('...finding areas to extract from ' + rg_surface + '...')
                                
                                # index maximum counts for mainland and Channel Isles
                                arr = arcpy.da.FeatureClassToNumPyArray(rg_surface, ('Value', 'Count'))
                                count = arr['Count']
                                value = arr['Value']
                                
                                index_to_extract_imax1 = numpy.argsort(count)[-1]           #;print('imax1 = ' + str(index_to_extract_imax1))
                                value_to_extract_imax1 = str(value[index_to_extract_imax1]) #;print('val1 = ' + str(value_to_extract_imax1))

                                index_to_extract_imax3 = numpy.argsort(count)[-3]           #;print('imax3 = ' + str(index_to_extract_imax3))
                                value_to_extract_imax3 = str(value[index_to_extract_imax3]) #;print('val3 = ' + str(value_to_extract_imax3))

                                print('   Result: first and third maximum pixel counts correspond to values ' + value_to_extract_imax1 + ' (mainland) & ' + value_to_extract_imax3 + ' (Channel Isles)')

                                # build SQL clause syntax
                                print('...building SQL string...')
                                inSQLClause = 'Value = ' + value_to_extract_imax1 + ' Or Value = ' + value_to_extract_imax3  #values from regional conditions

                                print('Extracting values ' + value_to_extract_imax1 + ' & ' + value_to_extract_imax3 + ' from ' + rg_surface)
                                print('...extracting values...')

                                # extract connected areas       west start 2:27 end 2:54 ~30 min
                                if len(arcpy.ListRasters(outname_extract))==0:
                                    attExtract = ExtractByAttributes(rg_surface, inSQLClause)
                                    attExtract.save(outname_extract)
                                    print('Extracted connected areas from ' + rg_surface)
                        
                        # ----------------------------------------------------------------------------------------               
                        elif area=='pnw_area':
                            print('Extracting west coast ' + str(area))

                            # Get region group raster
                            rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW' .format(flood_frequency, year, projection, region))[0]
                            rg_surface_path = gdb_in + '/' + rg_surface
                            print(area + ' file to extract is ' + rg_surface)

                            # Set outname
                            outname_extract = 'extract_' + rg_surface
                            
                            # Check that extract doesn't already exist
                            check_extract = arcpy.ListRasters('outname_extract')
                            if len(check_extract)!=0: print('Extraction already performed for ' + area)

                            elif len(check_extract)==0:
                                # Set groups for extraction - PNW ONLY
                                print('...finding groups to extract from ' + rg_surface + '...')

                                # Index maximum count
                                arr = arcpy.da.FeatureClassToNumPyArray(rg_surface, ('Value', 'Count'))
                                count = arr['Count']
                                value = arr['Value']
                                index_to_extract = numpy.argmax(count)
                                value_to_extract = str(value[index_to_extract])

                                # Build SQL clause syntax
                                print('...building SQL string...')
                                inSQLClause = 'Value = ' + value_to_extract

                                print('Extracting value=' + value_to_extract + ' from ' + rg_surface)
                                print('...extracting values...')

                                # Extract connected areas
                                if len(arcpy.ListRasters(outname_extract))==0:
                                    attExtract = ExtractByAttributes(rg_surface, inSQLClause)
                                    attExtract.save(outname_extract)
                                    print('Extracted connected areas from ' + rg_surface)
                    
                    print ('STATUS: extract connected regions run successfully for mainland pacific & PNW ' + str(region) + ' ' + str(projection) + ' ' + str(year))
                    
                elif region=='east_coast':
                    print('Extracting east coast: Gulf, Atlantic, and RI island area')

                    # ----------------------------------------------------------------------------------------
                    # Perform region group for each area: gulf coast, atlantic coast, and RI isle
                    for area in areas:
                        if area=='gulf_area' or area=='atlantic_area':
                            print('Extracting east coast ' + area)

                            # Get region group raster
                            rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area))[0]
                            rg_surface_path = gdb_in + '/' + rg_surface
                            print(area + ' file to extract is ' + rg_surface)

                            # Set outname
                            outname_extract = 'extract_' + rg_surface
                            
                            # ----------------------------------------------------------------------------------------
                            # Check that extract doesn't already exist
                            check_extract = arcpy.ListRasters(outname_extract)
                            if len(check_extract)!=0:
                                print('Extraction already performed for ' + area)
                            
                            elif len(check_extract)==0:
                                # ----------------------------------------------------------------------------------------
                                # Set groups for extraction
                                print('...finding groups to extract from ' + rg_surface + '...')

                                # Index maximum count
                                arr = arcpy.da.FeatureClassToNumPyArray(rg_surface, ('Value', 'Count'))
                                count = arr['Count']
                                value = arr['Value']
                                index_to_extract = numpy.argmax(count)
                                value_to_extract = str(value[index_to_extract])

                                # Build SQL clause syntax
                                print('...building SQL string...')
                                inSQLClause = 'Value = ' + value_to_extract

                                print('Extracting value=' + value_to_extract + ' from ' + rg_surface)
                                print('...extracting values...')

                                # Extract connected areas
                                attExtract = ExtractByAttributes(rg_surface, inSQLClause)
                                attExtract.save(outname_extract)
                                print('Extracted connected areas from ' + rg_surface)
                        
                        # ----------------------------------------------------------------------------------------                                       
                        elif area=='RI_block_isle':
                            print('Extracting east coast ' + area)

                            # Get region group raster
                            rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area))[0]
                            rg_surface_path = gdb_in + '/' + rg_surface
                            print(area + ' file to extract is ' + rg_surface)

                            # Set outname
                            outname_extract = 'extract_' + rg_surface
                        
                            # ----------------------------------------------------------------------------------------
                            # Check that extract doesn't already exist
                            check_extract = arcpy.ListRasters(outname_extract)
                            if len(check_extract)!=0:
                                print('Extraction already performed for ' + area)
                            
                            elif len(check_extract)==0:
                                # ----------------------------------------------------------------------------------------
                                # Set areas for extraction - RI Block Island (maximum[1]) in RI inundation layer
                                print('...finding areas to extract from ' + rg_surface + '...')
                                
                                # Index maximum counts
                                arr = arcpy.da.FeatureClassToNumPyArray(rg_surface, ('Value', 'Count'))
                                count = arr['Count']
                                value = arr['Value']
                                
                                index_to_extract_imax2 = numpy.argsort(count)[-2]           #;print('imax2 = ' + str(index_to_extract_imax2))
                                value_to_extract_imax2 = str(value[index_to_extract_imax2]) #;print('val2 = ' + str(value_to_extract_imax2))

                                print('   Result: second maximum pixel count corresponds to value ' + value_to_extract_imax2 + ' (RI island)')

                                # Build SQL clause syntax
                                print('...building SQL string...')
                                inSQLClause = 'Value = ' + value_to_extract_imax2  #values from regional conditions

                                print('Extracting value=' + value_to_extract_imax2 + ' from ' + rg_surface)
                                print('...extracting values...')

                                # Extract connected areas       RI isle ~1 min
                                attExtract = ExtractByAttributes(rg_surface, inSQLClause)
                                attExtract.save(outname_extract)
                                print('Extracted connected areas from ' + rg_surface)
                
                    print ('STATUS: extract connected regions run successfully for ' + str(region) + ' ' + str(projection) + ' ' + str(year))
                
                elif region=='pr_coast' or region=='vi_coast' or region=='gu_coast':

                    print('Extracting ' + str.lower(region))

                    # Get region group raster
                    rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region))[0]
                    rg_surface_path = gdb_in + '/' + rg_surface
                    print(region + ' file to extract is ' + rg_surface)

                    # Set outname
                    outname_extract = 'extract_' + rg_surface

                    # Check that extract doesn't already exist
                    check_extract = arcpy.ListRasters(outname_extract)
                    if len(check_extract)!=0:
                        print('Extraction already performed for ' + area)

                    elif len(check_extract)==0:
                        # Set mask
                        arcpy.env.workspace = gdb
                        mask_name = arcpy.ListFeatureClasses('aev_tidal_surface_{0}_mhhw_*_polygon' .format(region))[0]
                        mask_path = gdb + '/' + mask_name
                        arcpy.env.workspace = gdb_in

                        # ----------------------------------------------------------------------------------------
                        # Set conditions for extraction; islands each have one DEM with several pieces
                        print('...finding areas to extract from ' + rg_surface + '...')

                        arr = arcpy.da.FeatureClassToNumPyArray(rg_surface, ('Value', 'Count'))
                        count = arr['Count']
                        value = arr['Value']

                        # Index maximum counts for largest areas around islands
                        if region=='pr_coast':   # get 4 largest regions

                            index_to_extract_imax1 = numpy.argsort(count)[-1]           #;print('imax1 = ' + str(index_to_extract_imax1))
                            value_to_extract_imax1 = str(value[index_to_extract_imax1]) #;print('val1 = ' + str(value_to_extract_imax1))

                            index_to_extract_imax2 = numpy.argsort(count)[-2]           #;print('imax2 = ' + str(index_to_extract_imax2))
                            value_to_extract_imax2 = str(value[index_to_extract_imax2]) #;print('val2 = ' + str(value_to_extract_imax2))

                            index_to_extract_imax3 = numpy.argsort(count)[-3]           #;print('imax3 = ' + str(index_to_extract_imax3))
                            value_to_extract_imax3 = str(value[index_to_extract_imax3]) #;print('val3 = ' + str(value_to_extract_imax3))                        

                            index_to_extract_imax4 = numpy.argsort(count)[-4]           #;print('imax4 = ' + str(index_to_extract_imax4))
                            value_to_extract_imax4 = str(value[index_to_extract_imax4]) #;print('val4 = ' + str(value_to_extract_imax4))

                            print('   Result: first through fourth maximum pixel counts correspond to') 
                            print('   values ' + value_to_extract_imax1 + ' (Puerto Rico & Vieques) & ' + value_to_extract_imax2 + ' (Culebra) & ' + value_to_extract_imax3 + ' (Mona) & ' + value_to_extract_imax4 + ' (Desecheo)')

                            # Build SQL clause syntax
                            print('...building SQL string...')
                            inSQLClause = 'Value = ' + value_to_extract_imax1 + ' Or Value = ' + value_to_extract_imax2 + ' Or Value = ' + value_to_extract_imax3 + ' Or Value = ' + value_to_extract_imax4  #values from regional conditions

                        elif region=='vi_coast':   # get 2 largest regions

                            index_to_extract_imax1 = numpy.argsort(count)[-1]           #;print('imax1 = ' + str(index_to_extract_imax1))
                            value_to_extract_imax1 = str(value[index_to_extract_imax1]) #;print('val1 = ' + str(value_to_extract_imax1))

                            index_to_extract_imax2 = numpy.argsort(count)[-2]           #;print('imax2 = ' + str(index_to_extract_imax2))
                            value_to_extract_imax2 = str(value[index_to_extract_imax2]) #;print('val2 = ' + str(value_to_extract_imax2))

                            print('   Result: first and second maximum pixel counts correspond to values ' + value_to_extract_imax1 + ' (St. Thomas & St. John) & ' + value_to_extract_imax2 + ' (St. Croix)')

                            # Build SQL clause syntax
                            print('...building SQL string...')
                            inSQLClause = 'Value = ' + value_to_extract_imax1 + ' Or Value = ' + value_to_extract_imax2   #values from regional conditions

                        elif region=='gu_coast':   # get 1 largest region

                            index_to_extract_imax1 = numpy.argsort(count)[-1]           #;print('imax1 = ' + str(index_to_extract_imax1))
                            value_to_extract_imax1 = str(value[index_to_extract_imax1]) #;print('val1 = ' + str(value_to_extract_imax1))

                            print('   Result: first maximum pixel count corresponds to value ' + value_to_extract_imax1 + ' (Guam and surrounding minor islands)')

                            # Build SQL clause syntax
                            print('...building SQL string...')
                            inSQLClause = 'Value = ' + value_to_extract_imax1   #values from regional conditions

                        print('Extracting values from ' + rg_surface)
                        print('Outname is: ' + outname_extract)
                        print('SQL clause: ' + inSQLClause)
                        print('...extracting values...')

                        # Extract connected areas
                        if len(arcpy.ListRasters(outname_extract))==0:
                            arcpy.env.mask = mask_path ; #print(mask_path)
                            attExtract = ExtractByAttributes(rg_surface, inSQLClause)
                            attExtract.save(outname_extract)
                            print('Extracted connected areas from ' + rg_surface)

                    print ('STATUS: extract connected regions run successfully for ' + str(region) + ' ' + str(projection) + ' ' + str(year))

                else:
                    print('region not yet built') ; exit()
                    print ('STATUS: extract connected regions run successfully for _')

                if showtimes==1: tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))

    print('END EXTRACTING CONNECTED REGION GROUPS FUNCTION')

# ----------------------------------------------------------------------------------------
# CONVERT RASTER TO POLYGON     #1:30 - 2 hours east
def raster_to_polygon (regions,projections,years, flood_frequency=flood_frequency, areas='all'):
    print('RUNNING RASTER TO POLYGON SCRIPT')
    for region in regions:
        for projection in projections:
            for year in years:

                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                if showtimes==1: tf1 = datetime.datetime.now() ; #print('Time start: ' + str(tf1))

                # ----------------------------------------------------------------------------------------
                # Set gdb in and out
                gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

                # Set workspace
                arcpy.env.workspace = gdb_in

                # ----------------------------------------------------------------------------------------
                # Standardize area names
                if len(areas)==1 and areas[0]=='all' or (len(areas)==0) or 'all' in areas:
                    if region=='west_coast': areas = ['pnw_area','pacific_area']
                    if region=='east_coast': areas = ['gulf_area','atlantic_area','RI_block_isle']
                else:
                    set_areas = []
                    for area in areas:
                        if region=='west_coast':
                            if 'pnw' in area.lower():
                                set_areas.append('pnw_area')
                            elif 'pacific' in area.lower():
                                set_areas.append('pacific_area')
                        if region=='east_coast':
                            if 'gulf' in area.lower():
                                set_areas.append('gulf_area')
                            elif 'atlantic' in area.lower():
                                set_areas.append('atlantic_area')
                            elif 'ri' or 'block' or 'rhode' in area.lower():
                                set_areas.append('RI_block_isle')
                        # pr and vi don't use areas; can remain 'all'
                    if len(set_areas)>0:
                        areas = set_areas

                # ----------------------------------------------------------------------------------------
                # Raster to polygon
                if region == 'west_coast':
                    print('Building polygons for both mainland west coast and PNW area')

                    # Perform polygon conversion for each area: mainland and pnw
                    for area in areas:

                        # Get extracted rasters
                        if area=='pacific_area':
                            to_convert = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency,year,projection,region))[0]
                        elif area=='pnw_area':    
                            to_convert = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW' .format(flood_frequency,year,projection,region))[0]
                            
                        print(area + ' file to convert is ' + str(to_convert))

                        # Set path and file out
                        outname_polygon = 'polygon_' + str(to_convert)  #final_polygon_

                        # ----------------------------------------------------------------------------------------
                        # Check polygon doesn't already exist; create polygons       4 min total
                        check_poly = arcpy.ListFeatureClasses(outname_polygon)
                        if len(check_poly)!=0:
                            print('Polygons already performed for ' + area)

                        elif len(check_poly)==0:                            
                            print('...creating ' + str(area) + ' polygon...')
                            arcpy.RasterToPolygon_conversion(to_convert, outname_polygon, "SIMPLIFY", "VALUE")
                            print ('Converted ' + str(to_convert) + ' to polygon')  

                    print ('STATUS: polygon conversion run successfully for mainland Pacific and PNW west coast, ' + str(projection) + ' ' + str(year))
                    
                elif region=='east_coast':
                    print('Building polygons for east coast: Gulf, Atlantic, and RI island area')

                    # Perform polygon conversion for each area: gulf coast, atlantic coast, and RI isle
                    for area in areas:
                        if area=='gulf_area' or area=='atlantic_area' or area=='RI_block_isle':
                            print('Building polygons for east coast ' + area)

                            # Get extracted rasters
                            to_convert = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area))[0]
                            print(area + ' file to convert is ' + str(to_convert))

                            # Set path and file out
                            outname_polygon = 'polygon_' + str(to_convert)  #final_polygon_
                            print('Outname is: ' + outname_polygon)
                           
                            # ----------------------------------------------------------------------------------------
                            # Check polygon doesn't already exist; create polygons
                            check_poly = arcpy.ListFeatureClasses(outname_polygon)
                            if len(check_poly)!=0:
                                print('Polygons already performed for ' + area)

                            elif len(check_poly)==0:
                                print('...creating ' + area + ' polygon...')
                                arcpy.RasterToPolygon_conversion(to_convert, outname_polygon, "SIMPLIFY", "VALUE")
                                print ('Converted ' + str(to_convert) + ' to polygon')                    

                    print ('STATUS: polygon conversion run successfully for east coast areas ' + str(projection) + ' ' + str(year))

                elif region=='pr_coast' or region=='vi_coast' or region=='gu_coast' or region=='hi_coast':
                    print('Building polygons for ' + region)
                    
                    # Get extracted rasters
                    to_convert = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency,year,projection,region))[0]
                    print('File to convert is ' + str(to_convert))

                    # Set path and file out
                    outname_polygon = 'polygon_' + str(to_convert)  #final_polygon_
                    print('Outname is: ' + outname_polygon)
                    
                    # ----------------------------------------------------------------------------------------
                    # Check polygon doesn't already exist; create polygons
                    check_poly = arcpy.ListFeatureClasses(outname_polygon)
                    if len(check_poly)!=0:
                        print('Polygons already performed for ' + region)

                    elif len(check_poly)==0:
                        print('...creating polygon...')
                        arcpy.RasterToPolygon_conversion(to_convert, outname_polygon, "SIMPLIFY", "VALUE")
                        print ('Converted ' + str(to_convert) + ' to polygon')                    

                    print ('STATUS: polygon conversion run successfully for ' + str(region) + ' ' + str(projection) + ' ' + str(year))

                else:
                    print('region not yet built!')

                print('FINISHED ' + str(region) + ' ' + str(projection) + ' ' + str(year))
                if showtimes==1: tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))

    print('END RASTER TO POLYGON SCRIPT')


# ----------------------------------------------------------------------------------------
#SET FUNCTION INPUTS
regions     = ['gu_coast']      #['east_coast','west_coast','pr_coast','vi_coast','gu_coast'] #['hi_coast','as_coast']
projections = ['high','int_low']       #['high','int_low']
# years       = ['2030']          #['2022','2100','2060','2045','2035','2080','2050','2030','2065']
years       = ['2022','2030','2035','2045','2050','2060','2065','2080','2100']

# ----------------------------------------------------------------------------------------
#RUN FUNCTIONS

# dem_fix (regions)
# prep_gauge_data_for_interpolation (regions, projections)
# interpolate_and_create_water_level_surfaces (regions, projections, years)

# interpolate_and_create_water_level_surfaces(['east_coast'],['high'],['2025'])
# interpolate_and_create_water_level_surfaces(['pr_coast','vi_coast'] , ['high','int_high','int','int_low','low'] , ['2022','2030','2035','2045','2050','2060','2065','2080','2100'])


# interpolate_and_create_water_level_surfaces (regions, projections, years)
# inundated_areas (regions, projections, years)
# mosaic_chunks (regions, projections, years)
shift_rasters_from_scratch_gdb (regions, projections, years)
group_regions (regions, projections, years)
extract_connected (regions, projections, years)
raster_to_polygon (regions, projections, years)

# group_regions (['pr_coast'],['int_low'],years)
# extract_connected (['pr_coast'],projections,years, flood_frequency=flood_frequency, areas='all')
# raster_to_polygon ( ['pr_coast','vi_coast'], ['high','int_low'], years )

#### DONE
# inundated_areas_select_dems (['east_coast'], ['int_low'], ['2065'], ['NJ_Middle_GCS_3m_NAVDm_FY21Update'])
# mosaic_chunks         (['east_coast'], ['int_low'], ['2065'], ['atlantic_area'])
# group_regions         (['east_coast'], ['int_low'], ['2065'])
# extract_connected     (['east_coast'], ['int_low'], ['2065'])
# raster_to_polygon     (['east_coast'], ['int_low'], ['2065'])


#RUNNING BELOW
# inundated_areas       ('26', ['east_coast'], ['high'], ['2100']) #to scratch
# mosaic_chunks         ('26', ['east_coast'], ['high'], ['2100'], ['all']) #in scratch
# shift_rasters_from_scratch_gdb    ('26', ['east_coast'], ['high'], ['2100'])
# group_regions          ('26', ['east_coast'], ['high'], ['2100'])
# extract_connected     ('26', ['east_coast'], ['high'], ['2100'])
# raster_to_polygon     ('26', ['east_coast'], ['high'], ['2100'])



# inundated_areas                (flood_frequency,regions,projections,years)
# mosaic_chunks                  (flood_frequency,regions,projections,years,areas)  
# shift_rasters_from_scratch_gdb (flood_frequency,regions,projections,years)
## group_regions                 (flood_frequency,regions,projections,years)
# extract_connected              (flood_frequency,regions,projections,years)
# raster_to_polygon              (flood_frequency,regions,projections,years)

print('FINISHED ALL')
if showtimes==1: t2 = datetime.datetime.now() ; print('Time start: ' + str(t1) + '  &  Time end: ' + str(t2))
print('end')