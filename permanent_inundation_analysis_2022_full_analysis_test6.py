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

print ('STATUS 1: imports successful')

import datetime
t1 = datetime.datetime.now() ; #print('Time start: ' + str(t1))

# ----------------------------------------------------------------------------------------------------
# # SET INPUTS
# flood_frequency = '26'
# regions = ['west_coast']     #['east_coast','west_coast','pr_coast','vi_coast','gu_coast','hi_coast','as_coast']
# years = ['2060']     #['2022','2100','2060','2045','2035','2080','2050','2030','2065', ...]
# projections = ['int_low']  #['high','int_low']

# ----------------------------------------------------------------------------------------------------
# SET WORKSPACE
gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
arcpy.env.workspace = gdb

runpnw=1 #temp turn on/off WA pnw, rg onward
runscratch=0
runscratchout=1

# ----------------------------------------------------------------------------------------------------
# DEM ERROR FIX
def dem_fix(regions):
    exit() #turn off if re-fixing
    print('RUNNING SCRIPT FIXING DEM ERRORS')
    print('Setting all -9999 values in DEMS to NoData and resaving in new regional geodatabases')

    # set gdb to pull dems from, and set as workspace
    gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
    arcpy.env.workspace = gdb

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
    exit()
    for region in regions:
        print('RUNNING PREP MHHW LAYERS SCRIPT')
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

# **other : create points layers, generate near table, attach Station_ID field

# This requires the near table to be generated and Station_ID field to be added before running; should be in projected coordinate system
def prep_gauge_data_for_interpolation (region,projection):
            
    print('RUNNING PREP GAUGE DATA SCRIPT')
    print ('...running ' + str(region) + ' ' + str(projection) + '...')
    # tf1 = datetime.datetime.now() ; print('Time start: ' + str(tf1))

    # ----------------------------------------------------------------------------------------
    # set gdb, env
    gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
    arcpy.env.workspace = gdb

    # ----------------------------------------------------------------------------------------
    # Make layers from empty transect points file and file with gauges and water levels for each year/projection

    # gauges_layer = arcpy.MakeFeatureLayer_management("all_{0}_stations_26x_{1}_proj" .format(region, projections[0])) #old
    gauges_layer = arcpy.MakeFeatureLayer_management("all_{0}_stations_26x_{1}_proj" .format(region, projection))  #from csv from tg gauges analysis
    
    transect_points_layer = arcpy.MakeFeatureLayer_management('{0}_empty_points_proj_2022' .format(region),'transect_points_layer') #2022 refers to year of project, not analysis/projection
    
    # Join empty points and gauges file
    # arcpy.AddJoin_management('transect_points_layer','Near_FID', gauges_layer, 'OBJECTID','KEEP_COMMON')
    arcpy.AddJoin_management(transect_points_layer,'Station_ID', gauges_layer, 'Station_ID','KEEP_COMMON')

    print('Joined data to points for ' + str(projection) + ' projection')

    print ('STATUS: gauge data prep run successfully for ' + str(region) + ' ' + str(projection))
    # tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2)) 

    print('END PREP GAUGE DATA SCRIPT; RETURNING TO INTERPOLATE AND CREATE WLS SCRIPT')
    return transect_points_layer          

def interpolate_and_create_water_level_surfaces (flood_frequency,regions,projections,years):  #WORKING, double check in full run
    for region in regions:
        for projection in projections:

            print('RUNNING CREATE WATER LEVEL SURFACES SCRIPT')
            print ('...running ' + region + ' ' + projection + '...')
            tf1 = datetime.datetime.now() ; print('Time start: ' + str(tf1))

            # ----------------------------------------------------------------------------------------
            # set gdb in and out
            gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
            # set path for saving
            gdb_out = gdb

            # set workspace
            arcpy.env.workspace = gdb

            # ----------------------------------------------------------------------------------------
            # run prep functions - check this bit

            print('...running prep functions...')
            #prep_noaa_mhhw_layers(region)  #function not checked
            # transect_points_layer = prep_gauge_data_for_interpolation(region, projections) # add parameters or not?
            transect_points_layer = prep_gauge_data_for_interpolation(region, projection)

            runreadresults=0
            if runreadresults==1:
                # Comment out following lines:
                # number_of_points = arcpy.GetCount_management('transect_points_layer')
                # count = int(number_of_points.getOutput(0))
                # print(count)

                print('List of fields in transect points layer:')
                fields = arcpy.ListFields(transect_points_layer)
                for field in fields:
                    print(field.name)

            print('Prep functions completed')

            # ----------------------------------------------------------------------------------------
            # load regional mhhw data
            print('...loading regional MHHW water data...')
            if region=='west_coast' or region=='east_coast':
                mhhw_raster_string = 'all_noaa_mhhw_navd_mosaic'
                mhhw_poly_string = 'all_noaa_mhhw_navd_mosaic_polygon_dissolve_{0}' .format(region)
            else:
                print('ERROR: No MHHW layer for ' + region + ' region') ; exit()

            mhhw_data = arcpy.ListRasters(mhhw_raster_string)[0]
            # mhhw_data = arcpy.ListRasters("all_noaa_mhhw*")[0] # Put in config file for west and east
            # mhhw_layer = arcpy.MakeRasterLayer_management(mhhw_data,'mhhw_layer')
            mhhw_surface = Raster(mhhw_data)

            print('MHHW water surface data loaded')
            
            # Set interpolation limit as polygon created from regional coastal counties
            mhhw_mask = arcpy.ListFeatureClasses(mhhw_poly_string)[0]
            arcpy.env.mask = mhhw_mask

            #  Get cell size properties of mhhw_layer to use in interpolation between transect points
            cellsizes = []
            #for cellsize_property in cellsize_properties:
            # these lines should likely be a function
            size_tmp = arcpy.GetRasterProperties_management(mhhw_surface, "CELLSIZEX")
            size = size_tmp.getOutput(0)
            cellsizes.append(size)
            cellsizes = [float(item) for item in cellsizes]
            cellsize = cellsizes[0]
            print('Cell size of mhhw layer is: ' + str(cellsize)) ; #cellsize = cellsize*4
            
            # ----------------------------------------------------------------------------------------
            # For each year/proj horizon, interpolate empty points file based on field from joined table
            for year in years:

                print('Region is: ' + region + ' and projection is: ' + projection + ' and year is: ' + year)

                # Interpolate between points
                print('...creating interpolated surface...')
                # interpolated_surface_wrt_mhhw = NaturalNeighbor(transect_points_layer,str('all_{0}_stations_{1}x_{2}_proj.F{3}_{4}' .format(region,flood_frequency,projection, year,projection) ), cellsize) #old syntax
                interpolated_surface_wrt_mhhw = NaturalNeighbor(transect_points_layer,str('all_{0}_stations_{1}x_{2}_proj.{3}_{4}' .format(region,flood_frequency,projection,  projection,year) ), cellsize)
                print('Created interpolated surface')

                print('...saving interpolated surface...')
                outname_interpol = 'interpol_points_{0}_proj_{1}_{2}' .format(region,projection,year)
                interpolated_surface_wrt_mhhw.save(outname_interpol)
                print('Saved interpolated surface layer')

                # Add mhhw water surface to interpolated projection to get total water level wrt navd datum, in meters
                print('...creating water level surface...')
                water_level_surface =  Raster(mhhw_surface) + Raster(interpolated_surface_wrt_mhhw)  #GOOD
                print('Created water level surface')

                outname_wls = 'water_level_surface_wrt_navd_{0}x_{1}_{2}_{3}' .format(flood_frequency,year,projection,region)
                print('Outname is: ' + outname_wls)
                
                print('...saving water level surface...')
                water_level_surface.save(outname_wls)
                print('Saved water level surface layer')

                # outfind = arcpy.ListRasters(outname_wls)
                # print(len(outfind))
                
                print ('STATUS: interpolation and create water level surfaces run successfully for ' + region + ' ' + projection + ' ' + year)
                tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2)) 

    print('END CREATE WATER LEVEL SURFACES SCRIPT')


# **other

# ----------------------------------------------------------------------------------------------------
# INUNDATED AREAS               # 8 hours east
def inundated_areas (flood_frequency,regions,projections,years):
    for region in regions:
        for projection in projections:
            for year in years:
                
                print('RUNNING INUNDATED SURFACES SCRIPT')
                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                tf1 = datetime.datetime.now() ; print('Time start: ' + str(tf1))

                # ----------------------------------------------------------------------------------------
                # set gdb in and out
                gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                # set path for loading dems
                gdb_dem = 'C:/Users/annik/Documents/UCS/Data/2022/DEMs_fix/{0}/{0}_dem_fix.gdb' .format(region)  #
                # set path for saving inundated raster
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                if runscratch==1:
                    gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/scratch.gdb'

                # set workspace
                arcpy.env.workspace = gdb

                # ----------------------------------------------------------------------------------------
                # get the specific WLS from the gdb
                water_level = arcpy.ListRasters('water_level_surface_wrt_navd_26x_{0}_{1}_{2}'.format(year,projection,region))
                water_level_surface = Raster(water_level)
                print('Water level surface is: ' + str(water_level_surface))

                # ----------------------------------------------------------------------------------------
                # get dems
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

                # reset workspace
                arcpy.env.workspace = gdb_dem

                dems = arcpy.ListRasters('*')
                print(dems)

                num_elem = len(dems)
                print (str(num_elem) + ' inundation layers to create')

                # reset workspace
                arcpy.env.workspace = gdb

                # ----------------------------------------------------------------------------------------
                # compare dems to water level surface and output inundated raster
                for dem in dems:

                    print ('DEM is: ' + str(dem))

                    #set path for pulling in corrected dems
                    dem_path = str(gdb_dem) + '/' + str(dem) ; #print(dem_path)

                    #set path out
                    outname_inundated_area_surface = 'inundated_area_surface_{0}x_{1}_{2}_{3}'.format(flood_frequency, year, projection, region) + '_' + dem[:20]
                    outname_inundated_area_surface_full = str(gdb_out) + '/' + str(outname_inundated_area_surface)

                    #check if inundation file already exists, create if not
                    arcpy.env.workspace = gdb_out
                    inundated_check = arcpy.ListRasters(outname_inundated_area_surface) ; #print(len(inundated_check))
                    arcpy.env.workspace = gdb
                    if len(inundated_check)==0:  #if not already created

                        # set mask to dem extent
                        arcpy.env.mask = dem_path
                        arcpy.env.cellSize = "MINOF"

                        # Compare DEM and water level surface. Where WLS >= DEM, set value of 1 ; i.e. water above land is TRUE
                        print ('...creating inundated area surface...')
                        inundated_area_surface = Con(Raster(water_level_surface) >= Raster(dem_path), 1) # creates flat inundation area raster
                        # print ('Created inundated area surface for ' + str(dem))

                        # save raster
                        print('...saving layer...')
                        inundated_area_surface.save(outname_inundated_area_surface_full)
                        # print('Saved layer')

                    print ('Created inundated area surface: ' + str(outname_inundated_area_surface))

                print ('STATUS: inundated areas run successfully for ' + str(projection) + ' ' + str(year))
                tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))       
    
    print('END INUNDATED SURFACES SCRIPT')

def inundated_areas_select_dems (flood_frequency,regions,projections,years,dems):  #for individual dems
    for region in regions:
        for projection in projections:
            for year in years:
                    
                print('RUNNING INUNDATED SURFACES (SELECT DEMS) SCRIPT')
                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                tf1 = datetime.datetime.now() ; print('Time start: ' + str(tf1))

                # ----------------------------------------------------------------------------------------
                # set gdb in and out
                gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                # set path for loading dems
                gdb_dem = 'C:/Users/annik/Documents/UCS/Data/2022/DEMs_fix/{0}/{0}_dem_fix.gdb' .format(region)  #
                # set path for saving inundated raster
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

                # set workspace
                arcpy.env.workspace = gdb

                # ----------------------------------------------------------------------------------------
                # get the specific WLS from the gdb
                water_level = arcpy.ListRasters('water_level_surface_wrt_navd_26x_{0}_{1}_{2}'.format(year,projection,region))
                water_level_surface = Raster(water_level)
                print('Water level surface is: ' + str(water_level_surface))

                # ----------------------------------------------------------------------------------------
                # get dems
                # rundemmanual=0
                # if rundemmanual==1:
                #     west_dems = ['CA_ChannelIslands_GCS_5m_NAVDm','CA_EKA1_GCS_5m_NAVD88m','CA_EKA2_GCS_5m_NAVD88m','CA_LOX1_GCS_5m_NAVD88m','CA_LOX2_GCS_5m_NAVD88m','CA_MTR1_GCS_5m_NAVD88m','CA_MTR2_GCS_5m_NAVD88m','CA_MTR3_GCS_5m_NAVD88m','CA_SGX_GCS_5m_NAVD88m','OR_MFR_GCS_5m_NAVD88m','OR_PQR1_GCS_5m_NAVD88m_TillamookFix','OR_PQR2_GCS_5m_NAVD88m','WA_PQR_GCS_5m_NAVD88m','WA_PugetSoundNW_GCS_3m_NAVDm','WA_PugetSoundSW_GCS_3m_NAVDm','WA_SEW1_GCS_5m_NAVD88m','WA_SEW1b_GCS_5m_NAVD88m','WA_SEW2_GCS_5m_NAVDm']
                #     east_dems = ['AL_GCS_3m_NAVDm','CT_GCS_3m_NAVDm','DC_GCS_3m_NAVDm','DE_GCS_3m_NAVDm','FL_A3_GCS_5m_NAVDm','FL_JAX_1_GCS_5m_NAVD88m','FL_JAX_2_GCS_5m_NAVD88m','FL_MFL_1_GCS_5m_NAVD88m','FL_MFL_2_GCS_5m_NAVD88m','FL_MLB_1_GCS_5m_NAVD88m','FL_MLB_2_GCS_5m_NAVD88m','FL_Pan_East_GCS_3m_NAVDm','FL_Pan_West_GCS_3m_NAVDm','FL_TBW_1_GCS_5m_NAVD88m','FL_TBW_2_GCS_5m_NAVD88m','GA_CHS_GCS_5m_NAVDm','GA_JAX_GCS_5m_NAVDm','LA_Central_GCS_3m_NAVDm','LA_CentralEast_GCS_3m_NAVDm','LA_CentralNorth_GCS_3m_NAVDm','LA_Delta_GCS_3m_NAVDm','LA_LakePontchartrain_GCS_3m_NAVDm','LA_West_GCS_3m_NAVDm','MA_GCS_3m_NAVDm','MD_East_GCS_3m_NAVDm','MD_North_GCS_3m_NAVDm','MD_Southeast_GCS_3m_NAVDm','MD_Southwest_GCS_3m_NAVDm','MD_West_GCS_3m_NAVDm','ME_East_GCS_3m_NAVDm','ME_West_GCS_3m_NAVDm','MS_GCS_3m_NAVDm','NC_Middle1_GCS_3m_NAVDm','NC_Middle2_GCS_3m_NAVDm','NC_Northern_GCS_3m_NAVDm','NC_Southern1_GCS_3m_NAVDm','NC_Southern2_GCS_3m_NAVDm','NH_GCS_3m_NAVDm','NJ_Middle_GCS_3m_NAVDm_FY21Update','NJ_Northern_GCS_3m_NAVDm','NJ_Southern_GCS_3m_NAVDm_FY21Update','NY_Hudson_GCS_3m_NAVDm','NY_Metro_GCS_3m_NAVDm','NY_Suffolk_GCS_3m_NAVDm','PA_GCS_3m_NAVDm','RI_GCS_3m_NAVDm','SC_Central_GCS_3m_NAVDm','SC_North_GCS_3m_NAVDm','SC_South_GCS_3m_NAVDm','TX_Central_GCS_3m_NAVDm','TX_North1_GCS_3m_NAVDm','TX_North2_GCS_3m_NAVDm','TX_South1_GCS_3m_NAVDm','TX_South2_GCS_3m_NAVDm','VA_EasternShore_GCS_3m_NAVDm','VA_Middle_GCS_3m_NAVDm','VA_Northern_GCS_3m_NAVDm','VA_Southern_GCS_3m_NAVDm']
                #     hi_dems = ['HI_HFO_Hawaii_UTM_3m_LMSLm','HI_HFO_Kauai_UTM_3m_LMSLm','HI_HFO_Lanai_UTM_3m_LMSLm','HI_HFO_Maui_UTM_3m_LMSLm','HI_HFO_Molokai_UTM_3m_LMSLm','HI_Oahu_GCS_3m_LMSLm']
                #     gu_dems = ['Guam_GCS_3m_GUVD04m']
                #     as_dems = ['AS_OfuOlosega_DEM_v2_UTM','AS_Tau_DEM_v2_UTM','AS_Tutuila_DEM_v2_UTM']
                #     pr_dems = ['PR_GCS_3m_PRVD02m']
                #     vi_dems = ['USVI_GCS_3m_NAVDm']

                # reset workspace
                # arcpy.env.workspace = gdb_dem

                # dems = arcpy.ListRasters('*')
                # dems = ['NJ_Middle_GCS_3m_NAVDm_FY21Update'] #running single dem
                print(dems)

                num_elem = len(dems)
                print (str(num_elem) + ' inundation layers to create')

                # ----------------------------------------------------------------------------------------
                # compare dems to water level surface and output inundated raster
                for dem in dems:

                    print ('DEM is: ' + str(dem))

                    #set path for pulling in corrected dems
                    dem_path = str(gdb_dem) + '/' + str(dem) ; #print(dem_path)

                    #set path out
                    outname_inundated_area_surface = 'inundated_area_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region) + '_' + dem[:20]
                    outname_inundated_area_surface_full = str(gdb_out) + '/' + str(outname_inundated_area_surface)
                    print('Outname is: ' + str(outname_inundated_area_surface_full))

                    #delete bad inundation files
                    arcpy.env.workspace = gdb_out
                    delweird=1
                    if delweird==1:
                        inund_check_TODEL = arcpy.ListRasters(outname_inundated_area_surface)
                        if len(inund_check_TODEL)!=0:
                            for del_file in inund_check_TODEL:
                                print('...deleting old inundation file...')
                                arcpy.Delete_management(del_file)
                                print('Deleted old inundation file: ' + str(del_file))
                        elif len(inund_check_TODEL)==0:
                            print('No inundation file yet exists')

                    #check if inundation file already exists, create if not - condition turned off
                    # inundated_check = arcpy.ListRasters(outname_inundated_area_surface) ; print(len(inundated_check))
                    arcpy.env.workspace = gdb
                    if 1==1: 

                        # set mask to dem extent
                        arcpy.env.mask = dem_path
                        arcpy.env.cellSize = "MINOF"

                        # Compare DEM and water level surface. Where WLS >= DEM, set value of 1 ; i.e. water above land is TRUE
                        print ('...creating inundated area surface...')
                        inundated_area_surface = Con(Raster(water_level_surface) >= Raster(dem_path), 1) # creates flat inundation area raster
                        print ('Created inundated area surface for ' + str(dem))

                        # save raster
                        print('...saving layer...')
                        inundated_area_surface.save(outname_inundated_area_surface_full)
                        print('Saved layer')

                    print ('Created inundated area surface: ' + str(outname_inundated_area_surface))

                print ('STATUS: inundated areas run successfully for ' + str(projection) + ' ' + str(year))
                tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))       
    
    print('END INUNDATED SURFACES (SELECT DEMS) SCRIPT')

# ----------------------------------------------------------------------------------------
# MOSAIC                        #5 hours east -- 11 ?
def mosaic_chunks (flood_frequency,regions,projections,years,areas):
    for region in regions:
        for projection in projections:
            for year in years:
                
                print('RUNNING MOSAIC INUNDATED SURFACES SCRIPT')
                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                tf1 = datetime.datetime.now() ; #print('Time start: ' + str(tf1))

                # ----------------------------------------------------------------------------------------
                # set gdb in
                gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                if runscratch==1:
                    gdb_scratch = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/scratch.gdb'
                    gdb_in = gdb_scratch

                # set workspace
                arcpy.env.workspace = gdb_in

                # ----------------------------------------------------------------------------------------
                # set path out
                # gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                gdb_out = gdb_in
                if runscratchout==1:
                    gdb_scratch = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/scratch.gdb'
                    gdb_out = gdb_scratch

                # ----------------------------------------------------------------------------------------
                # resample in certain regions due to processing speed, then mosaic

                if region == 'east_coast':
                    ##NOTES: resample 3 min per element
                    ##NOTES: full resample and mosaic, east coast ~5 hrs to run (pre-area break)
                    
                    # set cellsize based on speed and quality testing
                    cellsize_out = "0.000045"
                    print('Cellsize out is ' + str(cellsize_out))

                    #name areas if set to 'all' or if parameter empty
                    if (len(areas)==1 and areas[0]=='all') or (len(areas)==0):
                        areas=['gulf_area','atlantic_area']
                    
                    #standardize area names
                    set_areas = []
                    for area in areas:
                        if area=='gulf':
                            set_areas.append('gulf_area')
                        elif area=='atlantic':
                            set_areas.append('atlantic_area')
                    if len(set_areas)>0:
                        areas = set_areas
                    
                    # ----------------------------------------------------------------------------------
                    for area in areas:
                        print('...running ' + area + '...')
                        tta = datetime.datetime.now()
                        
                        #set states to areas
                        if area=='gulf_area':
                            area_set = ['TX','LA','MS','AL','FL','GA_JAX'] #gulf
                        elif area=='atlantic_area':
                            area_set = ['FL_JAX_1','GA','SC','NC','VA','DC','MD','DE','PA','NJ','NY','CT','RI','MA','NH','ME'] #atlantic

                        # get inundation raster surfaces raster list
                        raster_inund_files = []
                        for state in area_set:
                            rasters = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}_{4}*' .format(flood_frequency,year,projection,region,state))
                            for file in rasters:
                                raster_inund_files.append(file)

                        print('Files to resample and mosaic:')        
                        for file in raster_inund_files:
                            print(file)
                        
                        #check area files
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
                        # resample inundation rasters and save        #NOTES 3-5 min each raster !!!
                        list_resampled_rast_paths = []
                        for raster_inund_file in raster_inund_files:
                            #set path for pulling in inundation raster
                            rast_path = str(gdb_in) + '/' + str(raster_inund_file)
                            
                            #set resampled raster file name
                            outname_resample = 'scratch_' + str(raster_inund_file) + '_resamp_temp'  #cellsz_pt000045
                            outname_resample_path = str(gdb_out) + '/' + str(outname_resample)  

                            #delete old resampled files
                            delweird=1
                            if delweird==1:
                                resample_check_TODEL = arcpy.ListRasters('scratch_' + str(raster_inund_file) + '_resamp*')
                                if len(resample_check_TODEL)!=0:
                                    for del_file in resample_check_TODEL:
                                        print('...deleting old resampled file...')
                                        arcpy.Delete_management(del_file)
                                        print('Deleted old resampled file: ' + str(del_file))

                            #resample file, if resampled file doesn't already exist
                            resample_check = arcpy.ListRasters(outname_resample)
                            if len(resample_check)==0:  #if not already resampled
                                print ('File to resample: ' + str(raster_inund_file))

                                print ('Resample outname is: ' + outname_resample)

                                # getcellsize_temp = arcpy.GetRasterProperties_management((rast_test),"CELLSIZEX")
                                # cellsize_temp = getcellsize_temp.getOutput(0)
                                # print("Original cellsize is: " + str(cellsize_temp))

                                #resample
                                print('...resampling...')
                                rast_resamp = arcpy.Resample_management(rast_path,outname_resample_path,cellsize_out,"NEAREST")  #NEAREST does not change cell values
                            print('Resampled ' + str(outname_resample))
                            
                            #make new set of resampled raster path names to mosaic
                            list_resampled_rast_paths.append(outname_resample_path)
                        
                        ttb = datetime.datetime.now() ; print('Area resampling ' + area + ' Time start: ' + str(tta) + '  &  Time end: ' + str(ttb))
                    
                        # ----------------------------------------------------------------------------------
                        # combine chunks and create mosaic, separated by area
                        print ('...resampling done, setting up mosaic...')

                        # for file in list_resampled_rast_paths:
                        #     print(file)

                        #check number of files to be mosaiced
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

                        tta = datetime.datetime.now() ; print('Start time: ' + str(tta))

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
                        
                        ttb = datetime.datetime.now() ; print('Area mosaic ' + area + ' Time start: ' + str(tta) + '  &  Time end: ' + str(ttb))
                        
                        print('Mosaicing ' + str(region) + ' ' + str(projection) + ' ' + str(year) + ' area ' + area + ' complete')
                    
                    print('Mosaicing complete for ' + str(region) + ' ' + str(projection) + ' ' + str(year) + ' ' + str(areas))

                elif region == 'west_coast':
                    runareas=0
                    if runareas==1:
                        #name areas if set to 'all' or if parameter empty
                        if len(areas)==1 and areas[0]=='all' or (len(areas)==0):
                            areas = ['pnw_area','pacific_area']

                        for area in areas:
                            print('...running ' + area + '...')
                                                    
                            #set states to areas
                            if area=='pnw_area':
                                area_set = ['WA_PugetSoundNW','WA_PugetSoundSW','WA_SEW2'] #pacific northwest, pnw
                            elif area=='pacific_area':
                                area_set = ['WA_SEW1','WA_SEW1b','WA_PQR','OR','CA'] #pacific

                            # get inundation raster surfaces raster list
                            raster_inund_files = []
                            for state in area_set:
                                rasters = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}_{4}*' .format(flood_frequency,year,projection,region,state))
                                for file in rasters:
                                    raster_inund_files.append(file)

                            print('Files to resample and mosaic:')        
                            for file in raster_inund_files:
                                print(file)
                            
                            # create list with full paths of inundation raster surfaces rasters
                            list_rast_paths = []
                            for raster_inund_file in raster_inund_files:
                                #set path for pulling in inundation rasters
                                rast_path = str(gdb_in) + '/' + str(raster_inund_file)
                                list_rast_paths.append(rast_path)
                            print(list_rast_paths)

                            #check number of files to be mosaiced
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
                        
                            tta = datetime.datetime.now() ; print('Start time: ' + str(tta))    

                            outname_merged = 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency, year, projection, region, area)  #change test nums
                            print ('Outname is: ' + outname_merged)

                            # combine chunks and create mosaic
                            print ('...creating ' + area + ' mosaic...')
                            arcpy.MosaicToNewRaster_management(list_rast_paths, gdb_out, outname_merged,"","","",1,"MAXIMUM") 
                            print('Created mosaic for ' + area)
                                                       
                            ttb = datetime.datetime.now() ; print('Area mosaic ' + area + ' Time start: ' + str(tta) + '  &  Time end: ' + str(ttb))
                            print('Mosaicing area ' + area + ' complete')
                        
                        print('Mosaicing complete for ' + str(region) + ' ' + str(projection) + ' ' + str(year) + ' ' + str(areas))
                    
                    elif runareas==0:
                        # create list with full paths of inundation raster surfaces rasters
                        list_rast_paths = []
                        for raster_inund_file in raster_inund_files:
                            #set path for pulling in inundation rasters
                            rast_path = str(gdb_in) + '/' + str(raster_inund_file)
                            list_rast_paths.append(rast_path)
                        print(list_rast_paths)

                        outname_merged = 'merged_raw_raster_surface_%sx_%s_%s_%s' %(flood_frequency, year, projection, region)  #change test nums
                        print ('Outname is: ' + outname_merged)

                        # combine chunks and create mosaic
                        print ('...creating mosaic...')                    
                        arcpy.MosaicToNewRaster_management(list_rast_paths, gdb_out, outname_merged,"","","",1,"MAXIMUM") 
                        print('Created mosaic')     

                else:
                    print('region not built for mosaicing!')
                    exit()

                print ('STATUS: mosaic run successfully for ' + str(projection) + ' ' + str(year))
                tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))
    ## NOTES: original script had cellSize=10, bandsNum=1.  This script and ArcGIS Pro tool give errors with cellSize=10.
    #   Cellsize "0.000045" selected based on cell size of existing layers, speed testing, and quality control

    print('END MOSAIC INUNDATED SURFACES SCRIPT')

def shift_rasters_from_scratch_gdb (flood_frequency,regions,projections,years):
    for region in regions:
        for projection in projections:
            for year in years:
                if runscratchout==1:
                    print('RUNNING SCRIPT SHIFTING INUNDATION AND MOSAIC RASTERS FROM SCRATCH TO MAIN GDB')
                    print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                    tf1 = datetime.datetime.now() ; #print('Time start: ' + str(tf1))

                    # ----------------------------------------------------------------------------------------
                    # set gdb in and out
                    gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                    gdb_scratch = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/scratch.gdb'
                    gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)                   

                    # set workspace
                    arcpy.env.workspace = gdb_scratch
                    
                    # ----------------------------------------------------------------------------------------
                    # grab rasters in scratch gdb and move to out gdb                
                    rasters_to_shift = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}*' .format(flood_frequency,year,projection,region))
                    
                    for raster in rasters_to_shift:
                        print('Raster to shift: ' + str(raster))

                        print('...shifting...')
                        arcpy.conversion.RasterToGeodatabase(raster, gdb_out, '')
                        #print('Complete')
                    print('Inundation files shifted from scratch gdb')

                    for del_file in rasters_to_shift:
                        arcpy.Delete_management(del_file)
                    print('Deleted inundation files from scratch gdb')

                    rasters_to_shift = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}*' .format(flood_frequency,year,projection,region))
                    
                    for raster in rasters_to_shift:
                        print('Raster to shift: ' + str(raster))

                        print('...shifting...')
                        arcpy.conversion.RasterToGeodatabase(raster, gdb_out, '')
                        #print('Complete')
                    print('Mosaic files shifted from scratch gdb')

                    for del_file in rasters_to_shift:
                        arcpy.Delete_management(del_file)
                    print('Deleted mosaic files from scratch gdb')

                    print('STATUS: rasters shifted successfully')
                    tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))
    print('END SHIFT RASTERS FROM SCRATCH GDB SCRIPT')


# ----------------------------------------------------------------------------------------
# REGION GROUP                  #3 - 4:30 hours east
def region_group (flood_frequency,regions,projections,years):
    for region in regions:
        for projection in projections:
            for year in years:

                print('RUNNING REGION GROUPING SCRIPT FOR CONNECTIVITY')
                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                tf1 = datetime.datetime.now() ; #print('Time start: ' + str(tf1))

                # ----------------------------------------------------------------------------------------
                # set gdb in and out
                gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

                # set workspace
                arcpy.env.workspace = gdb_in

                # ----------------------------------------------------------------------------------------
                # region group
                if region=='west_coast':
                    print('Region grouping for both mainland west coast and PNW area')

                    # Perform region group for each area: mainland and pnw
                    areas = ['mainland','pnw']
                    for area in areas:
                        if area=='mainland':
                            # get merged raster
                            merged_raw_surface = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region))[0]
                            ## print(str(merged_raw_surface))

                            # ----------------------------------------------------------------------------------------
                            # Perform region group for mainland
                            print('Mainland raw surface name is: ' + str(merged_raw_surface))
                            surface_path = gdb_in + '/' +  str(merged_raw_surface) ; #print(surface_path)
                            ##filename = os.path.basename(fullname)

                            fullname = str(merged_raw_surface)
                            runrg = 1
                            if runrg == 1:
                                print('...region grouping ' + fullname + '...')
                                outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                                outname_rg = 'rg_' + fullname
                                outRegionGrp.save(outname_rg)
                                print('Region grouped ' + merged_raw_surface)  
                        # ----------------------------------------------------------------------------------------
                        elif area=='pnw':
                            if runpnw==1:
                                    
                                # get merged raster
                                merged_raw_surface = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW_mosaicAll' .format(flood_frequency, year, projection, region))[0]
                                ## print(str(merged_raw_surface))

                                # ----------------------------------------------------------------------------------------
                                # Perform region group for PNW
                                print('PNW raw surface name is: ' + str(merged_raw_surface))
                                surface_path_pnw = gdb_in + '/' +  str(merged_raw_surface) ; #print(surface_path_pnw)
                                ##filename = os.path.basename(fullname)

                                fullname = str(merged_raw_surface)
                                runrg = 1                                                             ##TURN BACK ON
                                if runrg == 1:
                                    print('...region grouping ' + fullname + '...')
                                    outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                                    outname_rg = 'rg_' + 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW' .format(flood_frequency, year, projection, region)
                                    outRegionGrp.save(outname_rg)
                                    print('Region grouped ' + merged_raw_surface)  
                            
                    print ('STATUS: region group run successfully for mainland & PNW ' + str(region) + ' ' + str(projection) + ' ' + str(year))
                        
                elif region=='east_coast':
                    print('Region grouping for three east coast areas: Gulf Coast, Atlantic Coast, and RI Block Island')
                    # Perform region group for each area: mainland Gulf+Florida, Atlantic, and RI island
                    # Gulf and Atlantic have overlapping DEM area to allow for proper connectivity assessment at edge
                    areas = ['gulf_area','atlantic_area','RI_block_isle']
                    for area in areas:
                        if area=='gulf_area' or area=='atlantic_area':

                            # get merged raster
                            merged_raw_surface = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area))[0]
                            ## print(str(merged_raw_surface))

                            # ----------------------------------------------------------------------------------------
                            # Perform region group
                            print(area + ' raw surface name is: ' + str(merged_raw_surface))
                            surface_path = gdb_in + '/' +  str(merged_raw_surface) ; #print(surface_path)
                            ##filename = os.path.basename(fullname)
                            fullname = str(merged_raw_surface)

                            #check that rg doesn't already exist
                            check_rg = arcpy.ListRasters('rg_' + fullname)
                            if len(check_rg)>0:
                                print('Region group already performed for ' + area)
                            elif len(check_rg)==0:
                                runrg = 1
                                if runrg == 1:
                                    print('...region grouping ' + fullname + '...')
                                    outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                                    outname_rg = 'rg_' + fullname
                                    outRegionGrp.save(outname_rg)
                                    print('Region grouped ' + str(merged_raw_surface))  
                        # ----------------------------------------------------------------------------------------
                        elif area=='RI_block_isle':

                            # get RI inundation raster, NOT merged raster
                            raw_surface = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}_RI*' .format(flood_frequency,year,projection,region))[0]
                            ## print(str(raw_surface))

                            # ----------------------------------------------------------------------------------------
                            # Perform region group
                            print(area + ' raw surface name is: ' + str(raw_surface))
                            surface_path_riisle = gdb_in + '/' +  str(raw_surface) ; #print(surface_path_riisle)
                            ##filename = os.path.basename(fullname)
                            fullname = str(raw_surface)

                            #check that rg doesn't already exist
                            check_rg = arcpy.ListRasters('rg_' + fullname)
                            if len(check_rg)>0:
                                print('Region group already performed for ' + area)
                            elif len(check_rg)==0:
                                runrg = 1                                                             ##TURN ON
                                if runrg == 1:
                                    print('...region grouping ' + fullname + '...')
                                    outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                                    outname_rg = 'rg_' + 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area)
                                    outRegionGrp.save(outname_rg)
                                    print('Region grouped ' + str(raw_surface)) 
                    
                    print ('STATUS: region group run successfully for ' + str(projection) + ' ' + str(year))
                
                # other regions
                else:
                    print ('Region not built yet!') ; exit()
                    print ('STATUS: region group run successfully for ' + str(projection) + ' ' + str(year))

                tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))

    print('END REGION GROUPING SCRIPT FOR CONNECTIVITY')

# ----------------------------------------------------------------------------------------
# EXTRACT CONNECTED REGIONS     #3 hours east
def extract_connected (flood_frequency,regions,projections,years):
    for region in regions:
        for projection in projections:
            for year in years:

                print('RUNNING EXTRACTING CONNECTED REGION GROUPS SCRIPT')
                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                tf1 = datetime.datetime.now() ; #print('Time start: ' + str(tf1))

                # ----------------------------------------------------------------------------------------
                # set gdb in and out
                gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

                # set workspace
                arcpy.env.workspace = gdb_in

                # ----------------------------------------------------------------------------------------
                # extract
                if region=='west_coast':
                    print('Extracting both mainland west coast and PNW area')

                    # ----------------------------------------------------------------------------------------
                    # Perform region group for each area: mainland and pnw
                    areas = ['Mainland','PNW']
                    for area in areas:
                        if area=='Mainland':
                            print('Extracting west coast' + str.lower(area))

                            # get region group raster
                            rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region))[0]
                            rg_surface_path = gdb_in + '/' + rg_surface
                            print(area + ' file to extract is ' + rg_surface)

                            #set outname
                            outname_extract = 'extract_' + rg_surface

                            # ----------------------------------------------------------------------------------------
                            # set areas for extraction - mainland (maximum[0]) and Channel Isles (maximum[2])
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
                            runextract=1
                            if runextract==1:
                                attExtract = ExtractByAttributes(rg_surface, inSQLClause)
                                attExtract.save(outname_extract)
                                print('Extracted connected areas from ' + rg_surface)
                        # ----------------------------------------------------------------------------------------               
                        elif area=='PNW':
                            if runpnw==1:
                                print('Extracting west coast' + str(area))

                                # get region group raster
                                rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW' .format(flood_frequency, year, projection, region))[0]
                                rg_surface_path = gdb_in + '/' + rg_surface
                                print(area + ' file to extract is ' + rg_surface)

                                #set outname
                                outname_extract = 'extract_' + rg_surface
                                
                                # ----------------------------------------------------------------------------------------
                                # set groups for extraction - PNW ONLY
                                print('...finding groups to extract from ' + rg_surface + '...')

                                # index maximum count
                                arr = arcpy.da.FeatureClassToNumPyArray(rg_surface, ('Value', 'Count'))
                                count = arr['Count']
                                value = arr['Value']
                                index_to_extract = numpy.argmax(count)
                                value_to_extract = str(value[index_to_extract])

                                # build SQL clause syntax
                                print('...building SQL string...')
                                inSQLClause = 'Value = ' + value_to_extract

                                print('Extracting value=' + value_to_extract + ' from ' + rg_surface)
                                print('...extracting values...')

                                # extract connected areas
                                runextract=1
                                if runextract==1:
                                    attExtract = ExtractByAttributes(rg_surface, inSQLClause)
                                    attExtract.save(outname_extract)
                                    print('Extracted connected areas from ' + rg_surface)
                    print ('STATUS: extract connected regions run successfully for mainland & PNW ' + str(region) + ' ' + str(projection) + ' ' + str(year))
                    
                elif region=='east_coast':
                    print('Extracting east coast: Gulf, Atlantic, and RI island area')

                    # ----------------------------------------------------------------------------------------
                    # Perform region group for each area: gulf coast, atlantic coast, and RI isle
                    areas = ['gulf_area','atlantic_area','RI_block_isle']
                    for area in areas:
                        if area=='gulf_area' or area=='atlantic_area':
                            print('Extracting east coast ' + area)

                            # get region group raster
                            rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area))[0]
                            rg_surface_path = gdb_in + '/' + rg_surface
                            print(area + ' file to extract is ' + rg_surface)

                            #set outname
                            outname_extract = 'extract_' + rg_surface
                            
                            # ----------------------------------------------------------------------------------------
                            #check that extract doesn't already exist
                            check_extract = arcpy.ListRasters('outname_extract')
                            if len(check_extract)!=0:
                                print('Extraction already performed for ' + area)
                            elif len(check_extract)==0:
                                # ----------------------------------------------------------------------------------------
                                # set groups for extraction
                                print('...finding groups to extract from ' + rg_surface + '...')

                                # index maximum count
                                arr = arcpy.da.FeatureClassToNumPyArray(rg_surface, ('Value', 'Count'))
                                count = arr['Count']
                                value = arr['Value']
                                index_to_extract = numpy.argmax(count)
                                value_to_extract = str(value[index_to_extract])

                                # build SQL clause syntax
                                print('...building SQL string...')
                                inSQLClause = 'Value = ' + value_to_extract

                                print('Extracting value=' + value_to_extract + ' from ' + rg_surface)
                                print('...extracting values...')

                                # extract connected areas
                                runextract=1
                                if runextract==1:
                                    attExtract = ExtractByAttributes(rg_surface, inSQLClause)
                                    attExtract.save(outname_extract)
                                    print('Extracted connected areas from ' + rg_surface)
                        # ----------------------------------------------------------------------------------------                                       
                        elif area=='RI_block_isle':
                            print('Extracting east coast ' + area)

                            # get region group raster
                            rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area))[0]
                            rg_surface_path = gdb_in + '/' + rg_surface
                            print(area + ' file to extract is ' + rg_surface)

                            #set outname
                            outname_extract = 'extract_' + rg_surface
                        
                            # ----------------------------------------------------------------------------------------
                            #check that extract doesn't already exist
                            check_extract = arcpy.ListRasters('outname_extract')
                            if len(check_extract)!=0:
                                print('Extraction already performed for ' + area)
                            elif len(check_extract)==0:
                                # ----------------------------------------------------------------------------------------
                                # set areas for extraction - RI Block Island (maximum[1]) in RI inundation layer
                                print('...finding areas to extract from ' + rg_surface + '...')
                                
                                # index maximum counts for Atlantic mainland and RI Block Island
                                arr = arcpy.da.FeatureClassToNumPyArray(rg_surface, ('Value', 'Count'))
                                count = arr['Count']
                                value = arr['Value']
                                
                                index_to_extract_imax2 = numpy.argsort(count)[-2]           #;print('imax2 = ' + str(index_to_extract_imax2))
                                value_to_extract_imax2 = str(value[index_to_extract_imax2]) #;print('val2 = ' + str(value_to_extract_imax2))

                                print('   Result: second maximum pixel count corresponds to value ' + value_to_extract_imax2 + ' (RI island)')

                                # build SQL clause syntax
                                print('...building SQL string...')
                                inSQLClause = 'Value = ' + value_to_extract_imax2  #values from regional conditions

                                print('Extracting value=' + value_to_extract_imax2 + ' from ' + rg_surface)
                                print('...extracting values...')

                                # extract connected areas       RI isle ~1 min
                                runextract=1
                                if runextract==1:
                                    attExtract = ExtractByAttributes(rg_surface, inSQLClause)
                                    attExtract.save(outname_extract)
                                    print('Extracted connected areas from ' + rg_surface)
                
                    print ('STATUS: extract connected regions run successfully for ' + str(region) + ' ' + str(projection) + ' ' + str(year))
                
                else:
                    print('region not yet built') ; exit()
                    print ('STATUS: extract connected regions run successfully for _')

                tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))

    print('END EXTRACTING CONNECTED REGION GROUPS SCRIPT')

# ----------------------------------------------------------------------------------------
# CONVERT RASTER TO POLYGON     #1:30 - 2 hours east
def raster_to_polygon (flood_frequency,regions,projections,years):
    for region in regions:
        for projection in projections:
            for year in years:

                print('RUNNING RASTER TO POLYGON SCRIPT')
                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')
                tf1 = datetime.datetime.now() ; #print('Time start: ' + str(tf1))

                # ----------------------------------------------------------------------------------------
                # set gdb in and out
                gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

                # set workspace
                arcpy.env.workspace = gdb_in

                # ----------------------------------------------------------------------------------------
                # raster to polygon
                if region == 'west_coast':
                    print('Building polygons for both mainland west coast and PNW area')

                    # ----------------------------------------------------------------------------------------
                    # Perform polygon conversion for each area: mainland and pnw
                    areas = ['Mainland','PNW']
                    for area in areas:
                        if area=='Mainland':
                            # get extracted rasters
                            to_convert = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency,year,projection,region))[0]
                            print(area + ' file to convert is ' + str(to_convert))

                            # set path and file out
                            # gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                            outname_polygon = 'polygon_' + str(to_convert)  #final_polygon_

                            # ----------------------------------------------------------------------------------------
                            # create polygons       4 min total
                            runpolygons=1
                            if runpolygons==1:
                                print('...creating ' + str.lower(area) + ' polygon...')
                                arcpy.RasterToPolygon_conversion(to_convert, outname_polygon, "SIMPLIFY", "VALUE")
                                print ('Converted ' + str(to_convert) + ' to polygon')  

                        elif area=='PNW':
                            if runpnw==1:

                                # get extracted rasters
                                to_convert = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW' .format(flood_frequency,year,projection,region))[0]
                                print('PNW file to convert is ' + str(to_convert))

                                # set path and file out
                                # gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                                outname_polygon = 'polygon_' + str(to_convert)

                                # ----------------------------------------------------------------------------------------
                                # create polygons       4 min total
                                runpolygons=1
                                if runpolygons==1:
                                    print('...creating PNW polygon...')
                                    arcpy.RasterToPolygon_conversion(to_convert, outname_polygon, "SIMPLIFY", "VALUE")
                                    print ('Converted ' + str(to_convert) + ' to polygon')  

                    print ('STATUS: polygon conversion run successfully for mainland and PNW west coast, ' + str(projection) + ' ' + str(year))
                    
                elif region=='east_coast':
                    print('Building polygons for east coast: Gulf, Atlantic, and RI island area')

                    # ----------------------------------------------------------------------------------------
                    # Perform polygon conversion for each area: mainland and RI isle
                    areas = ['gulf_area','atlantic_area','RI_block_isle']
                    # areas = ['RI_block_isle']

                    for area in areas:
                        if area=='gulf_area' or area=='atlantic_area' or area=='RI_block_isle':
                            print('Building polygons for east coast ' + area)

                            # get extracted rasters
                            to_convert = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area))[0]
                            print(area + ' file to convert is ' + str(to_convert))

                            # set path and file out
                            # gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                            outname_polygon = 'polygon_' + str(to_convert)  #final_polygon_
                            print('Outname is: ' + outname_polygon)
                           
                            # ----------------------------------------------------------------------------------------
                            #check that extract doesn't already exist
                            check_poly = arcpy.ListFeatureClasses(outname_polygon)
                            if len(check_poly)!=0:
                                print('Polygons already performed for ' + area)
                            elif len(check_poly)==0:
                                # ----------------------------------------------------------------------------------------
                                # create polygons       4 min total ?
                                runpolygons=1
                                if runpolygons==1:
                                    print('...creating ' + area + ' polygon...')
                                    arcpy.RasterToPolygon_conversion(to_convert, outname_polygon, "SIMPLIFY", "VALUE")
                                    print ('Converted ' + str(to_convert) + ' to polygon')                    

                    print ('STATUS: polygon conversion run successfully for east coast areas ' + str(projection) + ' ' + str(year))

                else:
                    print('region not yet built!')

                print('FINISHED ' + str(region) + ' ' + str(projection) + ' ' + str(year))
                tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))

    print('END RASTER TO POLYGON SCRIPT')


# ----------------------------------------------------------------------------------------
#SET FUNCTION INPUTS
flood_frequency = '26'
regions = ['pr_coast']    #['east_coast','west_coast','pr_coast','vi_coast'] #['gu_coast','hi_coast','as_coast']
projections = ['high']       #['high','int_low']
years = ['2022']          #['2022','2100','2060','2045','2035','2080','2050','2030','2065']

# ----------------------------------------------------------------------------------------
#RUN FUNCTIONS

# dem_fix          (regions)
# prep_noaa_mhhw_layers  (regions)
#
# inundated_areas     ('26', ['east_coast'], ['int_low'], ['2065'])

### TESTING SOME STUFF SEPARATE FROM MOSAICS
# interpolate_and_create_water_level_surfaces('26',['east_coast'],['high'],['2025']) #testing


#### DONE
# inundated_areas_select_dems ('26', ['east_coast'], ['int_low'], ['2065'], ['NJ_Middle_GCS_3m_NAVDm_FY21Update']) #run&done
# mosaic_chunks         ('26', ['east_coast'], ['int_low'], ['2065'], ['atlantic_area']) #run
# region_group          ('26', ['east_coast'], ['int_low'], ['2065'])  #run all below
# extract_connected     ('26', ['east_coast'], ['int_low'], ['2065'])
# raster_to_polygon     ('26', ['east_coast'], ['int_low'], ['2065'])


#RUNNING BELOW
# inundated_areas       ('26', ['east_coast'], ['high'], ['2100']) #to scratch
# mosaic_chunks         ('26', ['east_coast'], ['high'], ['2100'], ['all']) #in scratch
# shift_rasters_from_scratch_gdb    ('26', ['east_coast'], ['high'], ['2100'])
# region_group          ('26', ['east_coast'], ['high'], ['2100'])
# extract_connected     ('26', ['east_coast'], ['high'], ['2100'])
# raster_to_polygon     ('26', ['east_coast'], ['high'], ['2100'])



# inundated_areas                (flood_frequency,regions,projections,years)
# mosaic_chunks                  (flood_frequency,regions,projections,years,areas)  
# shift_rasters_from_scratch_gdb (flood_frequency,regions,projections,years)
# region_group                   (flood_frequency,regions,projections,years)
# extract_connected              (flood_frequency,regions,projections,years)
# raster_to_polygon              (flood_frequency,regions,projections,years)

print('FINISHED ALL')

t2 = datetime.datetime.now() ; print('Time start: ' + str(t1) + '  &  Time end: ' + str(t2))

# print('TEST completed')
# print ('STATUS 3: region group run successfully for ' + str(projection) + ' ' + str(year))
print('end')