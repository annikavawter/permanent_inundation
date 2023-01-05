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
runRIisle=1  #temp turn on/off RI block island, rg onward

# ----------------------------------------------------------------------------------------------------
# DEM ERROR FIX
def dem_fix (regions):

    print('SCRIPT FIXING DEM ERRORS')
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

    # east_dems = ['SC_Central_GCS_3m_NAVDm','SC_North_GCS_3m_NAVDm','SC_South_GCS_3m_NAVDm','TX_Central_GCS_3m_NAVDm','TX_North1_GCS_3m_NAVDm','TX_North2_GCS_3m_NAVDm','TX_South1_GCS_3m_NAVDm','TX_South2_GCS_3m_NAVDm','VA_EasternShore_GCS_3m_NAVDm','VA_Middle_GCS_3m_NAVDm','VA_Northern_GCS_3m_NAVDm','VA_Southern_GCS_3m_NAVDm'] #testing

    for region in regions:
        print('Region is: ' + str(region))

        # set save location
        outfolder_str_dem = 'C:/Users/annik/Documents/UCS/Data/2022/DEMs_fix/{0}/{0}_dem_fix_0.gdb' .format(region)

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

            if dem=="RI_GCS_3m_NAVDm":
                mask_connect_RI_DEM_raster = arcpy.ListRasters("mask_connect_RI_DEM_raster")
                ## NOTES: use the following in the RasterCalculator tool in ArcGIS Pro program (syntax not working in this script)
                #Con (~IsNull ("RI_GCS_3m_NAVDm") , "RI_GCS_3m_NAVDm" , ( Con(IsNull("RI_GCS_3m_NAVDm") & IsNull("mask_connect_RI_DEM_raster"), "RI_GCS_3m_NAVDm" , (-99) ) ) )
                #FYI copied python command is as follows:
                # with arcpy.EnvManager(outputCoordinateSystem='GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]', snapRaster="RI_GCS_3m_NAVDm", extent="-71.908450714182 41.0949285433521 -71.0873815613665 42.0196992121698", mask=None):
                #     output_raster = arcpy.ia.RasterCalculator('Con (~IsNull ("RI_GCS_3m_NAVDm") , "RI_GCS_3m_NAVDm" , ( Con(IsNull("RI_GCS_3m_NAVDm") & IsNull("mask_connect_RI_DEM_raster"), "RI_GCS_3m_NAVDm" , (-99) ) ) )'); output_raster.save(r"C:\Users\annik\Documents\UCS\Data\2022\DEMs_fix\east_coast\east_coast_dem_fix.gdb\RI_GCS_3m_NAVDm")
                #print('Area connecting RI island to mainland added to DEM')
                print('Need to add area connecting RI island to mainland ("mask_connect_RI_DEM_raster"')

            outname_dem_fixnulls = str(outfolder_str_dem) + '/' + str(dem)

            print('...saving raster...')
            dem_fixnulls.save(outname_dem_fixnulls)
            print('New DEM created: ' + str(outname_dem_fixnulls))
            
        print ('All DEMs fixed for ' + region)
    print('STATUS: DEMs fixed successfully')

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
                gdb_dem = 'C:/Users/annik/Documents/UCS/Data/2022/DEMs_fix/{0}/{0}_dem_fix_0.gdb' .format(region)  #
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

                #     if region == 'west_coast':
                #         dems = west_dems #18 total
                #     elif region == 'east_coast':
                #         dems = east_dems  #58 total
                #     elif region == 'hi_coast':
                #         dems = hi_dems 
                #     elif region == 'gu_coast':
                #         dems = gu_dems
                #     elif region == 'as_coast':
                #         dems = as_dems
                #     elif region == 'pr_coast':
                #         dems = pr_dems
                #     elif region == 'vi_coast':
                #         dems = vi_dems

                # reset workspace
                arcpy.env.workspace = gdb_dem

                dems = arcpy.ListRasters('*')
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
                    outname_inundated_area_surface = 'inundated_area_surface_%sx_%s_%s_%s' %(flood_frequency, year, projection, region) + '_' + dem[:20]
                    outname_inundated_area_surface_full = str(gdb_out) + '/' + str(outname_inundated_area_surface)

                    #check if inundation file already exists, create if
                    arcpy.env.workspace = gdb_out
                    inundated_check = arcpy.ListRasters(outname_inundated_area_surface) ; #print(len(inundated_check))
                    arcpy.env.workspace = gdb_dem
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

# ----------------------------------------------------------------------------------------
# MOSAIC                        #5 hours east
def mosaic_chunks (flood_frequency,regions,projections,years):
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

                # set workspace
                arcpy.env.workspace = gdb_in

                # ----------------------------------------------------------------------------------------
                # get inundation raster surfaces raster list
                raster_inund_files = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}*' .format(flood_frequency, year, projection, region))
                # print(raster_inund_files)

                # ----------------------------------------------------------------------------------------
                # set path out
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

                # ----------------------------------------------------------------------------------------
                # reset workspace
                # arcpy.env.workspace = gdb

                # ----------------------------------------------------------------------------------------
                # resample in certain regions due to processing speed, then mosaic

                if region == 'east_coast':
                    ##NOTES: resample 3 min per element
                    ##NOTES: full resample and mosaic, east coast ~5 hrs to run
                    
                    # set cellsize based on speed and quality testing
                    cellsize_out = "0.000045"
                    print('Cellsize out is ' + str(cellsize_out))

                    # ----------------------------------------------------------------------------------
                    # resample inundation rasters and save        #NOTES 3-5 min each raster !!!
                    list_all_resampled_rast_paths = []
                    for raster_inund_file in raster_inund_files:

                        #set path for pulling in inundation raster
                        rast_path = str(gdb_in) + '/' + str(raster_inund_file)
                        
                        #set resampled raster file name
                        outname_resample = 'scratch_' + str(raster_inund_file) + '_resamp_temp'  #cellsz_pt000045
                        outname_resample_path = str(gdb_out) + '/' + str(outname_resample)  

                        #delete resampled files BAD
                        delweird=0
                        if delweird==1:
                            resample_check_TODEL = arcpy.ListRasters('scratch_' + str(raster_inund_file) + '_resampled_temp')
                            if len(resample_check_TODEL)==1:
                                for del_file in resample_check_TODEL:
                                    print(str(del_file))
                                    arcpy.Delete_management(del_file)
                                print('Deleted weird old resampled file')

                        #resample file, if resampled file doesn't already exist
                        resample_check = arcpy.ListRasters('scratch_' + str(raster_inund_file) + '_resamp_test')
                        
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
                        list_all_resampled_rast_paths.append(outname_resample_path)
                    
                    # ----------------------------------------------------------------------------------
                    # combine chunks and create mosaic
                    runmos=0  #turn on 
                    if runmos==1:
                        print ('...resampling done, setting up mosaic...')

                        num_elem = len(list_resampled_rast_paths)
                        print (str(num_elem) + ' elements total to mosaic')

                        outname_mosaic = 'merged_raw_raster_surface_%sx_%s_%s_%s' %(flood_frequency, year, projection, region)
                        print ('Mosaic outname is: ' + outname_mosaic)

                        print ('...creating mosaic...')
                        arcpy.MosaicToNewRaster_management(list_resampled_rast_paths, gdb_out, outname_mosaic,"","",cellsize_out,1,"MAXIMUM")   #cellsz_pt000045
                        # with arcpy.EnvManager(extent="-98.0633652753512 24.393838762558 -66.8803114143969 46.3976170001712"): #east coast extent
                        #     arcpy.MosaicToNewRaster_management(list_resampled_rast_paths, gdb_out, outname_mosaic,"","",cellsize_out,1,"MAXIMUM")   #cellsz_pt000045
                        print('Created mosaic')

                        #delete resampled files
                        for del_file in list_resampled_rast_paths:
                            arcpy.Delete_management(del_file)
                        print('Deleted temporary resampled files')

                    # ----------------------------------------------------------------------------------
                    # combine chunks and create mosaic FOR GULF AND EAST AREAS SEPARATELY, WITH OVERLAP FOR RG
                    runmos=1  #turn on 
                    if runmos==1:
                        print ('...resampling done, setting up mosaic...')

                        #make new set of resampled raster path names to mosaic based on area
                        area_set_gulf=['TX','LA','MS','AL','FL','GA_JAX']
                        area_set_atlantic=['FL_JAX_1','GA','SC','NC','VA','DC','MD','DE','PA','NJ','NY','CT','RI','MA','NH','ME']
                        resampled_gulf_paths = [] ; resampled_atlantic_paths = []

                        for state in area_set_gulf:
                            area_rast = arcpy.ListRasters('scratch_inundated_area_surface_{0}x_{1}_{2}_{3}_{4}*_resamp_temp' .format(flood_frequency,year,projection,region,state))
                            for lyr in area_rast:
                                resampled_gulf_paths.append(str(gdb_out) + '/' + str(lyr)) ; #print(lyr)
                        # print(resampled_gulf_paths)
                        # for el in resampled_gulf_paths: print(el)

                        for state in area_set_atlantic:
                            area_rast = arcpy.ListRasters('scratch_inundated_area_surface_{0}x_{1}_{2}_{3}_{4}*_resamp_temp' .format(flood_frequency,year,projection,region,state))
                            for lyr in area_rast:
                                resampled_atlantic_paths.append(str(gdb_out) + '/' + str(lyr)) ; #print(lyr)
                        # print(resampled_atlantic_paths)
                        # for el in resampled_atlantic_paths: print(el)
                        
                        # ----------------------------------------------------------------------------------
                        areas=['gulf_area','atlantic_area']
                        areas=['atlantic_area']

                        for area in areas:
                            tta = datetime.datetime.now()
                            #set list of resampled inundation raster paths in by region
                            if area=='gulf_area':
                                list_resampled_rast_paths = resampled_gulf_paths
                            elif area=='atlantic_area':
                                list_resampled_rast_paths = resampled_atlantic_paths
                            
                            print('Area to mosaic: ' + area)
                            for file in list_resampled_rast_paths:
                                print(file)

                            num_elem = len(list_resampled_rast_paths)
                            print (str(num_elem) + ' elements total to mosaic in ' + area)
                            if area=='gulf_area':
                                if num_elem!=25:
                                    print('Wrong number of resampled elements for ' + area + ' - need 25!'); exit()
                            elif area=='atlantic_area':
                                if num_elem!=35:
                                    print('Wrong number of resampled elements for ' + area + ' - need 35!'); exit()

                            outname_mosaic = 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency, year, projection, region, area)
                            print ('Mosaic outname is: ' + outname_mosaic)

                            print ('...creating ' + area + ' mosaic...')
                            arcpy.MosaicToNewRaster_management(list_resampled_rast_paths, gdb_out, outname_mosaic,"","",cellsize_out,1,"MAXIMUM")   #cellsz_pt000045
                            # with arcpy.EnvManager(extent="-98.0633652753512 24.393838762558 -66.8803114143969 46.3976170001712"): #east coast extent
                            #     arcpy.MosaicToNewRaster_management(list_resampled_rast_paths, gdb_out, outname_mosaic,"","",cellsize_out,1,"MAXIMUM")   #cellsz_pt000045
                            print('Created mosaic for ' + area)
                            moscreated=1
                            ttb = datetime.datetime.now() ; print('Area mosaic ' + area + ' Time start: ' + str(tta) + '  &  Time end: ' + str(ttb))

                        if moscreated==1:
                            #delete resampled files
                            for del_file in list_all_resampled_rast_paths:
                                arcpy.Delete_management(del_file)
                            print('Deleted temporary resampled files')
                
                elif region == 'west_coast':
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

                print ('STATUS 3: mosaic run successfully for ' + str(projection) + ' ' + str(year))
                tf2 = datetime.datetime.now() ; print('Time start: ' + str(tf1) + '  &  Time end: ' + str(tf2))
    ## NOTES: original script had cellSize=10, bandsNum=1.  This script and ArcGIS Pro tool gives errors with cellSize=10.
    #   Cellsize "0.000045" selected based on cell size of existing layers, speed testing, and quality control

    print('END MOSAIC INUNDATED SURFACES SCRIPT')
    
# ----------------------------------------------------------------------------------------
# REGION GROUP
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
                            surface_path = gdb_in + '/' +  str(merged_raw_surface) ; print(surface_path)
                            ##filename = os.path.basename(fullname)

                            fullname = str(merged_raw_surface)
                            test_rg = 1
                            if test_rg == 1:
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
                                surface_path_pnw = gdb_in + '/' +  str(merged_raw_surface) ; print(surface_path_pnw)
                                ##filename = os.path.basename(fullname)

                                fullname = str(merged_raw_surface)
                                test_rg = 1                                                             ##TURN BACK ON
                                if test_rg == 1:
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
                            surface_path = gdb_in + '/' +  str(merged_raw_surface) ; print(surface_path)
                            ##filename = os.path.basename(fullname)

                            fullname = str(merged_raw_surface)
                            test_rg = 1
                            if test_rg == 1:
                                print('...region grouping ' + fullname + '...')
                                outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                                outname_rg = 'rg_' + fullname
                                outRegionGrp.save(outname_rg)
                                print('Region grouped ' + str(merged_raw_surface))  
                        # ----------------------------------------------------------------------------------------
                        elif area=='RI_block_isle':
                            if runRIisle==1:

                                # get RI inundation raster, NOT merged raster
                                raw_surface = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}_RI*' .format(flood_frequency,year,projection,region))[0]
                                ## print(str(raw_surface))

                                # ----------------------------------------------------------------------------------------
                                # Perform region group
                                print(area + ' raw surface name is: ' + str(raw_surface))
                                surface_path_riisle = gdb_in + '/' +  str(raw_surface) ; print(surface_path_riisle)
                                ##filename = os.path.basename(fullname)

                                fullname = str(raw_surface)
                                test_rg = 1                                                             ##TURN BACK ON
                                if test_rg == 1:
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
# EXTRACT CONNECTED REGIONS
def extract_connected (flood_frequency,regions,projections,years):
    for region in regions:
        for projection in projections:
            for year in years:

                print('EXTRACTING CONNECTED REGION GROUPS SCRIPT')
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
                            if runRIisle==1:
                                print('Extracting east coast ' + area)

                                # get region group raster
                                rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area))[0]
                                rg_surface_path = gdb_in + '/' + rg_surface
                                print(area + ' file to extract is ' + rg_surface)

                                #set outname
                                outname_extract = 'extract_' + rg_surface

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
# CONVERT RASTER TO POLYGON
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
                    for area in areas:
                        if area=='gulf_area' or area=='atlantic_area':
                            print('Building polygons for east coast ' + area)

                            # get extracted rasters
                            to_convert = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area))[0]
                            print(area + ' file to convert is ' + str(to_convert))

                            # set path and file out
                            # gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                            outname_polygon = 'polygon_' + str(to_convert)  #final_polygon_

                            # ----------------------------------------------------------------------------------------
                            # create polygons       4 min total ?
                            runpolygons=1
                            if runpolygons==1:
                                print('...creating ' + str.lower(area) + ' polygon...')
                                arcpy.RasterToPolygon_conversion(to_convert, outname_polygon, "SIMPLIFY", "VALUE")
                                print ('Converted ' + str(to_convert) + ' to polygon') 
                        
                        elif area=='RI_block_isle':
                            if runRIisle==1:
                                print('Building polygons for east coast ' + area)

                                # get extracted rasters
                                to_convert = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_{4}' .format(flood_frequency,year,projection,region,area))[0]
                                print(area + ' file to convert is ' + str(to_convert))

                                # set path and file out
                                # gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                                outname_polygon = 'polygon_' + str(to_convert)

                                # ----------------------------------------------------------------------------------------
                                # create polygons       4 min total
                                runpolygons=1
                                if runpolygons==1:
                                    print('...creating ' + str.lower(area) + ' polygon...')
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
regions = ['east_coast']    #['east_coast','west_coast','pr_coast','vi_coast','gu_coast','hi_coast','as_coast']
projections = ['high']       #['high','int_low']
years = ['2050']          #['2022','2100','2060','2045','2035','2080','2050','2030','2065', ...]

# ----------------------------------------------------------------------------------------
#RUN FUNCTIONS

# dem_fix          (regions)
#
# inundated_areas     ('26', ['east_coast'], ['int_low'], ['2065'])

mosaic_chunks       ('26', ['east_coast'], ['int_low'], ['2060'])  #test area split  
# region_group          ('26', ['east_coast'], ['int_low'], ['2060'])  #test area split  rg gulf 1:55
# extract_connected     ('26', ['east_coast'], ['int_low'], ['2060'])  #test area split
# raster_to_polygon     ('26', ['east_coast'], ['int_low'], ['2060'])  #test area split

# inundated_areas     (flood_frequency,regions,projections,years)
# mosaic_chunks       (flood_frequency,regions,projections,years)  
# region_group        (flood_frequency,regions,projections,years)
# extract_connected   (flood_frequency,regions,projections,years)
# raster_to_polygon   (flood_frequency,regions,projections,years)


print('FINISHED ALL')

t2 = datetime.datetime.now() ; print('Time start: ' + str(t1) + '  &  Time end: ' + str(t2))

# print('TEST completed')
# print ('STATUS 3: region group run successfully for ' + str(projection) + ' ' + str(year))
print('end')