print('MOSAIC & REGION GROUP SCRIPT - WA PNW AREA FOR CONNECTIVITY')

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
# SET INPUTS
flood_frequency = '26'
regions = ['west_coast']
# regions = ['east_coast'] #not in this script
years = ['2022']     #['2022','2100','2060','2045','2035','2080','2050','2030','2065', ...]
projections = ['high']  #['high','int_low']

# ----------------------------------------------------------------------------------------------------
# SET WORKSPACE
gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
arcpy.env.workspace = gdb
#this is reset, then reset again in loop below

# ----------------------------------------------------------------------------------------------------
# MOSAIC THEN REGION GROUP
for region in regions:
    for projection in projections:
        for year in years:

            runmosaic = 0
            if runmosaic==1:
                
                mosversion = 2                  
                #mosversion 1 joins the mosaiced ocean raster with 3 inundated areas, but cell sizes are different and there are some shifts
                #mosversion 2 joins the 4 Puget Sound rasters first then does second mosaic to add SEW2 which has a different cell size - best

                # if mosversion==1:
                #     print ('...running PNW mosaic for ' + str(projection) + ' ' + str(year) + '...')

                #     # ----------------------------------------------------------------------------------------
                #     # set gdb in
                #     gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                #     arcpy.env.workspace = gdb_in
                #     # print('GDB is: ' + gdb_in)

                #     # ----------------------------------------------------------------------------------------
                #     # get inundation raster surfaces raster list: WA PNW only (PugetSoundNW, PugetSoundSW, SEW2)
                #     files_in_1 = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}_WA_PugetSound*' .format(flood_frequency, year, projection, region))
                #     files_in_2 = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}_WA_SEW2*' .format(flood_frequency, year, projection, region))
                    
                #     # ----------------------------------------------------------------------------------------
                #     # set path out
                #     gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

                #     #combine inundation rasters into one set
                #     raster_inund_files = []
                #     for file in files_in_1:
                #         raster_inund_files.append(file)
                #     for file in files_in_2:
                #         raster_inund_files.append(file)
                #     #for file in raster_ocean_area:
                #     #    raster_inund_files.append(file)
                #     print(raster_inund_files)

                #     # ----------------------------------------------------------------------------------------
                #     # reset workspace
                #     gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                #     arcpy.env.workspace = gdb

                #     # ----------------------------------------------------------------------------------------
                #     # get WA ocean raster
                #     raster_ocean_area = arcpy.ListRasters('DEM_WaterBuff_WA_MaskExtract_MosaicAll')
                #     print(raster_ocean_area)          

                #     # ----------------------------------------------------------------------------------------
                #     # create list with full paths of inundation raster surfaces rasters, including WA ocean raster
                #     list_rast_paths = []
                #     # test_path_list = []
                #     for raster_inund_file in raster_inund_files:
                #         #set path for pulling in inundation rasters
                #         rast_path = str(gdb_in) + '/' + str(raster_inund_file)
                #         list_rast_paths.append(rast_path)

                #     for ocean_file in raster_ocean_area:
                #         ocean_rast_path = gdb + '/' + str(ocean_file)  #path from project gdb, not in/out results gdb
                #         list_rast_paths.append(ocean_rast_path)

                #     print('Listing all raster paths to mosaic:')
                #     print(list_rast_paths)

                #     outname_merged = 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW_testmosord' .format(flood_frequency, year, projection, region)
                #     print ('Outname is: ' + outname_merged)

                #     runme = 1
                #     if runme==1:
                #         # ----------------------------------------------------------------------------------------
                #         # combine chunks and create mosaic      TESTING1 start 12:05, end <1 am
                #         print ('...creating mosaic...')

                #         ## SYNTAX: arcpy.MosaicToNewRaster_management(rastersIn,saveLocation,newRasterName,spatialReference,pixelType,cellSize,bandsNum,mosaicOperator,colormapMode)
                #         ## NOTES: original script had cellSize=10, bandsNum=1.  This script and ArcGIS Pro tool gives errors with cellSize=10.
                        
                #         # arcpy.MosaicToNewRaster_management(list_rast_paths, gdb_out, outname_merged,"","",10,1,"MAXIMUM")                       #no rasters out
                #         arcpy.MosaicToNewRaster_management(list_rast_paths, gdb_out, outname_merged,"","","",1,"MAXIMUM")                       #works without the cellsize=10 but takes ages to run

                #         print ('Created PNW mosaic for ' + str(year) + '_' + str(projection))

                if mosversion==2:
                    print ('...running PNW mosaics for ' + str(projection) + ' ' + str(year) + '...')
                    print('...prepping first mosaic...')

                    # ----------------------------------------------------------------------------------------
                    # set gdb in
                    gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                    arcpy.env.workspace = gdb_in
                    # print('GDB is: ' + gdb_in)

                    # ----------------------------------------------------------------------------------------
                    # get inundation raster surfaces raster list: WA PNW only (PugetSoundNW, PugetSoundSW, SEW2); gdb_in
                    file_in_inund_PSNW = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}_WA_PugetSoundNW*' .format(flood_frequency, year, projection, region))
                    file_in_inund_PSSW = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}_WA_PugetSoundSW*' .format(flood_frequency, year, projection, region))
                    file_in_inund_SEW = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}_WA_SEW2*' .format(flood_frequency, year, projection, region))
                    
                    # ----------------------------------------------------------------------------------------
                    # set path out
                    gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

                    # ----------------------------------------------------------------------------------------
                    # reset workspace
                    gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
                    arcpy.env.workspace = gdb

                    # ----------------------------------------------------------------------------------------
                    # get WA ocean rasters; project gdb
                    file_in_ocean_NW = arcpy.ListRasters('DEM_WaterBuff_WA_NW_MaskExtract')
                    file_in_ocean_SW = arcpy.ListRasters('DEM_WaterBuff_WA_SW_MaskExtract')
                    # print(raster_ocean_area)          

                    # ----------------------------------------------------------------------------------------
                    #combine inundation rasters into one set; also create list with full paths rasters
                    raster_inund_files = []
                    list_rast_paths = []
                    for file in file_in_inund_PSNW:     #from gdb_in
                        raster_inund_files.append(file)
                        rast_path = str(gdb_in) + '/' + str(file)
                        list_rast_paths.append(rast_path)
                    for file in file_in_inund_PSSW:
                        raster_inund_files.append(file)
                        rast_path = str(gdb_in) + '/' + str(file)
                        list_rast_paths.append(rast_path)
                    for file in file_in_ocean_NW:       #from gdb
                        raster_inund_files.append(file)
                        rast_path = str(gdb) + '/' + str(file)
                        list_rast_paths.append(rast_path)                        
                    for file in file_in_ocean_SW:
                        raster_inund_files.append(file)
                        rast_path = str(gdb) + '/' + str(file)
                        list_rast_paths.append(rast_path)  

                    print('Listing all raster file names for first mosaic:')
                    print(raster_inund_files)

                    # print('Listing all raster paths to first mosaic:')
                    # print(list_rast_paths)

                    # # ----------------------------------------------------------------------------------------
                    # # create list with full paths of inundation raster surfaces rasters, including WA ocean raster
                    # list_rast_paths = []
                    # for raster_inund_file in raster_inund_files:
                    #     #set path for pulling in inundation rasters
                    #     rast_path = str(gdb_in) + '/' + str(raster_inund_file)
                    #     list_rast_paths.append(rast_path)
                    # print('Listing all raster paths to first mosaic:')
                    # print(list_rast_paths)

                    outname_merged = 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW_mosaicPuget' .format(flood_frequency, year, projection, region)
                    print ('Outname is: ' + outname_merged)

                    runme = 1
                    if runme==1:
                        # ----------------------------------------------------------------------------------------
                        # combine chunks and create mosaic      TESTING1 start 12:05, end <1 am  TEST2V2 start 3:20, end <4:05 both mosaics
                        print ('...creating first mosaic...')
                        puget_sound_mosaic = arcpy.MosaicToNewRaster_management(list_rast_paths, gdb_out, outname_merged,"","","",1,"MAXIMUM") 
                        print ('Created Puget Sound mosaic for ' + str(year) + '_' + str(projection))
                    elif runme==0:
                        print('...skipped first mosaic...')
                    
                    # ----------------------------------------------------------------------------------------
                    # create new list with SEW2 and puget_sound_mosaic rasters; and full paths of rasters
                    print('...prepping second mosaic...')

                    arcpy.env.workspace = gdb_in
                    new_rasters = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW_mosaicPuget' .format(flood_frequency, year, projection, region))
                    # print(new_rasters)
                    
                    arcpy.env.workspace = gdb

                    raster_inund_files_2 = []
                    list_rast_paths_2 = []

                    for file in new_rasters:
                        raster_inund_files_2.append(file)
                        rast_path = str(gdb_in) + '/' + str(file)
                        list_rast_paths_2.append(rast_path)                        
                    for file in file_in_inund_SEW:
                        raster_inund_files_2.append(file)
                        rast_path = str(gdb_in) + '/' + str(file)
                        list_rast_paths_2.append(rast_path)                        

                    print('Listing all raster files for second mosaic:')
                    print(raster_inund_files_2)

                    # print('Listing all raster paths to second mosaic:')
                    # print(list_rast_paths_2)

                    outname_merged_2 = 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW_mosaicAll' .format(flood_frequency, year, projection, region)
                    print ('Outname is: ' + outname_merged_2)

                    runme = 1
                    if runme==1:
                        # ----------------------------------------------------------------------------------------
                        # combine chunks and create mosaic      TESTING1 start 12:05, end <1 am  TEST3PART2only start 4:30, 
                        print ('...creating second mosaic...')
                        puget_sound_mosaic = arcpy.MosaicToNewRaster_management(list_rast_paths_2, gdb_out, outname_merged_2,"","","",1,"MAXIMUM")   
                        print ('Created PNW mosaic for ' + str(year) + '_' + str(projection))
                    elif runme==0:
                        print('...skipped second mosaic...')
            else:
                print('...no mosaicking...')


            # REGION GROUPING
            runregiongrp = 1
            if runregiongrp==1:
                print ('...running region group for ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')

                # ----------------------------------------------------------------------------------------
                # set gdb in
                gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                arcpy.env.workspace = gdb_in

                # ----------------------------------------------------------------------------------------
                # set path out
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

                # ----------------------------------------------------------------------------------------
                # get merged raster for PNW area
                merged_raw_surfaces = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW_mosaicAll' .format(flood_frequency, year, projection, region))
                # for surface in merged_raw_surfaces: print(str(surface))

                # ----------------------------------------------------------------------------------------
                # set paths for inputs
                for surface in merged_raw_surfaces:  #only one per region,year,projection
                    print('Raw surface name is: ' + str(surface))
                    surface_path = gdb_in + '/' +  str(surface) ; print(surface_path)
                    #filename = os.path.basename(fullname)

                    # --------------------------------------------------------------------------------
                    # Perform region group              TIMETEST start 12:53,   only rg no mosaic
                    fullname = str(surface)
                    test_rg = 1
                    if test_rg == 1:
                        print('...region grouping ' + fullname + '...')
                        outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                        # outname_rg = 'rg_' + fullname
                        outname_rg = 'rg_' + 'merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW' .format(flood_frequency, year, projection, region)
                        outRegionGrp.save(outname_rg)

                        print("Region grouped successfully for " + surface)  
            else:
                print('...no region grouping...')
                            

# print('TEST completed')
print ('STATUS 3: region group run successfully for ' + str(projection) + ' ' + str(year) + ' PNW')
print('end')

exit()

#extra code from s4_region_group_test1
# REGION GROUP
for region in regions:
    for projection in projections:
        for year in years:
            print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')

            # ----------------------------------------------------------------------------------------
            # set gdb in
            gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
            arcpy.env.workspace = gdb_in

            # ----------------------------------------------------------------------------------------
            # get merged raster list
            merged_raw_surfaces = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region))
            # for surface in merged_raw_surfaces: print(str(surface))

            # ----------------------------------------------------------------------------------------
            # set path out
            gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)

            # ----------------------------------------------------------------------------------------
            # reset workspace
            #####gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
            ##### # arcpy.env.workspace = gdb

            # ----------------------------------------------------------------------------------------
            # 
            for surface in merged_raw_surfaces:  #only one per region,year,projection
                print('Raw surface name is: ' + str(surface))
                surface_path = gdb_in + '/' +  str(surface) ; print(surface_path)
                #filename = os.path.basename(fullname)

                # ------------------------------------------------------------------------------------            
                #mosaic extra rasters to allow proper connectivity
                if region == 'west_coast':
                    
                    withmosaic = 0  #set mosaicking option on or off
                    if withmosaic == 0:   #no mosaicking, only region grouping
                        #set files in
                        coastal_WA_mosaic = arcpy.ListRasters('merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW' .format(flood_frequency, year, projection, region))
                        for file in coastal_WA_mosaic:
                            print(file)

                        # --------------------------------------------------------------------------------
                        # Perform region group
                        fullname = str(coastal_WA_mosaic)
                        test_rg = 0
                        if test_rg == 1:
                            print('...region grouping ' + fullname + '...')
                            outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                            outname_rg = 'rg_' + fullname
                            outRegionGrp.save(outname_rg)

                            print("Region grouped " + hc_surface)                    
                    
                    elif withmosaic == 1:    #mosaicking then region grouping
                        print('...setting ' + str(region) + ' conditions...')

                        #mosaic the inundated raster surface to surface with filled in coastal water areas to allow complete connectivity for west region
                        #   west coast: areas with connectivity issues are  WA --> use DEM_WaterBuff_WA_MaskExtract_MosaicAll
                        #                                                   CA Channel Isles --> use west_coast_ocean_channel_isles_connector_raster
                        #       DEM_WaterBuff_WA_MaskExtract_MosaicAll: created manually from NoData DEM areas and drawn polygon masks
                        #       west_coast_ocean_channel_isles_connector_raster: created manually to connect Channel Isles DEM water area to mainland water area
                        #   mosaic surface and DEM_WaterBuff_WA_MaskExtract_MosaicAll and west_coast_ocean_channel_isles_connector_raster
                        #   *will need to extract two max cell numbers from the region group: mainland and CA Channel Isles area, and WA PNW area 
                        #       (separated by missing DEM on WA Olympic Peninsula)
                        
                        # --------------------------------------------------------------------------------
                        # set up the coastal water connectivity rasters
                        print('...grabbing coastal waters buffer rasters...')
                        arcpy.env.workspace = gdb
                        coastal_buffer_WA = arcpy.ListRasters('DEM_WaterBuff_WA_MaskExtract_MosaicAll')
                        # coastal_connect_CA = arcpy.ListRasters('west_coast_ocean_channel_isles_connector_raster')

                        for file in coastal_buffer_WA:
                            print('WA coastal buffer is: ' + str(file))
                            #set path for pulling in inundation rasters
                            coastal_buffer_WA_path = str(gdb) + '/' + str('DEM_WaterBuff_WA_MaskExtract_MosaicAll') ; print(coastal_buffer_WA_path)
                        
                        # for file in coastal_connect_CA:
                        #     print('CA coastal buffer is: ' + str(file))
                        #     #set path for pulling in inundation rasters
                        #     coastal_connect_CA_path = str(gdb) + '/' + str('west_coast_ocean_channel_isles_connector_raster') ; print(coastal_connect_CA_path)

                        # --------------------------------------------------------------------------------
                        # set up and mosaic inundation raster and connective buffers
                        mosaic_set = []
                        mosaic_set.append(surface_path)
                        mosaic_set.append(coastal_buffer_WA_path)
                        # mosaic_set.append(coastal_connect_CA_path)
                        print(mosaic_set)

                        # getcellsize_temp = arcpy.GetRasterProperties_management((surface_path),"CELLSIZEX")
                        # cellsize_temp = getcellsize.getOutput(0)
                        # print("Original cellsize is: " + str(cellsize_temp))
                        cellsize_out = "0.000045"

                        outname_hc = "hc_merged_raw_raster_surface_{0}x_{1}_{2}_{3}" .format(flood_frequency,year,projection,region)
                        print('Outname is: ' + outname_hc)

                        # #reset workspace
                        # arcpy.env.workspace = gdb_in

                        print('...mosaicking surface with coastal waters buffer...')

                        test_hc = 1
                        if test_hc == 1:   #start 12:30, still running 8:22... 
                            with arcpy.EnvManager(extent="-131.174129619082 31.736249574995 -112.735719523409 49.5923940886995"):  #west coast extent
                                # arcpy.management.MosaicToNewRaster(r"C:\Users\annik\Documents\UCS\Data\2022\permanent_inundation\west_coast\west_coast.gdb\merged_raw_raster_surface_26x_2022_high_west_coast;DEM_WaterBuff_WA_MaskExtract_MosaicAll_Resamp;west_coast_ocean_channel_isles_connector_raster", r"C:\Users\annik\Documents\UCS\Data\2022\permanent_inundation\west_coast\west_coast.gdb", "hc_merged_raw_raster_surface_26x_2022_high_west_coast", None, "8_BIT_UNSIGNED", "4.4994329757407E-05", 1, "MAXIMUM", "FIRST")
                                hc_surface = arcpy.management.MosaicToNewRaster(mosaic_set, gdb_out, outname_hc, None, "8_BIT_UNSIGNED", cellsize_out, 1, "MAXIMUM", "FIRST")

                            print ('Connectivity mosaic completed for ' + str(year) + ' ' + str(projection))

                        # --------------------------------------------------------------------------------
                        # Perform region group
                        fullname = str(hc_surface)
                        test_rg = 1
                        if test_rg == 1:
                            print('...region grouping ' + fullname + '...')
                            outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                            outname_rg = 'rg_' + fullname
                            outRegionGrp.save(outname_rg)

                            print("Region grouped " + hc_surface)


                # other regions
                else:
                    print ('Region not built yet!')