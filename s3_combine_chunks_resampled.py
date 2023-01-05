print('MOSAIC SCRIPT - WITH RESAMPLING')

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
# regions = ['west_coast']
regions = ['east_coast']
years = ['2100']     #['2022','2100','2060','2045','2035','2080','2050','2030','2065', ...]
projections = ['high']  #['high','int_low']

# ----------------------------------------------------------------------------------------------------
# SET WORKSPACE
gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
arcpy.env.workspace = gdb
#this is reset, then reset again in loop below

# ----------------------------------------------------------------------------------------------------
# MOSAIC

for region in regions:
    for projection in projections:
        for year in years:
            print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')

            # ----------------------------------------------------------------------------------------
            # set gdb in
            gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
            arcpy.env.workspace = gdb_in
            # print('GDB is: ' + gdb_in)

            # ----------------------------------------------------------------------------------------
            # get inundation raster surfaces raster list
            raster_inund_files = arcpy.ListRasters('inundated_area_surface_{0}x_{1}_{2}_{3}*' .format(flood_frequency, year, projection, region))
            # print(raster_inund_files)

            # ----------------------------------------------------------------------------------------
            # set path out
            gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
            outname_merged = 'merged_raw_raster_surface_%sx_%s_%s_%s' %(flood_frequency, year, projection, region)
            print ('Outname is: ' + outname_merged)

            # ----------------------------------------------------------------------------------------
            # reset workspace - don't need to
            gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
            # arcpy.env.workspace = gdb

            # ----------------------------------------------------------------------------------------
            # create list with full paths of inundation raster surfaces rasters - don't need
            # list_rast_paths = []
            # for raster_inund_file in raster_inund_files:
            #     #set path for pulling in inundation rasters
            #     rast_path = str(gdb_in) + '/' + str(raster_inund_file)
            #     list_rast_paths.append(rast_path)
            # print(list_rast_paths)

            # NOTES: raster_inund_files is the set of short file names; list_rast_paths is the set of long path names
            
            # ----------------------------------------------------------------------------------------
            if region == 'east_coast':
                # print ('...testing...')     #NOTES: two pieces start 6:32, end resample 6:33, end mosaic 6:38  --> 6 minutes / 2 elements  = 3 min per element
                                            #NOTES: full east coast start 6:45, end resample _sc south 7:04, end mosaic 1.34   (should be ~5 hrs)
                

                # ----------------------------------------------------------------------------------
                # resample inundation rasters and save
                # NOTES: cell size 10 is too large and won't run due to error. (arcgisscripting.ExecuteError: ERROR 000995: Resampling failed. The specified cell size is too large or too small.)
                #        Tested 0.01, visual inspection shows it does not appropriately represent inundated area.
                #        Tested 0.0001, which has short run (<1 min each) time while covering inundated area.  *best of these!
                #        Tested 0.00001, which has longer run time (3-5 min each) for rasters with very small original cells (4.5E-05)
                #        Tested 0.000045, which has short run time (<1 min each) and matches the max cell size of east coast rasters.  *BEST OVERALL
                
                cellsize_out = "0.000045"
                print('Cellsize out is 0.000045')

                # # test two first
                # raster_inund_files = raster_inund_files[0:2]
                # print(raster_inund_files)

                list_resampled_rast_paths = []
                for raster_inund_file in raster_inund_files:

                    print ("...resampling...")
                    print (raster_inund_file)

                    #set path for pulling in inundation raster
                    rast_path = str(gdb_in) + '/' + str(raster_inund_file)

                    #set resampled raster file name
                    outname_resample = str(gdb_out) + "/" + str(raster_inund_file) + "_resample_cellsz_pt000045"  #cellsz_pt000045
                    print ("Resample outname is: " + outname_resample)

                    # getcellsize_temp = arcpy.GetRasterProperties_management((rast_test),"CELLSIZEX")
                    # cellsize_temp = getcellsize.getOutput(0)
                    # print("Original cellsize is: " + str(cellsize_temp))

                    #resample
                    arcpy.Resample_management(rast_path,outname_resample,cellsize_out,"NEAREST")  #NEAREST does not change cell values

                    #make new set of resampled raster path names to mosaic
                    list_resampled_rast_paths.append(outname_resample)
 
                # print (list_resampled_rast_paths)

                # ----------------------------------------------------------------------------------
                # combine chunks and create mosaic
                # NOTES: start 5:37 end 5:44 -- still 3-5 min each raster !!!

                print ('...resampling done, setting up mosaic...')

                num_elem = len(list_resampled_rast_paths)
                print (str(num_elem) + ' elements total')

                outname_mosaic = 'merged_raw_raster_surface_%sx_%s_%s_%s' %(flood_frequency, year, projection, region)
                print ('Mosaic outname is: ' + outname_mosaic)

                print ('...creating mosaic...')

                arcpy.MosaicToNewRaster_management(list_resampled_rast_paths, gdb_out, outname_mosaic,"","",cellsize_out,1,"MAXIMUM")   #cellsz_pt000045

                print ('Created mosaic for ' + str(region) + '_' + str(year) + '_' + str(projection))

                
 

  

                    

                    


                ## WORKING ABOVE HERE, TESTING RESAMPLING AND CELL SIZES

            elif region == 'west_coast':
                # ----------------------------------------------------------------------------------
                # resample inundation rasters and save
                
                cellsize_out = "0.000045"
                print('Cellsize out is 0.000045')

                list_resampled_rast_paths = []
                for raster_inund_file in raster_inund_files:

                    print ("...resampling...")
                    print (raster_inund_file)

                    #set path for pulling in inundation raster
                    rast_path = str(gdb_in) + '/' + str(raster_inund_file)

                    #set resampled raster file name
                    outname_resample = str(gdb_out) + "/" + str(raster_inund_file) + "_resample_cellsz_pt000045"  #cellsz_pt000045
                    print ("Resample outname is: " + outname_resample)

                    # getcellsize_temp = arcpy.GetRasterProperties_management((rast_test),"CELLSIZEX")
                    # cellsize_temp = getcellsize.getOutput(0)
                    # print("Original cellsize is: " + str(cellsize_temp))

                    #resample
                    arcpy.Resample_management(rast_path,outname_resample,cellsize_out,"NEAREST")  #NEAREST does not change cell values

                    #make new set of resampled raster path names to mosaic
                    list_resampled_rast_paths.append(outname_resample)
 
                # print (list_resampled_rast_paths)

                # ----------------------------------------------------------------------------------
                # combine chunks and create mosaic
                # NOTES: start 5:37 end 5:44 -- still 3-5 min each raster !!! FOR EAST NOT WEST

                print ('...resampling done, setting up mosaic...')

                num_elem = len(list_resampled_rast_paths)
                print (str(num_elem) + ' elements total')

                outname_mosaic = 'merged_raw_raster_surface_%sx_%s_%s_%s' %(flood_frequency, year, projection, region)
                print ('Mosaic outname is: ' + outname_mosaic)

                print ('...creating mosaic...')

                arcpy.MosaicToNewRaster_management(list_resampled_rast_paths, gdb_out, outname_mosaic,"","",cellsize_out,1,"MAXIMUM")   #cellsz_pt000045

                print ('Created mosaic for ' + str(region) + '_' + str(year) + '_' + str(projection))

            # # ----------------------------------------------------------------------------------------
            # elif region == 'west_coast':
            #     # combine chunks and create mosaic
            #     print ('...creating mosaic...')

            #     ## SYNTAX: arcpy.MosaicToNewRaster_management(rastersIn,saveLocation,newRasterName,spatialReference,pixelType,cellSize,bandsNum,mosaicOperator,colormapMode)
            #     ## NOTES: original script had cellSize=10, bandsNum=1.  This script and ArcGIS Pro tool gives errors with cellSize=10.
                
            #     # arcpy.MosaicToNewRaster_management(list_rast_paths, gdb_out, outname_merged,"","",10,1,"MAXIMUM")                       #no rasters out
            #     arcpy.MosaicToNewRaster_management(list_rast_paths, gdb_out, outname_merged,"","","",1,"MAXIMUM")                       #works without the cellsize=10 but takes ages to run

            #     print ('Created mosaic for ' + str(year) + '_' + str(projection))

            # ----------------------------------------------------------------------------------------
            else:
                print('wrong region!')
            
            print ('Created mosaic for ' + str(year) + '_' + str(projection))



# ----------------------------------------------------------------------------------------------------
##from ArcGIS Pro geoprocessing tool, testing only five rasters (TX):
# arcpy.management.MosaicToNewRaster("inundated_area_surface_26x_2022_high_east_coast_TX_Central_GCS_3m_NA;inundated_area_surface_26x_2022_high_east_coast_TX_Central_GCS_3m_NA;inundated_area_surface_26x_2022_high_east_coast_TX_North1_GCS_3m_NAV;inundated_area_surface_26x_2022_high_east_coast_TX_North2_GCS_3m_NAV;inundated_area_surface_26x_2022_high_east_coast_TX_South1_GCS_3m_NAV;inundated_area_surface_26x_2022_high_east_coast_TX_South2_GCS_3m_NAV", r"C:\Users\annik\Documents\ArcGIS\Projects\UCS_Flood_2022\UCS_Flood_2022.gdb", "TEST31_merged_raw_raster_surface_26x_2022_high_east_coast_TX", 'GEOGCS["GCS_NAD_1983_2011",DATUM["D_NAD_1983_2011",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]', "8_BIT_UNSIGNED", None, 1, "MAXIMUM", "FIRST")
# ----------------------------------------------------------------------------------------------------

# print('TEST completed')
print ('STATUS 3: mosaic run successfully for ' + str(projection) + ' ' + str(year))
print('end')