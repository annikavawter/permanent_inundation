print('MOSAIC SCRIPT')

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
# regions = ['east_coast']
years = ['2045','2060','2065','2080']     #['2022','2100','2060','2045','2035','2080','2050','2030','2065', ...]
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
            print ('...running ' + str(projection) + ' ' + str(year) + '...')

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
            outname_merged = 'merged_raw_raster_surface_%sx_%s_%s_%s' %(flood_frequency, year, projection, region)  #change test nums
            print ('Outname is: ' + outname_merged)

            # ----------------------------------------------------------------------------------------
            # reset workspace
            gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
            # arcpy.env.workspace = gdb

            # ----------------------------------------------------------------------------------------
            # create list with full paths of inundation raster surfaces rasters
            list_rast_paths = []
            test_path_list = []
            for raster_inund_file in raster_inund_files:
                #set path for pulling in inundation rasters
                rast_path = str(gdb_in) + '/' + str(raster_inund_file)
                list_rast_paths.append(rast_path)
            print(list_rast_paths)

            # ----------------------------------------------------------------------------------------
            # combine chunks and create mosaic
            print ('...creating mosaic...')

            ## SYNTAX: arcpy.MosaicToNewRaster_management(rastersIn,saveLocation,newRasterName,spatialReference,pixelType,cellSize,bandsNum,mosaicOperator,colormapMode)
            ## NOTES: original script had cellSize=10, bandsNum=1.  This script and ArcGIS Pro tool gives errors with cellSize=10.
            
            # arcpy.MosaicToNewRaster_management(list_rast_paths, gdb_out, outname_merged,"","",10,1,"MAXIMUM")                       #no rasters out
            arcpy.MosaicToNewRaster_management(list_rast_paths, gdb_out, outname_merged,"","","",1,"MAXIMUM")                       #works without the cellsize=10 but takes ages to run

            print ('Created mosaic for ' + str(year) + '_' + str(projection))


# ----------------------------------------------------------------------------------------------------
##from ArcGIS Pro geoprocessing tool, testing only five rasters (TX):
# arcpy.management.MosaicToNewRaster("inundated_area_surface_26x_2022_high_east_coast_TX_Central_GCS_3m_NA;inundated_area_surface_26x_2022_high_east_coast_TX_Central_GCS_3m_NA;inundated_area_surface_26x_2022_high_east_coast_TX_North1_GCS_3m_NAV;inundated_area_surface_26x_2022_high_east_coast_TX_North2_GCS_3m_NAV;inundated_area_surface_26x_2022_high_east_coast_TX_South1_GCS_3m_NAV;inundated_area_surface_26x_2022_high_east_coast_TX_South2_GCS_3m_NAV", r"C:\Users\annik\Documents\ArcGIS\Projects\UCS_Flood_2022\UCS_Flood_2022.gdb", "TEST31_merged_raw_raster_surface_26x_2022_high_east_coast_TX", 'GEOGCS["GCS_NAD_1983_2011",DATUM["D_NAD_1983_2011",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]', "8_BIT_UNSIGNED", None, 1, "MAXIMUM", "FIRST")
# ----------------------------------------------------------------------------------------------------

# print('TEST completed')
print ('STATUS 3: mosaic run successfully for ' + str(projection) + ' ' + str(year))
print('end')