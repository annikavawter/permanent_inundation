print('RASTER TO POLYGON SCRIPT')

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
# SET INPUTS
flood_frequency = '26'
regions = ['west_coast']
# regions = ['east_coast']
years = ['2060']     #['2022','2100','2060','2045','2035','2080','2050','2030','2065', ...]
projections = ['int_low']  #['high','int_low']

# ----------------------------------------------------------------------------------------------------
# SET WORKSPACE
gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
arcpy.env.workspace = gdb

# ----------------------------------------------------------------------------------------------------
# CONVERT RASTER TO POLYGON
for region in regions:
    for projection in projections:
        for year in years:
            print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')

            # ----------------------------------------------------------------------------------------
            # set gdb in
            gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
            arcpy.env.workspace = gdb_in

            if region == 'west_coast':

                # get extracted rasters, mainland and PNW
                to_convert_main = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region))[0]
                print('Mainland file to convert is ' + str(to_convert_main))

                to_convert_pnw = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW' .format(flood_frequency, year, projection, region))[0]
                print('PNW file to convert is ' + str(to_convert_pnw))

                # ----------------------------------------------------------------------------------------
                # set path and file out
                gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                outname_polygon_main = 'polygon_' + str(to_convert_main)  #final_polygon_
                outname_polygon_pnw = 'polygon_' + str(to_convert_pnw)

                # ----------------------------------------------------------------------------------------
                # create polygons       4 min total

                runpolygonsmain=1
                if runpolygonsmain==1:
                    print('...creating mainland polygon...')
                    arcpy.RasterToPolygon_conversion(to_convert_main, outname_polygon_main, "SIMPLIFY", "VALUE")
                    print ('Converted ' + str(to_convert_main) + ' to polygon')  

                runpolygonspnw=1
                if runpolygonspnw==1:
                    print('...creating PNW polygon...')
                    arcpy.RasterToPolygon_conversion(to_convert_pnw, outname_polygon_pnw, "SIMPLIFY", "VALUE")
                    print ('Converted ' + str(to_convert_pnw) + ' to polygon')  
                                    
                print ('STATUS: polygon conversion run successfully for mainland and PNW west coast, ' + str(projection) + ' ' + str(year))
                

            elif region=='east_coast':
                print('east coast not yet built!')

t2 = datetime.datetime.now() ; print('Time start: ' + str(t1) + '  &  Time end: ' + str(t2))

# print('TEST completed')
print ('STATUS 3: polygons created successfully')
print('end')
