print('REGION GROUPING SCRIPT FOR CONNECTIVITY')

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
years = ['2100']     #['2022','2100','2060','2045','2035','2080','2050','2030','2065', ...]
projections = ['int_low']  #['high','int_low']

# ----------------------------------------------------------------------------------------------------
# SET WORKSPACE
gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
arcpy.env.workspace = gdb

# ----------------------------------------------------------------------------------------------------
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
            ##### # arcpy.env.workspace = gdb

            # ----------------------------------------------------------------------------------------
            # 
            for surface in merged_raw_surfaces:  #only one per region,year,projection
                print('Raw surface name is: ' + str(surface))
                surface_path = gdb_in + '/' +  str(surface) ; print(surface_path)
                #filename = os.path.basename(fullname)

                if region == 'west_coast':
                    
                    # --------------------------------------------------------------------------------
                    # Perform region group
                    fullname = str(surface)
                    test_rg = 1
                    if test_rg == 1:
                        print('...region grouping ' + fullname + '...')
                        outRegionGrp = RegionGroup(fullname, "EIGHT", "WITHIN",'NO_LINK')
                        outname_rg = 'rg_' + fullname
                        outRegionGrp.save(outname_rg)

                        print("Region grouped " + surface)                    
                    
                # other regions
                else:
                    print ('Region not built yet!')

t2 = datetime.datetime.now() ; print('Time start: ' + str(t1) + '  &  Time end: ' + str(t2))

# print('TEST completed')
print ('STATUS 3: region group run successfully for ' + str(projection) + ' ' + str(year))
print('end')