# exit()
print('SCRIPT ALL RG THRU POLYGONS')

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
# regions = ['west_coast']
# # regions = ['east_coast']
# years = ['2060']     #['2022','2100','2060','2045','2035','2080','2050','2030','2065', ...]
# projections = ['int_low']  #['high','int_low']

# ----------------------------------------------------------------------------------------------------
# SET WORKSPACE
gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
arcpy.env.workspace = gdb

runpnw=1 #temp turn on/off pnw

# ----------------------------------------------------------------------------------------------------
# REGION GROUP
def region_group (flood_frequency,regions,projections,years):
    for region in regions:
        for projection in projections:
            for year in years:

                print('REGION GROUPING SCRIPT FOR CONNECTIVITY')
                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')

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

                                print("Region grouped " + merged_raw_surface)  
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

                                    print("Region grouped " + merged_raw_surface)  
                            
                    print ('STATUS: region group run successfully for mainland & PNW ' + str(region) + ' ' + str(projection) + ' ' + str(year))
                        
                # other regions
                else:
                    print ('Region not built yet!') ; exit()
                    print ('STATUS: region group run successfully for ' + str(projection) + ' ' + str(year))


# ----------------------------------------------------------------------------------------
# EXTRACT CONNECTED REGIONS
def extract_connected (flood_frequency,regions,projections,years):
    for region in regions:
        for projection in projections:
            for year in years:

                print('EXTRACT CONNECTED REGION GROUPS SCRIPT')
                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')

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
                    print('Extracting for both mainland west coast and PNW area')

                    # ----------------------------------------------------------------------------------------
                    # Perform region group for each area: mainland and pnw
                    areas = ['mainland','pnw']
                    for area in areas:
                        if area=='mainland':
                            # get region group raster
                            rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region))[0]
                            rg_surface_path = gdb_in + '/' + rg_surface
                            print('Mainland file to extract is ' + rg_surface)

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
                        elif area=='pnw':
                            if runpnw==1:

                                # get region group raster
                                rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW' .format(flood_frequency, year, projection, region))[0]
                                rg_surface_path = gdb_in + '/' + rg_surface
                                print('PNW file to extract is ' + rg_surface)

                                #set outname
                                outname_extract = 'extract_' + rg_surface
                                
                                # ----------------------------------------------------------------------------------------
                                # set areas for extraction - PNW ONLY
                                print('...finding areas to extract from ' + rg_surface + '...')

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
                    print('region not yet built')
                    exit()
                    print ('STATUS: extract connected regions run successfully for _')


# ----------------------------------------------------------------------------------------
# CONVERT RASTER TO POLYGON
def raster_to_polygon (flood_frequency,regions,projections,years):
    for region in regions:
        for projection in projections:
            for year in years:

                print('RASTER TO POLYGON SCRIPT')
                print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')

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
                    print('Extracting for both mainland west coast and PNW area')

                    # ----------------------------------------------------------------------------------------
                    # Perform polygon conversion for each area: mainland and pnw
                    areas = ['mainland','pnw']
                    for area in areas:
                        if area=='mainland':
                            # get extracted rasters
                            to_convert = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region))[0]
                            print('Mainland file to convert is ' + str(to_convert))

                            # set path and file out
                            # gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
                            outname_polygon = 'polygon_' + str(to_convert)  #final_polygon_

                            # ----------------------------------------------------------------------------------------
                            # create polygons       4 min total
                            runpolygons=1
                            if runpolygons==1:
                                print('...creating mainland polygon...')
                                arcpy.RasterToPolygon_conversion(to_convert, outname_polygon, "SIMPLIFY", "VALUE")
                                print ('Converted ' + str(to_convert) + ' to polygon')  

                        elif area=='pnw':
                            if runpnw==1:

                                # get extracted rasters
                                to_convert = arcpy.ListRasters('extract_rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW' .format(flood_frequency, year, projection, region))[0]
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
                    print('east coast not yet built!')

                print('FINISHED ' + str(region) + ' ' + str(projection) + ' ' + str(year))


# ----------------------------------------------------------------------------------------
#SET FUNCTION INPUTS
flood_frequency = '26'
regions = ['west_coast']
projections = ['high']    #['high','int_low']
years = ['2050','2030','2065']      #['2022','2100','2060','2045','2035','2080','2050','2030','2065', ...]

# ----------------------------------------------------------------------------------------
#RUN FUNCTIONS
region_group        (flood_frequency,regions,projections,years)
extract_connected   (flood_frequency,regions,projections,years)
raster_to_polygon   (flood_frequency,regions,projections,years)

print('FINISHED ALL')



t2 = datetime.datetime.now() ; print('Time start: ' + str(t1) + '  &  Time end: ' + str(t2))

# print('TEST completed')
# print ('STATUS 3: region group run successfully for ' + str(projection) + ' ' + str(year))
print('end')