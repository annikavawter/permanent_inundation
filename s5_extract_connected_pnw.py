print('EXTRACT CONNECTED REGION GROUPS SCRIPT - PNW')

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
years = ['2022']     #['2022','2100','2060','2045','2035','2080','2050','2030','2065', ...]
projections = ['high']  #['high','int_low']

# ----------------------------------------------------------------------------------------------------
# SET WORKSPACE
gdb = 'C:/Users/annik/Documents/ArcGIS/Projects/UCS_Flood_2022/UCS_Flood_2022.gdb'
arcpy.env.workspace = gdb

# ----------------------------------------------------------------------------------------------------
# EXTRACT CONNECTED REGIONS
for region in regions:
    for projection in projections:
        for year in years:
            print ('...running ' + str(region) + ' ' + str(projection) + ' ' + str(year) + '...')

            # ----------------------------------------------------------------------------------------
            # set gdb in
            gdb_in = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
            arcpy.env.workspace = gdb_in

            # ----------------------------------------------------------------------------------------
            # get region group raster
            # rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region))[0]
            rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}_WA_PNW' .format(flood_frequency, year, projection, region))[0]
            rg_surface_path = gdb_in + '/' + rg_surface
            print('File to extract is ' + rg_surface)

            # ----------------------------------------------------------------------------------------
            # set path and file out
            gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
            outname_extract = 'extract_' + rg_surface
            print('Outname is: ' + outname_extract)
            # outname_extract = str(gdb_out) + '/' + 'extract_' + rg_surface

            # ----------------------------------------------------------------------------------------
            # set areas for extraction - PNW ONLY
            print('...finding areas to extract from ' + rg_surface + '...')
            # arcpy.env.workspace = gdb

            arr = arcpy.da.FeatureClassToNumPyArray(rg_surface, ('Value', 'Count'))
            count = arr['Count']
            value = arr['Value']
            index_to_extract = numpy.argmax(count)
            value_to_extract = str(value[index_to_extract])

            print('...building SQL string...')
            inSQLClause = 'Value = ' + value_to_extract

            print('Extracting value=' + value_to_extract + ' from ' + rg_surface)
            print('...extracting values...')

            # ----------------------------------------------------------------------------------------
            # extract connected areas
            runextract=1
            if runextract==1:
                attExtract = ExtractByAttributes(rg_surface, inSQLClause)
                attExtract.save(outname_extract)

                print('Extracted connected areas from ' + rg_surface)


# print('TEST completed')
print ('STATUS 3: extract connected regions run successfully for PNW')
print('end')


exit()
####
#EXTENTS
                # with arcpy.EnvManager(extent="-124.013841065705 46.9204239753356 -121.254802230646 49.1018394753356", cellSize="MINOF"): #WA PNW extent
                # out_raster = arcpy.sa.ExtractByAttributes(r"C:\Users\annik\Documents\UCS\Data\2022\permanent_inundation\west_coast\west_coast.gdb\rg_merged_raw_raster_surface_26x_2022_high_west_coast", "Value = 1 Or Value = 2"); out_raster.save(r"C:\Users\annik\Documents\ArcGIS\Projects\UCS_Flood_2022\UCS_Flood_2022.gdb\Extract_rg_m1")

                # with arcpy.EnvManager(extent="-120.582289787904 33.4127402244847 -118.96721832204 34.1605909795068", cellSize="MINOF", mask=r"C:\Users\annik\Documents\UCS\Data\2022\DEMs_regional\west_coast\west_coast_dem.gdb\CA_ChannelIslands_GCS_5m_NAVDm"):  #CA Channel Isles extent and mask
