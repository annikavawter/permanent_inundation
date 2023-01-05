print('EXTRACT CONNECTED REGION GROUPS SCRIPT')

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
            rg_surface = arcpy.ListRasters('rg_merged_raw_raster_surface_{0}x_{1}_{2}_{3}' .format(flood_frequency, year, projection, region))[0]
            print('File to extract is ' + rg_surface)

            # ----------------------------------------------------------------------------------------
            # set path and file out
            gdb_out = 'C:/Users/annik/Documents/UCS/Data/2022/permanent_inundation/{0}/{0}.gdb' .format(region)
            outname_extract = 'extract_' + rg_surface

            # ----------------------------------------------------------------------------------------
            # set areas for extraction - mainland (maximum[0]) and Channel Isles (maximum[2])
            print('...finding areas to extract from ' + rg_surface + '...')

            if region == 'west_coast':
                
                # Index maximum counts for mainland and Channel Isles
                arr = arcpy.da.FeatureClassToNumPyArray(rg_surface, ('Value', 'Count'))
                count = arr['Count']
                value = arr['Value']
                
                index_to_extract_imax1 = numpy.argsort(count)[-1]           #;print('imax1 = ' + str(index_to_extract_imax1))
                value_to_extract_imax1 = str(value[index_to_extract_imax1]) #;print('val1 = ' + str(value_to_extract_imax1))

                index_to_extract_imax3 = numpy.argsort(count)[-3]           #;print('imax3 = ' + str(index_to_extract_imax3))
                value_to_extract_imax3 = str(value[index_to_extract_imax3]) #;print('val3 = ' + str(value_to_extract_imax3))

                print('   Result: first and third maximum pixel counts correspond to values ' + value_to_extract_imax1 + ' (mainland) & ' + value_to_extract_imax3 + ' (Channel Isles)')

                #build SQL clause syntax
                print('...building SQL string...')
                inSQLClause = 'Value = ' + value_to_extract_imax1 + ' Or Value = ' + value_to_extract_imax3  #new code line, CHECKING, values from regional conditions

                print('Extracting values ' + value_to_extract_imax1 + ' & ' + value_to_extract_imax3 + ' from ' + rg_surface)
                print('...extracting values...')

            # ----------------------------------------------------------------------------------------
            # extract connected areas       west start 2:27 end 2:54 ~30 min
            runextract=1
            if runextract==1:
                attExtract = ExtractByAttributes(rg_surface, inSQLClause)
                attExtract.save(outname_extract)

                print('Extracted connected areas from ' + rg_surface)
            

            if region=='west_coast':
                print ('STATUS 3: extract connected regions run successfully for mainland west coast')
            elif region=='east_coast':
                exit()
                # print ('STATUS 3: extract connected regions run successfully for _')

t2 = datetime.datetime.now() ; print('Time start: ' + str(t1) + '  &  Time end: ' + str(t2))

# print('TEST completed')
# print ('STATUS 3: extract connected regions run successfully)
print('end')


exit()

#INDEXING NOTES
# #use numpy to find mainland and Channel Isles count indices and values, test using west coast 2022 high
# # values: 46760 is the mainland max, 27198 is the partial PNW area, 85746 is the Channel Isles area
#
# temp3 = numpy.argsort(count)
# imax1 = temp3[-1] ; cmax1 = count[imax1]  ; vmax1 = value[imax1]    ;print('imax1 = ' + str(imax1) + ' where count = ' + str(cmax1) + ' and value = ' + str(vmax1) )   #mainland
# imax3 = temp3[-3] ; cmax3 = count[imax3]  ; vmax3 = value[imax3]    ;print('imax3 = ' + str(imax3) + ' where count = ' + str(cmax3) + ' and value = ' + str(vmax1) )   #channel isles
#
# testcount_maxmainland = count[imax1]                          ; print(testcount_maxmainland)    #gets 922552459.0 which is first max count
# testcount_maxmainland2 = count[numpy.argsort(count)][-1]      ; print(testcount_maxmainland2)   #gets 922552459.0
# testcount_maxpnwpartial = count[temp3[-2]]                    ; print(testcount_maxpnwpartial)  #gets 134488393.0 which is second max count
# testcount_maxchannel = count[imax3]                           ; print(testcount_maxchannel)     #gets 82602506.0 which is third max count
#
# testvalue_maxmainland = value[imax1]                          ; print(testvalue_maxmainland)    #gets 46760 which is value of first max count
# testvalue_maxmainland2 = value[numpy.argsort(count)][-1]      ; print(testvalue_maxmainland2)   #gets 46760
# testvalue_maxpnwpartial = value[temp3[-2]]                    ; print(testvalue_maxpnwpartial)  #gets 27198 which is value of second max count
# testvalue_maxchannel = value[imax3]                           ; print(testvalue_maxchannel)     #gets 85746 which is value of third max count            