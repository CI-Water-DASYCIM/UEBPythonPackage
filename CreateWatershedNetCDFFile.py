#-------------------------------------------------------------------------------
# Name:         CreateWatershedNetCDFFile.py
# Purpose:      To generate a watershed netcdf file from a given watershed shape file
#
#
# Author:      Pabitra
#
# Created:     21/02/2013
# Copyright:   (c) Pabitra 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# set desktop license used to ArcView
# ref: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//002z0000000z000000
import arcview

import arcpy
import arcgisscripting
import os
import sys
import traceback

# check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# local variables
bufferedWatershedRasterFile = None
netCDFWatershedFile = None
gp = None

# settings for runnning this code locally not part of the workflow. To run this code on remote app server as part of the workflow
# comment out the following 6 lines
# to run locally not part of a workflow, uncomment the following 6 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\Scratch\TestwatershedDEMRasterFileRowsCols\watershed_buffered.tif') # r'E:\CIWaterData\Temp\ws_buffered.tif'
####argumentList.append(r'E:\Scratch\TestwatershedDEMRasterFileRowsCols\ws_dem.tif') # r'E:\CIWaterData\Temp\ws_dem.tif'
##argumentList.append(r'E:\Scratch\TestwatershedDEMRasterFileRowsCols\watershed.nc') # r'E:\CIWaterData\Temp\watershed.nc'
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 2 more arguments, so total of 3
if (len(sys.argv) < 3):
    print('Invalid arguments:')
    print('1st argument: Input buffered watershed raster file name with path')
    print('2nd argument: Output watershed netcdf file name with path')
    raise Exception("There has to be 2 arguments to convert a shape file to netCDF file.")
    exit()

# retrieve passed arguments
bufferedWatershedRasterFile = sys.argv[1]
netCDFWatershedFile = sys.argv[2]

try:
    # check if provided buffered watershed raster file exists
    if(os.path.isfile(bufferedWatershedRasterFile) == False):
        sys.exit("Specified buffered watershed raster file was not found:" + bufferedWatershedRasterFile)
    else:
        filePath = os.path.dirname(bufferedWatershedRasterFile)

    # if the output watershed netcdf file already exist, then delete it
    if(os.path.isfile(netCDFWatershedFile) == True):
        os.unlink(netCDFWatershedFile)

    # run raster to NetCDF conversion
    variable = "watershed"
    units = "meter"
    XDimension = "x"
    YDimension = "y"
    bandDimension = ""
    arcpy.RasterToNetCDF_md(bufferedWatershedRasterFile, netCDFWatershedFile, variable, units, XDimension, YDimension, bandDimension)

    print('NetCDF watershed domain file was created:' + netCDFWatershedFile)
    print('Done...')

except:
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pyErrMsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
    print(pyErrMsg)
    print('>>>Done...with exception')
    raise Exception(pyErrMsg)

finally:
    # check in any necessary licenses
    arcpy.CheckInExtension("spatial")
    if(gp != None):
        gp.CheckInExtension("spatial")
        del gp




