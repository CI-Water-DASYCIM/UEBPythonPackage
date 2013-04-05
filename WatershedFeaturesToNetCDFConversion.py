#-------------------------------------------------------------------------------
# Name:         WatershedFeaturesToNetCDFConversion.py
# Purpose:      Generate a watershed netcdf file from a given watershed shape file
#               as well genrates the buferred watershed shape file
#
# Author:      Pabitra
#
# Created:     21/02/2013
# Copyright:   (c) Pabitra 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
# Import of arcview is necessary prior to importing arcpy
import arcview
import arcpy
import os
import glob
import shutil
import sys
import traceback

# local variables
referenceDEMFile = None
watershedOriginalShapeFile = None
projectedWatershedShapeFile = None
bufferedWatreshedFile = None
rasterWatershedFilePath = None
netCDFWatershedFile = None
bufferSize = 500
variable = "watershed"
units = "meter"
XDimension = "x"
YDimension = "y"
bandDimension = ""

# settings for runnning this code locally. To run this code on remote app server comment out the following 6 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\DEM\gsl100.tif')
##argumentList.append(r'E:\CIWaterData\Temp\Watershed.shp')
##argumentList.append(r'E:\CIWaterData\Temp\ws_buffered.shp')
##argumentList.append(r'E:\CIWaterData\Temp\watershed.nc')
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 4 more arguments, so total of 5
if (len(sys.argv) < 5):
    print('Invalid arguments:')
    print('1st argument: Input reference DEM raster file name with path.')
    print('2nd argument: Input watershed shape file name with path')
    print('3rd argument: Output watershed buffered shape file name with path')
    print('4th argument: Output watershed netcdf file name with path')
    raise Exception("There has to be 4 arguments to convert a shape file to netCDF file.")
    exit()

# retrieve passed arguments
referenceDEMFile = sys.argv[1]
watershedOriginalShapeFile = sys.argv[2]
bufferedWatreshedFile = sys.argv[3]
netCDFWatershedFile = sys.argv[4]

# check if provided watershed file exists
if(os.path.isfile(watershedOriginalShapeFile) == True):
    # extract the filename from the path of the original shape file and prefix with 'projected'
    # for the projected shape file name
    projectedWatershedShapeFile = 'projected' + os.path.basename(watershedOriginalShapeFile)
##    bufferedWatreshedFile = 'buffered' + os.path.basename(watershedOriginalShapeFile)
    filePath = os.path.dirname(watershedOriginalShapeFile)
    projectedWatershedShapeFile = os.path.join(filePath, projectedWatershedShapeFile)
##    bufferedWatreshedFile = os.path.join(filePath, bufferedWatreshedFile)
    rasterWatershedFilePath = os.path.join(filePath, 'raster')
else:
    print('Exception')
    raise Exception("Exception: Specified original watershed shape file was not found:" + watershedOriginalShapeFile)
    exit()

try:
    # check if provided DEM file exists
    if(os.path.isfile(referenceDEMFile) == False):
        raise Exception("Specified reference DEM file was not found:" + referenceDEMFile)
        exit()

    if(os.path.isfile(netCDFWatershedFile) == True):
        os.unlink(netCDFWatershedFile)

    # delete any pre-existing projected shape files and buffered shape files
    projectedShapeFileDirPath = os.path.dirname(projectedWatershedShapeFile)
    os.chdir(projectedShapeFileDirPath)
    projectedShapeFiles = glob.glob('projected*.*')

    for filename in projectedShapeFiles:
        os.unlink(filename)

    bufferedShapeFiles = glob.glob('ws_buffered.*')
    for filename in bufferedShapeFiles:
        os.unlink(filename)

    # delete the raster folder if exists
    if(os.path.isdir(rasterWatershedFilePath) == True):
        shutil.rmtree(rasterWatershedFilePath)

    #1. project the original shape file to the coordinate system of the DEM file
    # get spatial reference (cooordinate system) from the DEM file
    # ref: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//000v000000p6000000
    DEMdesc = arcpy.Describe(referenceDEMFile)
    demCS = ""
    demCS = DEMdesc.spatialReference.name
    demCS = demCS.replace("_", " ")
    if(not demCS):
        # set output coordinate system if DEM file does not have coordinate system
        outCS = arcpy.SpatialReference('NAD 1983 UTM Zone 12N')
    else:
        outCS = arcpy.SpatialReference(demCS)

    #set transformation method
    transformMethod = 'WGS_1984_(ITRF00)_To_NAD_1983'

    # run the project tool
    arcpy.Project_management(watershedOriginalShapeFile, projectedWatershedShapeFile, outCS, transformMethod)
    print('projected watershed file created:' + projectedWatershedShapeFile)

    # run buffer tool
    bufferSize = str(bufferSize) + ' Meters'
    arcpy.Buffer_analysis(projectedWatershedShapeFile, bufferedWatreshedFile, bufferSize, "FULL", "ROUND", "NONE", "")
    print('Buffered watershed file created:' + bufferedWatreshedFile)

    # run Feature to Raster conversion
    # Note: The grid size of the output raster file(rasterWatershedFilePath) will be same as the referenceDEMFile
    arcpy.FeatureToRaster_conversion(bufferedWatreshedFile, "Id", rasterWatershedFilePath, referenceDEMFile)
    print('Raster watershed file created:' + rasterWatershedFilePath)

    # run raster to NetCDF conversion
    arcpy.RasterToNetCDF_md(rasterWatershedFilePath, netCDFWatershedFile, variable, units, XDimension, YDimension, bandDimension)

    print('NetCDF watershed domain file created:' + netCDFWatershedFile)
    print('Done...')

except:
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pyErrMsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
    print(pyErrMsg)
    print('>>>done...with exception')
    raise Exception(pyErrMsg)




