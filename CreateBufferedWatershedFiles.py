#-------------------------------------------------------------------------------
# Name:         CreateBufferedWatershedFiles.py
# Purpose:      To generate a watershed netcdf file from a given watershed shape file
#               as well generate the buffered watershed shape file
#Note:          The input watershed shape file must need to have one of the following
#               geographic coordinate systems:
#               NAD_1927, NAD_1983, WGS_1984
#               If the input DEM file has a geographic coordinate system then it must be 'NAD 1983'
#
# Author:      Pabitra
#
# Created:     21/02/2013
# Updates:     4/15/2014
# Copyright:   (c) Pabitra 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# set desktop license used to ArcView
# ref: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//002z0000000z000000
import arcview

import arcpy
import os
import glob
import sys
import traceback

# check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# local variables
referenceDEMFile = None
watershedOriginalShapeFile = None
projectedWatershedShapeFile = None
bufferedWatershedShapeFile = None
bufferedWatershedRasterFile = None
bufferSize = None

# settings for running this code locally not part of the workflow. To run this code on remote app server as part of the workflow
# comment out the following 8 lines
# to run locally not part of a workflow, uncomment the following 8 lines
# argumentList = []
# argumentList.append('') #this argument is reserved for the name of this script file
# argumentList.append(r'E:\CIWaterData\DEM\gsl100.tif')
# argumentList.append(r'E:\CIWaterData\Temp\watershed.shp')
# argumentList.append(r'E:\CIWaterData\Temp\watershed_buffered.shp') # don't change this outfile name when testing
# argumentList.append(r'E:\CIWaterData\Temp\watershed_buffered.tif') # don't change this outfile name when testing
# argumentList.append('500') # watershed buffer size
# sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 5 more arguments, so total of 6
if len(sys.argv) < 6:
    print('Invalid arguments:')
    print('1st argument: Input reference DEM raster file name with path')
    print('2nd argument: Input watershed shape file name with path')
    print('3rd argument: Output watershed buffered shape file name with path')
    print('4th argument: Output watershed buffered raster file name with path')
    print('5th argument: Watershed buffer size (in meters)')
    raise Exception("There has to be 5 arguments to create buffered watershed shape and raster files.")
    exit()

# retrieve passed arguments
referenceDEMFile = sys.argv[1]
watershedOriginalShapeFile = sys.argv[2]
bufferedWatershedShapeFile = sys.argv[3]
bufferedWatershedRasterFile = sys.argv[4]
bufferSize = sys.argv[5]

# check if provided watershed shape file exists
if not os.path.isfile(watershedOriginalShapeFile):
    # Look for a file with .shp extension. If no .shp file can be found then raise exception.
    shape_file_lookup_dir = os.path.dirname(watershedOriginalShapeFile)
    # find the first file with extension '.shp' in the lookup_dir
    shape_file_found = False
    for _file in os.listdir(shape_file_lookup_dir):
        if _file.endswith(".shp"):
            watershedOriginalShapeFile = os.path.join(shape_file_lookup_dir, _file)
            shape_file_found = True
            break

    if not shape_file_found:
        print('Exception')
        sys.exit("Exception: Specified original watershed shape file ({0}) was not "
                        "found.".format(watershedOriginalShapeFile))

# extract the filename from the path of the original shape file and prefix with 'projected'
# for the projected shape file name
projectedWatershedShapeFile = 'projected' + os.path.basename(watershedOriginalShapeFile)
filePath = os.path.dirname(watershedOriginalShapeFile)
projectedWatershedShapeFile = os.path.join(filePath, projectedWatershedShapeFile)

if not bufferSize.isdigit():
    print('Exception')
    sys.exit("Exception: Specified buffer size ({0}) is invalid. Buffer size needs to be a "
                    "number.".format(bufferSize))

if float(bufferSize) <= 0:
    print('Exception')
    sys.exit("Exception: Specified buffer size ({0}) is invalid. Buffer size needs to be a "
                    "number greater than zero.".format(bufferSize))

try:
    # check if provided DEM file exists
    if not os.path.isfile(referenceDEMFile):
        sys.exit("Specified reference DEM file ({0}) was not found.".format(referenceDEMFile))

    # delete any pre-existing projected shape files and buffered shape files
    projectedShapeFileDirPath = os.path.dirname(projectedWatershedShapeFile)
    os.chdir(projectedShapeFileDirPath)
    projectedShapeFiles = glob.glob('projected*.*')

    for filename in projectedShapeFiles:
        os.unlink(filename)

    bufferedShapeFiles = glob.glob('watershed_buffered.*')
    for filename in bufferedShapeFiles:
        os.unlink(filename)

    #1. project the original shape file to the coordinate system of the DEM file
    # get spatial reference (cooordinate system) from the DEM file
    # ref: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//000v000000p6000000
    DEMdesc = arcpy.Describe(referenceDEMFile)
    demCS = ""
    demCS = DEMdesc.spatialReference.name
    demCS = demCS.replace("_", " ")
    if not demCS:
        # set output coordinate system if DEM file does not have coordinate system
        outCS = arcpy.SpatialReference('NAD 1983 UTM Zone 12N')
    else:
        outCS = arcpy.SpatialReference(demCS)

    # get the spatial reference of the watershed
    WSdesc = arcpy.Describe(watershedOriginalShapeFile)
    wsCS = ""
    wsCS = WSdesc.spatialReference.name
    wsCS = wsCS.replace("_", " ")
    print(wsCS)

    # if the watershed file does not have a geo-coordinate raise exception
    if wsCS == "Unknown":
        raise Exception('Input watershed has unknown geographic coordinate system.')

    try:
        # run the project tool
        arcpy.Project_management(watershedOriginalShapeFile, projectedWatershedShapeFile, outCS)
        print('projected watershed file created:' + projectedWatershedShapeFile)
    except arcpy.ExecuteError:
        raise Exception("Failed to project the provided watershed file." + arcpy.GetMessage(2))

    # run buffer tool
    bufferSize = str(bufferSize) + ' Meters'
    try:
        arcpy.Buffer_analysis(projectedWatershedShapeFile, bufferedWatershedShapeFile,
                          bufferSize, "FULL", "ROUND", "NONE", "")
        print('Buffered watershed shape file was created:' + bufferedWatershedShapeFile)
    except arcpy.ExecuteError:
        raise Exception("Failed to add buffer to the projected watershed file." + arcpy.GetMessage(2))

    # run Feature to Raster conversion
    # Note: The grid size of the output raster file(bufferedWatershedRasterFile) will be same as the referenceDEMFile
    arcpy.FeatureToRaster_conversion(bufferedWatershedShapeFile, "Id", bufferedWatershedRasterFile, referenceDEMFile)
    print('Buffered watershed raster file created:' + bufferedWatershedRasterFile)

    print('Done...')

except arcpy.ExecuteError:
    print('>>>Done...with exception')
    raise Exception(arcpy.GetMessage(2))
except:
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pyErrMsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + str(sys.exc_type)+ ": " + \
               str(sys.exc_value) + "\n"
    print(pyErrMsg)
    print('>>>Done...with exception')
    raise Exception(pyErrMsg)

finally:
    # check in any necessary licenses
    arcpy.CheckInExtension("spatial")




