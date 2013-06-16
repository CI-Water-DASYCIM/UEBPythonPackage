#-------------------------------------------------------------------------------
# Name:         WatershedFeaturesToNetCDFConversion.py
# Purpose:      Generate a watershed netcdf file from a given watershed shape file
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
##watershedDEMRasterFile = None
bufferedWatershedRasterFile = None
##clippedWatershedRasterFileName = 'ws_buffered_clip.tif' # temporary file name
##clippedWatershedRasterFile = None
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
##    print('2nd argument: Input watershed DEM raster file name with path')
    print('3rd argument: Output watershed netcdf file name with path')
    raise Exception("There has to be 2 arguments to convert a shape file to netCDF file.")
    exit()

# retrieve passed arguments
bufferedWatershedRasterFile = sys.argv[1]
##watershedDEMRasterFile = sys.argv[2]
netCDFWatershedFile = sys.argv[2]

try:
##    # check if provided reference DEM file exists
##    if(os.path.isfile(referenceDEMFile) == False):
##        sys.exit("Specified reference DEM file was not found:" + referenceDEMFile)

    # check if provided watershed DEM raster file exists
##    if(os.path.isfile(watershedDEMRasterFile) == False):
##        sys.exit("Specified reference DEM file was not found:" + watershedDEMRasterFile)

    # check if provided buffered watershed raster file exists
    if(os.path.isfile(bufferedWatershedRasterFile) == False):
        sys.exit("Specified buffered watershed raster file was not found:" + bufferedWatershedRasterFile)
    else:
        filePath = os.path.dirname(bufferedWatershedRasterFile)
##        clippedWatershedRasterFile = os.path.join(filePath, clippedWatershedRasterFileName)

    # if the temporary output watershed clipped raster file already exist, then delete it
##    if(os.path.isfile(clippedWatershedRasterFile) == True):
##        os.unlink(clippedWatershedRasterFile)

    # if the output watershed netcdf file already exist, then delete it
    if(os.path.isfile(netCDFWatershedFile) == True):
        os.unlink(netCDFWatershedFile)

    # get spatial reference (cooordinate system) from the watershed DEM raster file
    # ref: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//000v000000p6000000
##    DEMdesc = arcpy.Describe(watershedDEMRasterFile)
##    demCS = ""
##    demCS = DEMdesc.spatialReference.name
##    demCS = demCS.replace("_", " ")
##    if(not demCS):
##        # set output coordinate system if DEM file does not have coordinate system
##        outCS = arcpy.SpatialReference('NAD 1983 UTM Zone 12N')
##    else:
##        outCS = arcpy.SpatialReference(demCS)
##
##    gp = arcgisscripting.create()
##
##    # check out any necessary licenses
##    gp.CheckOutExtension("spatial")
##
##    gp.SnapRaster = watershedDEMRasterFile
##
##    # reclip the buffered watershed raster file to match with the size of the watershed dem raster file
##    wsBoundingBox = str(repr(DEMdesc.extent.XMin)) + " " + str(repr(DEMdesc.extent.YMin)) + " " + str(repr(DEMdesc.extent.XMax)) + " " + str(repr(DEMdesc.extent.YMax))
##
##    # ref:http://webhelp.esri.com/arcgisdesktop/9.3/index.cfm?TopicName=Clip_(Data_Management)
##    gp.clip_management(bufferedWatershedRasterFile, wsBoundingBox, clippedWatershedRasterFile, clipping_geometry ="NONE")
##
##    # delete the input buffered watershed raster file
##    os.unlink(bufferedWatershedRasterFile)
##
##    # rename the clipped watershed raster file to the original name (input watershed raster file name)
##    os.rename(clippedWatershedRasterFile, bufferedWatershedRasterFile)

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




