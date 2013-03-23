# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# CreateWatershedNLCDDataSet.py
# Created on: 2013-01-10 10:23:37.00000
#
# Description:Creates a NLCD dataset file (.img) for the given whatershed shape file based on
# on an original NLCD dataset file that covers at least the target watershed and saved in the same
# floder as the ws shape file
# ---------------------------------------------------------------------------


# this says all code at indent level 0 is part of the main() function
def main():
    pass

# this says run  all code at indent level 0 when this script file is ran as standalone script
# meaning not imported in another script
if __name__ == '__main__':
    main()

# Import arcpy module
import arcpy
import arcgisscripting
import os
import sys
import traceback
from arcpy import env
from arcpy.sa import *
# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# Local variables:
ProjectedNLCDDataSetFile = None
##ExtractedNLCDFileName = "ExtractedWSNLCD.tif" #temporary file
outWSNLCDFileName = None
##ResampledNLCDFileName = "ResampledWSNLCD.img"
##ElevationLayerNetCDF = None
BufferedWSFile = None
WSDEMFile = None
##ExtractedNLCDFile = None
outWSNLCDFile = None
##ResampledNLCDFile = None
##ClippedWSDEMFile = None

# settings for runnning this code locally. To run this code on remote app server comment out the following 7 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\NLCDDataSetUSA\ProjNLCD2006_LC_N36W096_v1.img')
##argumentList.append(r'E:\CIWaterData\Temp\ws_buffered.shp')
##argumentList.append('E:\CIWaterData\Temp\ws_dem.tif')
##argumentList.append('ws_nlcd.img')
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 4 more arguments, so total of 5
if (len(sys.argv) < 5):
    print('Invalid arguments:')
    print('1st argument: Input NLCD datatset file name with path.')
    print('2nd argument: Input buffered watershed shape file name with path.')
    print('3rd argument: Input watershed DEM file name with path.')
    print('4th argument: Output watershed NLCD data set file name.')
    raise Exception("Exception: There has to be 4 arguments to generate land cover dataset for the watershed." + str(len(sys.argv)))
    exit()

# retrieve passed arguments
ProjectedNLCDDataSetFile = sys.argv[1]
BufferedWSFile = sys.argv[2]
WSDEMFile = sys.argv[3]
outWSNLCDFileName = sys.argv[4]
##CellSize = None
##if(len(sys.argv) > 5):
##    CellSize = sys.argv[5]

# check if provided projected lncd dataset file exists
if(os.path.isfile(ProjectedNLCDDataSetFile) == False):
    print('Exception')
    raise Exception("Exception: Specified projected NLCD dataset file was not found:" + ProjectedNLCDDataSetFile)
    exit()

# check if provided ws DEM file exists
if(os.path.isfile(WSDEMFile) == False):
    print('Exception')
    raise Exception("Exception: Specified watershed DEM file was not found:" + WSDEMFile)
    exit()

# check if provided buffered watershed file exists
if(os.path.isfile(BufferedWSFile) == True):
    filePath = os.path.dirname(BufferedWSFile)
    # set the path for the temporary extracted ws DEM file
##    ExtractedNLCDFile = os.path.join(filePath, ExtractedNLCDFileName)
    # set the path for the output watershed nlcd dataset file
    outWSNLCDFile = os.path.join(filePath, outWSNLCDFileName)
##    ResampledNLCDFile = os.path.join(filePath, ResampledNLCDFileName)
##    WSDEMFile = os.path.join(filePath, WSDEMFileName)
else:
    print('Exception')
    raise Exception("Exception: Specified buffered watershed shape file was not found:" + BufferedWSFile)
    exit()

try:
    # if there exists a previously extracted DEM file delete it
##    if(os.path.isfile(ExtractedNLCDFile) == True):
##        os.unlink(ExtractedNLCDFile)

    # if there exists a previously created watershed sepecfic nlcd datatset delete it
    if(os.path.isfile(outWSNLCDFile) == True):
        os.unlink(outWSNLCDFile)

    # if there exists a previously resampled DEM file delete it
##    if(os.path.isfile(ResampledNLCDFile) == True):
##        os.unlink(ResampledNLCDFile)

    # if there exists a previously projected nlcd dataset file delete it
##    if(os.path.isfile(ProjectedNLCDDataSetFile) == True):
##        os.unlink(ProjectedNLCDDataSetFile)

    #1. project the original shape file to the coordinate system of the DEM file
    # get spatial reference (cooordinate system) from the DEM file
    # ref: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//000v000000p6000000
    WSdesc = arcpy.Describe(BufferedWSFile)
    wsCS = ""
    wsCS = WSdesc.spatialReference.name
    wsCS = wsCS.replace("_", " ")
    if(not wsCS):
        # set output coordinate system if DEM file does not have coordinate system
        outCS = arcpy.SpatialReference('NAD 1983 UTM Zone 12N')
    else:
        outCS = arcpy.SpatialReference(wsCS)

    # Create the Geoprocessor object
    gp = arcgisscripting.create()

     # Check out any necessary licenses
    gp.CheckOutExtension("spatial")

    gp.SnapRaster = WSDEMFile

    # Process: Extract by Rectangle (NOT needed)
    # Ref: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//009z0000002r000000.htm
    # get the rectangular boundary of the WS shape file
    shapeFileDesc = arcpy.Describe(BufferedWSFile)

    # Process: Clip the extracted DEM file to size of the buffered WS shape file
    wsBoundingBox = str(shapeFileDesc.extent.XMin) + " " + str(shapeFileDesc.extent.YMin) + " " + str(shapeFileDesc.extent.XMax) + " " + str(shapeFileDesc.extent.YMax)

    # ref:http://webhelp.esri.com/arcgisdesktop/9.3/index.cfm?TopicName=Clip_(Data_Management)
    gp.clip_management(ProjectedNLCDDataSetFile, wsBoundingBox, outWSNLCDFile, clipping_geometry ="ClippingGeometry")

    # Process: Resample
    # Resampling is necessary if we give the option for the user to provide a different
    # cell size than the cell size of the original dem file (OriginalDEMCoveringWSFile)
    # Ref: http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=resample_(data_management)
##    if(CellSize != None):
##        gp.Resample_management(outWSNLCDFile, ResampledNLCDFile, str(CellSize), "NEAREST")

    print('>>>done..')

except:
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pyErrMsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
    print(pyErrMsg)
    print('>>>done...with exception')
    raise Exception(pyErrMsg)
