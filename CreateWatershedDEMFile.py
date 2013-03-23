# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# CreateWatershedDEMFile.py
# Created on: 2013-01-10 10:23:37.00000
#
# Description:Creates a DEM file for the given whatershed shape file based on
# on an original DEM file that covers at least the target watershed
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
import sys
import os
import traceback
from arcpy import env
# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")


# Local variables:
OriginalDEMCoveringWSFile = None
ExtractedWSDEMFileName = "ExtractedWSDEM.tif" # temporary file
ClippedWSDEMFileName = None
ResampledWSDEMFileName = "ws_resample_dem.tif" # temporary file
BufferedWSFile = None

# settings for runnning this code locally. To run this code on remote app server comment out the following 5 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\DEM\gsl100.tif')
##argumentList.append(r'E:\CIWaterData\Temp\ws_buffered.shp')
##argumentList.append('ws_dem.tif')
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 3 more arguments, so total of 4
if (len(sys.argv) < 4):
    print('Invalid arguments:')
    print('1st argument: Input DEM raster file name with path from which the watershed DEM raster file to be created.')
    print('2nd argument: Input buffered watershed shape file name with file path')
    print('3rd argument: Output watershed DEM raster file name')
    raise Exception("There has to be 3 arguments to generate DEM file for the watershed.")
    exit()

# retrieve passed arguments
OriginalDEMCoveringWSFile = sys.argv[1]
BufferedWSFile = sys.argv[2]
ClippedWSDEMFileName = sys.argv[3]
CellSize = None
if(len(sys.argv) > 4):
    CellSize = sys.argv[4]

# check if provided DEM file exists
if(os.path.isfile(OriginalDEMCoveringWSFile) == False):
    print('Exception')
    raise Exception("Specified original watershed shape file was not found:" + OriginalDEMCoveringWSFile)
    exit()

ExtractedWSDEMFile = None
ClippedWSDEMFile = None
ResampledWSDEMFile = None

# check if provided buffered watershed file exists
if(os.path.isfile(BufferedWSFile) == True):
    filePath = os.path.dirname(BufferedWSFile)
    # set the path for the temporary extracted ws DEM file
    ExtractedWSDEMFile = os.path.join(filePath, ExtractedWSDEMFileName)
    # set the path for the temporary clipped ws DEM file
    ClippedWSDEMFile = os.path.join(filePath, ClippedWSDEMFileName)
    ResampledWSDEMFile = os.path.join(filePath, ResampledWSDEMFileName)
else:
    print('Exception')
    raise Exception("Specified buffered watershed shape file was not found:" + BufferedWSFile)
    exit()

try:
    # if there exists a previously extracted DEM file delete it
    if(os.path.isfile(ExtractedWSDEMFile) == True):
        os.unlink(ExtractedWSDEMFile)

    # if there exists a previously clipped DEM file delete it
    if(os.path.isfile(ClippedWSDEMFile) == True):
        os.unlink(ClippedWSDEMFile)

    # if there exists a previously resampled DEM file delete it
    if(os.path.isfile(ResampledWSDEMFile) == True):
        os.unlink(ResampledWSDEMFile)

    gp = arcgisscripting.create()

    gp.SnapRaster = OriginalDEMCoveringWSFile

    # get the rectangular boundary of the shape file
    shapeFileDesc = arcpy.Describe(BufferedWSFile)

    # Process: Extract by Rectangle
    gp.ExtractByRectangle_sa(OriginalDEMCoveringWSFile, BufferedWSFile, ExtractedWSDEMFile, "INSIDE")

    # Process: Clip the extracted DEM file to size of the buffered shape file
    wsBoundingBox = str(shapeFileDesc.extent.XMin) + " " + str(shapeFileDesc.extent.YMin) + " " + str(shapeFileDesc.extent.XMax) + " " + str(shapeFileDesc.extent.YMax)

    # ref:http://webhelp.esri.com/arcgisdesktop/9.3/index.cfm?TopicName=Clip_(Data_Management)
    gp.clip_management(ExtractedWSDEMFile,wsBoundingBox,ClippedWSDEMFile, clipping_geometry ="NONE")

    # Process: Resample
    # Resampling is necessary if we give the option for the user to provide a different
    # cell size than the cell size of the original dem file (OriginalDEMCoveringWSFile)
    # Ref: http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=resample_(data_management)
    if(CellSize != None):
        gp.Resample_management(ClippedWSDEMFile, ResampledWSDEMFile, str(CellSize), "NEAREST")
        # delete the clipped file
        os.unlink(ClippedWSDEMFile)
        # rename the resample dem file as the clipped file
        os.rename(ResampledWSDEMFile, ClippedWSDEMFile)

    print('>>>Done...')
except:
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pyErrMsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
    print(pyErrMsg)
    print('>>>done...with exception')
    raise Exception(pyErrMsg)