# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Name:     CreateWatershedDEMFile.py
# Author:   Pabitra Dash (pabitra.dash@usu.edu)
# Created on: 2013-01-10 10:23:37.00000
#
# Description:
#   Creates a DEM raster file for the given whatershed shape file based on
#   on a source DEM raster file that covers at least the target watershed
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
ClippedWSDEMRasterFileName = None
BufferedWSRasterFile = None
ResampledWSDEMFileName = "ws_resample_dem.tif" # temporary file
gp = None

# settings for runnning this code locally not part of the workflow. To run this code on remote app server as part of the workflow
# comment out the following 6 lines
# to run locally not part of a workflow, uncomment the following 6 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\DEM\gsl100.tif')
##argumentList.append(r'E:\CIWaterData\Temp\ws_buffered.tif')
##argumentList.append('ws_dem.tif')
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 3 more arguments, so total of 4
if (len(sys.argv) < 4):
    print('Invalid arguments:')
    print('1st argument: Input DEM raster file name with path from which the watershed DEM raster file to be created.')
    print('2nd argument: Input buffered watershed raster file name with file path')
    print('3rd argument: Output watershed DEM raster file name')
    raise Exception("There has to be 3 arguments to generate DEM file for the watershed.")
    exit()

# retrieve passed arguments
OriginalDEMCoveringWSFile = sys.argv[1]
BufferedWSRasterFile = sys.argv[2]
ClippedWSDEMRasterFileName = sys.argv[3]
CellSize = None
if(len(sys.argv) > 4):
    CellSize = sys.argv[4]

# check if provided DEM file exists
if(os.path.isfile(OriginalDEMCoveringWSFile) == False):
    print('Exception')
    raise Exception("Specified source DEM file ({0}) was not found.".format(OriginalDEMCoveringWSFile))
    exit()

ExtractedWSDEMFile = None
ClippedWSDEMRasterFile = None
ResampledWSDEMFile = None

# check if provided buffered watershed file exists
if(os.path.isfile(BufferedWSRasterFile) == True):
    filePath = os.path.dirname(BufferedWSRasterFile)

    # set the path for the temporary extracted ws DEM file
    ExtractedWSDEMFile = os.path.join(filePath, ExtractedWSDEMFileName)

    # set the path for the temporary clipped ws DEM file
    ClippedWSDEMRasterFile = os.path.join(filePath, ClippedWSDEMRasterFileName)
    ResampledWSDEMFile = os.path.join(filePath, ResampledWSDEMFileName)
else:
    print('Exception')
    raise Exception("Specified buffered watershed shape file ({0}) was not found.".format(BufferedWSRasterFile))
    exit()

try:
    # if there exists a previously extracted DEM file delete it
    if(os.path.isfile(ExtractedWSDEMFile) == True):
        os.unlink(ExtractedWSDEMFile)

    # if there exists a previously clipped DEM file delete it
    if(os.path.isfile(ClippedWSDEMRasterFile) == True):
        os.unlink(ClippedWSDEMRasterFile)

    # if there exists a previously resampled DEM file delete it
    if(os.path.isfile(ResampledWSDEMFile) == True):
        os.unlink(ResampledWSDEMFile)

    gp = arcgisscripting.create()

    # check out any necessary licenses
    gp.CheckOutExtension("spatial")

    gp.SnapRaster = BufferedWSRasterFile #OriginalDEMCoveringWSFile

    # get the rectangular boundary of the shape file
    rasterFileDesc = arcpy.Describe(BufferedWSRasterFile)

    # Process: Extract by rectangle
    gp.ExtractByRectangle_sa(OriginalDEMCoveringWSFile, BufferedWSRasterFile, ExtractedWSDEMFile, "INSIDE")

    # Process: Clip the extracted DEM file to size of the buffered shape file
    # repr function is used to preserve the floating value precision
    wsBoundingBox = str(repr(rasterFileDesc.extent.XMin)) + " " + str(repr(rasterFileDesc.extent.YMin)) + " " + str(repr(rasterFileDesc.extent.XMax)) + " " + str(repr(rasterFileDesc.extent.YMax))

    # ref:http://webhelp.esri.com/arcgisdesktop/9.3/index.cfm?TopicName=Clip_(Data_Management)
    gp.clip_management(ExtractedWSDEMFile, wsBoundingBox, ClippedWSDEMRasterFile, clipping_geometry ="NONE")

    # Process: Resample if cell size has been provided
    # Resampling is necessary if we give the option for the user to provide a different
    # cell size than the cell size of the original dem file (OriginalDEMCoveringWSFile)
    # Ref: http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=resample_(data_management)
    if(CellSize != None):
        gp.Resample_management(ClippedWSDEMRasterFile, ResampledWSDEMFile, str(CellSize), "NEAREST")
        # delete the clipped file
        os.unlink(ClippedWSDEMRasterFile)
        # rename the resample dem file as the clipped file
        os.rename(ResampledWSDEMFile, ClippedWSDEMRasterFile)

    print('>>>Done...')
except:
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pyErrMsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
    print(pyErrMsg)
    print('>>>done...with exception')
    raise Exception(pyErrMsg)

finally:
    # check in any necessary licenses
    arcpy.CheckInExtension("spatial")
    # check in any necessary licenses
    if(gp != None):
        gp.CheckInExtension("spatial")
        del gp