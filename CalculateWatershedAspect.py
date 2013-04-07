# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Name:     CalculateWatershedAspect.py
# Author:   Pabitra Dash (pabitra.dash@usu.edu)
# Created on: 2013-01-14 09:17:55.00000
#
# Description:
#   Generates wateshed aspect netcdf file from the watershed  DEM raster file
# ---------------------------------------------------------------------------

# this says all code at indent level 0 is part of the main() function
def main():
    pass

# this says run  all code at indent level 0 when this script file is ran as standalone script
# meaning not imported in another script
if __name__ == '__main__':
    main()

import arcgisscripting
import os
import sys
import traceback
import glob
import shutil

# Local variables:
InWS_DEM_File = None
OutWSDEMAspectFileName = "ws_aspect.tif" # temporary file
OutAspectNetCDFFileName = None

# settings for runnning this code locally. To run this code on remote app server comment out the following 5 lines
# To run this code locally, uncomment the following 5 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\Temp\ws_dem.tif')
##argumentList.append('ws_aspect.nc')
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 2 more argument, so total of 3
if (len(sys.argv) < 3):
    print('Invalid arguments:')
    print('1st argument: Input watershed DEM file name with file path.')
    print('2nd argument: Output watershed aspect netcdf file name.')
    raise Exception("There has to be 2 arguments to calculate aspect at each grid of the watershed.")
    exit()

InWS_DEM_File = sys.argv[1]
OutAspectNetCDFFileName = sys.argv[2]

# check if provided DEM file exists
if(os.path.isfile(InWS_DEM_File) == False):
    raise Exception("Specified watershed DEM file ({0}) was not found.".format(InWS_DEM_File))
    exit()

filePath = os.path.dirname(InWS_DEM_File)
OutAspectNetCDFFile = os.path.join(filePath, OutAspectNetCDFFileName)
OutWSDEMAspectFile = os.path.join(filePath, OutWSDEMAspectFileName)

# delete any pre-existing slope raster files
aspectRasterFilePath = os.path.dirname(OutWSDEMAspectFile)
os.chdir(aspectRasterFilePath)
aspectRasterFiles = glob.glob('ws_aspect.*')

for filename in aspectRasterFiles:
    os.unlink(filename)

# create the Geoprocessor object
gp = arcgisscripting.create()

# check out any necessary licenses
gp.CheckOutExtension("spatial")

try:
    # Process: Aspect
    gp.Aspect_sa(InWS_DEM_File, OutWSDEMAspectFile)

    # generate netcdf aspect file
    InWSDEMAspectFile = OutWSDEMAspectFile
    variable = "aspect"
    units = "degree"
    XDimension = "x"
    YDimension = "y"
    bandDimension = ""

    # Process: RasterToNetCDF
    gp.RasterToNetCDF_md(InWSDEMAspectFile, OutAspectNetCDFFile, variable, units,
                        XDimension, YDimension, bandDimension)
    print('Done..')

except:
    tbinfo = traceback.format_tb(tb)[0]
    pyErrMsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
    print(pyErrMsg)
    print('>>>Done...with exception')
    raise Exception(pyErrMsg)
