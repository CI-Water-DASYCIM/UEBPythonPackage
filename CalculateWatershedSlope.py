#-------------------------------------------------------------------------------
# Name:     CalculateWatershedSlope/py
# Author:   Pabitra Dash (pabitra.dash@usu.edu)
# Purpose:
#   Generates watershed slope netcdf file and saves to the same
#   directory where the ws dem file exists
#
# Created:     13/01/2013
# Copyright:   (c) 2013
# Licence:     <your licence>
# Ref:http://webhelp.esri.com/arcgisdesktop/9.3/index.cfm?TopicName=Slope
#-------------------------------------------------------------------------------

# this says all code at indent level 0 is part of the main() function
def main():
    pass

# this says run  all code at indent level 0 when this script file is ran as standalone script
# meaning not imported in another script
if __name__ == '__main__':
    main()

# set desktop license used to ArcView
# ref: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//002z0000000z000000
import arcview

import arcgisscripting
import os
import sys
import traceback
import glob
import shutil

# local variables
InWSDEMRasterFile = None # this is the file generated after regriding the DEM file to match with buffered watershed
OutSlopeNetCDFFileName = None
OutWSDEMSlopeFileName = "ws_slope.tif" # temporary file
gp = None

# settings for runnning this code locally not part of the workflow. To run this code on remote app server as part of the workflow
# comment out the following 5 lines
# to run locally not part of a workflow, uncomment the following 5 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\Temp\ws_dem.tif')
##argumentList.append('ws_slope.nc')
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 2 more argument, so total of 3
if (len(sys.argv) < 3):
    print('Invalid arguments:')
    print('1st argument: Input watershed DEM raster file')
    print('2nd argument: Output watershed slope netcdf file name.')
    raise Exception("There has to be 2 argument to calculate aspect at each grid of the watershed.")
    exit()

# retrieve passed argument
InWSDEMRasterFile = sys.argv[1]
OutSlopeNetCDFFileName = sys.argv[2]

# check if provided DEM file exists
if(os.path.isfile(InWSDEMRasterFile) == False):
    print('Exception')
    raise Exception("Specified watershed DEM file ({0}) was not found.".format(InWSDEMRasterFile))
    exit()

filePath = os.path.dirname(InWSDEMRasterFile)
OutSlopeNetCDFFile = os.path.join(filePath, OutSlopeNetCDFFileName)
OutWSDEMSlopeFile = os.path.join(filePath, OutWSDEMSlopeFileName)

# delete any pre-existing slope raster files
slopeRasterFilePath = os.path.dirname(OutWSDEMSlopeFile)
os.chdir(slopeRasterFilePath)
slopeRasterFiles = glob.glob('ws_slope.*')
for filename in slopeRasterFiles:
    os.unlink(filename)

# create the Geoprocessor object
gp = arcgisscripting.create()

# check out any necessary licenses
gp.CheckOutExtension("spatial")

try:
    InMeasurementType = "DEGREE"
    ZFactor = "1"

    # Process: Slope
    gp.Slope_sa(InWSDEMRasterFile, OutWSDEMSlopeFile, InMeasurementType, ZFactor)

    # generate netcdf slope file
    InWSDEMSlopeFile = OutWSDEMSlopeFile
    variable = "slope"
    units = "degree"
    XDimension = "x"
    YDimension = "y"
    bandDimension = ""

    # Process: RasterToNetCDF
    gp.RasterToNetCDF_md(InWSDEMSlopeFile, OutSlopeNetCDFFile, variable, units,
                        XDimension, YDimension, bandDimension)

    print('>>Done...')
except:
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pyErrMsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
    print(pyErrMsg)
    print('>>>done...with exception')
    raise Exception(pyErrMsg)
finally:
    # check in any necessary licenses
    if(gp != None):
        gp.CheckInExtension("spatial")
        del gp