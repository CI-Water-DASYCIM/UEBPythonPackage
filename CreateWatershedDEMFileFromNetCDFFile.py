# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Name: CreateWatershedDEMFileFromNetCDFFile.py
# Author: Pabitra Dash (pabitra.dash@usu.edu)
# Created on: 2013-07-07
#
# Description:
#   Creates watershed DEM raster file from an input watershed NetCDF file
#   Uses 'NAD 1983 UTM Zone 12N' as the coordinate system for the generated
#   raster file
# ---------------------------------------------------------------------------

# this says all code at indent level 0 is part of the main() function
def main():
    pass

# this says run  all code at indent level 0 when this script file is ran as standalone script
# meaning not imported in another script
if __name__ == '__main__':
    main()


import arcpy
import sys
import os
import traceback

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

#local variables
inWatershedNetCDFFile = None
outWatershedDEMRasterFile = None
outWatershedLayer = 'watershed_layer' # temporary in memory file

# settings for runnning this code locally not part of the workflow. To run this code on remote app server as part of the workflow
# comment out the following 5 lines
# to run locally not part of a workflow, uncomment the following 5 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\Scratch\watershed.nc')
##argumentList.append(r'E:\Scratch\ws_dem.tif')
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 2 more arguments, so total of 3
if (len(sys.argv) < 3):
    print('Invalid arguments:')
    print('1st argument: Input watershed NetCDF file name with path from which the watershed DEM raster file to be created.')
    print('2nd argument: Output watershed DEM raster file name with path')
    raise Exception("There has to be 2 arguments to generate DEM raster file for the watershed.")
    exit()

# retrieve passed arguments
inWatershedNetCDFFile = sys.argv[1]
outWatershedDEMRasterFile = sys.argv[2]

# check that the provided input watershed netcdf file has extension '.nc'
filename, extension = os.path.splitext(inWatershedNetCDFFile)
if extension != '.nc':
    raise Exception("Input watershed file ({0}) is not a netCDF file.".format(inWatershedNetCDFFile))
    exit()

# check that the provided output watershed raster file has extension '.tif'
filename, extension = os.path.splitext(outWatershedDEMRasterFile)
if extension != '.tif':
    raise Exception("Output watershed raster file ({0}) is not a tif file.".format(outWatershedDEMRasterFile))
    exit()

# check if provided watershed NetCDF file exists
if(os.path.isfile(inWatershedNetCDFFile) == False):
    raise Exception("Input watershed NetCDF file ({0}) was not found.".format(inWatershedNetCDFFile))
    exit()


try:
    # if there exists a previously extracted DEM file delete it
    if(os.path.isfile(outWatershedDEMRasterFile) == True):
        os.unlink(outWatershedDEMRasterFile)

    # delete if the in-memory raster layer still exists
    if outWatershedLayer != None:
        arcpy.Delete_management(outWatershedLayer)

    # Process: Make NetCDF Raster Layer
    dataVariableName = 'watershed'
    dataXdimension = 'x'
    dataYdimension = 'y'
    dataBandDimension = ''
    dataDimensionValues = ''
    dataValueSelectionMethod = 'BY_VALUE'
    arcpy.MakeNetCDFRasterLayer_md(inWatershedNetCDFFile, dataVariableName, dataXdimension, dataYdimension, outWatershedLayer, dataBandDimension, dataDimensionValues, dataValueSelectionMethod)

    # Process: Project Raster
    # Ref: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//00170000007q000000
    inWatershedLayer = outWatershedLayer
    # set output coordinate system  (the input watershed necdf file does not have the coordinate system information)
    outCS = arcpy.SpatialReference('NAD 1983 UTM Zone 12N')
    # resampling type
    resamplingType = 'BILINEAR' # Bilinear interpolation
    #set transformation method
    transformMethod = 'WGS_1984_(ITRF00)_To_NAD_1983'
    arcpy.ProjectRaster_management(inWatershedLayer, outWatershedDEMRasterFile, outCS, resamplingType, transformMethod, "", "")
    print('Done')
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