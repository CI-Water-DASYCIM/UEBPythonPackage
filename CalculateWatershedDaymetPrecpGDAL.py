# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# CalculateWatershedDaymetPrecpGDAL.py
# Created on: 2013-01-28 14:54:03.00000
# Author: Pabitra Dash (pabitra.dash@usu.edu)
#
# Description:
#   Generates a netcdf precipitation data file for the
#   domain watershed from a number of
#   netcdf data files obtanied from Daymet website and stored locally

# Assumptions:
#   There will be one or more input daymet netcdf files
#   Each of these data files will have 365 time bands
# Dependencies:
#   netCDF4 python library
#   GDAL python library
# ---------------------------------------------------------------------------

# this says all code at indent level 0 is part of the main() function
def main():
    pass

# this says run  all code at indent level 0 when this script file is ran as standalone script
# meaning not imported in another script
if __name__ == '__main__':
    main()

# Import modules
import arcpy
import os
import sys
import traceback
import glob
import arcgisscripting
import datetime
import time
from netCDF4 import Dataset
import gdal_merge
import threading


# Local variables:
sourcePrecpNetCDFilePath = None
outNetCDFFilePath = None
outNetCDFFileName = None
outRasterFilePath = None
clippedWSDEMFile = None
inPrecpDataFileNames = None

# TODO: read these 2 magic strings for dir path from a config file
localPyScriptsPath = r'E:\SoftwareProjects\CIWaterPythonScripts'
remotePyScriptsPath = r'C:\CIWaterPythonScripts'
if(os.path.isdir(remotePyScriptsPath) == True):
    sys.path.append(remotePyScriptsPath)
elif(os.path.isdir(localPyScriptsPath) == True):
    sys.path.append(localPyScriptsPath)
else:
     raise Exception("Script file (DaymetNumberOfDaysToProcess.py) was not found.")
     exit()

import DaymetNumberOfDaysToProcess

# settings for runnning this code locally. To run this code on remote app server comment out the following 9 lines
# To run this code locally uncomment the following 9 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\precdatasets')
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\precdatasets\OutNetCDF')
##argumentList.append('precp_daily_one_data.nc')
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\precdatasets\Raster')
##argumentList.append(r"E:\CIWaterData\Temp\ws_dem.tif")
##argumentList.append('prcp1.nc;prcp2.nc')
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 6 more argument, so total of 7
if (len(sys.argv) < 7):
    print('Invalid arguments:')
    print('1st argument: Input precp netcdf file path')
    print('2nd argument: Output netcdf file path (must be different from input netcdf file path)')
    print('3rd argument: Output netcdf file name')
    print('4th argument: Output raster file path')
    print('5th argument: Input watershed DEM file name with file path')
    print('6th argument: List of input precp netcdf file names')
    raise Exception("There has to be 6 arguments to calculate Dayment precpitation at each grid of the watershed.")
    exit()

# retrieve passed arguments
sourcePrecpNetCDFilePath = sys.argv[1]
outNetCDFFilePath = sys.argv[2]
outNetCDFFileName = sys.argv[3]
outRasterFilePath = sys.argv[4]
clippedWSDEMFile = sys.argv[5]
inPrecpDataFileNames = sys.argv[6]

# array to store all input Daymet netcdf file names
precpNetCDFFilesArray = []
precpNetCDFFilesArray = inPrecpDataFileNames.split(';')

if(len(precpNetCDFFilesArray) == 0):
    raise Exception("Netcdf input Dayment data file names are missing.")
    exit()


# array to store raster layers created from each of the input netcdf prep data file
inMemoryRasterLayersArray = []

# param: directoryPath dir path from which all files need to be deleted
def deleteAllFilesFromDirectory(directoryPath):
    os.chdir(directoryPath)
    prcpRasterFiles = glob.glob('*.*')
    for filename in prcpRasterFiles:
        os.unlink(filename)

# param: inputNetCDF_File netcdef file from which the value for the
#        global attribute start_year need to be read
def getInputDataFileStartYear(inputNetCDF_File):
    rootGrp = Dataset(inputNetCDF_File, 'r', format='NETCDF3_CLASSIC')
    dataYear = None
    for gAttribute in rootGrp.ncattrs():
        if(gAttribute == 'start_year'):
            dataYear = getattr(rootGrp,gAttribute)
            return dataYear

    return dataYear

def deleteAllInMemoryRasterLayers():
    # delete all in memory raster layes
    global inMemoryRasterLayersArray
    for inMemeoryRasterLayer in inMemoryRasterLayersArray:
        try:
            arcpy.Delete_management(inMemeoryRasterLayer)
        except:
            pass

# param: inRasterFilePath the file path where the raster files need to be merged exist
# param: outMergedRasterFile merged output file name and path
# param: inRasterFileNamePatternToMatch file name matching pattern to be used for selecting a
#        and merging raster files
def mergeRasters(inRasterFilePath, outMergedRasterFile, inRasterFileNamePatternToMatch):
    argv = []
    argv.append('')
    argv.append('-q')
    argv.append('-v')
    argv.append('-separate')    #this would add each raster as a separate band
    argv.append('-o')
    argv.append(outMergedRasterFile)

    # select all raster files that matches the filename pattern
    inRasterFile = os.path.join(inRasterFilePath, inRasterFileNamePatternToMatch)
    argv.append(inRasterFile)
    gdal_merge.main(argv)

# param: inputRasterFile raster file from which a netcdf file needs to be generated
# param: outNetCDFFile name and file path for the output netcdf file
def convertRasterToNetCDF(inputRasterFile, outNetCDFFile):
    variable = "Prec"
    units = "" #actual unit for precp data is 'mm/day' but this unit does not work with netcdf create function - correct unit is set after the file is created
    XDimension = "x"
    YDimension = "y"
    bandDimension = "time"

    # create the Geoprocessor object
    gp = arcgisscripting.create()

    # Check out any necessary licenses
    gp.CheckOutExtension("spatial")

    gp.RasterToNetCDF_md(inputRasterFile, outNetCDFFile, variable, units,
                            XDimension, YDimension, bandDimension)

# param: netCDF_File netcdf file for which data variable units to be set
# param: dataSetVariableName name of the variable for which units to be set
# param: units the units to be set for the data variable
def setNetCDFDDataUnits(netCDF_File, dataSetVariableName, units):
    # open the file with read/write mode (r+)
    rootGrp = Dataset(netCDF_File, 'r+', format='NETCDF3_CLASSIC')

    # get the data variable
    taVariable = rootGrp.variables[dataSetVariableName]

    # set the units of temperature variable to kelivin
    taVariable.units = units

    # close the file
    rootGrp.close()

# param: netCDF_File the netcdf file for which the start_year global attribute to be set
# param: dataYear the value for the global attribute start_year
def setNetCDFDataYearAttribute(netCDF_File, dataYear):
    # open the file with read/write mode (r+)
    rootGrp = Dataset(netCDF_File, 'r+', format='NETCDF3_CLASSIC')

    rootGrp.start_year = dataYear

    # close the file
    rootGrp.close()

# param: netCDF_File the netcdf file for which the global attribute data_varaible_name needs to be set
# pram: variableName name of the variable to be assigned to the global attribute data_variable_name
def setNetCDFDataVariableNameAttribute(netCDF_File, variableName):
    # open the file with read/write mode (r+)
    rootGrp = Dataset(netCDF_File, 'r+', format='NETCDF3_CLASSIC')

    rootGrp.data_variable_name = variableName

    #close the file
    rootGrp.close()

try:
    deleteAllInMemoryRasterLayers()

    # check that each of the source daymet netcdf data files exit
    for netcdfFileName in precpNetCDFFilesArray:
        netcdfFile = os.path.join(sourcePrecpNetCDFilePath, netcdfFileName)
        if(os.path.isfile(netcdfFile) == False):
            raise Exception("Input netcdf file ({0}) was not found.".format(netcdfFile))
            exit()

    if(os.path.isfile(clippedWSDEMFile) == False):
        raise Exception("Input watershed DEM file ({0}) was not found.".format(clippedWSDEMFile))
        exit()

    # check the raster output directory exists
    if(os.path.isdir(outRasterFilePath) == False):
        raise Exception("Raster output directory ({0}) was not found.".format(outRasterFilePath))
        exit()

    # check the netcdf output directory exists
    if(os.path.isdir(outNetCDFFilePath) == False):
        raise Exception("Netcdf output directory ({0}) was not found.".format(outNetCDFFilePath))
        exit()

    # check the netcdf output directory exists
    if(os.path.isdir(sourcePrecpNetCDFilePath) == False):
        raise Exception("Netcdf input directory ({0}) was not found.".format(sourcePrecpNetCDFilePath))
        exit()

    # check that the netcdf input file directory and the output directory are not the same directory
    if(sourcePrecpNetCDFilePath == outNetCDFFilePath):
        raise Exception("Netcdf input directory and output directory found to be the same directory.")
        exit()

    # delete any pre-existing output data raster files
    deleteAllFilesFromDirectory(outRasterFilePath)

    # delete any pre-existing output netcdf files
    deleteAllFilesFromDirectory(outNetCDFFilePath)

    # find the data year from one of the input Daymet netcdf files
    # open the first netCDF file in the array
    inputNetCDFDataFile = os.path.join(sourcePrecpNetCDFilePath, precpNetCDFFilesArray[0] )
    dataYear = getInputDataFileStartYear(inputNetCDFDataFile)

    if(dataYear == None):
        raise Exception("Data year was not found in the specified Daymet input data files.")
        exit()

    WSdesc = arcpy.Describe(clippedWSDEMFile)
    wsCS = ""
    wsCS = WSdesc.spatialReference.name
    wsCS = wsCS.replace("_", " ")
    if(not wsCS):
        # set output coordinate system if DEM file does not have coordinate system
        outCS = arcpy.SpatialReference('NAD 1983 UTM Zone 12N')
    else:
        outCS = arcpy.SpatialReference(wsCS)

    wsGridMeanCellHeight = WSdesc.MeanCellHeight

    # get the rectangular boundary of the WS dem file for clipping the
    # precp data raster file
    wsBoundingBox = str(WSdesc.extent.XMin) + " " + str(WSdesc.extent.YMin) + " " + str(WSdesc.extent.XMax) + " " + str(WSdesc.extent.YMax)

    # generate raster layer names and store in the above array
    # as well as create rsater layer from each of the input netcdf data files
    dataVariableName = 'prcp'
    dataXdimension = 'x'
    dataYdimension = 'y'
    dataBandDimension = ''
    dataDimensionValues = ''
    dataValueSelectionMethod = 'BY_VALUE'

    for fileNo in range(0, len(precpNetCDFFilesArray)):
        inMemoryRasterLayersArray.append('rasterLayer' + str(fileNo + 1))

        # Process: Make NetCDF Raster Layer (in memory)
        inputNetcdfFile = os.path.join(sourcePrecpNetCDFilePath, precpNetCDFFilesArray[fileNo])
        outRasterLayer = inMemoryRasterLayersArray[fileNo]
        arcpy.MakeNetCDFRasterLayer_md(inputNetcdfFile, dataVariableName, dataXdimension,
                                        dataYdimension, outRasterLayer, dataBandDimension, dataDimensionValues, dataValueSelectionMethod)

    # create the Geoprocessor object
    gp = arcgisscripting.create()

    # check out any necessary licenses
    gp.CheckOutExtension("spatial")

    # snap raster is neeed for projection of the data raster file to the watershed projection
    gp.SnapRaster = clippedWSDEMFile

    startDay = 1
    startMonth = 1
    startDate = datetime.date(dataYear, startMonth, startDay)
    timeValueFormat = "%m/%d/%Y %H:%M:%S %p"
    rasterFileNameDateFormat = "%Y_%m_%d"
    valueSelectionMethod = "BY_VALUE"

    # create the list of input raster layers for mosaic operation
    inputRastersToMosaic = ''
    for fileNo in range(0, len(inMemoryRasterLayersArray)):
        if((fileNo + 1) < len(inMemoryRasterLayersArray)):
                inputRastersToMosaic = inputRastersToMosaic + inMemoryRasterLayersArray[fileNo] + ';'
        else:
            inputRastersToMosaic = inputRastersToMosaic + inMemoryRasterLayersArray[fileNo]


    print('>>>raster start time:' + str(datetime.datetime.now()))

    # select data for each time band from each of the in memory raster layers and then
    # mosiac the raster layer to save to disk for a given day
    start_time = time.time()
    daysInYear = DaymetNumberOfDaysToProcess.getNumberOfDays();
    for day in range(0, daysInYear):
        newDay = startDate + datetime.timedelta(days= day)
        outRasterFileName = 'prcp' + newDay.strftime(rasterFileNameDateFormat) + '.tif'
        outProjRasterFileName = 'wsProjprcp' + newDay.strftime(rasterFileNameDateFormat) + '.tif'
        outClippedRasterFileName = 'wsClippedprcp' + newDay.strftime(rasterFileNameDateFormat) + '.tif'

        valueSelect = ["time", newDay.strftime(timeValueFormat)]
        print('Day:' + str(day + 1))
        dimensionValues = [valueSelect]

        # Execute SelectByDimension to select the current time band from the in memory raster layers
        for fileNo in range(0, len(inMemoryRasterLayersArray)):
            arcpy.SelectByDimension_md(inMemoryRasterLayersArray[fileNo], dimensionValues, valueSelectionMethod)

        # mosaic all the raster layers for the current time selection
        #ref: http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=mosaic_to_new_raster_(data_management)
        numberOfBands = '1'
        coordinate_system_for_the_output_raster = '' # this setting will use the cordinate system for the input rasters for the output raster
        pixel_type_of_output_raster = '32_BIT_FLOAT'
        cell_size_of_output_raster = '' # this setting will use the cell size of the input rasters for the ouput raster
        mosiac_method = 'MEAN' # the output cell value of the overlapping areas will be the mean of the overlapped cells
        mosaic_colormap_mode = 'FIRST' # this setting will use the colormap of the first input raster
        arcpy.MosaicToNewRaster_management(inputRastersToMosaic, outRasterFilePath, outRasterFileName,
                                            coordinate_system_for_the_output_raster, pixel_type_of_output_raster,
                                            cell_size_of_output_raster, numberOfBands, mosiac_method, mosaic_colormap_mode)

        # project the saved raster to projection of the watershed  dem
        inRasterFile = os.path.join(outRasterFilePath, outRasterFileName)
        outProjRasterFile = os.path.join(outRasterFilePath, outProjRasterFileName)
        resampling_type = "BILINEAR" # Bilinear interpolation

        lockGP = threading.Lock()
        lockGP.acquire()

        gp.ProjectRaster_management(inRasterFile, outProjRasterFile, outCS, resampling_type, wsGridMeanCellHeight)

        # clip the raster data file to the extent of the buffered watershed
        # ref:http://webhelp.esri.com/arcgisdesktop/9.3/index.cfm?TopicName=Clip_(Data_Management)
        outClippedRasterFile = os.path.join(outRasterFilePath, outClippedRasterFileName)
        gp.clip_management(outProjRasterFile, wsBoundingBox, outClippedRasterFile, clipping_geometry ="ClippingGeometry")

        lockGP.release()

    end_time = time.time()
    raster_elapsed_time = end_time - start_time
    print('>>>>Elapsed time in creating raster files:' + str(raster_elapsed_time))

    # delete all in memory raster layers
    deleteAllInMemoryRasterLayers()

    start_time = time.time()
    # merge all the watershed data raster files using the gdal_merge function
    outMergedRasterFileName = 'mergedPrcp.tif'
    outMergedRasterFile = os.path.join(outRasterFilePath, outMergedRasterFileName)
    inRasterFileNamePatternToMatch = 'wsClippedprcp' + str(dataYear) + '_??_??.tif'

    mergeRasters(outRasterFilePath, outMergedRasterFile, inRasterFileNamePatternToMatch)

    print('>>>>Elapsed time in creating raster files:' + str(raster_elapsed_time))

    end_time = time.time()
    elapsed_time = end_time - start_time
    print('>>>>Elapsed time in merging the raster files:' + str(elapsed_time))

    start_time = time.time()
    # create the final netcdf data file for the watershed
    # Process: RasterToNetCDF
    # parameters for creating netcdf file from the raster data file
    variable = "Prec"
    units = "" #actual unit for precp data is 'mm/day' but this unit does not work with netcdf create function - correct unit is set after the file is created
    XDimension = "x"
    YDimension = "y"
    bandDimension = "time"

    lockGP = threading.Lock()
    lockGP.acquire()

    outNetCDFFile = os.path.join(outNetCDFFilePath, outNetCDFFileName)
    gp.RasterToNetCDF_md(outMergedRasterFile, outNetCDFFile, variable, units,
                            XDimension, YDimension, bandDimension)
    lockGP.release()

    del gp
    end_time = time.time()
    elapsed_time = end_time - start_time
    print('>>>>Elapsed time in creating the final netcdf file:' + str(elapsed_time))

    # set the units of the Precp variable in the output netcdf file
    setNetCDFDDataUnits(outNetCDFFile, variable, 'mm/day')

    # set the start_year global attribute to dataYear
    setNetCDFDataYearAttribute(outNetCDFFile, dataYear)

    # set the data_varaible_name global attribute to as set in variable
    setNetCDFDataVariableNameAttribute(outNetCDFFile, variable)

    # delete all files from the outRasterFilePath folder
    deleteAllFilesFromDirectory(outRasterFilePath)

    print('done...')

except:
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pyErrMsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
    print(pyErrMsg)
    print('>>>done...with exception')
    raise Exception(pyErrMsg)
