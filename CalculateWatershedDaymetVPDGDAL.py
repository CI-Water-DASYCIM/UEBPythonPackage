# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# CalculateWatershedDayemtVPDGDAL.py
# Created on: 2013-01-28 14:54:03.00000
# Author: Pabitra Dash (pabitra.dash@usu.edu)
#
# Description:
#   Generates a netcdf vapor pressure deficit (VPD)data file for the
#   domain watershed from a number of
#   netcdf data files obtanied from Daymet website

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
# set desktop license used to ArcView
# ref: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//002z0000000z000000
import arcview

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

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# Local variables:
sourceVpNetCDFilePath = None
outNetCDFFilePath = None
outNetCDFFileName = None
outRasterFilePath = None
clippedWSDEMRasterFile = None
inVpDataFileNames = None
gp = None
inSimulationStartDate = None
inSimulationEndDate = None

# find the dir path of this python script location
thisScriptPath = os.path.dirname(sys.argv[0])

# append the found path to the python search path
sys.path.append(thisScriptPath)

# check if the script file to read the number of days to simulate exists
##scriptFileToReadNumberOfDaysToSimulate = os.path.join(thisScriptPath,'DaymetNumberOfDaysToProcess.py')
##if(os.path.isfile(scriptFileToReadNumberOfDaysToSimulate) == False):
##    raise Exception("Script file ({0}) was not found.".format(scriptFileToReadNumberOfDaysToSimulate))
##    exit()


##import DaymetNumberOfDaysToProcess

# settings for runnning this script locally not part of the workflow.
# to run this code locally NOT as part of the workflow uncomment the following 12 lines
##thisScriptFullFilePath = os.path.join(thisScriptPath,'CalculateWatershedDaymetVPGDAL.py')
##argumentList = []
##argumentList.append(thisScriptFullFilePath) #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\vpdatasets')
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\vpdatasets\OutNetCDF')
##argumentList.append('vp_daily_one_data.nc')
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\vpdatasets\Raster')
##argumentList.append(r"E:\CIWaterData\Temp\ws_dem.tif")
##argumentList.append('vp1.nc;vp2.nc')
##argumentList.append('2011/01/01') # simulation start date
##argumentList.append('2011/01/05') # simulation end date
##sys.argv = argumentList


# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 9 more argument, so total of 9
if (len(sys.argv) < 9):
    print('Invalid arguments:')
    print('1st argument: Input vp netcdf file path')
    print('2nd argument: Output vp netcdf file path')
    print('3rd argument: Output vp netcdf file name')
    print('4th argument: Output vp temporary raster file path')
    print('5th argument: Input watershed DEM raster file name with file path')
    print('6th argument: List of input vp netcdf file names')
    print('7th argument: Simulation start date')
    print('8th argument: Simulation end date')
    raise Exception("There has to be 8 arguments to calculate Dayment vapor pressure at each grid of the watershed.")
    exit()

# retrieve passed argument values
sourceVpNetCDFilePath = sys.argv[1]
outNetCDFFilePath = sys.argv[2]
outNetCDFFileName = sys.argv[3]
outRasterFilePath = sys.argv[4]
clippedWSDEMRasterFile = sys.argv[5]
inVpDataFileNames = sys.argv[6]
inSimulationStartDate = sys.argv[7]
inSimulationEndDate = sys.argv[8]

# array to store all input Daymet netcdf file names
vpNetCDFFilesArray = []
vpNetCDFFilesArray = inVpDataFileNames.split(';')

if(len(vpNetCDFFilesArray) == 0):
    raise Exception("Netcdf input Daymet data file names are missing.")
    exit()

# validate the start and end simulation dates by converting to date types
try:
    startDate = datetime.datetime.strptime(inSimulationStartDate, '%Y/%m/%d')
    endDate = datetime.datetime.strptime(inSimulationEndDate, '%Y/%m/%d')
except ValueError as ex:
    raise Exception(str(ex))
    exit()

if(startDate >= endDate):
    raise Exception("Simulation end date must be date after the simulation start date")
    exit()

# array to store raster layers created from each of the input netcdf prep data file
inMemoryRasterLayersArray = []

def deleteAllFilesFromDirectory(directoryPath):
    os.chdir(directoryPath)
    rasterFiles = glob.glob('*.*')
    for filename in rasterFiles:
        os.unlink(filename)

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

def convertRasterToNetCDF(inputRasterFile, outNetCDFFile):
    variable = "Vp"
    units = "" #actual unit for RH data is 'kPa' but this unit does not work with netcdf create function - correct unit is set after the file is created
    XDimension = "x"
    YDimension = "y"
    bandDimension = "time"
    # Create the Geoprocessor object
    gp = arcgisscripting.create()

    # Check out any necessary licenses
    gp.CheckOutExtension("spatial")

    gp.RasterToNetCDF_md(inputRasterFile, outNetCDFFile, variable, units,
                            XDimension, YDimension, bandDimension)

    # check in any necessary licenses
    gp.CheckInExtension("spatial")

def setNetCDFDataUnits(vpNetCDF_File, dataSetVariableName, variableUnit):
    #open the file with read/write mode (r+)
    rootGrp = Dataset(vpNetCDF_File, 'r+', format='NETCDF3_CLASSIC')

    #get the vp variable
    precpVariable = rootGrp.variables[dataSetVariableName]

    #set the units of vp variable to killopascal
    precpVariable.units = variableUnit

    #close the file
    rootGrp.close()

def setNetCDFDataYearAttribute(netCDF_File, dataYear):
    #open the file with read/write mode (r+)
    rootGrp = Dataset(netCDF_File, 'r+', format='NETCDF3_CLASSIC')

    rootGrp.start_year = dataYear

    #close the file
    rootGrp.close()

def setNetCDFDataVariableNameAttribute(netCDF_File, variableName):
    #open the file with read/write mode (r+)
    rootGrp = Dataset(netCDF_File, 'r+', format='NETCDF3_CLASSIC')

    rootGrp.data_variable_name = variableName

    #close the file
    rootGrp.close()

try:
    deleteAllInMemoryRasterLayers()

    # check that each of the source netcdf data files exit
    for netcdfFileName in vpNetCDFFilesArray:
        netcdfFile = os.path.join(sourceVpNetCDFilePath, netcdfFileName)
        if(os.path.isfile(netcdfFile) == False):
            raise Exception("Input netcdf file ({0}) was not found.".format(netcdfFile))
            exit()

    if(os.path.isfile(clippedWSDEMRasterFile) == False):
        raise Exception("Input watershed DEM file ({0}) was not found.".format(clippedWSDEMRasterFile))
        exit()

    # check the raster output directory exists
    if(os.path.isdir(outRasterFilePath) == False):
        raise Exception("Raster output directory ({0}) was not found.".format(outRasterFilePath))
        exit()

    # check the netcdf output directory exists
    if(os.path.isdir(outNetCDFFilePath) == False):
        raise Exception("Netcdf output directory ({0}) outNetCDFFilePath was not found.".format(outNetCDFFilePath))
        exit()

    # check the netcdf output directory exists
    if(os.path.isdir(sourceVpNetCDFilePath) == False):
        raise Exception("Netcdf input directory ({0}) was not found.".format(sourceVpNetCDFilePath))
        exit()

    # check that the netcdf input file directory and the output directory are not the same directory
    if(sourceVpNetCDFilePath == outNetCDFFilePath):
        raise Exception("Netcdf input directory and output directory found to be the same directory.")
        exit()

    # delete any pre-existing output data raster files
    deleteAllFilesFromDirectory(outRasterFilePath)

    # delete any pre-existing output netcdf files
    deleteAllFilesFromDirectory(outNetCDFFilePath)

    # Find the data year from one of the input Daymet netcdf files
    # Open the first netCDF file in the array
    inputNetCDFDataFile = os.path.join(sourceVpNetCDFilePath, vpNetCDFFilesArray[0] )
    dataYear = getInputDataFileStartYear(inputNetCDFDataFile)

    if(dataYear == None):
        raise Exception("Data year was not found in the specified Daymet input data files.")
        exit()

    if (dataYear != startDate.year):
        raise Exception("Specified simulation start year ({0}) does not match with Daymet data year.".format(startDate.year))
        exit()

    if (dataYear != endDate.year):
        raise Exception("Specified simulation end year ({0}) does not match with Daymet data year.".format(endDate.year))
        exit()

    # grab the start time
    start_time = time.clock()

    WSdesc = arcpy.Describe(clippedWSDEMRasterFile)
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
    # vp data raster file
    wsBoundingBox = str(repr(WSdesc.extent.XMin)) + " " + str(repr(WSdesc.extent.YMin)) + " " + str(repr(WSdesc.extent.XMax)) + " " + str(repr(WSdesc.extent.YMax))

    # generate raster layer names and store in the above array
    # as well as create rsater layer from each of the input netcdf data files
    dataVariableName = 'vp'
    dataXdimension = 'x'
    dataYdimension = 'y'
    dataBandDimension = ''
    dataDimensionValues = ''
    dataValueSelectionMethod = 'BY_VALUE'

    for fileNo in range(0, len(vpNetCDFFilesArray)):
        inMemoryRasterLayersArray.append('rasterLayer' + str(fileNo + 1))

        # Process: Make NetCDF Raster Layer (in memory)
        inputNetcdfFile = os.path.join(sourceVpNetCDFilePath, vpNetCDFFilesArray[fileNo])
        outRasterLayer = inMemoryRasterLayersArray[fileNo]
        arcpy.MakeNetCDFRasterLayer_md(inputNetcdfFile, dataVariableName, dataXdimension, dataYdimension,
                                        outRasterLayer, dataBandDimension, dataDimensionValues, dataValueSelectionMethod)

    # Create the Geoprocessor object
    gp = arcgisscripting.create()

    # Check out any necessary licenses
    gp.CheckOutExtension("spatial")

    # snap raster is neeed for projection of the data raster file to the watershed projection
    gp.SnapRaster = clippedWSDEMRasterFile

##    startDay = 1
##    startMonth = 1
##    startDate = datetime.date(dataYear, startMonth, startDay)
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
    # get the number of days to simulate from simulation start and end dates
    #daysInYear = DaymetNumberOfDaysToProcess.getNumberOfDays();
    numberOfSimulationDays = (endDate - startDate).days + 1
    for day in range(0, numberOfSimulationDays):
        newDay = startDate + datetime.timedelta(days= day)
        outRasterFileName = 'vp' + newDay.strftime(rasterFileNameDateFormat) + '.tif'
        outProjRasterFileName = 'wsProjVp' + newDay.strftime(rasterFileNameDateFormat) + '.tif'
        outClippedRasterFileName = 'wsClippedVp' + newDay.strftime(rasterFileNameDateFormat) + '.tif'

        valueSelect = ["time", newDay.strftime(timeValueFormat)]
        print('Day:' + str(day + 1))
        dimensionValues = [valueSelect]

        # execute SelectByDimension to select the current time band from the in memory raster layers
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

    end_time = time.clock()
    raster_elapsed_time = end_time - start_time


    print('>>>>Time took to generate rasters:' + str(raster_elapsed_time) + ' seconds')

    # delete all in memory raster layers
    deleteAllInMemoryRasterLayers()

    start_time = time.clock()

    # merge all the watershed data raster files using the gdal_merge function
    outMergedRasterFileName = 'mergedVp.tif'
    outMergedRasterFile = os.path.join(outRasterFilePath, outMergedRasterFileName)
    inRasterFileNamePatternToMatch = 'wsClippedVp' + str(dataYear) + '_??_??.tif'

    mergeRasters(outRasterFilePath, outMergedRasterFile, inRasterFileNamePatternToMatch)

    end_time = time.clock()
    elapsed_time = end_time - start_time

    print('>>>>Time took to generate rasters:' + str(raster_elapsed_time) + ' seconds')
    print('>>>>Time took to merge rasters:' + str(elapsed_time) + ' seconds')


    start_time = time.clock()
    # create the final netcdf data file for the watershed
    # Process: RasterToNetCDF
    # parameters for creating netcdf file from the raster data file
    variable = "vp"
    units = "" #actual unit for vp data is 'kPa' but this unit does not work with netcdf create function - correct unit is set after the file is created
    XDimension = "x"
    YDimension = "y"
    bandDimension = "time"

    lockGP = threading.Lock()
    lockGP.acquire()

    outNetCDFFile = os.path.join(outNetCDFFilePath, outNetCDFFileName)
    gp.RasterToNetCDF_md(outMergedRasterFile, outNetCDFFile, variable, units,
                            XDimension, YDimension, bandDimension)
    lockGP.release()
    print('>>>>netcdf end time:' + str(datetime.datetime.now()))

    # set the units of the vp variable in the output netcdf file
    setNetCDFDataUnits(outNetCDFFile, variable, 'Pa')

    #set the start_year global attribute to dataYear
    setNetCDFDataYearAttribute(outNetCDFFile, dataYear)

    #set the data_varaible_name global attribute to as set in variable
    setNetCDFDataVariableNameAttribute(outNetCDFFile, variable)

    end_time = time.clock()
    elapsed_time = end_time - start_time

    print('>>>>Time took to generate netcdf file:' + str(elapsed_time) + ' seconds')

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
finally:
    # check in any necessary licenses
    arcpy.CheckInExtension("spatial")
    # check in any necessary licenses
    if(gp != None):
        gp.CheckInExtension("spatial")
        del gp