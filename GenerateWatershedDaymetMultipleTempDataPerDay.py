#-------------------------------------------------------------------------------
# Name:     GenerateWatershedDaymetMultipleTempDataPerDay.py
# Author:   Pabitra Dash (pabitra.dash@usu.edu)
# Purpose:
#   Generates multiple air temp data points per day and writes to a netcdf file using a single data point per day
#   from an input netcdf file
#
# Created:     25/02/2013
# Copyright:   (c) 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from netCDF4 import Dataset
import time
import datetime
import math
import numpy
import shutil
import os
import sys
import traceback
from netCDF4 import num2date, date2num

# this says all code at indent level 0 is part of the main() function
def main():
    pass

# this says run  all code at indent level 0 when this script file is ran as standalone script
# meaning not imported in another script
if __name__ == '__main__':
    main()


# Local variables:
sourceTaNetCDFilePath = None
outNetCDFFilePath = None
inTminNetCDF_FileName = None
inTmaxNetCDF_FileName = None
outNetCDF_FileName = None
destNetCDF_FilePath = None
outDataVariableName = None
inTminNetCDF_File = None
inTmaxNetCDF_File = None
inTimeStep = None
inRootGrpTmin = None
inRootGrpTmax = None
outRootGrp = None
inSimulationStartDate = None

# settings for runnning this script locally not part of the workflow.
# to run this code locally NOT as part of the workflow uncomment the following 12 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\tempdatasets\OutNetCDF')
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\tempdatasets\OutNetCDF')
##argumentList.append(r'E:\CIWaterData\Temp')
##argumentList.append('tmin_daily_one_data.nc')
##argumentList.append('tmax_daily_one_data.nc')
##argumentList.append('ta_daily_multiple_data.nc')
##argumentList.append('T')
##argumentList.append(6)
##argumentList.append('2011/01/01') # simulation start date
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 9 more argument, so total of 10
if (len(sys.argv) < 10):
    print('Invalid arguments:')
    print('1st argument: Source directory path for watershed specific daymet temperature netcdf files')
    print('2nd argument: Temporary directory path for output netcdf file')
    print('3rd argument: Final directory path for output netcdf file')
    print('4th argument: Name of the netcdf file that has watreshed daily minimum temperature data')
    print('5th argument: Name of the netcdf file that has watreshed daily maximum temperature data')
    print('6th argument: Name of of the output netcdf file' )
    print('7th argument: Name of the temperature data variable for the output netcdf file')
    print('8th argument: Time step value (in hours) allowed values area: 1, 2, 3, 4,, 6')
    print('9th argument: Simulation start date')
    raise Exception("There has to be 9 arguments to calculate multiple temperature data point per day.")
    exit()


# retrieve passed arguments
sourceTaNetCDFilePath = sys.argv[1]
outNetCDFFilePath = sys.argv[2]
destNetCDF_FilePath = sys.argv[3]
inTminNetCDF_FileName = sys.argv[4]
inTmaxNetCDF_FileName = sys.argv[5]
outNetCDF_FileName = sys.argv[6]
outDataVariableName = sys.argv[7]
inTimeStep = int(sys.argv[8])
inSimulationStartDate = sys.argv[9]

# param: hour - hour in 24-hour cycle at which the factor is needed
# param: interval - time interval at which factor is calculated
def getTminFactor(hour, interval):
    piAtGivenTime = math.pi/2 + (math.pi/12) * hour
    piRangeStartValue = piAtGivenTime - (math.pi/12) * interval/2
    piRangeEndVale = piAtGivenTime + (math.pi/12) * interval/2
    tMinFactor = 2 *(math.cos(piRangeStartValue) - math.cos(piRangeEndVale))/math.pi

    # after taking the absolute value of the factor and then rounding it to 3 decimal places if the factor is zero
    if(round(math.fabs(tMinFactor), 3) == 0):
        tMinFactor = 0.5
    elif(tMinFactor < 0):
        tMinFactor = 1 + tMinFactor

    return tMinFactor

# validate the simulation start date by converting to date type
try:
    startDate = datetime.datetime.strptime(inSimulationStartDate, '%Y/%m/%d')
except ValueError as ex:
    raise Exception(str(ex))
    exit()

try:
    # check the netcdf temporary output directory exists
    if(os.path.isdir(outNetCDFFilePath) == False):
        sys.exit("Netcdf output directory ({0}) was not found.".format(outNetCDFFilePath))

    # check the netcdf final output directory exists
    if(os.path.isdir(destNetCDF_FilePath) == False):
        sys.exit("Netcdf final output directory ({0}) was not found.".format(destNetCDF_FilePath))

    # check the netcdf input directory exists
    if(os.path.isdir(sourceTaNetCDFilePath) == False):
        sys.exit("Netcdf input directory ({0}) was not found.".format(sourceTaNetCDFilePath))

    inTminNetCDF_File = os.path.join(sourceTaNetCDFilePath, inTminNetCDF_FileName)
    if(os.path.isfile(inTminNetCDF_File) == False):
        sys.exit("Input netcdf file ({0}) was not found.".format(inTminNetCDF_File))

    inTmaxNetCDF_File = os.path.join(sourceTaNetCDFilePath, inTmaxNetCDF_FileName)
    if(os.path.isfile(inTmaxNetCDF_File) == False):
        sys.exit("Input netcdf file ({0}) was not found.".format(inTmaxNetCDF_File))

    if(inTimeStep != 1 and inTimeStep != 2 and inTimeStep != 3 and inTimeStep != 4 and inTimeStep != 6):
        errMsg = "Provided time step value ("  + inTimeStep +  ") is not a valid time step value.\n"
        errMsg += "Valid values are: 1, 2, 3, 4, and 6."
        sys.exit(errMsg)

    # open the netCDF file that has minimum daily temp data. Open in readonly mode based on which we will be creating a new netcdf file
    inRootGrpTmin = Dataset(inTminNetCDF_File, 'r', format='NETCDF3_CLASSIC')

    # open the netCDF file that has maximum daily temp data. Open in readonly mode based on which we will be creating a new netcdf file
    inRootGrpTmax = Dataset(inTmaxNetCDF_File, 'r', format='NETCDF3_CLASSIC')

    # check that both input temp data files have the same data year
    if (inRootGrpTmin.start_year != inRootGrpTmax.start_year):
        raise Exception("Specified temp data files seem to have different data years.")
        exit()

    # check that the simulation start year is same as the data year of one of the input temp data files
    if (inRootGrpTmin.start_year != startDate.year):
        raise Exception("Specified simulation start year ({0}) does not match with input temp data files data year.".format(startDate.year))
        exit()

    inputTminVar = inRootGrpTmin.variables['tmin']
    inputTmaxVar = inRootGrpTmax.variables['tmax']
    inputXminVar = inRootGrpTmin.variables['x']
    inputYminVar = inRootGrpTmin.variables['y']

    # open a new blank netcdf file to which we will be writting data
    outNetCDF_File = os.path.join(outNetCDFFilePath, outNetCDF_FileName)
    outRootGrp = Dataset(outNetCDF_File, 'w', format='NETCDF3_CLASSIC')

    # add global file level attributes to the new netcdf file
    outRootGrp.start_year = inRootGrpTmin.start_year
    outRootGrp.data_variable_name = outDataVariableName
    outRootGrp.data_time_step = inTimeStep
    outRootGrp.orginal_data_source = 'Daymet Software Version 2.0'
    outRootGrp.conventions = 'CF-1.0'
    outRootGrp.modified_data_source = 'CI Water System'
    outRootGrp.spatial_reference = 'NAD83_UTM_Zone_12N'
    outRootGrp.datum = 'D_North_America_1983'

    print(inputTminVar.shape)

    # Create 3 dimensions for the 3 variables: time, x and y
    dataPointsPerDayNeeded = 24/inTimeStep
    outTimeDimensionSize = inputTminVar.shape[0] * dataPointsPerDayNeeded
    outYvarDimensionSize = inputTminVar.shape[1]
    outXvarDimensionSize = inputTminVar.shape[2]
    outRootGrp.createDimension('time', outTimeDimensionSize)
    outRootGrp.createDimension('x', outXvarDimensionSize)
    outRootGrp.createDimension('y', outYvarDimensionSize)

    #print each dimension name, dimension length
##    for dimName, dimObj in outRootGrp.dimensions.iteritems():
##        print dimName, len(dimObj)

    # create a Temp variable of data type f4 (32-bit floating point data type) that has data in all three dimensions
    vTemp= outRootGrp.createVariable(outDataVariableName, 'f4',('time', 'y', 'x'), least_significant_digit=3)
    # create a variable for each dimension to hold data for that specific dimension
    vTime = outRootGrp.createVariable('time', 'f8', ('time'))   #f8 is 64-bit floating point
    vX = outRootGrp.createVariable('x', 'f4', ('x'))            #f4 is 32-bit floating point
    vY = outRootGrp.createVariable('y', 'f4', ('y'))

    #DEBUG: print dimensions of the variable temp and tmin
    print(vTemp.shape)
    print(inputTminVar.shape)

    # add attributes to time variable
    vTime.units = 'hours since 0001-01-01 00:00:00.0'
    vTime.calendar = 'gregorian'
    vTime.long_name = 'Band'

    # add attributes to Prec variable
    vTemp.long_name = outDataVariableName
    vTemp.esri_pe_string = inputTminVar.esri_pe_string
    vTemp.coordinates  = inputTminVar.coordinates
    vTemp.grid_mapping  = inputTminVar.grid_mapping
    vTemp. missing_value  = inputTminVar. missing_value
    vTemp.units  = 'deg C' #Celcious

    # add attributes to x variable
    vX.long_name = inputXminVar.long_name
    vX.standard_name = inputXminVar.standard_name
    vX.units = inputXminVar.units

    # add attributes to y variable
    vY.long_name = inputYminVar.long_name
    vY.standard_name = inputYminVar.standard_name
    vY.units = inputYminVar.units

    # assign data to x variable
    vX[:] = inputXminVar[:]

    # assign data to y variable
    vY[:] = inputYminVar[:]

    # assign time data to time variable
    #dataYear = outRootGrp.start_year
    #startDate = datetime.date(dataYear, 1, 1)

    #generate input time step based date time values starting with the simulation start date
    #dates =[datetime.datetime(dataYear,1,1,0,0,0) + n*datetime.timedelta(hours=inTimeStep) for n in range(vTime.shape[0])]
    dates =[datetime.datetime(startDate.year, startDate.month, startDate.day, 0, 0, 0) + n*datetime.timedelta(hours=inTimeStep) for n in range(vTime.shape[0])]
    vTime[:] = date2num(dates, units=vTime.units, calendar=vTime.calendar)

    times, cols, rows = inputTminVar.shape
    startTimeIndex = 0
    outTempValue = 0.0
    start_time = time.clock()
    tMinFactorsArray = []
    for dataPoint in range(0, dataPointsPerDayNeeded):
        tMinFactor = getTminFactor(dataPoint*inTimeStep, inTimeStep)
        tMinFactorsArray.append(tMinFactor)

    for time_step in range(0, times):
        startTimeIndex = time_step * dataPointsPerDayNeeded
        outTempDataArray = numpy.empty((dataPointsPerDayNeeded,cols, rows), dtype=numpy.float32)
        #create an empty 3D array to hold slice of the input tmin data array
        inTminDataSlice = numpy.empty((1, cols, rows), dtype=numpy.float32)
        inTminDataSlice[:] = inputTminVar[time_step:time_step+1, 0:cols, 0:rows]

        #create an empty 3D array to hold slice of the input tmax data array
        inTmaxDataSlice = numpy.empty((1, cols, rows), dtype=numpy.float32)
        inTmaxDataSlice[:] = inputTmaxVar[time_step:time_step+1, 0:cols, 0:rows]

        for row in range(0, rows):
            for col in range(0, cols):
                tMin = inTminDataSlice[0][col][row]
                tMax = inTmaxDataSlice[0][col][row]
                for dataPoint in range(0, dataPointsPerDayNeeded):
                    tMinFactor = tMinFactorsArray[dataPoint]
                    tMaxFactor = 1 - tMinFactor
                    outTempValue = tMinFactor * tMin + tMaxFactor * tMax
                    outTempDataArray[dataPoint][col][row] = outTempValue

        # write the Temp data to the output netcdf file
        vTemp[startTimeIndex:startTimeIndex+dataPointsPerDayNeeded, 0:cols, 0:rows] = outTempDataArray[0:dataPointsPerDayNeeded, 0:cols, 0:rows]

    # closing of the output netcdf file is necessary before we can move this to the destination folder
    outRootGrp.close()
    outRootGrp = None

    end_time = time.clock()
    elapsed_time = end_time - start_time
    print(str(elapsed_time) + ' seconds')

    # if the output netcdf file already exists at the destination folder, delete it before we can move the file there
    destNetCDFFile = os.path.join(destNetCDF_FilePath, outNetCDF_FileName)
    if(os.path.isfile(destNetCDFFile) == True):
        os.unlink(destNetCDFFile)

    # move the generated netcdf file to the destination folder
    shutil.move(outNetCDF_File, destNetCDF_FilePath)
    print(">>>Done...")

except:
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pyErrMsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
    print(pyErrMsg)
    print('>>>done...with exception')
    raise Exception(pyErrMsg)

finally:
    #close the netcdf files
    if(inRootGrpTmin != None):
        inRootGrpTmin.close()
    if(inRootGrpTmax != None):
        inRootGrpTmax.close()
    if(outRootGrp != None):
        outRootGrp.close()