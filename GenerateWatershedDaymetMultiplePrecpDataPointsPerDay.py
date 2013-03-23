#-------------------------------------------------------------------------------
# Name:         GenerateMultiplePrecpDataPointsPerDay.py
# Purpose:      Generates multiple prcp data points per day from a single data point per day
#
# Author:      Pabitra
#
# Created:     21/02/2013
# Copyright:   (c) Pabitra 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# this says all code at indent level 0 is part of the main() function
def main():
    pass

# this says run  all code at indent level 0 when this script file is ran as standalone script
# meaning not imported in another script
if __name__ == '__main__':
    main()

from netCDF4 import Dataset
import time
import datetime
import numpy
import shutil
import os
import sys
import traceback
from netCDF4 import num2date, date2num

#local variables
inNetCDF_File = None
outNetCDF_File = None
inTimeStep = None
destNetCDF_FilePath = None
outRootGrp = None
inRootGrp = None

# settings for runnning this code locally. To run this code on remote app server comment out the following 5 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\precdatasets\OutNetCDF\precp_daily_one_data.nc')
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\precdatasets\OutNetCDF\precp_daily_multiple_data.nc')
##argumentList.append(r'E:\CIWaterData\Temp')
##argumentList.append(6)
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 4 more argument, so total of 5
if (len(sys.argv) < 4):
    print('Invalid arguments:')
    print('1st argument: Input watershed specific daymet precp netcdf file name with directory path')
    print('2nd argument: Output watershed specific daymet precp netcdf file name with temporary directory path')
    print('3rd argument: Destination path for the output precp netcdf file')
    print('4th argument: Time step value (in hours) allowed values area: 1, 2, 3, 4,, 6')
    raise Exception("There has to be 4 arguments to calculate multiple precipitation data points per day.")
    exit()


# retrieve the passed arguments
inNetCDF_File = sys.argv[1]
outNetCDF_File = sys.argv[2]
destNetCDF_FilePath = sys.argv[3]
inTimeStep = int(sys.argv[4])


try:

    #check if the input netcdf file exists
    if(os.path.isfile(inNetCDF_File) == False):
        raise Exception("Input netcdf file " + inNetCDF_File +  " was not found.")
        exit()

    #check if the output netcdf file temporary directory exists
    filePath = os.path.dirname(outNetCDF_File)
    if(os.path.isdir(filePath) == False):
        raise Exception("Netcdf output temporary directory " + filePath +  " was not found.")
        exit()

    #check if the output netcdf file destination directoryt exists
    if(os.path.isdir(destNetCDF_FilePath) == False):
        raise Exception("Netcdf output destination directory " + destNetCDF_FilePath +  " was not found.")
        exit()

     # validate input time step value
    if(inTimeStep != 1 and inTimeStep != 2 and inTimeStep != 3 and inTimeStep != 4 and inTimeStep != 6):
        errMsg = "Provided time step value ("  + inTimeStep +  ") is not a valid time step value.\n"
        errMsg += "Valid values are: 1, 2, 3, 4, and 6."
        raise Exception(errMsg)
        exit()

    #open the netCDF file in readonly mode based on which we will be creating a new netcdf file
    inRootGrp = Dataset(inNetCDF_File, 'r', format='NETCDF4')

    #open a new blank netcdf file to which we will be writting data
    outRootGrp = Dataset(outNetCDF_File, 'w', format='NETCDF4')

    # DEBUG: print global attributes of the original file:
##    for gAttribute in inRootGrp.ncattrs():
##        print 'Global attribute name:', gAttribute, '=', getattr(inRootGrp,gAttribute)

    # add global file level attributes to the new netcdf file
    outRootGrp.start_year = inRootGrp.start_year #2010 # TODO: avoid this hard codeed year value by using 'inRootGrp.start_year' once  we run the new code in CalculateWatershedDayemetPrecpGDAL.py
    outRootGrp.data_variable_name = 'Prec'
    outRootGrp.data_time_step = inTimeStep
    outRootGrp.orginal_data_source = 'Daymet Software Version 2.0'
    outRootGrp.conventions = 'CF-1.0'
    outRootGrp.modified_data_source = 'CI Water System'
    outRootGrp.spatial_reference = 'NAD83_UTM_Zone_12N'
    outRootGrp.datum = 'D_North_America_1983'

    #  DEBUG: print global attributes of the new file:
##    for gAttribute in outRootGrp.ncattrs():
##        print 'Global attribute name:', gAttribute, '=', getattr(outRootGrp,gAttribute)


    # get dimension values from the original netcdf file
    inputTimeVar = inRootGrp.variables['time']
    inputXvar = inRootGrp.variables['x']
    inputYvar = inRootGrp.variables['y']
    inputPrecVar = inRootGrp.variables['Prec']

    print(inputTimeVar.shape[0])
    print(inputXvar.shape[0])
    print(inputYvar.shape[0])
    print(inputPrecVar.shape)

    #DEBUG:
##    print(inputXvar[:])

##    for vAttribute in inputPrecVar.ncattrs():
##        print 'Prec attribute', vAttribute, '=', getattr(inputPrecVar, vAttribute)
##
##    for vAttribute in inputTimeVar.ncattrs():
##        print 'time attribute', vAttribute, '=', getattr(inputTimeVar, vAttribute)
##
##    for vAttribute in inputXvar.ncattrs():
##        print 'x attribute', vAttribute, '=', getattr(inputXvar, vAttribute)
##
##    for vAttribute in inputYvar.ncattrs():
##        print 'y attribute', vAttribute, '=', getattr(inputYvar, vAttribute)

    # Create 3 dimensions for the 3 variables: time, x and y
    dataPointsPerDayNeeded = 24/inTimeStep
    outTimeDimensionSize = inputTimeVar.shape[0] * dataPointsPerDayNeeded
    outXvarDimensionSize = inputXvar.shape[0]
    outYvarDimensionSize = inputYvar.shape[0]

    outRootGrp.createDimension('time', outTimeDimensionSize)
    outRootGrp.createDimension('x', outXvarDimensionSize)
    outRootGrp.createDimension('y', outYvarDimensionSize)

    #DEBUG: print each dimension name, dimension length
##    for dimName, dimObj in outRootGrp.dimensions.iteritems():
##        print dimName, len(dimObj)


    # create a Prec variable of data type f4 that has data in all three dimensions
    vPrec= outRootGrp.createVariable('Prec', 'f4',('time', 'y', 'x'))
    # create a variable for each dimension to hold data for that specific dimension
    vTime = outRootGrp.createVariable('time', 'f8', ('time'))   #f8 is 64-bit floating point
    vX = outRootGrp.createVariable('x', 'f4', ('x'))            #f4 is 32-bit floating point
    vY = outRootGrp.createVariable('y', 'f4', ('y'))

    print(vPrec.shape)
    print(inputPrecVar.shape)

    #DEBUG:
##    for varKVP in outRootGrp.variables.iteritems():
##        print varKVP[0]


    # add attributes to time variable
    vTime.units = 'hours since 0001-01-01 00:00:00.0'
    vTime.calendar = 'gregorian'
    vTime.long_name = 'Band'

    # add attributes to Prec variable
    vPrec.long_name = inputPrecVar.long_name
    vPrec.esri_pe_string = inputPrecVar.esri_pe_string
    vPrec.coordinates  = inputPrecVar.coordinates
    vPrec.grid_mapping  = inputPrecVar.grid_mapping
    vPrec. missing_value  = inputPrecVar. missing_value
    vPrec.units  = 'mm'

    # add attributes to x variable
    vX.long_name = inputXvar.long_name
    vX.standard_name = inputXvar.standard_name
    vX.units = inputXvar.units

    # add attributes to y variable
    vY.long_name = inputYvar.long_name
    vY.standard_name = inputYvar.standard_name
    vY.units = inputYvar.units

    # DEBUG:
##    for vAttribute in vPrec.ncattrs():
##        print 'Prec attribute', vAttribute, '=', getattr(vPrec, vAttribute)
##
##    for vAttribute in vTime.ncattrs():
##        print 'time attribute', vAttribute, '=', getattr(vTime, vAttribute)
##
##    for vAttribute in vX.ncattrs():
##        print 'x attribute', vAttribute, '=', getattr(vX, vAttribute)
##
##    for vAttribute in vY.ncattrs():
##        print 'y attribute', vAttribute, '=', getattr(vY, vAttribute)

    # assign data to x variable
    vX[:] = inputXvar[:]

    # DEBUG:
##    print(vX[:])

    # assign data to y variable
    vY[:] = inputYvar[:]

    # DEBUG:
##    print(vY[:])

    # assign time data to time variable
    dataYear = outRootGrp.start_year
    startDate = datetime.date(dataYear, 1, 1)
    dates =[datetime.datetime(dataYear,1,1,0,0,0) + n*datetime.timedelta(hours=inTimeStep) for n in range(vTime.shape[0])]
    vTime[:] = date2num(dates, units=vTime.units, calendar=vTime.calendar)
##    print(vTime[:])

    times, cols, rows = inputPrecVar.shape
    startTimeIndex = 0
    outPrecValue = 0.0

    start_time = time.clock()

    for time_step in range(0, times):
        startTimeIndex = time_step * dataPointsPerDayNeeded
        outTempDataArray = numpy.empty((dataPointsPerDayNeeded,cols, rows), dtype=numpy.float32)
        inDataSlice = numpy.empty((dataPointsPerDayNeeded, cols, rows), dtype=numpy.float32)
        inDataSlice[:] = inputPrecVar[time_step:time_step+1, 0:cols, 0:rows]
        for row in range(0, rows):
            for col in range(0, cols):
                outPrecValue = inDataSlice[0][col][row]/dataPointsPerDayNeeded
                for dataPoint in range(0, dataPointsPerDayNeeded):
                    outTempDataArray[dataPoint][col][row] = outPrecValue

        # write the Prec data to the output netcdf file
        vPrec[startTimeIndex:startTimeIndex + dataPointsPerDayNeeded, 0:cols, 0:rows] = outTempDataArray[0:dataPointsPerDayNeeded, 0:cols, 0:rows]

    end_time = time.clock()
    elapsed_time = end_time - start_time
    print('Time taken for the script to finish: ' + str(elapsed_time) + ' seconds')

    #close the netcdf files
    outRootGrp.close()
    inRootGrp.close()

    # if the output netcdf file already exists at the destination folder, delete it before we can move the file there
    outNetCDF_FileName = os.path.basename(outNetCDF_File)
    destNetCDFFile = os.path.join(destNetCDF_FilePath, outNetCDF_FileName)
    if(os.path.isfile(destNetCDFFile) == True):
        os.unlink(destNetCDFFile)
    # move the generated netcdf file to the destination folder
    shutil.move(outNetCDF_File, destNetCDF_FilePath)
    print('>>>done...')

except:
    if(outRootGrp != None):
        outRootGrp.close()
    if(inRootGrp != None):
        inRootGrp.close()
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pyErrMsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
    print(pyErrMsg)
    print('>>>done...with exception')
    raise Exception(pyErrMsg)