#-------------------------------------------------------------------------------
# Name:     GenerateWatershedDaymetMultipleWindDataPerDay.py
# Author:   Pabitra Dash (pabitra.dash@usu.edu)
# Purpose:
#   To generate a gridded netcdf file containing constant wind speed data based on the
#   data for grid varaiable x, y and time from the prep netcdf file
#
# Created:     21/02/2013
# Copyright:   (c) 2013
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

#local variables
inPrecpNetCDF_File = None
outWindNetCDF_File = None
destNetCDF_FilePath = None
inWindSpeed = None
outRootGrp = None
inRootGrp = None

# settings for runnning this code locally. To run this code on remote app server comment out the following 7 lines
# To run locally, uncomment the following 7 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\Temp\precp_daily_multiple_data.nc')
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\winddatasets\OutNetCDF\wind_daily_multiple_data.nc')
##argumentList.append(r'E:\CIWaterData\Temp')
##argumentList.append(2)
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 4 more argument, so total of 5
if (len(sys.argv) < 5):
    print('Invalid arguments:')
    print('1st argument: Input watershed specific daymet precp netcdf file name with directory path')
    print('2nd argument: Output watershed specific daymet wind speed netcdf file name with directory path')
    print('3rd argument: Destination path for the output precp netcdf file')
    print('4th argument: Constant wind speed (in m/sec) that will be written to the ouput netcdf file')
    raise Exception("There has to be 4 arguments to calculate multiple wind speed data points per day.")
    exit()


inPrecpNetCDF_File = sys.argv[1]
outWindNetCDF_File = sys.argv[2]
destNetCDF_FilePath = sys.argv[3]
inWindSpeed = int(sys.argv[4])

try:

    # check if the input netcdf file exists
    if(os.path.isfile(inPrecpNetCDF_File) == False):
        raise Exception("Input netcdf file ({0}) was not found.".format(inPrecpNetCDF_File))
        exit()

    # check if the output netcdf file directoryt exists
    filePath = os.path.dirname(outWindNetCDF_File)
    if(os.path.isdir(filePath) == False):
        raise Exception("Netcdf output directory ({0}) was not found.".format(filePath))
        exit()

    # check if the output netcdf file destination directoryt exists
    if(os.path.isdir(destNetCDF_FilePath) == False):
        raise Exception("Netcdf output destination directory ({0}) was not found.".format(destNetCDF_FilePath))
        exit()

     # validate input wind speed value
    if(inWindSpeed < 0 and inWindSpeed > 20):
        errMsg = "Provided wind speed value ("  + inWindSpeed +  ") is not a valid value.\n"
        errMsg += "Valid values are in the range of 0 to 20."
        raise Exception(errMsg)
        exit()

    # open the netCDF file in readonly mode based on which we will be creating a new netcdf file
    inRootGrp = Dataset(inPrecpNetCDF_File, 'r', format='NETCDF3_CLASSIC')

    # open a new blank netcdf file to which we will be writting data
    outRootGrp = Dataset(outWindNetCDF_File, 'w', format='NETCDF3_CLASSIC')

    # DEBUG: print global attributes of the original file:
##    for gAttribute in inRootGrp.ncattrs():
##        print 'Global attribute name:', gAttribute, '=', getattr(inRootGrp,gAttribute)

    # add global file level attributes to the new netcdf file
    outRootGrp.start_year = inRootGrp.start_year
    outRootGrp.data_variable_name = 'V'
    outRootGrp.data_time_step = inRootGrp.data_time_step
    outRootGrp.orginal_data_source = 'CI Water Generated Data Source'
    outRootGrp.conventions = 'CF-1.0'
    outRootGrp.modified_data_source = 'CI Water System'
    outRootGrp.spatial_reference = 'NAD83_UTM_Zone_12N'
    outRootGrp.datum = 'D_North_America_1983'

    # get dimension values from the original netcdf file
    inputTimeVar = inRootGrp.variables['time']
    inputXvar = inRootGrp.variables['x']
    inputYvar = inRootGrp.variables['y']
    inputPrecVar = inRootGrp.variables['Prec']

    print(inputTimeVar.shape[0])
    print(inputXvar.shape[0])
    print(inputYvar.shape[0])
    print(inputPrecVar.shape)

    # Create 3 dimensions for the 3 variables: time, x and y based on the dimension of the
    # corresponding variables from the precp netcdf file

    outTimeDimensionSize = inputTimeVar.shape[0]
    outXvarDimensionSize = inputXvar.shape[0]
    outYvarDimensionSize = inputYvar.shape[0]

    outRootGrp.createDimension('time', outTimeDimensionSize)
    outRootGrp.createDimension('x', outXvarDimensionSize)
    outRootGrp.createDimension('y', outYvarDimensionSize)

    #DEBUG: print each dimension name, dimension length
##    for dimName, dimObj in outRootGrp.dimensions.iteritems():
##        print dimName, len(dimObj)


    # create a Prec variable of data type f4 (32-bit floating point data type) that has data in all three dimensions
    vWind= outRootGrp.createVariable('V', 'f4',('time', 'y', 'x'))
    # create a variable for each dimension to hold data for that specific dimension
    vTime = outRootGrp.createVariable('time', 'f8', ('time'))   #f8 is 64-bit floating point
    vX = outRootGrp.createVariable('x', 'f4', ('x'))            #f4 is 32-bit floating point
    vY = outRootGrp.createVariable('y', 'f4', ('y'))

    # add attributes to time variable
    vTime.units = 'hours since 0001-01-01 00:00:00.0'
    vTime.calendar = 'gregorian'
    vTime.long_name = 'Band'

    # add attributes to V (wind speed) variable
    vWind.long_name = 'Wind speed'
    vWind.esri_pe_string = inputPrecVar.esri_pe_string
    vWind.coordinates  = inputPrecVar.coordinates
    vWind.grid_mapping  = inputPrecVar.grid_mapping
    vWind. missing_value  = inputPrecVar. missing_value
    vWind.units  = 'm/sec'

    # add attributes to x variable
    vX.long_name = inputXvar.long_name
    vX.standard_name = inputXvar.standard_name
    vX.units = inputXvar.units

    # add attributes to y variable
    vY.long_name = inputYvar.long_name
    vY.standard_name = inputYvar.standard_name
    vY.units = inputYvar.units

    # assign data to output x variable same as the input x variable
    vX[:] = inputXvar[:]

    # assign data to output y variable same as the input y variable
    vY[:] = inputYvar[:]

    # assign data to output time variable same as the input time variable
    vTime[:] = inputTimeVar[:]

    # get the dimension of the prec varaiable array shape
    times, cols, rows = inputPrecVar.shape

    # assign wind data to the output wind variable array by each time step
    # in order to avoid memory error. Assigning data to the full wind varaible array in one
    # step causes memory error
    for time_step in range(0, times):
        shape = (1, cols, rows)
        vWind[time_step: time_step +1, 0:cols, 0:rows] = numpy.ones(shape, dtype=numpy.float32) * inWindSpeed

    print(vWind.shape)

    # close the netcdf files
    outRootGrp.close()
    inRootGrp.close()

    # if the output netcdf file already exists at the destination folder, delete it before we can move the file there
    outNetCDF_FileName = os.path.basename(outWindNetCDF_File)
    destNetCDFFile = os.path.join(destNetCDF_FilePath, outNetCDF_FileName)
    if(os.path.isfile(destNetCDFFile) == True):
        os.unlink(destNetCDFFile)

    # move the generated netcdf file to the destination folder
    shutil.move(outWindNetCDF_File, destNetCDF_FilePath)

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
