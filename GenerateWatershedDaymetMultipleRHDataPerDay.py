#-------------------------------------------------------------------------------
# Name:     GenerateWatershedDaymetMultipleRHDataPerDay.py
# Author:   Pabitra Dash (pabitra.dash@usu.edu)
#
# Purpose:
#   Generates multiple rh data points per day and writes to a netcdf file using
#   a single data point per day from a input netcdf file
#
# Created:     25/02/2013
# Copyright:   (c) 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from netCDF4 import Dataset
import time
import datetime
import numpy
import shutil
import os
import sys
import traceback
from netCDF4 import num2date, date2num
import math

# this says all code at indent level 0 is part of the main() function
def main():
    pass

# this says run  all code at indent level 0 when this script file is ran as standalone script
# meaning not imported in another script
if __name__ == '__main__':
    main()


# Local variables:
sourceTaNetCDFile = None
sourceVpdNetCDFFile = None
outRH_NetCDFFile = None
destNetCDF_FilePath = None
inTimeStep = None
outNetCDFDataVariableName = None
outRootGrp = None
inRootGrpTa = None
inRootGrpVp = None

# settings for runnning this code locally not part of the workflow. To run this code on remote app server as part of the workflow
# comment out the following 9 lines
# to run locally not part of a workflow, uncomment the following 9 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\Temp\ta_daily_multiple_data.nc')
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\vpdatasets\OutNetCDF\vp_daily_multiple_data.nc')
##argumentList.append(r'E:\CIWaterData\DaymetTimeSeriesData\Logan\vpdatasets\OutNetCDF\rh_daily_multiple_data.nc')
##argumentList.append(r'E:\CIWaterData\Temp')
##argumentList.append('rh')
##argumentList.append(6)
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 6 more argument, so total of 7
if (len(sys.argv) < 7):
    print('Invalid arguments:')
    print('1st argument: Input watershed specific daymet temperature netcdf file name with file path')
    print('2nd argument: Input watershed specific daymet vapor pressure deficit netcdf file name with file path')
    print('3rd argument: Output watershed specific daymet RH netcdf data file name with temporary file path')
    print('4th argument: Destination path for the output RH netcdf file')
    print('5th argument: Name of the relative humidity data variable for the output netcdf file')
    print('6th argument: Time step value (in hours) allowed values area: 1, 2, 3, 4,, 6')
    raise Exception("There has to be 6 arguments to calculate multiple RH data points per day.")
    exit()

# retrieve the passed arguments
sourceTaNetCDFile = sys.argv[1]
sourceVpdNetCDFFile = sys.argv[2]
outRH_NetCDFFile = sys.argv[3]
destNetCDF_FilePath = sys.argv[4]
outNetCDFDataVariableName = sys.argv[5]
inTimeStep = int(sys.argv[6])

try:

    # check the netcdf vpd data file exists
    if(os.path.isfile(sourceVpdNetCDFFile) == False):
        sys.exit("Watershed vapor presssure netcdf input file ({0}) was not found.".format(sourceVpdNetCDFFile))

    # check the netcdf temp data file exists
    if(os.path.isfile(sourceTaNetCDFile) == False):
        sys.exit("Watershed temperature netcdf input file ({0}) was not found.".format(sourceTaNetCDFile))

    # check if the output netcdf file temporary directory exists
    filePath = os.path.dirname(outRH_NetCDFFile)
    if(os.path.isdir(filePath) == False):
        sys.exit("Netcdf output temporary directory ({0}) was not found.".format(filePath))

    # check the netcdf final output directory exists
    if(os.path.isdir(destNetCDF_FilePath) == False):
        sys.exit("Destination Netcdf output directory ({0}) was not found.".format(destNetCDF_FilePath))

    # validate input time step value
    if(inTimeStep != 1 and inTimeStep != 2 and inTimeStep != 3 and inTimeStep != 4 and inTimeStep != 6):
        errMsg = "Provided time step value ("  + inTimeStep +  ") is not a valid time step value.\n"
        errMsg += "Valid values are: 1, 2, 3, 4, and 6."
        sys.exit(errMsg)

    # check that the netcdf input file directory and the output directory are not the same directory
    inTaFilePath = os.path.dirname(sourceTaNetCDFile)
    inVpFilePath = os.path.dirname(sourceVpdNetCDFFile)

    # open the netCDF file that has temp data. Open in readonly mode based on which we will be creating a new netcdf file
    inRootGrpTa= Dataset(sourceTaNetCDFile, 'r', format='NETCDF3_CLASSIC')

    # open the netCDF file that has vp data. Open in readonly mode based on which we will be creating a new netcdf file
    inRootGrpVp = Dataset(sourceVpdNetCDFFile, 'r', format='NETCDF3_CLASSIC')

    # check the input time step for the output rh netcdf file is same as the
    # time step used in the input vp and temp files
    if(inTimeStep != inRootGrpTa.data_time_step or inTimeStep != inRootGrpVp.data_time_step):
        errMsg = "Provided time step value ("  + inTimeStep +  ") for the ouput RH netcdf is not a valid time step value.\n"
        errMsg += "It must match with the time step used in the temperature and vapor pressure files."
        sys.exit(errMsg)

    inputTaVar = inRootGrpTa.variables[inRootGrpTa.data_variable_name]
    inputTaXVar = inRootGrpTa.variables['x']
    inputTaYVar = inRootGrpTa.variables['y']
    inputTaTimeVar = inRootGrpTa.variables['time']

    inputVpVar = inRootGrpVp.variables[inRootGrpVp.data_variable_name]
    inputVpTimeVar = inRootGrpVp.variables['time']

    #DEBUG: print netcdf variable diemensions
    print('Dimension of Temp var: ' + str(inputTaVar.shape))
    print('Dimension of Vp var: ' + str(inputVpVar.shape))

    #open a new blank netcdf file to which we will be writting data
    outRootGrp = Dataset(outRH_NetCDFFile, 'w', format='NETCDF3_CLASSIC')

    # add global file level attributes to the new netcdf file
    outRootGrp.start_year = inRootGrpTa.start_year
    outRootGrp.data_variable_name = outNetCDFDataVariableName
    outRootGrp.data_time_step = inRootGrpTa.data_time_step
    outRootGrp.orginal_data_source = 'Daymet Software Version 2.0'
    outRootGrp.conventions = 'CF-1.0'
    outRootGrp.modified_data_source = 'CI Water System'
    outRootGrp.spatial_reference = 'NAD83_UTM_Zone_12N'
    outRootGrp.datum = 'D_North_America_1983'

    # create 3 dimensions for the output netcdf file
    outTimeDimensionSize = inputVpVar.shape[0]
    outYvarDimensionSize = inputVpVar.shape[1]
    outXvarDimensionSize = inputVpVar.shape[2]

    outRootGrp.createDimension('time', outTimeDimensionSize)
    outRootGrp.createDimension('x', outXvarDimensionSize)
    outRootGrp.createDimension('y', outYvarDimensionSize)

    #DEBUG: print each dimension name, dimension length
    for dimName, dimObj in outRootGrp.dimensions.iteritems():
        print dimName, len(dimObj)


    # create a RH variable of data type f4 (32-bit floating point) that has data in all three dimensions
    vRH= outRootGrp.createVariable(outRootGrp.data_variable_name, 'f4',('time', 'y', 'x'))
    # create a variable for each dimension to hold data for that specific dimension
    vTime = outRootGrp.createVariable('time', 'f8', ('time'))   #f8 64-bit floating point
    vX = outRootGrp.createVariable('x', 'f4', ('x'))            #f4: 32-bit floating point
    vY = outRootGrp.createVariable('y', 'f4', ('y'))

    #DEBUG: print demensions of output variable rh and input variable vp
    print(vRH.shape)
    print(inputVpVar.shape)

    # add attributes to time variable
    vTime.units = 'hours since 0001-01-01 00:00:00.0'
    vTime.calendar = 'gregorian'
    vTime.long_name = 'Band'

    # add attributes to Prec variable
    vRH.long_name = outRootGrp.data_variable_name
    vRH.esri_pe_string = inputTaVar.esri_pe_string
    vRH.coordinates  = inputTaVar.coordinates
    vRH.grid_mapping  = inputTaVar.grid_mapping
    vRH. missing_value  = inputTaVar.missing_value
    vRH.units  = ''

    # add attributes to x variable
    vX.long_name = inputTaXVar.long_name
    vX.standard_name = inputTaXVar.standard_name
    vX.units = inputTaXVar.units

    # add attributes to y variable
    vY.long_name = inputTaYVar.long_name
    vY.standard_name = inputTaYVar.standard_name
    vY.units = inputTaYVar.units

    # assign data to output x variable same as the input x data
    vX[:] = inputTaXVar[:]

    # assign data to output y variable same as the input y data
    vY[:] = inputTaYVar[:]

    # assign data to output time variable same as the input time data
    vTime[:]  = inputVpTimeVar
    times, cols, rows = inputVpVar.shape

    start_time = time.clock()

    # assign data to Rh netcdf variable
    for time_step in range(0, times):
        outRhDataArray = numpy.empty((1,cols, rows), dtype=numpy.float32)

        # create an empty 3D array to hold slice of the input temperature data array
        inTaDataSlice = numpy.empty((1, cols, rows), dtype=numpy.float32)
        inTaDataSlice[:] = inputTaVar[time_step:time_step+1, 0:cols, 0:rows]

        # create an empty 3D array to hold slice of the input vpd data array
        inVpDataSlice = numpy.empty((1, cols, rows), dtype=numpy.float32)
        inVpDataSlice[:] = inputVpVar[time_step:time_step+1, 0:cols, 0:rows]

        for row in range(0, rows):
            for col in range(0, cols):
                ta = inTaDataSlice[0][col][row]
                vpd = inVpDataSlice[0][col][row]

                # calculate saturated vapor pressure (units Pa)
                sVp = 611 * math.exp((17.3 * ta) /(237.3 + ta))
                rh = 1 - (vpd / sVp)

                outRhDataArray[0][col][row] = rh

        # write the RH data to the output netcdf file
        vRH[time_step:time_step+1, 0:cols, 0:rows] = outRhDataArray[0:1, 0:cols, 0:rows]

    # closing of the output netcdf file is necessary here before we can move this file to the destination folder
    outRootGrp.close()
    outRootGrp = None

    end_time = time.clock()
    elapsed_time = end_time - start_time
    print('Time taken for the script to finish: ' + str(elapsed_time) + ' seconds')

    # if the output netcdf file already exists at the destination folder, delete it before we can move the file there
    outNetCDF_FileName = os.path.basename(outRH_NetCDFFile)
    destNetCDFFile = os.path.join(destNetCDF_FilePath, outNetCDF_FileName)
    if(os.path.isfile(destNetCDFFile) == True):
        os.unlink(destNetCDFFile)

    # move the generated netcdf file to the destination folder
    shutil.move(outRH_NetCDFFile, destNetCDF_FilePath)
    print("Done...")

except:
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pyErrMsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
    print(pyErrMsg)
    print('>>>done...with exception')
    raise Exception(pyErrMsg)

finally:
    #close the netcdf files
    if(inRootGrpTa != None):
        inRootGrpTa.close()
    if(inRootGrpVp != None):
        inRootGrpVp.close()
    if(outRootGrp != None):
        outRootGrp.close()
