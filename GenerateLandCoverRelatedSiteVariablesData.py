#-------------------------------------------------------------------------------
# Name:     GenerateLandCoverRelatedSiteVariablesData.py
# Author:   Pabitra Dash (pabitra.dash@usu.edu)
# Purpose:
#   To generate netcdf file for each of the UEB land cover related
#   site variables (cc, hcan, lay, ycage) for a given watershed DEM file. The
#   generated netcdf files are stored in the same directory where the
#   waterhsed NLCD dataset file exists
#
# Created:      22/01/2013
# Copyright:    (c) 2013
# Licence:      <your licence>
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

import arcpy
import numpy
import os
import sys
import traceback

# check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# local variables
wsNLCDRasterFile = None
ccValueRasterFileName = "WSCCValues.img"        #temporary file
hcanValueRasterFileName = "WSHCANValues.img"    #temporary file
laiValueRasterFileName = "WSLAIValues.img"      #temporary file
ycageValueRasterFileName = "WSYCAGEValues.img"  #temporary file
netCDFCCValueFileName = None
netCDFHCANValueFileName = None
netCDFLAIValueFileName = None
netCDFYCAGEValueFileName = None

# settings for runnning this code locally not part of the workflow. To run this code on remote app server as part of the workflow
# comment out the following 8 lines
# to run locally not part of a workflow, uncomment the following 8 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\Temp\ws_nlcd.img')
##argumentList.append('ws_cc_nlcd.nc')
##argumentList.append('ws_hc_nlcd.nc')
##argumentList.append('ws_lai_nlcd.nc')
##argumentList.append('ws_ycage_nlcd.nc')
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 5 more arguments, so total of 6
if (len(sys.argv) < 6):
    print('Invalid arguments:')
    print('1st argument: Input watershed NLCD dataset raster file name with file path.')
    print('2nd argument: Output watershed netcd file name for canopy cover data.')
    print('3rd argument: Output watershed netcd file name for height of canopy data.')
    print('4th argument: Output watershed netcd file name for LAI data.')
    print('5th argument: Output watershed netcd file name for ycage data.')
    raise Exception("Exception:There has to be 5 argument to generate landcover related site variables data files.")
    exit()

# retrieve passed arguments
wsNLCDRasterFile = sys.argv[1]
netCDFCCValueFileName = sys.argv[2]
netCDFHCANValueFileName = sys.argv[3]
netCDFLAIValueFileName = sys.argv[4]
netCDFYCAGEValueFileName =  sys.argv[5]

ccValueRasterFile = None
hcanValueRasterFile = None
laiValueRasterFile = None
ycageValueRasterFile = None
netCDFCCValueFile = None
netCDFHCANValueFile = None
netCDFLAIValueFile = None
netCDFYCAGEValueFile = None

# ref: http://www.mrlc.gov/nlcd06_leg.php
# create the land classification lookup dictionary for CC values
# the 1st column (key) is the classification code values and the 2nd column
# represents cc values
landClassificationCCValueLookup = {
        '11':   '0',
        '12':   '0',
        '21':   '0',
        '22':   '0',
        '23':   '0',
        '24':   '0',
        '31':   '0',
        '41':   '0.5',
        '42':   '0.7',
        '43':   '0.8',
        '51':   '0',
        '52':   '0.5',
        '71':   '0',
        '72':   '0',
        '73':   '0',
        '74':   '0',
        '81':   '0',
        '82':   '0',
        '90':   '0.5',
        '95':   '0'
    }

# create the land classification lookup dictionary for hcan values
# the 1st column (key) is the classification code values and the 2nd column
# represents hcan values
landClassificationHCANValueLookup = {
        '11':   '0',
        '12':   '0',
        '21':   '0',
        '22':   '0',
        '23':   '0',
        '24':   '0',
        '31':   '0',
        '41':   '8',
        '42':   '15',
        '43':   '10',
        '51':   '0',
        '52':   '3',
        '71':   '0',
        '72':   '0',
        '73':   '0',
        '74':   '0',
        '81':   '0',
        '82':   '0',
        '90':   '8',
        '95':   '0'
    }

# create the land classification lookup dictionary for lai values
# the 1st column (key) is the classification code values and the 2nd column
# represents lai values
landClassificationLAIValueLookup = {
        '11':   '0',
        '12':   '0',
        '21':   '0',
        '22':   '0',
        '23':   '0',
        '24':   '0',
        '31':   '0',
        '41':   '1',
        '42':   '4.5',
        '43':   '4',
        '51':   '0',
        '52':   '1',
        '71':   '0',
        '72':   '0',
        '73':   '0',
        '74':   '0',
        '81':   '0',
        '82':   '0',
        '90':   '1',
        '95':   '0'
    }

# create the land classification lookup dictionary for lai values
# the 1st column (key) is the classification code values and the 2nd column
# represents ycage values
landClassificationYCAGEValueLookup = {
        '11':   '2',
        '12':   '2',
        '21':   '2',
        '22':   '2',
        '23':   '2',
        '24':   '2',
        '31':   '2',
        '41':   '2',
        '42':   '3',
        '43':   '2',
        '51':   '2',
        '52':   '2',
        '71':   '2',
        '72':   '2',
        '73':   '2',
        '74':   '2',
        '81':   '2',
        '82':   '2',
        '90':   '2',
        '95':   '2'
    }

def generateSiteVariableData(siteVariableDataArray, dataRasterFileToCreate, netCDFDataFileToCreate, variableName):
        # write cc values array to a raster file
        siteVariableDataRaster = arcpy.NumPyArrayToRaster(siteVariableDataArray, lower_left_corner= point, x_cell_size=rasMeanCellWidth, y_cell_size=rasMeanCellHeight)
        siteVariableDataRaster.save(dataRasterFileToCreate)

        # define coordinate system for the ccValueRasterFile same as the watershed NLCD raster file
        arcpy.DefineProjection_management(dataRasterFileToCreate, outCS)

        # write the landocover cc site variable raster file to nedCDF file
        if(variableName == 'hcan'):
            units = "meter"
        else:
            units = ""  # no units

        variable = variableName
        XDimension = "x"
        YDimension = "y"
        bandDimension = ""

        arcpy.RasterToNetCDF_md(dataRasterFileToCreate, netCDFDataFileToCreate, variable, units, XDimension, YDimension, bandDimension)
try:
    # check if provided watershed NLCD dataset raster file exists
    if(os.path.isfile(wsNLCDRasterFile) == True):
        filePath = os.path.dirname(wsNLCDRasterFile)

        # set the path for the temporary ws cc raster file
        ccValueRasterFile = os.path.join(filePath, ccValueRasterFileName)

        # set the path for the temporary ws cc raster file
        hcanValueRasterFile = os.path.join(filePath, hcanValueRasterFileName)

        # set the path for the temporary ws lai raster file
        laiValueRasterFile = os.path.join(filePath, laiValueRasterFileName)

        # set the path for the temporary ws ycage raster file
        ycageValueRasterFile = os.path.join(filePath, ycageValueRasterFileName)

        # set the path for the output cc value netcdf file
        netCDFCCValueFile = os.path.join(filePath, netCDFCCValueFileName)

        # set the path for the output hcan value netcdf file
        netCDFHCANValueFile = os.path.join(filePath, netCDFHCANValueFileName)

        # set the path for the output lai value netcdf file
        netCDFLAIValueFile = os.path.join(filePath, netCDFLAIValueFileName)

        # set the path for the output ycage value netcdf file
        netCDFYCAGEValueFile = os.path.join(filePath, netCDFYCAGEValueFileName)
    else:
        errMsg = "Specified watershed NLCD raster file ({0}) was not found.".format(wsNLCDRasterFile)
        print(errMsg)
        sys.exit(errMsg)

    # if there exists a previously cc netcdf file delete it
    if(os.path.isfile(netCDFCCValueFile) == True):
        os.unlink(netCDFCCValueFile)

    # if there exists a previously hcan netcdf file delete it
    if(os.path.isfile(netCDFHCANValueFile) == True):
        os.unlink(netCDFHCANValueFile)

    # if there exists a previously lai netcdf file delete it
    if(os.path.isfile(netCDFLAIValueFile) == True):
        os.unlink(netCDFLAIValueFile)

     # if there exists a previously ycage netcdf file delete it
    if(os.path.isfile(netCDFYCAGEValueFile) == True):
        os.unlink(netCDFYCAGEValueFile)

    # if there exists a previously cc raster file delete it
    if(os.path.isfile(ccValueRasterFile) == True):
        os.unlink(ccValueRasterFile)

    # if there exists a previously hcan raster file delete it
    if(os.path.isfile(hcanValueRasterFile) == True):
        os.unlink(hcanValueRasterFile)

    # if there exists a previously lai raster file delete it
    if(os.path.isfile(laiValueRasterFile) == True):
        os.unlink(laiValueRasterFile)

     # if there exists a previously ycage raster file delete it
    if(os.path.isfile(ycageValueRasterFile) == True):
        os.unlink(ycageValueRasterFile)

    rstArray = arcpy.RasterToNumPyArray(wsNLCDRasterFile) # Change rasterFile to numpy array

    rasRows, rasCols = rstArray.shape                     # Return the rows, columns

    # create an empty array to store cc values
    wsCCValueArray = numpy.zeros(shape=(rasRows,rasCols))

    # create an empty array to store hcan values
    wsHCANValueArray = numpy.zeros(shape=(rasRows,rasCols))

    # create an empty array to store lai values
    wsLAIValueArray = numpy.zeros(shape=(rasRows,rasCols))

    # create an empty array to store ycage values
    wsYCAGEValueArray = numpy.zeros(shape=(rasRows,rasCols))

    # read each classification code from the ws nlcd raster array
    # and fill the wsCCValueArray with corresponding cc values
    for rowNum in xrange(rasRows):                     # Loop through the rows
        for colNum in xrange(rasCols):                 # Loop through the row's columns
            classificationCode = rstArray.item(rowNum, colNum)
            print(classificationCode)

            if(classificationCode == -128): # no data value
                classificationCode = 31

            ccValue = landClassificationCCValueLookup[str(classificationCode)]
            hcanValue = landClassificationHCANValueLookup[str(classificationCode)]
            laiValue = landClassificationLAIValueLookup[str(classificationCode)]
            ycageValue = landClassificationYCAGEValueLookup[str(classificationCode)]
            wsCCValueArray[rowNum, colNum] = ccValue
            wsHCANValueArray[rowNum, colNum] = hcanValue
            wsLAIValueArray[rowNum, colNum] = laiValue
            wsYCAGEValueArray[rowNum, colNum] = ycageValue

    # get properties of the input watershed nlcd raster
    wsNLCDRasterDesc = arcpy.Describe(wsNLCDRasterFile)

    # coordinates of the lower left corner of watershed nlcd raster file
    rasXmin = wsNLCDRasterDesc.Extent.XMin
    rasYmin = wsNLCDRasterDesc.Extent.YMin

    # cell size, raster size of the watershed nlcd raster file
    rasMeanCellHeight = wsNLCDRasterDesc.MeanCellHeight
    rasMeanCellWidth = wsNLCDRasterDesc.MeanCellWidth

    # find the coordinate system used in the ws NLCD raster file
    wsCS = ""
    wsCS = wsNLCDRasterDesc.spatialReference.name
    wsCS = wsCS.replace("_", " ")

    if(not wsCS):
        # set output coordinate system if DEM file does not have coordinate system
        outCS = arcpy.SpatialReference('NAD 1983 UTM Zone 12N')
    else:
        outCS = arcpy.SpatialReference(wsCS)

    point = arcpy.Point(rasXmin, rasYmin)

    # generate cc  netcdf data file
    generateSiteVariableData(wsCCValueArray, ccValueRasterFile, netCDFCCValueFile, 'cc')

    # generate hcan netcdf data file
    generateSiteVariableData(wsHCANValueArray, hcanValueRasterFile, netCDFHCANValueFile, 'hcan')

    # generate lai netcdf data  file
    generateSiteVariableData(wsLAIValueArray, laiValueRasterFile, netCDFLAIValueFile, 'lai')

    # generate ycage netcdf data  file
    generateSiteVariableData(wsYCAGEValueArray, ycageValueRasterFile, netCDFYCAGEValueFile, 'ycage')

    print('Done..')
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