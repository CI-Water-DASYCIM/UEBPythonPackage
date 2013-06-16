# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# CalculateWatershedAtmosphericPressure.py
# Author: Pabitra Dash (pabitra.dash@usu.edu)

# Description:
#   Calculates the atmospheric pressure for a watershed based on
#   provided watershed DEM file and writes a text file and saves it to the same directory
#   where the ws DEM file exists
# ---------------------------------------------------------------------------

# ref: http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?id=1663&pid=1644&topicname=Get%20Raster%20Properties%20(Data%20Management)&

# set desktop license used to ArcView
# ref: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//002z0000000z000000
import arcview

import arcgisscripting
import sys
import os
import traceback

# Local variables:
WSDEMFile = None # this is the file generated after regriding the DEM file to match with buffered watershed
OutputWsAtomPresTextFileName = None
gp = None

# settings for runnning this code locally not part of the workflow. To run this code on remote app server as part of the workflow
# comment out the following 5 lines
# to run locally not part of a workflow, uncomment the following 5 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\Temp\ws_dem.tif')
##argumentList.append('ws_atom_pres.txt')
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 2 more argument, so total of 3
if (len(sys.argv) < 3):
    print('Invalid arguments:')
    print('1st argument: Input watershed DEM raster file name with file path')
    print('2nd argument: Output watershed atmospherice pressure text file name')
    raise Exception("There has to be 2 argument to calculate atmospheric pressure for the watershed.")
    exit()

# retrieve the passed argument
WSDEMFile = sys.argv[1]
OutputWsAtomPresTextFileName = sys.argv[2]

# check if provided DEM file exists
if(os.path.isfile(WSDEMFile) == False):
    raise Exception("Exception: Specified watershed DEM file was not found:" + WSDEMFile)
    exit()

gp = arcgisscripting.create()

# Check out any necessary licenses
gp.CheckOutExtension("spatial")

# find average elevation of watershed DEM
try:
    avgElevation = gp.GetRasterProperties (WSDEMFile, "MEAN")
    print(avgElevation)

    # calculate autmospheric pressure
    Po = 101325.0 # in Pascal unit for standard atmospheric pressure at sea level
    To = 288.15 # in Kelvin unit for standard temperature at sea level
    watershedAtmoPressure = eval("Po * (1- 0.0065 * avgElevation / To) ** 5.257")

    filePath = os.path.dirname(WSDEMFile)
    outputFilePath = os.path.join(filePath, OutputWsAtomPresTextFileName)
    textFileWriter = None
    try:
        textFileWriter = open(outputFilePath, "w")
        textFileWriter.write(str(watershedAtmoPressure))
    finally:
        textFileWriter.close()

    print ('Done..')

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