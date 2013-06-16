#-------------------------------------------------------------------------------
# Name:     CalculateWatershedLatLonValues.py
# Author:   Pabitra Dash (pabitra.dash@usu.edu)

# Purpose:
#   To calculate lat/lon values for the watershed domain. For smaller
#   watershed only the lat/lon value for the domain midpoint will be
#   alculated and written to a text file. For larger watershed, lat/lon
#   values for each of the grid cell will be calculated and then be written to
#   seprarately two different NetCDF files (one for lat values and the other for
#   lon values. the genrated files will be saved in the same floder where
#   watershed dem file exists
#
#
# Created:     17/01/2013
# Copyright:   (c) Pabitra 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# this says all code at indent level 0 is part of the main() function
def main():
    pass

# this says run  all code at indent level 0 when this script file is ran as a standalone script
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

#REF: http://gis.stackexchange.com/questions/20783/how-to-get-x-y-coordinates-and-cell-value-of-each-pixel-in-a-raster-using-python

# local variables
wsDEMRasterFile = None
latRasterFileName = "WSDEMLat.tif" #temporary file
lonRasterFileName = "WSDEMLon.tif" #temporary file
latFileName = None
lonFileName= None
latRasterFile = None
lonRasterFile = None
netCDFLatFile = None
netCDFLonFile = None

# settings for runnning this code locally not part of the workflow. To run this code on remote app server as part of the workflow
# comment out the following 6 lines
# to run locally not part of a workflow, uncomment the following 6 lines
##argumentList = []
##argumentList.append('') #this argument is reserved for the name of this script file
##argumentList.append(r'E:\CIWaterData\Temp\ws_dem.tif')
##argumentList.append('ws_lat')
##argumentList.append('ws_lon')
##sys.argv = argumentList

# the first argument sys.argv[0] is the name of this script file
# excluding the name of the script file we need 3 more arguments, so total of 4
if (len(sys.argv) < 4):
    print('Invalid arguments:')
    print('1st argument: Input watershed DEM raster file name with path')
    print('2nd argument: Output watershed lat values file name without file extension.')
    print('3rd argument: Output watershed lon values file name without file extension.')
    raise Exception("There has to be 3 argument to generate lat/lon values for the watershed.")
    exit()

# retrieve the passed argument
wsDEMRasterFile = sys.argv[1]
latFileName = sys.argv[2]
lonFileName = sys.argv[3]

# check if provided ws DEM file exists
if(os.path.isfile(wsDEMRasterFile) == True):
    filePath = os.path.dirname(wsDEMRasterFile)

    # set the path for the temporary lat raster file
    latRasterFile = os.path.join(filePath, latRasterFileName)

    # set the path for the temporary lon raster file
    lonRasterFile = os.path.join(filePath, lonRasterFileName)

    # set the path for the lat netcdf file
    netCDFLatFileName = latFileName + '.nc'
    netCDFLatFile = os.path.join(filePath, netCDFLatFileName)

    # set the path for the lon netcdf file
    netCDFLonFileName = lonFileName + '.nc'
    netCDFLonFile = os.path.join(filePath, netCDFLonFileName)
else:
    print('Exception')
    raise Exception("Specified watershed DEM file ({0}) was not found.:".format(wsDEMRasterFile))
    exit()

#Get properties of the input raster
inRasterDesc = arcpy.Describe(wsDEMRasterFile)

try:
    # coordinates of the lower left corner
    rasXmin = inRasterDesc.Extent.XMin
    rasYmin = inRasterDesc.Extent.YMin
    rasXMax= inRasterDesc.Extent.XMax
    rasYMax = inRasterDesc.Extent.YMax

    # Cell size, raster size
    rasMeanCellHeight = inRasterDesc.MeanCellHeight
    rasMeanCellWidth = inRasterDesc.MeanCellWidth
    rasHeight = inRasterDesc.Height
    rasWidth = inRasterDesc.Width

    # lat at the center of the raster
    rasCenterX = rasXmin + ((0.5 + rasWidth) * rasMeanCellWidth)

    # lon at the center of the raster
    rasCenterY = rasYmin + ((0.5 + rasHeight) * rasMeanCellHeight)

    isCreateGriddedLatFile = True
    isCreateGriddedLonFile= True

    # calculate east-west extent of the wsdemfile (XMax-XMin)
    # If (XMax-XMin) > 0.5 then create the lat netcdf file
    # otherwise calculate the mid point lon value for the dem file
    # and write to a text file Lon.txt
    if(rasYMax - rasYmin < 0.5):
        filePath = os.path.dirname(wsDEMRasterFile)
        outputFilePath = os.path.join(filePath,  lonFileName + '.txt')
        textFileWriter = None
        try:
            textFileWriter = open(outputFilePath, "w")
            textFileWriter.write(str(rasCenterY))
        finally:
            textFileWriter.close()

        isCreateGriddedLonFile = False

    # calculate north-soutth extent of the wsdemfile (YMax-YMin)
    # If (YMax-YMin) > 0.5 then create the lat netcdf file
    # otherwise calculate the mid point lat value for the dem file
    # and write to a text file Lat.txt
    if(rasYMax - rasYmin < 0.5):
        filePath = os.path.dirname(wsDEMRasterFile)
        outputFilePath = os.path.join(filePath, latFileName + '.txt')
        textFileWriter = None
        try:
            textFileWriter = open(outputFilePath, "w")
##            textFileWriter.write("Lat: %s\n"%rasCenterX)
            textFileWriter.write(str(rasCenterX))
        finally:
            textFileWriter.close()

        isCreateGriddedLatFile = False

    # if there exists a previously lat raster file delete it
    if(os.path.isfile(latRasterFile) == True):
        os.unlink(latRasterFile)

    # if there exists a previously lon raster file delete it
    if(os.path.isfile(lonRasterFile) == True):
        os.unlink(lonRasterFile)

    # if there exists a previously lat netcdf file delete it
    if(os.path.isfile(netCDFLatFile) == True):
        os.unlink(netCDFLatFile)

    # if there exists a previously lon netcdf file delete it
    if(os.path.isfile(netCDFLonFile) == True):
        os.unlink(netCDFLonFile)

    # Define projection system for the latRasterFile
    wsCS = ""
    wsCS = inRasterDesc.spatialReference.name
    wsCS = wsCS.replace("_", " ")

    if(not wsCS):
        # set output coordinate system if DEM file does not have coordinate system
        outCS = arcpy.SpatialReference('NAD 1983 UTM Zone 12N')
    else:
        outCS = arcpy.SpatialReference(wsCS)

    if(isCreateGriddedLatFile == True):
        # create numpy array of coordinates of cell centroids
        def rasCentrX(rasHeight, rasWidth):
            coordX = rasXmin + ((0.5 + rasWidth) * rasMeanCellWidth)
            return coordX

        inRasterCoordX = numpy.fromfunction(rasCentrX, (rasHeight,rasWidth)) #numpy array of X coord

        # write lat values array to a raster file
        point = arcpy.Point(rasXmin, rasYmin)
        latRaster = arcpy.NumPyArrayToRaster(inRasterCoordX, lower_left_corner= point, x_cell_size=rasMeanCellWidth, y_cell_size=rasMeanCellHeight)
        latRaster.save(latRasterFile)

        # define coordinate system for the latRasterFile same as the watershed DEM file
        arcpy.DefineProjection_management(latRasterFile, outCS)

        # write the lat raster file to nedCDF file
        variable = "latitude"
        units = "degree"
        XDimension = "x"
        YDimension = "y"
        bandDimension = ""

        arcpy.RasterToNetCDF_md(latRasterFile, netCDFLatFile, variable, units, XDimension, YDimension, bandDimension)

    if(isCreateGriddedLonFile == True):
        def rasCentrY(rasHeight, rasWidth):
            coordY = rasYmin + ((0.5 + rasHeight) * rasMeanCellHeight)
            return coordY

        inRasterCoordY = numpy.fromfunction(rasCentrY, (rasHeight,rasWidth)) #numpy array of Y coord

        # write lon values array to a raster file
        point = arcpy.Point(rasXmin, rasYmin)
        lonRaster = arcpy.NumPyArrayToRaster(inRasterCoordY, lower_left_corner= point, x_cell_size=rasMeanCellWidth, y_cell_size=rasMeanCellHeight)
        lonRaster.save(lonRasterFile)

        # define coordinate system for the LonRasterFile same as the watershed DEM file
        arcpy.DefineProjection_management(lonRasterFile, outCS)

        # write the lon raster file to nedCDF file
        variable = "longitude"
        units = "degree"
        XDimension = "x"
        YDimension = "y"
        bandDimension = ""

        arcpy.RasterToNetCDF_md(lonRasterFile, netCDFLonFile, variable, units, XDimension, YDimension, bandDimension)

    print ("Done...")

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