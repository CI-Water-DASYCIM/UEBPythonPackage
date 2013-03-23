#-------------------------------------------------------------------------------
# Name:         DaymetNumberOfDaysToProcess.py
# Purpose:      Reads a number from a text file
#               as number of days to process daymet data. Scripts related to Dayment
#               data then can call the getNumberOfDays() in this file to generate data
#               for any number of days by changing the number of days in the text file
#
# Author:      Pabitra
#
# Created:     21/03/2013
# Copyright:   (c) Pabitra 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os

def main():
    pass

if __name__ == '__main__':
    main()

# test reading a text file

def getNumberOfDays():
    fileToOpenOnLocal = r"E:\SoftwareProjects\CIWaterPythonScripts\NumberOfDaysToProcessDaymentData.txt"
    fileToOpenOnRemote = r"C:\CIWaterPythonScripts\NumberOfDaysToProcessDaymentData.txt"
    fileToOpen = None

    if(os.path.isfile(fileToOpenOnRemote) == True):
        fileToOpen = fileToOpenOnRemote
    elif(os.path.isfile(fileToOpenOnLocal) == True):
        fileToOpen = fileToOpenOnLocal

##    print(fileToOpen)

    numberOfDays = 365
    with open(fileToOpen) as file:
            numberOfDays = file.readline()

    try:
        numberOfDays = int(numberOfDays)
        if(numberOfDays < 0 or numberOfDays > 365):
            numberOfDays = 365
    except ValueError:
        numberOfDays = 365

    return numberOfDays
##    print("Done...")

