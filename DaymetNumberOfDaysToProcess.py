#-------------------------------------------------------------------------------
# Name:     DaymetNumberOfDaysToProcess.py
# Author:   Pabitra Dash (pabitra.dash@usu.edu)
# Purpose:
#   Reads a number from a text file (NumberOfDaysToProcessDaymetData.txt)
#   as number of days to process daymet data. Scripts related to Dayment
#   data then can call the getNumberOfDays() in this file to generate data
#   for any number of days by changing the number of days in the text file
#
# Created:     21/03/2013
# Copyright:   (c) 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os
import sys

def main():
    pass

if __name__ == '__main__':
    main()

def getNumberOfDays():
    # find the dir path of this python script location
    thisScriptPath = os.path.dirname(sys.argv[0])
    fileToOpen = os.path.join(thisScriptPath, 'NumberOfDaysToProcessDaymetData.txt')

    # default number of days to simulate
    numberOfDays = 365

    # if the text file can't be found return 365 days as number of days to simulate
    if(os.path.isfile(fileToOpen) == False):
        return numberOfDays

    with open(fileToOpen) as file:
            numberOfDays = file.readline()

    try:
        numberOfDays = int(numberOfDays)
        if(numberOfDays < 0 or numberOfDays > 365):
            numberOfDays = 365
    except ValueError:
        numberOfDays = 365

    return numberOfDays

