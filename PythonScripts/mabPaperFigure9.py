'''   
Copyright (c) 2010-2013, Delft University of Technology
Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
All rights reserved.
See http://bit.ly/12SHPLR for license details.

This script can be used to generate a sample plot of the time-history of Mab's (action) orbital
elements.
'''

################################################################################################### 
# Set up input deck
###################################################################################################

# Set Monte Carlo run IDs as a list.
monteCarloRunId =

# Set absolute path to directory with data files.
datafilesPath   =

# Set absolute path to SQLite database with simulation data.
databasePath    = 

# Set case name.
caseName        = 

# Set absolute path to output directory.
outputPath      = 

# Set figure dpi.
figureDPI       = 

###################################################################################################

'''
    
                            DO NOT EDIT PARAMETERS BEYOND THIS POINT!!!

'''


################################################################################################### 
# Import Python packages and modules
###################################################################################################

# # Import necessary external packages.
# import math
import matplotlib.pyplot as plt
import numpy
import numpy.lib.recfunctions
import os
import sqlite3
import time

# # Import user-defined modules.
import constants

###################################################################################################


################################################################################################### 
# Start timer
###################################################################################################

startTime = time.time()

###################################################################################################


################################################################################################### 
# Check output directory
###################################################################################################

# Check if output directory exists; if not, create it.
if not os.path.exists(outputPath):
    os.makedirs(outputPath)

###################################################################################################


###################################################################################################
# Execute database operations
###################################################################################################

# Connect to SQLite database.
database = sqlite3.connect(databasePath)

# Enter database and retrieve desired data.
with database:
    # Set cursor to scan through database and execute queries.
    cursor = database.cursor()
    
    # Select the case ID corresponding to the case name provided.
    cursor.execute("SELECT caseId FROM test_particle_case WHERE caseName == \"" + caseName + "\";")
    caseId = cursor.fetchall()[0][0]  
    
    # Select all the case data associated with the case ID.
    cursor.execute("SELECT * FROM test_particle_case WHERE caseId == " + str(caseId) + ";")
    rawCaseData = cursor.fetchall()
    caseDataColumnNameList = [column_name[0] for column_name in cursor.description]  
    
# Close database connection if it is open.    
if database:
    database.close()    
    
# Store case data in dictionary, using column names from database.    
caseData = {} 
for i,name in enumerate(caseDataColumnNameList):
    caseData[name] = rawCaseData[0][i]    

###################################################################################################


###################################################################################################
# Read and store datafiles
###################################################################################################

# Change the working directory to the directory containing the data files.
os.chdir(datafilesPath)

# Declare empty array to store Keplerian element arrays.
keplerianElements = []

# Set up base for filenames.
filenameBase = 'monteCarloRun' + str(monteCarloRunId) + "_"
keplerianElementsFilename = filenameBase + 'keplerianActionElements.csv'

# Read in Keplerian elements history for test particle.
keplerianElements.append(numpy.genfromtxt(keplerianElementsFilename, \
                         delimiter = ',', comments = '#', names = True) )

###################################################################################################


###################################################################################################
# Generate figures
###################################################################################################

# Set output path and case-prefix for files generated.
outputPathAndCasePrefix = outputPath + "/case" + str(caseData['caseId']) + "_"

plt.figure()
plt.xlabel("Epoch [Julian years]")
plt.ylabel("Mab semi-major axis wrt T0 [km]")
plt.xlim(xmin = 0.0, xmax = 50.0)
for keplerData in keplerianElements:
    plt.plot(keplerData['epoch']/constants.julianYear, \
             (keplerData['semiMajorAxis'] - caseData['perturbedBodySemiMajorAxisAtT0']) \
             * constants.meterInKilometers, 'k')
plt.savefig(outputPathAndCasePrefix + "semiMajorAxisHistory.pdf", \
            dpi = figureDPI)    
plt.close()

plt.figure()
plt.xlabel("Epoch [Julian years]")
plt.ylabel("Mab eccentricity wrt T0 [km]")
plt.xlim(xmin = 0.0, xmax = 50.0)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
for keplerData in keplerianElements:
    plt.plot(keplerData['epoch']/constants.julianYear, \
             keplerData['eccentricity'] - caseData['perturbedBodyEccentricityAtT0'], 'k')
plt.savefig(outputPathAndCasePrefix + "eccentricityHistory.pdf", \
            dpi = figureDPI)    
plt.close()


plt.figure()
plt.xlabel("Epoch [Julian years]")
plt.ylabel("Mab inclination [km]")
plt.xlim(xmin = 0.0, xmax = 50.0)
for keplerData in keplerianElements:
    plt.plot(keplerData['epoch']/constants.julianYear,keplerData['inclination'], 'k')
plt.savefig(outputPathAndCasePrefix + "inclinationHistory.pdf", \
            dpi = figureDPI)    
plt.close()


###################################################################################################


###################################################################################################
# Finalize timer and print elapsed time.
###################################################################################################

# Finalize timer.
endTime = time.time()

# Print elapsed time for script [s].
print "This script took " + str(endTime - startTime) + "s"

###################################################################################################
