'''   
Copyright (c) 2010-2013, Delft University of Technology
Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
All rights reserved.
See http://bit.ly/12SHPLR for license details.

This script can be used to generate histograms of Mab's orbital element variations due to a 
perturber ring.
'''

################################################################################################### 
# Set up input deck
###################################################################################################

# Set absolute path to SQLite database with simulation data.
databasePath        = "/Users/kartikkumar/Documents/Education/PhD/Simulations/" \
                      + "Tudat/Workspace/tudatApplications/stochasticMigration/" \
                      + "stochasticMigrationResults.sqlite"

# Set case name.
caseName        = "circular_equatorial_nominal"

# Set absolute path to output directory.
outputPath      = "/Users/kartikkumar/Desktop/Figure10"

# Set figure dpi.
figureDPI       = 600

###################################################################################################

'''
    
                            DO NOT EDIT PARAMETERS BEYOND THIS POINT!!!

'''

################################################################################################### 
# Import Python packages and modules
###################################################################################################

# # Import necessary external packages.
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy
import os
import sqlite3
import time

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
    
    # Select the case ID corresponding to the random walk case name provided.
    cursor.execute("SELECT caseId FROM random_walk_case WHERE caseName == \"" + caseName + "\";")
    caseId = cursor.fetchall()[0][0]  
    
    # Select all the random walk case data associated with the case ID.
    cursor.execute("SELECT * FROM random_walk_case WHERE caseId == " + str(caseId) + ";")
    rawCaseData = cursor.fetchall()
    caseDataColumnNameList = [column_name[0] for column_name in cursor.description]  

    # Select all the random walk input data associated with the case ID.
    cursor.execute("SELECT * FROM random_walk_input WHERE randomWalkCaseId == " \
                    + str(caseId) + ";")
    rawInputData = cursor.fetchall()
    inputDataColumnNameList = [column_name[0] for column_name in cursor.description]     

    # Select all the random walk output data associated with the case ID.
    cursor.execute("SELECT * FROM random_walk_output WHERE random_walk_output.monteCarloRunId \
                    IN (SELECT monteCarloRunId FROM random_walk_input \
                        WHERE randomWalkCaseId == " + str(caseId) + ");")
    rawOutputData = cursor.fetchall()
    outputDataColumnNameList = [column_name[0] for column_name in cursor.description]     
    
# Close database connection if it is open.    
if database:
    database.close()    
    
# Store case data in dictionary, using column names from database.    
caseData = {} 
for i,name in enumerate(caseDataColumnNameList):
    caseData[name] = rawCaseData[0][i]    

# Store input data in structured array, using column names from database.
inputDataTypeList = [] 
for name in inputDataColumnNameList:
   inputDataTypeList.append((name, '<f8'))

inputData = numpy.array(rawInputData, dtype = numpy.dtype(inputDataTypeList))          

# Store output data in structured array, using column names from database.
outputDataTypeList = [] 
for name in outputDataColumnNameList:
   outputDataTypeList.append((name, '<f8'))
   
outputData = numpy.array(rawOutputData, dtype = numpy.dtype(outputDataTypeList))  

###################################################################################################


###################################################################################################
# Plot histograms of random walk output data.
###################################################################################################

matplotlib.rcParams.update({'font.size': 18})

# Set output path and case-prefix for files generated.
outputPathAndCasePrefix = outputPath + "/randomWalkCase" + str(caseData['caseId'])

# Plot histogram of obsevation period start epochs used for random walk simulations.
plt.figure()
plt.hist(inputData['observationPeriodStartEpoch'], facecolor='w', edgecolor='k')
plt.xlabel("Epoch at start of observation period [s]")
plt.ylabel("Frequency [-]")
plt.savefig(outputPathAndCasePrefix + "HistogramObservationPeriodStartEpoch.pdf", dpi = figureDPI)
plt.close()

# Plot histograms of maximum changes of perturbed body's orbital elements, based on reduction of 
# random walk simulation data.
figure = plt.figure()
plt.hist(numpy.rad2deg(outputData['maximumLongitudeResidualChange']), \
         facecolor='w', edgecolor='k', bins=20)
plt.xlabel("Maximum longitude residual change [deg]")
plt.ylabel("Frequency [-]")
plt.yticks([0.0,50.0,100.0,150.0,200.0])
plt.savefig(outputPathAndCasePrefix + "HistogramMaximumLongitudeResidualChange.pdf", \
            dpi = figureDPI)
plt.close()

figure = plt.figure()
figureAxes = plt.gca()        
plt.hist(outputData['maximumEccentricityChange'], facecolor='w', edgecolor='k', bins=20)
figureAxes.xaxis.major.formatter.set_powerlimits((0,0)) 
plt.xlabel("Maximum eccentricity change [-]")
plt.ylabel("Frequency [-]")
plt.yticks([0.0,50.0,100.0,150.0,200.0,250.0,300.0])
plt.savefig(outputPathAndCasePrefix + "HistogramMaximumEccentricityChange.pdf", \
            dpi = figureDPI)
plt.close()

figure = plt.figure()
figureAxes = plt.gca()        
plt.hist(numpy.rad2deg(outputData['maximumInclinationChange']), \
         facecolor='w', edgecolor='k', bins=20)
figureAxes.xaxis.major.formatter.set_powerlimits((0,0)) 
plt.xlabel("Maximum inclination change [deg]")
plt.ylabel("Frequency [-]")
plt.yticks([0.0,50.0,100.0,150.0])
plt.savefig(outputPathAndCasePrefix + "HistogramMaximumInclinationChange.pdf", \
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
