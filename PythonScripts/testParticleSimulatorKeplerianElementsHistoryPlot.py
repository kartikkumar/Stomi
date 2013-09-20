'''   
Copyright (c) 2010-2013, Delft University of Technology
Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
All rights reserved.
See http://bit.ly/12SHPLR for license details.

This script can be used to generate plots based on the output from a single simulation executed
with the test particle simulator.
'''

################################################################################################### 
# Set up input deck
###################################################################################################

# Set simulation IDs as a list.
simulationIds   =

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

# Import necessary external packages.
import itertools
import matplotlib.pyplot as plt
import numpy
import os
import sqlite3
import time

# Import user-defined modules.
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

# Loop through simulation IDs, read Keplerian elements datafiles, store in array of arrays.
for (counter, simulationId) in enumerate(simulationIds):
    # Set up base for filenames.
    filenameBase = 'simulation' + str(simulationId) + "_"

    keplerianElementsFilename = filenameBase + 'keplerianElements.csv'

    # Read in Keplerian elements history for test particle.
    keplerianElements.append(numpy.genfromtxt(keplerianElementsFilename, \
                             delimiter = ',', comments = '#', names = True) )

###################################################################################################


###################################################################################################
# Generate figures
###################################################################################################

# Set output path and case-prefix for files generated.
outputPathAndCasePrefix = outputPath + "/case" + str(caseData['caseId']) + "_"

# Set plot markers to cycle through.
linestyles = itertools.cycle(['solid','dashed'])

plt.figure()
figureAxes = plt.gca()
figureAxes.spines['left'].set_position(('data',0))
figureAxes.spines['right'].set_color('none')
figureAxes.spines['bottom'].set_position(('data',0))
figureAxes.spines['top'].set_color('none')
figureAxes.spines['left'].set_smart_bounds(True)
figureAxes.spines['bottom'].set_smart_bounds(True)
figureAxes.xaxis.set_ticks_position('bottom')
figureAxes.yaxis.set_ticks_position('left')
plt.xlabel("Epoch [Julian years]")
plt.ylabel("Semi-major axis wrt Mab [km]")
figureAxes.xaxis.set_label_coords(0.85, 0.6)
plt.xlim(xmin = 0.0, xmax = 50.0)
for keplerData in keplerianElements:
    plt.plot(keplerData['epoch']/constants.julianYear, \
             (keplerData['semiMajorAxis'] - caseData['perturbedBodySemiMajorAxisAtT0']) \
             * constants.meterInKilometers, 'k', linestyle = next(linestyles))
plt.savefig(outputPathAndCasePrefix + "exampleRelativeSemiMajorAxisHistories.pdf", \
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