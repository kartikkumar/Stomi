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
monteCarloRunIds    = (1,1)

# Set absolute path to directory with data files.
datafilesPath       = "/Users/kartikkumar/Documents/Education/PhD/PhD Thesis/2_MabsOrbitalMotion/MabIcarusJournalPaper/Figures/Figure9"

# Set absolute path to SQLite database with simulation data.
databasePath        = "/Users/kartikkumar/Documents/Education/PhD/Simulations/" \
                      + "Tudat/Workspace/tudatApplications/stochasticMigration/" \
                      + "stochasticMigrationResults.sqlite"

# Set case name.
caseName            = "circular_equatorial_nominal"

# Set absolute path to output directory.
outputPath          = "/Users/kartikkumar/Documents/Education/PhD/PhD Thesis/2_MabsOrbitalMotion/MabIcarusJournalPaper/Figures/Figure9"

# Set figure dpi.
figureDPI           = 600

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
import matplotlib
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
    
    # Select test particle case ID associated with random walk case.
    cursor.execute("SELECT testParticleCaseId FROM random_walk_case " \
                   "WHERE caseName == \"" + caseName + "\";")
    testParticleCaseId = cursor.fetchall()[0][0]

    # Select observation period associated with random walk case.
    cursor.execute("SELECT observationPeriod FROM random_walk_case " \
                   "WHERE caseName == \"" + caseName + "\";")
    observationPeriod = cursor.fetchall()[0][0]    

    # Select all the test particle case data associated with the case ID.
    cursor.execute("SELECT * FROM test_particle_case WHERE caseId == " \
                   + str(testParticleCaseId) + ";")
    rawTestParticleCaseData = cursor.fetchall()
    caseDataColumnNameList = [column_name[0] for column_name in cursor.description]      

# Store random walk case data in dictionary, using column names from database.    
testParticleCaseData = {} 
for i,name in enumerate(caseDataColumnNameList):
    testParticleCaseData[name] = rawTestParticleCaseData[0][i]   

# Close database connection if it is open.    
if database:
    database.close()    

###################################################################################################


###################################################################################################
# Read and store datafiles
###################################################################################################

# Change the working directory to the directory containing the data files.
os.chdir(datafilesPath)

# Declare empty array to store Keplerian action element arrays.
keplerianActionElements = []

# Declare empty array to store longitude residual arrays.
longitudeResiduals = []

for monteCarloRunId in monteCarloRunIds:
    # Set up base for filenames.
    filenameBase = 'monteCarloRun' + str(monteCarloRunId) + "_"

    keplerianActionElementsFilename = filenameBase + 'keplerianActionElements.csv'

    # Read in Keplerian action elements history for Monte Carlo run.
    keplerianActionElements.append(numpy.genfromtxt(keplerianActionElementsFilename, \
                                   delimiter = ',', comments = '#', names = True) )

    longitudeResidualsFilename = filenameBase + 'longitudeResiduals.csv'

    # Read in longitude residual history for Monte Carlo run.
    longitudeResiduals.append(numpy.genfromtxt(longitudeResidualsFilename, \
                              delimiter = ',', comments = '#', names = True) )    

###################################################################################################


###################################################################################################
# Generate figures
###################################################################################################

matplotlib.rcParams.update({'font.size': 18})

for i,keplerData in enumerate(keplerianActionElements):
    # Set output path and case-prefix for files generated.
    outputPathAndCasePrefix = outputPath + "/monteCarloRun" + str(monteCarloRunIds[i]) + "_"

    # Plot time-histories of Keplerian elements.
    fig = plt.figure()
    plt.xlabel("Epoch [Julian years]")
    plt.ylabel("$\Delta a_{Mab}$ [km]")
    plt.xlim(xmin = 0.0, xmax = 50.0)
    plt.plot(keplerData['epoch']/constants.julianYear, \
             (keplerData['semiMajorAxis'] - testParticleCaseData['perturbedBodySemiMajorAxisAtT0']) \
             * constants.meterInKilometers, 'k')
    plt.tight_layout(True)    
    plt.savefig(outputPathAndCasePrefix + "semiMajorAxisHistory.pdf", dpi = figureDPI)    
    plt.close()

    fig = plt.figure()
    plt.xlabel("Epoch [Julian years]")
    plt.ylabel("$\Delta e_{Mab}$ [-]")
    plt.xlim(xmin = 0.0, xmax = 50.0)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(keplerData['epoch']/constants.julianYear, \
             keplerData['eccentricity'] - testParticleCaseData['perturbedBodyEccentricityAtT0'], 'k')
    plt.tight_layout(True)        
    plt.savefig(outputPathAndCasePrefix + "eccentricityHistory.pdf", dpi = figureDPI)    
    plt.close()    

    fig = plt.figure()
    plt.xlabel("Epoch [Julian years]")
    plt.ylabel("$\Delta i_{Mab}$ [deg]")
    plt.xlim(xmin = 0.0, xmax = 50.0)
    plt.plot(keplerData['epoch']/constants.julianYear, \
             ( keplerData['inclination'] - testParticleCaseData['perturbedBodyInclinationAtT0'] ) \
             * constants.radiansInDegrees, 'k')
    plt.tight_layout(True)        
    plt.savefig(outputPathAndCasePrefix + "inclinationHistory.pdf", dpi = figureDPI)    
    plt.close()        

    fig = plt.figure()
    plt.xlabel("Epoch [Julian years]")
    plt.ylabel("$\Delta L_{Mab}$ [deg]")
    plt.xlim(xmin = 0.0, xmax = observationPeriod/constants.julianYear)
    plt.plot(longitudeResiduals[i]['epoch']/constants.julianYear, \
             longitudeResiduals[i]['longitudeResidual'] * constants.radiansInDegrees, 'k')
    plt.tight_layout(True)        
    plt.savefig(outputPathAndCasePrefix + "longitudeResidualHistory.pdf", dpi = figureDPI)    
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
