'''   
Copyright (c) 2010-2013, Delft University of Technology
Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
All rights reserved.
See http://bit.ly/12SHPLR for license details.

This script can be used to generate plots to illustrated how conjunction and oppositions are 
detected by the testParticleSimulator application.
'''

################################################################################################### 
# Set up input deck
###################################################################################################

# Set simulation ID.
simulationId    = 17860

# Set absolute path to SQLite database with simulation data.
databasePath    = "/Users/kartikkumar/Documents/Education/PhD/Simulations/" \
                  + "Tudat/Workspace/tudatApplications/stochasticMigration/" \
                  + "stochasticMigrationResults.sqlite"

# Set absolute path to directory with data files.
datafilesPath   = "/Users/kartikkumar/Documents/Education/PhD/PhD Thesis/2_MabsOrbitalMotion/MabIcarusJournalPaper/Figures/Figure15/17860/"

# Set case name.
caseName        = "circular_equatorial"

# Set absolute path to output directory.
outputPath      = "/Users/kartikkumar/Documents/Education/PhD/PhD Thesis/2_MabsOrbitalMotion/MabIcarusJournalPaper/Figures/Figure15/"

# Set figure dpi.
figureDPI       = 600

# Set font size for axes labels.
fontSize        = 24

###################################################################################################

'''
    
                            DO NOT EDIT PARAMETERS BEYOND THIS POINT!!!

'''

################################################################################################### 
# Import Python packages and modules
###################################################################################################

# Import necessary external packages.
import constants
import matplotlib
import matplotlib.pyplot as plt
import numpy
import os
import sqlite3
import time

# ###################################################################################################


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
    
    # Select the case ID corresponding to the test particle case name provided.
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

# Set up base for filenames.
filenameBase = 'simulation' + str(simulationId) + "_"

distanceFilename = filenameBase + 'mutualDistance.csv'

# Read in distance history data and store in array of arrays.
distanceHistory = numpy.genfromtxt(distanceFilename, delimiter = ',', comments = '#', names = True)

oppositionFilename = filenameBase + 'oppositionEvents.csv'

# Read in opposition events data and store in array of arrays.
oppositionEvents = numpy.genfromtxt(oppositionFilename, delimiter = ',', \
                                    comments = '#', names = True)

conjunctionFilename = filenameBase + 'conjunctionEvents.csv'

# Read in conjunction events data and store in array of arrays.
conjunctionEvents = numpy.genfromtxt(conjunctionFilename, delimiter = ',', \
                                     comments = '#', names = True)

# Compute conjunction/opposition detection line crossing points.
rawConjunctionCrossings = []
rawOppositionCrossings = []
switch = False
conjunctionLine = caseData['conjunctionEventDetectionDistance']
oppositionLine = caseData['oppositionEventDetectionDistance']

for i in range(len(distanceHistory) - 2):
    if (switch == False 
        and distanceHistory['mutualDistance'][i] > conjunctionLine
        and distanceHistory['mutualDistance'][i+1] < conjunctionLine ):
        rawConjunctionCrossings.append((distanceHistory['epoch'][i], \
                                        distanceHistory['mutualDistance'][i]))
        switch = True

    elif (switch == True 
          and distanceHistory['mutualDistance'][i] < oppositionLine
          and distanceHistory['mutualDistance'][i+1] > oppositionLine ):
        rawOppositionCrossings.append((distanceHistory['epoch'][i], \
                                       distanceHistory['mutualDistance'][i]))
        switch = False

conjunctionCrossings = numpy.array(rawConjunctionCrossings, dtype=conjunctionEvents.dtype)
oppositionCrossings = numpy.array(rawOppositionCrossings, dtype=oppositionEvents.dtype)

###################################################################################################


###################################################################################################
# Generate figures
###################################################################################################

matplotlib.rcParams.update({'font.size': fontSize})

# Set output path and case-prefix for files generated.
outputPathAndCasePrefix = outputPath + "/case" + str(caseData['caseId']) + "_" + filenameBase   

# Plot conjunction and opposition detection schema.

# Figure #1
figure = plt.figure()
figureAxes = plt.gca()      
figureAxes.yaxis.major.formatter.set_powerlimits((0,0)) 
plt.xlabel("Epoch [Julian years]")
plt.ylabel("Distance wrt. Mab [km]")
plt.ylim(ymin=-1.0e4,ymax=2.1e5)
plt.plot(distanceHistory['epoch']/constants.julianYear, \
         distanceHistory['mutualDistance'] * constants.meterInKilometers, 'k')
plt.plot(oppositionEvents['epoch']/constants.julianYear, \
         oppositionEvents['mutualDistance'] * constants.meterInKilometers, \
         marker='o', markersize='15.0', color='w', linestyle='none')
plt.plot(conjunctionEvents['epoch']/constants.julianYear, \
         conjunctionEvents['mutualDistance'] * constants.meterInKilometers, \
         marker='o', markersize='15.0', color='0.75', linestyle='none')   
xmin, xmax = plt.xlim()
plt.plot([xmin, xmax], \
         [conjunctionLine * constants.meterInKilometers, \
          conjunctionLine * constants.meterInKilometers], \
         linestyle='dashed', color='k')
plt.plot([xmin, xmax], \
         [oppositionLine * constants.meterInKilometers, \
          oppositionLine * constants.meterInKilometers], \
         linestyle='dashed', color='k')
plt.arrow(distanceHistory['epoch'][0]/constants.julianYear - 3.0, \
          distanceHistory['mutualDistance'][0] * constants.meterInKilometers, \
          0.0, -distanceHistory['mutualDistance'][0] * 0.2 * constants.meterInKilometers, \
          head_width=3.0, head_length=1.0e4, fc='k', ec='k', linewidth=5.0)
distanceSegment = distanceHistory[(distanceHistory['epoch'] > distanceHistory['epoch'][0]) \
                                  & (distanceHistory['epoch'] < oppositionCrossings['epoch'][0])]
plt.plot(distanceSegment['epoch']/constants.julianYear, \
         distanceSegment['mutualDistance'] * constants.meterInKilometers, color='k', linewidth=3.0)
plt.plot(distanceHistory['epoch'][0]/constants.julianYear, \
         distanceHistory['mutualDistance'][0] * constants.meterInKilometers, \
         marker='^', markersize='15.0', color='k', linestyle='none')
plt.plot(oppositionCrossings['epoch'][0]/constants.julianYear, \
         oppositionCrossings['mutualDistance'][0] * constants.meterInKilometers, \
         marker='v', markersize='15.0', color='k', linestyle='none')
plt.tight_layout(True)        
plt.savefig(outputPathAndCasePrefix + "eventDetectionSchema1.pdf", \
            dpi = figureDPI)    
plt.close() 

# Figure #2
figure = plt.figure()
figureAxes = plt.gca()      
figureAxes.yaxis.major.formatter.set_powerlimits((0,0)) 
plt.xlabel("Epoch [Julian years]")
plt.ylabel("Distance wrt. Mab [km]")
plt.ylim(ymin=-1.0e4,ymax=2.1e5)
plt.plot(distanceHistory['epoch']/constants.julianYear, \
         distanceHistory['mutualDistance'] * constants.meterInKilometers, 'k')
plt.plot(oppositionEvents['epoch']/constants.julianYear, \
         oppositionEvents['mutualDistance'] * constants.meterInKilometers, \
         marker='o', markersize='15.0', color='w', linestyle='none')
plt.plot(conjunctionEvents['epoch']/constants.julianYear, \
         conjunctionEvents['mutualDistance'] * constants.meterInKilometers, \
         marker='o', markersize='15.0', color='0.75', linestyle='none')   
xmin, xmax = plt.xlim()
plt.plot([xmin, xmax], \
         [conjunctionLine * constants.meterInKilometers, \
          conjunctionLine * constants.meterInKilometers], \
         linestyle='dashed', color='k')
plt.plot([xmin, xmax], \
         [oppositionLine * constants.meterInKilometers, \
          oppositionLine * constants.meterInKilometers], \
         linestyle='dashed', color='k')
distanceSegment = distanceHistory[(distanceHistory['epoch'] > oppositionCrossings['epoch'][0]) \
                                  & (distanceHistory['epoch'] < conjunctionCrossings['epoch'][1])]
plt.plot(distanceSegment['epoch']/constants.julianYear, \
         distanceSegment['mutualDistance'] * constants.meterInKilometers, color='k', linewidth=3.0)
plt.plot(oppositionCrossings['epoch'][0]/constants.julianYear, \
         oppositionCrossings['mutualDistance'][0] * constants.meterInKilometers, \
         marker='^', markersize='15.0', color='k', linestyle='none')
plt.plot(conjunctionCrossings['epoch'][1]/constants.julianYear, \
         conjunctionCrossings['mutualDistance'][1] * constants.meterInKilometers, \
         marker='v', markersize='15.0', color='k', linestyle='none')
plt.plot(oppositionEvents['epoch'][0]/constants.julianYear, \
         oppositionEvents['mutualDistance'][0] * constants.meterInKilometers, \
         marker='o', markersize='15.0', color='k', linestyle='none')
plt.tight_layout(True)        
plt.savefig(outputPathAndCasePrefix + "eventDetectionSchema2.pdf", \
            dpi = figureDPI)    
plt.close() 

# Figure #3
figure = plt.figure()
figureAxes = plt.gca()      
figureAxes.yaxis.major.formatter.set_powerlimits((0,0)) 
plt.xlabel("Epoch [Julian years]")
plt.ylabel("Distance wrt. Mab [km]")
plt.ylim(ymin=-1.0e4,ymax=2.1e5)
plt.plot(distanceHistory['epoch']/constants.julianYear, \
         distanceHistory['mutualDistance'] * constants.meterInKilometers, 'k')
plt.plot(oppositionEvents['epoch']/constants.julianYear, \
         oppositionEvents['mutualDistance'] * constants.meterInKilometers, \
         marker='o', markersize='15.0', color='w', linestyle='none')
plt.plot(conjunctionEvents['epoch']/constants.julianYear, \
         conjunctionEvents['mutualDistance'] * constants.meterInKilometers, \
         marker='o', markersize='15.0', color='0.75', linestyle='none')   
xmin, xmax = plt.xlim()
plt.plot([xmin, xmax], \
         [conjunctionLine * constants.meterInKilometers, \
          conjunctionLine * constants.meterInKilometers], \
         linestyle='dashed', color='k')
plt.plot([xmin, xmax], \
         [oppositionLine * constants.meterInKilometers, \
          oppositionLine * constants.meterInKilometers], \
         linestyle='dashed', color='k')
distanceSegment = distanceHistory[(distanceHistory['epoch'] > conjunctionCrossings['epoch'][1]) \
                                  & (distanceHistory['epoch'] < oppositionCrossings['epoch'][1])]
plt.plot(distanceSegment['epoch']/constants.julianYear, \
         distanceSegment['mutualDistance'] * constants.meterInKilometers, color='k', linewidth=3.0)
plt.plot(conjunctionCrossings['epoch'][1]/constants.julianYear, \
         conjunctionCrossings['mutualDistance'][1] * constants.meterInKilometers, \
         marker='^', markersize='15.0', color='k', linestyle='none')
plt.plot(oppositionCrossings['epoch'][1]/constants.julianYear, \
         oppositionCrossings['mutualDistance'][1] * constants.meterInKilometers, \
         marker='v', markersize='15.0', color='k', linestyle='none')
plt.plot(conjunctionEvents['epoch'][0]/constants.julianYear, \
         conjunctionEvents['mutualDistance'][0] * constants.meterInKilometers, \
         marker='o', markersize='15.0', color='k', linestyle='none')
plt.tight_layout(True)        
plt.savefig(outputPathAndCasePrefix + "eventDetectionSchema3.pdf", \
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
