'''   
Copyright (c) 2010-2013, Delft University of Technology
Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
All rights reserved.
See http://bit.ly/12SHPLR for license details.

This script can be used to generate a plot of a test particle's mutual distance history.
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
datafilesPath   = "/Users/kartikkumar/Documents/Education/PhD/PhD Thesis/2_MabsOrbitalMotion/MabIcarusJournalPaper/Figures/Figure14/17860/"

# Set case name.
caseName        = "circular_equatorial"

# Set absolute path to output directory.
outputPath      = "/Users/kartikkumar/Documents/Education/PhD/PhD Thesis/2_MabsOrbitalMotion/MabIcarusJournalPaper/Figures/Figure14/"

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
from matplotlib.path import Path
import matplotlib.patches as patches
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

###################################################################################################


###################################################################################################
# Generate figures
###################################################################################################

matplotlib.rcParams.update({'font.size': fontSize})

# Set output path and case-prefix for files generated.
outputPathAndCasePrefix = outputPath + "/case" + str(caseData['caseId']) + "_" + filenameBase

# Set up vertices for zoom box.
vertices = [ (6.0, 1.85e5), # left, bottom
             (6.0, 2.05e5), # left, top
             (11.55, 2.05e5), # right, top
             (11.55, 1.85e5), # right, bottom
             (0., 0.), # ignored
           ]

# Set up instructions to plot zoom box.
codes = [ Path.MOVETO,
          Path.LINETO,
          Path.LINETO,
          Path.LINETO,
          Path.CLOSEPOLY,
        ]

# Set up path for plotting zoom box.
path = Path(vertices, codes)        

# Plot distance history with zoom box, conjunction events, and opposition events included.
figure = plt.figure()
axes = figure.add_subplot(111)
patch = patches.PathPatch(path, facecolor='none', linestyle='solid', linewidth=2)
axes.add_patch(patch)
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
         [caseData['conjunctionEventDetectionDistance'] * constants.meterInKilometers, \
          caseData['conjunctionEventDetectionDistance'] * constants.meterInKilometers], \
         linestyle='dashed', color='k')
plt.plot([xmin, xmax], \
         [caseData['oppositionEventDetectionDistance'] * constants.meterInKilometers, \
          caseData['oppositionEventDetectionDistance'] * constants.meterInKilometers], \
         linestyle='dashed', color='k')
plt.plot()
plt.tight_layout(True)        
plt.savefig(outputPathAndCasePrefix + "distanceHistory.pdf", \
            dpi = figureDPI)    
plt.close()

# Plot zoomed in view of distance history.
plt.figure()
figureAxes = plt.gca()        
figureAxes.locator_params(tight=True, nbins=4)
figureAxes.yaxis.major.formatter.set_powerlimits((0,0)) 
plt.xlabel("Epoch [Julian years]")
plt.ylabel("Distance wrt. Mab [km]")
plt.xlim(xmin=8.7,xmax=8.85)    
plt.ylim(ymin=1.95e5,ymax=1.96e5)
plt.plot(distanceHistory['epoch']/constants.julianYear, \
         distanceHistory['mutualDistance'] * constants.meterInKilometers, 'k')
plt.plot(oppositionEvents['epoch']/constants.julianYear, \
         oppositionEvents['mutualDistance'] * constants.meterInKilometers, \
         marker='o', markersize='20.0', color='w')
plt.tight_layout(True)                       
plt.savefig(outputPathAndCasePrefix + "distanceHistoryZoom.pdf", \
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
