'''   
Copyright (c) 2010-2013, Delft University of Technology
Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
All rights reserved.
See http://bit.ly/12SHPLR for license details.

This script can be used to generate sample plots showing the data reduction process to derive the
effect of a ring of low-mass perturbers on Mab's orbit.
'''

################################################################################################### 
# Set up input deck
###################################################################################################

# Set Monte Carlo run ID.
monteCarloRunId     = 1

# Set absolute path to directory with data files.
datafilesPath       = "/Users/kartikkumar/Documents/Education/PhD/PhD Thesis/" \
                      + "2_MabsOrbitalMotion/MabIcarusJournalPaper/Figures/Figure16"

# Set absolute path to SQLite database with simulation data.
databasePath        = "/Users/kartikkumar/Desktop/" \
                      + "mabPaperResults.sqlite"

# Set case name.
caseName            = "circular_equatorial_nominal"

# Set absolute path to output directory.
outputPath          = "/Users/kartikkumar/Documents/Education/PhD/PhD Thesis/" \
                      + "2_MabsOrbitalMotion/MabIcarusJournalPaper/Figures/Figure16"

# Set figure dpi.
figureDPI           = 600

# Set font size for axes labels.
fontSize            = 24

###################################################################################################

'''
    
                            DO NOT EDIT PARAMETERS BEYOND THIS POINT!!!

'''


################################################################################################### 
# Import Python packages and modules
###################################################################################################

# Import necessary external packages.
import matplotlib.pyplot as plt
import matplotlib
import numpy
import os
import sqlite3
import time

# Import user-defined modules.
import constants
import epochWindowStepFunctionAverage

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

    # Select observation period start epoch associated with random walk Monte Carlo run.
    cursor.execute("SELECT observationPeriodStartEpoch FROM random_walk_input " \
                   "WHERE monteCarloRunId == " + str(monteCarloRunId))
    observationPeriodStartEpoch = cursor.fetchall()[0][0]

    # Select observation period associated with random walk case.
    cursor.execute("SELECT * FROM random_walk_case " \
                   "WHERE caseName == \"" + caseName + "\";")
    rawRandomWalkCaseData = cursor.fetchall()
    caseDataColumnNameList = [column_name[0] for column_name in cursor.description]   

# Store random walk case data in dictionary, using column names from database.    
randomWalkCaseData = {} 
for i,name in enumerate(caseDataColumnNameList):
    randomWalkCaseData[name] = rawRandomWalkCaseData[0][i]   

# Close database connection if it is open.    
if database:
    database.close()    

###################################################################################################


###################################################################################################
# Read and store datafiles
###################################################################################################

# Change the working directory to the directory containing the data files.
os.chdir(datafilesPath)

# Set up filenames.
filenameBase = 'monteCarloRun' + str(monteCarloRunId) + "_"
keplerianActionElementsFilename = filenameBase + 'keplerianActionElements.csv'
longitudesFilename = filenameBase + 'longitudeHistory.csv'
reducedLongitudesFilename = filenameBase + 'reducedLongitudeHistory.csv'
longitudeResidualsFilename = filenameBase + 'longitudeResiduals.csv'
epochWindowAveragesFilename = filenameBase + "epochWindowAverages.csv"

# Store Keplerian action element array.
keplerianActionElements = numpy.genfromtxt(keplerianActionElementsFilename, delimiter = ',', \
                                           comments = '#', names = True)

# Store longitude array.
longitudes = numpy.genfromtxt(longitudesFilename, delimiter = ',', comments = '#', names = True)

# Store reduced longitude array.
reducedLongitudes = numpy.genfromtxt(reducedLongitudesFilename,\
                                     delimiter = ',', comments = '#', names = True)

# Store longitude residual array.
longitudeResiduals = numpy.genfromtxt(longitudeResidualsFilename, delimiter = ',', \
                                      comments = '#', names = True)

# Store epoch window averages array.
epochWindowAverages = numpy.genfromtxt(epochWindowAveragesFilename, delimiter = ',', \
                                       comments = '#', names = True)

###################################################################################################


###################################################################################################
# Compute derived parameters
###################################################################################################

# Compute epoch window spacing.
epochWindowSpacing = randomWalkCaseData['observationPeriod'] \
    / ( randomWalkCaseData['numberOfEpochWindows'] - 1 )

# Compute central epochs for epoch windows.
centralEpochs = []
centralEpochIndices = []
for i in xrange(0, randomWalkCaseData['numberOfEpochWindows']):
    centralEpochs.append((observationPeriodStartEpoch + i*epochWindowSpacing))
    centralEpochIndices.append( \
        numpy.abs(keplerianActionElements['epoch'] - centralEpochs[i]).argmin())

# Set font size for figures.
matplotlib.rcParams.update({'font.size': fontSize})

# Set figure labels.
figures = ('16', '17', '18')

# Set subfigure labels.
subfigures = ('a', 'b', 'c', 'd')

# Set filename base.
filenameBase  = outputPath + "/figure{0}{1}_monteCarloRun" + str(monteCarloRunId)

# Set up plotting data.
plotCentralEpoch = [centralEpoch/constants.julianYear for centralEpoch in centralEpochs]

plotEpoch = keplerianActionElements['epoch']/constants.julianYear
plotSemiMajorAxis = keplerianActionElements['semiMajorAxis']*constants.meterInKilometers

plotLongitudeEpoch = longitudes['epoch']/constants.julianYear
plotLongitude = longitudes['longitude']*constants.radiansInDegrees

plotReducedLongitudeEpoch = reducedLongitudes['epoch']/constants.julianYear
plotReducedLongitude = reducedLongitudes['longitude']*constants.radiansInDegrees

plotLongitudeResidualEpoch = longitudeResiduals['epoch']/constants.julianYear
plotLongitudeResidual = longitudeResiduals['longitudeResidual']*constants.radiansInDegrees

plotEccentricity = keplerianActionElements['eccentricity']
plotInclination = keplerianActionElements['inclination']*constants.radiansInDegrees

plotObservationPeriodStartEpoch = observationPeriodStartEpoch/constants.julianYear
plotObservationPeriodEndEpoch = (observationPeriodStartEpoch \
                                 + randomWalkCaseData['observationPeriod'])/constants.julianYear

plotEpochWindowAveragesEpoch = epochWindowAverages['epoch']/constants.julianYear
plotEpochWindowLongitudeAverages = epochWindowAverages['longitudeResidual'] \
                                   *constants.radiansInDegrees
plotEpochWindowEccentricityAverages = epochWindowAverages['eccentricity']
plotEpochWindowInclinationAverages = epochWindowAverages['inclination']*constants.radiansInDegrees

plotHalfEpochWindowSize = 0.5*randomWalkCaseData['epochWindowSize']/constants.julianYear

###################################################################################################


###################################################################################################
# Generate longitude residual figures
###################################################################################################

# Plot time-history of Mab's semi-major axis.
fig = plt.figure()
plt.xlabel("Epoch [Julian years]")
plt.ylabel("$a_{Mab}$ [km]")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# Plot semi-major axis data.
plt.plot(plotEpoch, plotSemiMajorAxis, 'k')
# Plot central epochs for the epoch windows.
plt.plot(plotCentralEpoch, plotSemiMajorAxis[centralEpochIndices], \
         marker='o', color='0.80', linestyle='none', markersize = '8.0')
# Get y-axis limits.
ymin, ymax = plt.ylim()
# Plot vertical lines to denote the observation period.
plt.plot([plotObservationPeriodStartEpoch, plotObservationPeriodStartEpoch], [ymin, ymax], \
         linestyle='dashed', color='k')
plt.plot([plotObservationPeriodEndEpoch, plotObservationPeriodEndEpoch], [ymin, ymax], \
         linestyle='dashed', color='k')
# Save to file and close figure.
plt.tight_layout(True)       
outputFilenameBase = filenameBase.format(figures[0],subfigures[0])
plt.savefig(outputFilenameBase + "SemiMajorAxisHistory.pdf", dpi = figureDPI)    
plt.close()

# Plot time-history of Mab's longitude.
fig = plt.figure()
plt.xlabel("Epoch [Julian years]")
plt.ylabel("$L_{Mab}$ [deg]")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# Plot longitude data.
plt.plot(plotLongitudeEpoch, plotLongitude, 'k')
# Plot reduced longitudes for the epoch windows.
plt.plot(plotReducedLongitudeEpoch, plotReducedLongitude, \
         marker='.', linestyle='none', color='k')
# Get y-axis limits.
ymin, ymax = plt.ylim()
# Plot vertical lines to denote the observation period.
plt.plot([plotObservationPeriodStartEpoch, plotObservationPeriodStartEpoch], [ymin, ymax], \
         linestyle='dashed', color='k')
plt.plot([plotObservationPeriodEndEpoch, plotObservationPeriodEndEpoch], [ymin, ymax], \
         linestyle='dashed', color='k')
# Save to file and close figure.
plt.tight_layout(True)       
outputFilenameBase = filenameBase.format(figures[0], subfigures[1]) 
plt.savefig(outputFilenameBase + "LongitudeHistory.pdf", dpi = figureDPI)    
plt.close()  

# Plot epoch window data for longitude residuals.
fig = plt.figure()
axes = fig.gca()
plt.xlabel("Epoch [Julian years]")
plt.ylabel("$\Delta L_{Mab}$ [deg]")
# Plot longitude residual data points in epoch windows.
plt.plot(plotLongitudeResidualEpoch, plotLongitudeResidual, \
         marker='.', color='k', linestyle='none')
# Plot longitude residual averages computed per epoch window.
plt.plot(plotEpochWindowAveragesEpoch, plotEpochWindowLongitudeAverages,\
         marker='o', color='0.25', markersize=15.0, linestyle='none')
# Get y-axis limits.
ymin, ymax = plt.ylim()
# Plot shaded areas to indicate epoch window bounds.
for i in xrange(0,randomWalkCaseData['numberOfEpochWindows']):
  epochWindowLowerBound = plotCentralEpoch[i] - plotHalfEpochWindowSize
  epochWindowUpperBound = plotCentralEpoch[i] + plotHalfEpochWindowSize                         
  axes.fill_between([epochWindowLowerBound, epochWindowUpperBound], \
                     ymin, ymax, facecolor='0.75', alpha=0.5)
# Save to file and close figure. 
plt.ylim(ymin, ymax)
plt.tight_layout(True)       
outputFilenameBase = filenameBase.format(figures[0], subfigures[2])
plt.savefig(outputFilenameBase + "LongitudeResidualEpochWindows.pdf", dpi = figureDPI)    
plt.close()    

# Plot longitude residuals in first epoch window.
fig = plt.figure()
axes = fig.gca() 
plt.xlabel("Epoch [Julian years]")
plt.ylabel("$\Delta L_{Mab}$ [deg]")
# Extract data for first epoch window.
firstEpochWindowLowerBound = plotCentralEpoch[0] - plotHalfEpochWindowSize
firstEpochWindowUpperBound = plotCentralEpoch[0] + plotHalfEpochWindowSize   
firstWindowIndices = (plotLongitudeResidualEpoch > firstEpochWindowLowerBound) \
                     & (plotLongitudeResidualEpoch < firstEpochWindowUpperBound)
for i in xrange(0,len(firstWindowIndices)-1):
  if firstWindowIndices[i] == False and firstWindowIndices[i+1] == True:
    firstWindowIndices[i] = True
    break
firstWindowLongitudeResidualEpoch = plotLongitudeResidualEpoch[firstWindowIndices]
firstWindowLongitudeResidual = plotLongitudeResidual[firstWindowIndices]    
# Manually create step-function.
stepEpoch = []
stepLongitudeResidual = []
for i in xrange(0,len(firstWindowLongitudeResidualEpoch)-1):
  stepEpoch.append(firstWindowLongitudeResidualEpoch[i])
  stepLongitudeResidual.append(firstWindowLongitudeResidual[i])
  stepEpoch.append(firstWindowLongitudeResidualEpoch[i+1])
  stepLongitudeResidual.append(firstWindowLongitudeResidual[i])
# Plot longitude residual step-function.
plt.step(firstWindowLongitudeResidualEpoch, firstWindowLongitudeResidual, \
         marker='o', color='k', linestyle='dashed', linewidth = 2.0, where='post')
plt.step([firstWindowLongitudeResidualEpoch[-1], firstEpochWindowUpperBound], \
         [firstWindowLongitudeResidual[-1], firstWindowLongitudeResidual[-1]], \
         color='k', linestyle='dashed', linewidth = 2.0, where='post')
# Plot epoch window average.
plt.plot(plotEpochWindowAveragesEpoch[0], plotEpochWindowLongitudeAverages[0],\
         marker='o', color='0.25', markersize=15.0, linestyle='none')
# Plot shaded area for epoch window.
ymin, ymax = plt.ylim()
stepEpoch[0] = firstEpochWindowLowerBound
stepEpoch[-1] = firstEpochWindowUpperBound
axes.fill_between(stepEpoch, ymin, stepLongitudeResidual, facecolor='0.75', alpha=0.5)
# Save to file and close figure. 
plt.ylim(ymin, ymax)
plt.tight_layout(True)       
outputFilenameBase = filenameBase.format(figures[0], subfigures[3])
plt.savefig(outputFilenameBase + "LongitudeResidualFirstEpochWindow.pdf", dpi = figureDPI)    
plt.close()    

# Compute epoch window averages to compare with averages in data file.
# for i in xrange(0,randomWalkCaseData['numberOfEpochWindows']):  
#   print epochWindowStepFunctionAverage.computeStepFunctionAverage( \
#     longitudeResiduals['epoch'], longitudeResiduals['longitudeResidual'], \
#     centralEpochs[i] - 0.5*randomWalkCaseData['epochWindowSize'], \
#     centralEpochs[i] + 0.5*randomWalkCaseData['epochWindowSize'])*constants.radiansInDegrees
#   print plotEpochWindowLongitudeAverages[i]

###################################################################################################


###################################################################################################
# Generate eccentricity figures
###################################################################################################

# Plot time-history of Mab's eccentricity.
fig = plt.figure()
plt.xlabel("Epoch [Julian years]")
plt.ylabel("$e_{Mab}$ [-]")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# Plot eccentricity data.
plt.plot(plotEpoch, plotEccentricity, 'k')
# Plot central epochs for the epoch windows.
plt.plot(plotCentralEpoch, plotEccentricity[centralEpochIndices], \
         marker='o', color='0.80', linestyle='none', markersize = '8.0')
# Get y-axis limits.
ymin, ymax = plt.ylim()
# Plot vertical lines to denote the observation period.
plt.plot([plotObservationPeriodStartEpoch, plotObservationPeriodStartEpoch], [ymin, ymax], \
         linestyle='dashed', color='k')
plt.plot([plotObservationPeriodEndEpoch, plotObservationPeriodEndEpoch], [ymin, ymax], \
         linestyle='dashed', color='k')
# Save to file and close figure.
plt.tight_layout(True)       
outputFilenameBase = filenameBase.format(figures[1], subfigures[0])
plt.savefig(outputFilenameBase + "EccentricityHistory.pdf", dpi = figureDPI)    
plt.close()  

# Plot epoch window data for eccentricity.
fig = plt.figure()
axes = fig.gca()
plt.xlabel("Epoch [Julian years]")
plt.ylabel("$e_{Mab}$ [-]")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# Extract data for epoch windows.
epochWindowIndices = []
for j in xrange(0,len(plotCentralEpoch)):
  epochWindowLowerBound = plotCentralEpoch[j] - plotHalfEpochWindowSize
  epochWindowUpperBound = plotCentralEpoch[j] + plotHalfEpochWindowSize   
  singleEpochWindowIndices = (plotEpoch > epochWindowLowerBound) \
                              & (plotEpoch < epochWindowUpperBound)
  for i in xrange(0,len(singleEpochWindowIndices)-1):
    if singleEpochWindowIndices[i] == False and singleEpochWindowIndices[i+1] == True:
      singleEpochWindowIndices[i] = True
      break                              
  epochWindowIndices.append(singleEpochWindowIndices)

# Plot eccentricity data points in epoch windows.
for i in xrange(0,len(epochWindowIndices)):  
  epochWindowEccentricityEpoch = plotEpoch[epochWindowIndices[i]]
  epochWindowEccentricity = plotEccentricity[epochWindowIndices[i]]    
  plt.plot(epochWindowEccentricityEpoch, epochWindowEccentricity, \
           marker='.', color='k', linestyle='none')
# Plot eccentricity averages computed per epoch window.
plt.plot(plotEpochWindowAveragesEpoch, plotEpochWindowEccentricityAverages,\
         marker='o', color='0.25', markersize=15.0, linestyle='none')
# Get y-axis limits.
ymin, ymax = plt.ylim()
# Plot shaded areas to indicate epoch window bounds.
for i in xrange(0,randomWalkCaseData['numberOfEpochWindows']):
  epochWindowLowerBound = plotCentralEpoch[i] - plotHalfEpochWindowSize
  epochWindowUpperBound = plotCentralEpoch[i] + plotHalfEpochWindowSize                         
  axes.fill_between([epochWindowLowerBound, epochWindowUpperBound], \
                     ymin, ymax, facecolor='0.75', alpha=0.5)
# Save to file and close figure. 
plt.ylim(ymin, ymax)
plt.tight_layout(True)       
outputFilenameBase = filenameBase.format(figures[1], subfigures[1]) 
plt.savefig(outputFilenameBase + "EccentricityEpochWindows.pdf", dpi = figureDPI)    
plt.close()  

# Plot eccentricity in first epoch window.
fig = plt.figure()
axes = fig.gca() 
plt.xlabel("Epoch [Julian years]")
plt.ylabel("$e_{Mab}$ [-]")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# Extract data for first epoch window.
firstEpochWindowLowerBound = plotCentralEpoch[0] - plotHalfEpochWindowSize
firstEpochWindowUpperBound = plotCentralEpoch[0] + plotHalfEpochWindowSize   
firstWindowIndices = (plotEpoch > firstEpochWindowLowerBound) \
                     & (plotEpoch < firstEpochWindowUpperBound)
for i in xrange(0,len(firstWindowIndices)-1):
  if firstWindowIndices[i] == False and firstWindowIndices[i+1] == True:
    firstWindowIndices[i] = True
    break
firstWindowEccentricityEpoch = plotEpoch[firstWindowIndices]
firstWindowEccentricity = plotEccentricity[firstWindowIndices]    
# Manually create step-function.
stepEpoch = []
stepEccentricity = []
for i in xrange(0,len(firstWindowEccentricityEpoch)-1):
  stepEpoch.append(firstWindowEccentricityEpoch[i])
  stepEccentricity.append(firstWindowEccentricity[i])
  stepEpoch.append(firstWindowEccentricityEpoch[i+1])
  stepEccentricity.append(firstWindowEccentricity[i])
# Plot longitude residual step-function.
plt.step(firstWindowEccentricityEpoch, firstWindowEccentricity, \
         marker='o', color='k', linestyle='dashed', linewidth = 2.0, where='post')
plt.step([firstWindowEccentricityEpoch[-1], firstEpochWindowUpperBound], \
         [firstWindowEccentricity[-1], firstWindowEccentricity[-1]], \
         color='k', linestyle='dashed', linewidth = 2.0, where='post')
# Plot epoch window average.
plt.plot(plotEpochWindowAveragesEpoch[0], plotEpochWindowEccentricityAverages[0],\
         marker='o', color='0.25', markersize=15.0, linestyle='none')
# Plot shaded area for epoch window.
ymin, ymax = plt.ylim()
stepEpoch[0] = firstEpochWindowLowerBound
stepEpoch[-1] = firstEpochWindowUpperBound
axes.fill_between(stepEpoch, ymin, stepEccentricity, facecolor='0.75', alpha=0.5)
# Save to file and close figure. 
plt.ylim(ymin, ymax)
plt.tight_layout(True)       
outputFilenameBase = filenameBase.format(figures[1], subfigures[2])
plt.savefig(outputFilenameBase + "EccentricityFirstEpochWindow.pdf", dpi = figureDPI)    
plt.close()     

# Compute epoch window averages to compare with averages in data file.
# for i in xrange(0,randomWalkCaseData['numberOfEpochWindows']):  
#   print epochWindowStepFunctionAverage.computeStepFunctionAverage( \
#     keplerianActionElements['epoch'], keplerianActionElements['eccentricity'], \
#     centralEpochs[i] - 0.5*randomWalkCaseData['epochWindowSize'], \
#     centralEpochs[i] + 0.5*randomWalkCaseData['epochWindowSize'])
#   print plotEpochWindowEccentricityAverages[i]

###################################################################################################


###################################################################################################
# Generate inclination figures
###################################################################################################

# Plot time-history of Mab's inclination.
fig = plt.figure()
plt.xlabel("Epoch [Julian years]")
plt.ylabel("$i_{Mab}$ [deg]")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# Plot inclination data.
plt.plot(plotEpoch, plotInclination, 'k')
# Plot central epochs for the epoch windows.
plt.plot(plotCentralEpoch, plotInclination[centralEpochIndices], \
         marker='o', color='0.80', linestyle='none', markersize = '8.0')
# Get y-axis limits.
ymin, ymax = plt.ylim()
# Plot vertical lines to denote the observation period.
plt.plot([plotObservationPeriodStartEpoch, plotObservationPeriodStartEpoch], [ymin, ymax], \
         linestyle='dashed', color='k')
plt.plot([plotObservationPeriodEndEpoch, plotObservationPeriodEndEpoch], [ymin, ymax], \
         linestyle='dashed', color='k')
# Save to file and close figure.
plt.tight_layout(True)       
outputFilenameBase = filenameBase.format(figures[2], subfigures[0]) 
plt.savefig(outputFilenameBase + "InclinationHistory.pdf", dpi = figureDPI)    
plt.close()  

# Plot epoch window data for eccentricity.
fig = plt.figure()
axes = fig.gca()
plt.xlabel("Epoch [Julian years]")
plt.ylabel("$i_{Mab}$ [deg]")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# Extract data for epoch windows.
epochWindowIndices = []
for j in xrange(0,len(plotCentralEpoch)):
  epochWindowLowerBound = plotCentralEpoch[j] - plotHalfEpochWindowSize
  epochWindowUpperBound = plotCentralEpoch[j] + plotHalfEpochWindowSize   
  singleEpochWindowIndices = (plotEpoch > epochWindowLowerBound) \
                              & (plotEpoch < epochWindowUpperBound)
  for i in xrange(0,len(singleEpochWindowIndices)-1):
    if singleEpochWindowIndices[i] == False and singleEpochWindowIndices[i+1] == True:
      singleEpochWindowIndices[i] = True
      break                              
  epochWindowIndices.append(singleEpochWindowIndices)

# Plot inclination data points in epoch windows.
for i in xrange(0,len(epochWindowIndices)):  
  epochWindowInclinationEpoch = plotEpoch[epochWindowIndices[i]]
  epochWindowInclination = plotInclination[epochWindowIndices[i]]    
  plt.plot(epochWindowInclinationEpoch, epochWindowInclination, \
           marker='.', color='k', linestyle='none')
# Plot inclination averages computed per epoch window.
plt.plot(plotEpochWindowAveragesEpoch, plotEpochWindowInclinationAverages,\
         marker='o', color='0.25', markersize=15.0, linestyle='none')
# Get y-axis limits.
ymin, ymax = plt.ylim()
# Plot shaded areas to indicate epoch window bounds.
for i in xrange(0,randomWalkCaseData['numberOfEpochWindows']):
  epochWindowLowerBound = plotCentralEpoch[i] - plotHalfEpochWindowSize
  epochWindowUpperBound = plotCentralEpoch[i] + plotHalfEpochWindowSize                         
  axes.fill_between([epochWindowLowerBound, epochWindowUpperBound], \
                     ymin, ymax, facecolor='0.75', alpha=0.5)
# Save to file and close figure. 
plt.ylim(ymin, ymax)
plt.tight_layout(True)       
outputFilenameBase = filenameBase.format(figures[2], subfigures[1]) 
plt.savefig(outputFilenameBase + "InclinationEpochWindows.pdf", dpi = figureDPI)    
plt.close()  

# Plot inclination in first epoch window.
fig = plt.figure()
axes = fig.gca() 
plt.xlabel("Epoch [Julian years]")
plt.ylabel("$i_{Mab}$ [deg]")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# Extract data for first epoch window.
firstEpochWindowLowerBound = plotCentralEpoch[0] - plotHalfEpochWindowSize
firstEpochWindowUpperBound = plotCentralEpoch[0] + plotHalfEpochWindowSize   
firstWindowIndices = (plotEpoch > firstEpochWindowLowerBound) \
                     & (plotEpoch < firstEpochWindowUpperBound)
for i in xrange(0,len(firstWindowIndices)-1):
  if firstWindowIndices[i] == False and firstWindowIndices[i+1] == True:
    firstWindowIndices[i] = True
    break
firstWindowInclinationEpoch = plotEpoch[firstWindowIndices]
firstWindowInclination = plotInclination[firstWindowIndices]    
# Manually create step-function.
stepEpoch = []
stepInclination = []
for i in xrange(0,len(firstWindowInclinationEpoch)-1):
  stepEpoch.append(firstWindowInclinationEpoch[i])
  stepInclination.append(firstWindowInclination[i])
  stepEpoch.append(firstWindowInclinationEpoch[i+1])
  stepInclination.append(firstWindowInclination[i])
# Plot longitude residual step-function.
plt.step(firstWindowInclinationEpoch, firstWindowInclination, \
         marker='o', color='k', linestyle='dashed', linewidth = 2.0, where='post')
plt.step([firstWindowInclinationEpoch[-1], firstEpochWindowUpperBound], \
         [firstWindowInclination[-1], firstWindowInclination[-1]], \
         color='k', linestyle='dashed', linewidth = 2.0, where='post')
# Plot epoch window average.
plt.plot(plotEpochWindowAveragesEpoch[0], plotEpochWindowInclinationAverages[0],\
         marker='o', color='0.25', markersize=15.0, linestyle='none')
# Plot shaded area for epoch window.
ymin, ymax = plt.ylim()
stepEpoch[0] = firstEpochWindowLowerBound
stepEpoch[-1] = firstEpochWindowUpperBound
axes.fill_between(stepEpoch, ymin, stepInclination, facecolor='0.75', alpha=0.5)
# Save to file and close figure. 
plt.ylim(ymin, ymax)
plt.tight_layout(True)       
outputFilenameBase = filenameBase.format(figures[2], subfigures[2]) 
plt.savefig(outputFilenameBase + "InclinationFirstEpochWindow.pdf", dpi = figureDPI)    
plt.close()   

# Compute epoch window averages to compare with averages in data file.
# for i in xrange(0,randomWalkCaseData['numberOfEpochWindows']):  
#   print epochWindowStepFunctionAverage.computeStepFunctionAverage( \
#     keplerianActionElements['epoch'], keplerianActionElements['inclination'], \
#     centralEpochs[i] - 0.5*randomWalkCaseData['epochWindowSize'], \
#     centralEpochs[i] + 0.5*randomWalkCaseData['epochWindowSize'])*constants.radiansInDegrees
#   print plotEpochWindowInclinationAverages[i]

###################################################################################################


###################################################################################################
# Finalize timer and print elapsed time.
###################################################################################################

# Finalize timer.
endTime = time.time()

# Print elapsed time for script [s].
print "This script took " + str(endTime - startTime) + "s"

###################################################################################################
