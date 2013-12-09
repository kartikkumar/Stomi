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
databasePath        = "/Users/kartikkumar/Documents/Education/PhD/Simulations/" \
                      + "Tudat/Workspace/tudatApplications/stochasticMigration/" \
                      + "stochasticMigrationResults.sqlite"

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

# Shift longitude residual epochs.
longitudeResiduals['epoch'] += observationPeriodStartEpoch

# Shift epoch window average epochs.
epochWindowAverages['epoch'] += observationPeriodStartEpoch

# Set font size for figures.
matplotlib.rcParams.update({'font.size': fontSize})

# Set subfigure labels.
subfigures = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i')

# Set filename base.
filenameBase  = outputPath + "/figure16%s_monteCarloRun" + str(monteCarloRunId)

# Set up plotting data.
plotCentralEpoch = [centralEpoch/constants.julianYear for centralEpoch in centralEpochs]

plotEpoch = keplerianActionElements['epoch']/constants.julianYear
plotSemiMajorAxis = keplerianActionElements['semiMajorAxis']*constants.meterInKilometers

plotLongitudeResidualEpoch = longitudeResiduals['epoch']/constants.julianYear
plotLongitudeResidual = longitudeResiduals['longitudeResidual']*constants.radiansInDegrees

plotObservationPeriodStartEpoch = observationPeriodStartEpoch/constants.julianYear
plotObservationPeriodEndEpoch = (observationPeriodStartEpoch \
                                 + randomWalkCaseData['observationPeriod'])/constants.julianYear

plotEpochWindowAveragesEpoch = epochWindowAverages['epoch']/constants.julianYear
plotEpochWindowLongitudeAverages = epochWindowAverages['longitudeResidual'] \
                                   *constants.radiansInDegrees

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
outputFilenameBase = filenameBase % subfigures[0] 
plt.savefig(outputFilenameBase + "SemiMajorAxisHistory.pdf", dpi = figureDPI)    
plt.close()  

# # Plot epoch window data for longitude residuals.
# fig = plt.figure()
# axes = fig.gca()
# plt.xlabel("Epoch [Julian years]")
# plt.ylabel("$\Delta L_{Mab}$ [deg]")
# # Plot longitude residual data points in epoch windows.
# plt.plot(plotLongitudeResidualEpoch, plotLongitudeResidual, \
#          marker='.', color='k', linestyle='none')
# # Plot longitude residual averages computed per epoch window.
# plt.plot(plotEpochWindowAveragesEpoch, plotEpochWindowLongitudeAverages,\
#          marker='o', color='0.25', markersize=15.0, linestyle='none')
# # Get y-axis limits.
# ymin, ymax = plt.ylim()
# # Plot shaded areas to indicate epoch window bounds.
# for i in xrange(0,randomWalkCaseData['numberOfEpochWindows']):
#   epochWindowLowerBound = plotCentralEpoch[i] - plotHalfEpochWindowSize
#   epochWindowUpperBound = plotCentralEpoch[i] + plotHalfEpochWindowSize                         
#   axes.fill_between([epochWindowLowerBound, epochWindowUpperBound], \
#                      ymin, ymax, facecolor='0.75', alpha=0.5)
# # Save to file and close figure. 
# plt.ylim(ymin, ymax)
# plt.tight_layout(True)       
# outputFilenameBase = filenameBase % subfigures[1] 
# plt.savefig(outputFilenameBase + "LongitudeResidualEpochWindows.pdf", dpi = figureDPI)    
# plt.close()    

# # Plot longitude residuals in first epoch window.
# fig = plt.figure()
# axes = fig.gca() 
# plt.xlabel("Epoch [Julian years]")
# plt.ylabel("$\Delta L_{Mab}$ [deg]")
# # Extract data for first epoch window.
# firstEpochWindowLowerBound = plotCentralEpoch[0] - plotHalfEpochWindowSize
# firstEpochWindowUpperBound = plotCentralEpoch[0] + plotHalfEpochWindowSize   
# firstWindowIndices = (plotLongitudeResidualEpoch > firstEpochWindowLowerBound) \
#                      & (plotLongitudeResidualEpoch < firstEpochWindowUpperBound)
# for i in xrange(0,len(firstWindowIndices)-1):
#   if firstWindowIndices[i] == False and firstWindowIndices[i+1] == True:
#     firstWindowIndices[i] = True
#     break
# for i in xrange(0,len(firstWindowIndices)-1):
#   if firstWindowIndices[i] == True and firstWindowIndices[i+1] == False:
#     firstWindowIndices[i+1] = True
#     break    
# firstWindowLongitudeResidualEpoch = plotLongitudeResidualEpoch[firstWindowIndices]
# firstWindowLongitudeResidual = plotLongitudeResidual[firstWindowIndices]    
# print firstEpochWindowLowerBound
# print plotLongitudeResidualEpoch
# # Manually create step-function.
# stepEpoch = []
# stepLongitudeResidual = []
# for i in xrange(0,len(firstWindowLongitudeResidualEpoch)-1):
#   stepEpoch.append(firstWindowLongitudeResidualEpoch[i])
#   stepLongitudeResidual.append(firstWindowLongitudeResidual[i])
#   stepEpoch.append(firstWindowLongitudeResidualEpoch[i+1])
#   stepLongitudeResidual.append(firstWindowLongitudeResidual[i])
# # Plot longitude residual step-function.
# plt.step(firstWindowLongitudeResidualEpoch, firstWindowLongitudeResidual, \
#          marker='o', color='k', linestyle='dashed', linewidth = 2.0, where='post')
# # Plot epoch window average.
# plt.plot(plotEpochWindowAveragesEpoch[0], plotEpochWindowLongitudeAverages[0],\
#          marker='o', color='0.25', markersize=15.0, linestyle='none')
# # Plot shaded area for epoch window.
# ymin, ymax = plt.ylim()
# stepEpoch[0] = firstEpochWindowLowerBound
# stepEpoch[-1] = firstEpochWindowUpperBound
# axes.fill_between(stepEpoch, ymin, stepLongitudeResidual, facecolor='0.75', alpha=0.5)
# # Save to file and close figure. 
# plt.ylim(ymin, ymax)
# plt.tight_layout(True)       
# outputFilenameBase = filenameBase % subfigures[2] 
# plt.savefig(outputFilenameBase + "LongitudeResidualFirstEpochWindow.pdf", dpi = figureDPI)    
# plt.close()    

# ###################################################################################################


# # ###################################################################################################
# # # Generate eccentricity figures
# # ###################################################################################################

# # # Plot time-histories of Keplerian elements.
# # output = outputFilename % subfigures[3]
# # fig = plt.figure()
# # plt.xlabel("Epoch [Julian years]")
# # plt.ylabel("$e_{Mab}$ [-]")
# # plt.xlim(xmin=0.0, xmax=50.0)
# # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# # plt.plot(keplerianActionElements['epoch']/constants.julianYear, \
# #          keplerianActionElements['eccentricity'], 'k')
# # plt.plot([centralEpoch / constants.julianYear for centralEpoch in centralEpochs], \
# #          keplerianActionElements['eccentricity'][centralEpochIndices], \
# #          marker='o', color='0.80', linestyle='none', markersize = '8.0')
# # ymin, ymax = plt.ylim()
# # plt.plot([observationPeriodStartEpoch/constants.julianYear, \
# #           observationPeriodStartEpoch/constants.julianYear],
# #          [ymin, ymax],
# #          linestyle='dashed', color='k')
# # plt.plot([(observationPeriodStartEpoch+randomWalkCaseData['observationPeriod'])/constants.julianYear,
# #           (observationPeriodStartEpoch+randomWalkCaseData['observationPeriod'])/constants.julianYear],
# #          [ymin, ymax],
# #          linestyle='dashed', color='k')
# # plt.tight_layout(True)        
# # plt.savefig(output + "EccentricityHistory.pdf", dpi = figureDPI)    
# # plt.close()  

# # output = outputFilename % subfigures[4]
# # fig = plt.figure()
# # axes = fig.gca()
# # plt.xlabel("Epoch [Julian years]")
# # plt.ylabel("$e_{Mab}$ [-]")
# # xmin = (observationPeriodStartEpoch - 0.5*randomWalkCaseData['epochWindowSize'])/constants.julianYear
# # xmax = (observationPeriodStartEpoch + randomWalkCaseData['observationPeriod'] + \
# #         0.5*randomWalkCaseData['epochWindowSize'])/constants.julianYear
# # plt.xlim(xmin, xmax) 
# # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# # eccentricityWindowIndices = []
# # for i in xrange(0, randomWalkCaseData['numberOfEpochWindows']):
# #   eccentricityWindowIndices.append( \
# #     (keplerianActionElements['epoch'] > (observationPeriodStartEpoch \
# #                                          - 0.5*randomWalkCaseData['epochWindowSize'] \
# #                                          + i*epochWindowSpacing)) \
# #     & (keplerianActionElements['epoch'] < (observationPeriodStartEpoch \
# #                                            + 0.5*randomWalkCaseData['epochWindowSize'] \
# #                                            + i*epochWindowSpacing)))

# # for i in xrange(0, len(eccentricityWindowIndices)):
# #   eccentricityWindowEpochs = keplerianActionElements['epoch'][eccentricityWindowIndices[i]]
# #   eccentricityWindows = keplerianActionElements['eccentricity'][eccentricityWindowIndices[i]]
# #   plt.plot(eccentricityWindowEpochs/constants.julianYear, eccentricityWindows, \
# #            marker='.', color='k', linestyle='none')

# # plt.plot([centralEpoch / constants.julianYear for centralEpoch in centralEpochs], \
# #          epochWindowAverages['eccentricity'], \
# #          marker='o', color='0.25', markersize=15.0, linestyle='none')
# # # ymin, ymax = plt.ylim()
# # # for i in xrange(0,randomWalkCaseData['numberOfEpochWindows']):
# # #     epochData = longitudeResiduals['epoch'][\
# # #                         (longitudeResiduals['epoch'] \
# # #                           > -randomWalkCaseData['epochWindowSize']*0.5 + i*epochWindowSpacing) \
# # #                         & (longitudeResiduals['epoch'] \
# # #                             < randomWalkCaseData['epochWindowSize']*0.5 + i*epochWindowSpacing)]
# # #     axes.fill_between(epochData/constants.julianYear,ymin, ymax, facecolor='0.75', alpha=0.5)
# # # plt.ylim(ymin, ymax)
# # plt.tight_layout(True)        
# # plt.savefig(output + "EccentricityEpochWindows.pdf", dpi = figureDPI)    
# # plt.close()    


# # output = outputFilename % subfigures[5]
# # fig = plt.figure()
# # axes = fig.gca() 
# # plt.xlabel("Epoch [Julian years]")
# # plt.ylabel("$e_{Mab}$ [-]")
# # xmin = (observationPeriodStartEpoch - 0.5*randomWalkCaseData['epochWindowSize'])/constants.julianYear
# # xmax = (observationPeriodStartEpoch + 0.5*randomWalkCaseData['epochWindowSize'])/constants.julianYear
# # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# # plt.xlim(xmin*0.995, xmax*1.005)
# # eccentricityFirstWindowIndices = (keplerianActionElements['epoch']/constants.julianYear > xmin) \
# #                                  & (keplerianActionElements['epoch']/constants.julianYear < xmax)

# # for i in xrange(0,len(keplerianActionElements['epoch'])-1):
# #   if (eccentricityFirstWindowIndices[i] == False) and (eccentricityFirstWindowIndices[i+1] == True):
# #      eccentricityFirstWindowIndices[i] = True
# #      break 

# # for i in xrange(0,len(keplerianActionElements['epoch'])-1):
# #   if (eccentricityFirstWindowIndices[i] == True) and (eccentricityFirstWindowIndices[i+1] == False):
# #      eccentricityFirstWindowIndices[i+1] = True
# #      break      

# # eccentricityFirstWindowEpochs = keplerianActionElements['epoch'][ \
# #                                   eccentricityFirstWindowIndices]/constants.julianYear
# # eccentricityFirstWindow = keplerianActionElements['eccentricity'][eccentricityFirstWindowIndices]

# # stepEpoch = []
# # stepEccentricities = []
# # for i in xrange(0,len(longitudeResidualsFirstWindowEpochs)-1):
# #   stepEpoch.append(eccentricityFirstWindowEpochs[i])
# #   stepEccentricities.append(eccentricityFirstWindow[i])
# #   stepEpoch.append(eccentricityFirstWindowEpochs[i+1])
# #   stepEccentricities.append(eccentricityFirstWindow[i])

# # plt.step(eccentricityFirstWindowEpochs, eccentricityFirstWindow, \
# #          marker='o', color='k', linestyle='dashed', linewidth = 2.0, where='post')
# # plt.plot(centralEpochs[0]/constants.julianYear, \
# #          epochWindowAverages['eccentricity'][0],\
# #          marker='o', color='0.25', markersize=15.0, linestyle='none')
# # ymin, ymax = plt.ylim()
# # # axes.fill_between(stepEpoch, ymin, stepEccentricities, facecolor='0.75', alpha=0.5)
# # plt.ylim(ymin, ymax)
# # plt.tight_layout(True)        
# # plt.savefig(output + "EccentricityFirstEpochWindow.pdf", dpi = figureDPI)    
# # plt.close()  

# # ###################################################################################################


# # ###################################################################################################
# # # Generate inclination figures
# # ###################################################################################################

# # # Set output path and case-prefix for files generated.
# # outputFilename  = outputPath + "/figure16%s_monteCarloRun" + str(monteCarloRunId)

# # # Plot time-histories of Keplerian elements.
# # output = outputFilename % subfigures[6]
# # fig = plt.figure()
# # plt.xlabel("Epoch [Julian years]")
# # plt.ylabel("$i_{Mab}$ [deg]")
# # plt.xlim(xmin=0.0, xmax=50.0)
# # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# # plt.plot(keplerianActionElements['epoch']/constants.julianYear, \
# #          keplerianActionElements['inclination']*constants.radiansInDegrees, 'k')
# # plt.plot([centralEpoch / constants.julianYear for centralEpoch in centralEpochs], \
# #          keplerianActionElements['inclination'][centralEpochIndices]*constants.radiansInDegrees, \
# #          marker='o', color='0.75', linestyle='none')
# # ymin, ymax = plt.ylim()
# # plt.plot([observationPeriodStartEpoch/constants.julianYear, \
# #           observationPeriodStartEpoch/constants.julianYear],
# #          [ymin, ymax],
# #          linestyle='dashed', color='k')
# # plt.plot([(observationPeriodStartEpoch+randomWalkCaseData['observationPeriod'])/constants.julianYear,
# #           (observationPeriodStartEpoch+randomWalkCaseData['observationPeriod'])/constants.julianYear],
# #          [ymin, ymax],
# #          linestyle='dashed', color='k')
# # plt.tight_layout(True)        
# # plt.savefig(output + "InclinationHistory.pdf", dpi = figureDPI)    
# # plt.close()  

# # output = outputFilename % subfigures[7]
# # fig = plt.figure()
# # axes = fig.gca()
# # plt.xlabel("Epoch [Julian years]")
# # plt.ylabel("$\Delta i_{Mab}$ [deg]")
# # plt.xlim(xmin = -randomWalkCaseData['epochWindowSize']/(constants.julianYear*2.0), \
# #          xmax = randomWalkCaseData['observationPeriod']/constants.julianYear \
# #                 + randomWalkCaseData['epochWindowSize']/(constants.julianYear*2.0))
# # plt.plot(keplerianActionElements['epoch']/constants.julianYear, \
# #          keplerianActionElements['inclination']*constants.radiansInDegrees, \
# #          marker='.', color='k', linestyle='none')
# # plt.plot(epochWindowAverages['epoch']/constants.julianYear, \
# #          epochWindowAverages['inclination']*constants.radiansInDegrees,\
# #          marker='o', color='0.75', markersize=15.0, linestyle='none')
# # ymin, ymax = plt.ylim()
# # for i in xrange(0,randomWalkCaseData['numberOfEpochWindows']):
# #     epochData = keplerianActionElements['epoch'][\
# #                         (keplerianActionElements['epoch'] \
# #                             > (-randomWalkCaseData['epochWindowSize']*0.5 \
# #                                + i*epochWindowSpacing)/constants.julianYear) \
# #                         & (keplerianActionElements['epoch'] \
# #                             < (randomWalkCaseData['epochWindowSize']*0.5 \
# #                                + i*epochWindowSpacing)/constants.julianYear)]
# #     axes.fill_between(epochData,ymin, ymax, facecolor='0.75', alpha=0.5)
# # plt.ylim(ymin, ymax)
# # plt.tight_layout(True)        
# # plt.savefig(output + "InclinationEpochWindows.pdf", dpi = figureDPI)    
# # plt.close()     

# # ###################################################################################################


###################################################################################################
# Finalize timer and print elapsed time.
###################################################################################################

# Finalize timer.
endTime = time.time()

# Print elapsed time for script [s].
print "This script took " + str(endTime - startTime) + "s"

###################################################################################################
