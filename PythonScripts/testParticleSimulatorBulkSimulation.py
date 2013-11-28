'''   
Copyright (c) 2010-2013, Delft University of Technology
Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
All rights reserved.
See http://bit.ly/12SHPLR for license details.

This script can be used to generate plots based on the bulk output from the test particle simulator
for a specified case.
'''

################################################################################################### 
# Set up input deck
###################################################################################################

# Set absolute path to SQLite database with simulation data.
databasePath    = "/Users/kartikkumar/Documents/Education/PhD/Simulations/" \
                      + "Tudat/Workspace/tudatApplications/stochasticMigration/" \
                      + "stochasticMigrationResults.sqlite"

# Set case name.
caseName        = "circular_equatorial"

# Set absolute path to output directory.
outputPath      = "/Users/kartikkumar/Desktop"

# Set cut-off for big kicks included in plots.
bigKicksCutOff  = 5000

# Show figures in interactive matplotlib window?
isShowFigures   = 1

# Set figure dpi.
figureDPI       = 600

###################################################################################################


'''
    
                            DO NOT EDIT PARAMETERS BEYOND THIS POINT!!!

'''


################################################################################################### 
# Import Python packages and modules
###################################################################################################

# Import necessary external packages.
import math
import matplotlib.pyplot as plt
import numpy
import numpy.lib.recfunctions
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

    # Select all the input data associated with the case ID.
    cursor.execute("SELECT * FROM test_particle_input WHERE testParticleCaseId == " + str(caseId) + ";")
    rawInputData = cursor.fetchall()
    inputDataColumnNameList = [column_name[0] for column_name in cursor.description]  

    # Select all the output data associated with the case ID.
    cursor.execute("SELECT * FROM test_particle_kicks WHERE test_particle_kicks.testParticleSimulationId \
                    IN (SELECT simulationId FROM test_particle_input \
                        WHERE testParticleCaseId == " + str(caseId) + ");")
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
# Compute derived parameters
###################################################################################################

# Compute mass of perturbed body [kg].
perturbedBodyMass = 4.0/3.0 * math.pi * caseData['perturbedBodyRadius']**3.0 \
                   * caseData['perturbedBodyBulkDensity']

# Compute perturbed body gravitational parameter [m^3 s^-2].
perturbedBodyGravitationalParameter = constants.universalGravitationalConstant * perturbedBodyMass

# Compute Hill sphere radius [m].
hillRadius = caseData['perturbedBodySemiMajorAxisAtT0'] \
            * (perturbedBodyGravitationalParameter \
               / (3.0 * caseData['centralBodyGravitationalParameter']))**(1.0 / 3.0)

# Compute Hill velocity using Hill radius [m s^-1].
hillVelocity = hillRadius * math.sqrt(caseData['centralBodyGravitationalParameter'] \
              / caseData['perturbedBodySemiMajorAxisAtT0']**3.0)      

# Compute semi-major axis distribution limit in Hill radii [Hill radii].
semiMajorAxisDistributionLimitInHillRadii = caseData['semiMajorAxisDistributionLimit'] / hillRadius       

###################################################################################################


###################################################################################################
# Execute data processing
###################################################################################################

# Declare structured array containing data using to generate plots.   
plottingData = numpy.array(numpy.zeros(len(outputData)),
                          dtype = numpy.dtype([('semiMajorAxisKick', '<f8'),
                                               ('eccentricityKick', '<f8'),
                                               ('inclinationKick', '<f8'),
                                               ('semiMajorAxisKickMagnitude', '<f8'),
                                               ('eccentricityKickMagnitude', '<f8'),
                                               ('inclinationKickMagnitude', '<f8'),
                                               ('simulationId', '<i4'),
                                               ('initialSemiMajorAxis', '<f8'),
                                               ('preConjunctionSemiMajorAxis', '<f8'),
                                               ('preConjunctionEccentricity', '<f8'),
                                               ('preConjunctionInclination', '<f8')]))

# Store semi-major axis kicks [km].
plottingData['semiMajorAxisKick'] = (outputData['postConjunctionSemiMajorAxis'] \
                                    - outputData['preConjunctionSemiMajorAxis']) \
                                   * constants.meterInKilometers
    
# Store eccentricity kicks [-].                                    
plottingData['eccentricityKick'] = (outputData['postConjunctionEccentricity'] \
                                   - outputData['preConjunctionEccentricity'])

# Store inclination kicks [deg].
plottingData['inclinationKick'] = (outputData['postConjunctionInclination'] \
                                  - outputData['preConjunctionInclination']) \
                                 * constants.radiansInDegrees                                

# Store (magnitude) semi-major axis kicks [km].
plottingData['semiMajorAxisKickMagnitude'] = numpy.abs(plottingData['semiMajorAxisKick'])
    
# Store (magnitude) eccentricity kicks [-].                                    
plottingData['eccentricityKickMagnitude'] = numpy.abs(plottingData['eccentricityKick']) 

# Store (magnitude) inclination kicks [deg].
plottingData['inclinationKickMagnitude'] = numpy.abs(plottingData['inclinationKick'])    

# Create a mapping from simulation IDs to array indices in input data array.
simulationIdMapping = dict(zip(inputData['simulationId'], range(len(inputData))))

# Store initial semi-major axes based on simulation ID mapping.
plottingData['initialSemiMajorAxis'] \
    = (numpy.array([inputData['semiMajorAxis'][simulationIdMapping[key]] \
                   for key in outputData['testParticleSimulationId']])
       - caseData['perturbedBodySemiMajorAxisAtT0']) * constants.meterInKilometers

# Store simulation numbers.
plottingData['simulationId'] = outputData['testParticleSimulationId']

# Sort plotting data in descending order based on semi-major axis kick (magnitude).
semiMajorAxisKickSortedData = plottingData[numpy.argsort(plottingData, \
                                           order='semiMajorAxisKickMagnitude')][::-1]

# Sort plotting data in descending order based on eccentricity kick (magnitude).
eccentricitySortedData = plottingData[numpy.argsort(plottingData, \
                                      order='eccentricityKickMagnitude')][::-1]                                                                

# Store pre-conjunction semi-major axes [km].                                               
plottingData['preConjunctionSemiMajorAxis'] = (outputData['preConjunctionSemiMajorAxis'] \
                                               - caseData['perturbedBodySemiMajorAxisAtT0']) \
                                              * constants.meterInKilometers

# Store pre-conjunction eccentricity [-].
plottingData['preConjunctionEccentricity'] = outputData['preConjunctionEccentricity'] 

# Store pre-conjunction eccentricity [km].
plottingData['preConjunctionInclination'] = outputData['preConjunctionInclination'] \
                                            * constants.radiansInDegrees
    
# # Extract indices for horseshoe orbits.
# horseshoeOrbitIndices = (numpy.abs(plottingData['initialSemiMajorAxis']) < 40.0) \
#                         & (plottingData['semiMajorAxisKickMagnitude'] > 1.0)
                        
# # Extract unique simulation numbers for horseshoe orbits.
# horseshoeOrbitSimulationNumbers \
#     = numpy.unique(plottingData['simulationId'][horseshoeOrbitIndices])     

###################################################################################################


###################################################################################################
# Generate figures
###################################################################################################

# Set output path and case-prefix for files generated.
outputPathAndCasePrefix = outputPath + "/case" + str(caseData['caseId']) + "_"

# Compute semi-major axis plot limit in Hill radii.
semiMajorAxisPlotLimitInHillRadii = int(math.floor(semiMajorAxisDistributionLimitInHillRadii \
                                                   / 5.0)) * 5

# Compute semi-major axis plot tick marks in Hill radii.
semiMajorAxisPlotTickMarks = numpy.arange(-semiMajorAxisPlotLimitInHillRadii, \
                                          semiMajorAxisPlotLimitInHillRadii+5, 5)

# Plot semi-major axis kicks (magnitude) vs. initial semi-major axes with respect to perturbed 
# body.
figure = plt.figure()  
axis1 = figure.add_subplot(111)
axis1.plot(plottingData['initialSemiMajorAxis'], plottingData['semiMajorAxisKickMagnitude'], \
           '.k', rasterized = True)
# axis1.plot(plottingData['initialSemiMajorAxis'][horseshoeOrbitIndices], \
#            plottingData['semiMajorAxisKickMagnitude'][horseshoeOrbitIndices], \
#            marker = '.', color = '0.5', linestyle = 'None')
axis1.set_yscale('log')
axis1.set_xlabel("Initial semi-major axis relative to Mab [km]")
axis1.set_ylabel("Semi-major axis kick magnitude [km]")
xlimits = axis1.set_xlim(
    -caseData['semiMajorAxisDistributionLimit'] * constants.meterInKilometers, \
    caseData['semiMajorAxisDistributionLimit'] * constants.meterInKilometers)
axis1.set_ylim(ymin = 0.0)
axis2 = axis1.twiny()
axis2.set_xlim(xlimits[0],xlimits[1])
axis2.set_xticks(semiMajorAxisPlotTickMarks * hillRadius * constants.meterInKilometers)
axis2.set_xticklabels(semiMajorAxisPlotTickMarks)
axis2.set_xlabel("Initial semi-major axis relative to Mab [Hill radii]")
plt.savefig(outputPathAndCasePrefix + "semiMajorAxisKickMagnitudeVsInitialSemiMajorAxes.pdf", \
            dpi = figureDPI)
plt.close()

# Plot eccentricity vs. semi-major axis kicks (magnitudes) (big kicks only, cut-off set by 
# bigKicksCutOff). 
plt.figure() 
figureAxes = plt.gca()        
plt.plot(semiMajorAxisKickSortedData['semiMajorAxisKickMagnitude'][:bigKicksCutOff], \
         semiMajorAxisKickSortedData['eccentricityKickMagnitude'][:bigKicksCutOff], '.k')
plt.plot(eccentricitySortedData['semiMajorAxisKickMagnitude'][:bigKicksCutOff], \
         eccentricitySortedData['eccentricityKickMagnitude'][:bigKicksCutOff], '.k')
figureAxes.yaxis.major.formatter.set_powerlimits((0,0)) 
plt.xlabel("Semi-major axis kick magnitude [km]")
plt.ylabel("Eccentricity kick magnitude [-]")
plt.xlim(xmin = 0.0)
plt.ylim(ymin = 0.0)
plt.savefig(outputPathAndCasePrefix + "eccentricityVsSemiMajorAxisBigKicksMagnitude.pdf", \
            dpi = figureDPI)    
plt.close()

# Plot semi-major axis kicks (magnitude) vs. initial semi-major axes with respect to Mab 
# (big kicks interactive plot that can be used to associate output data with simulation IDs). 
if (isShowFigures == 1):
    labels = ['{0}'.format(simulationId) for simulationId in \
              semiMajorAxisKickSortedData['simulationId'][:bigKicksCutOff]]  
    plt.figure()  
    plt.plot(semiMajorAxisKickSortedData['initialSemiMajorAxis'][:bigKicksCutOff], \
             semiMajorAxisKickSortedData['semiMajorAxisKickMagnitude'][:bigKicksCutOff], '.k')
    plt.plot(eccentricitySortedData['initialSemiMajorAxis'][:bigKicksCutOff], \
             eccentricitySortedData['semiMajorAxisKickMagnitude'][:bigKicksCutOff], '.k')
    for label, initialSemiMajorAxis, semiMajorAxisKickMagnitude in \
        zip(labels, semiMajorAxisKickSortedData['initialSemiMajorAxis'], \
            semiMajorAxisKickSortedData['semiMajorAxisKickMagnitude']):
        plt.annotate(label, xy = (initialSemiMajorAxis, semiMajorAxisKickMagnitude), \
                     xytext = (-5, 5), textcoords = 'offset points', ha = 'right', \
                     va = 'bottom', size = 'small')
    plt.yscale('log')
    plt.xlabel("Initial semi-major axis with respect to Mab [km]")
    plt.ylabel("Magnitude of semi-major axis kick [km]")
    plt.xlim(xmin = -caseData['semiMajorAxisDistributionLimit'] * constants.meterInKilometers, \
              xmax =caseData['semiMajorAxisDistributionLimit'] * constants.meterInKilometers)
    plt.ylim(ymin = 0.0)
    plt.show()

''' 
Plot semi-major axis kicks vs. pre-conjunction semi-major axes.
'''

# Plot semi-major axis kicks (magnitude) vs. pre-conjunction semi-major axes with respect to 
# perturbed body.   
figure = plt.figure() 
axis1 = figure.add_subplot(111)
plt.plot(plottingData['preConjunctionSemiMajorAxis'], \
         plottingData['semiMajorAxisKickMagnitude'], '.k', rasterized = True)
# axis1.plot(plottingData['preconjunctionSemiMajorAxis'][horseshoeOrbitIndices], \
#            plottingData['semiMajorAxisKickMagnitude'][horseshoeOrbitIndices], \
#            marker = '.', color = '0.5', linestyle = 'None')
axis1.set_yscale('log')
axis1.set_xlabel("Pre-conjunction semi-major axis relative to Mab [km]")
axis1.set_ylabel("Semi-major axis kick magnitude [km]")
xlimits = axis1.set_xlim(
    -caseData['semiMajorAxisDistributionLimit'] * constants.meterInKilometers, \
    caseData['semiMajorAxisDistributionLimit'] * constants.meterInKilometers)
axis1.set_ylim(ymin = 0.0)
axis2 = axis1.twiny()
axis2.set_xlim(xlimits[0],xlimits[1])
axis2.set_xticks(semiMajorAxisPlotTickMarks * hillRadius * constants.meterInKilometers)
axis2.set_xticklabels(semiMajorAxisPlotTickMarks)
axis2.set_xlabel("Pre-conjunction semi-major axis relative to Mab [Hill radii]")
plt.savefig(outputPathAndCasePrefix + "semiMajorAxisKickVsPreConjunctionSemiMajorAxes.pdf", \
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
