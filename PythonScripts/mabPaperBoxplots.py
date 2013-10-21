'''   
Copyright (c) 2010-2013, Delft University of Technology
Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
All rights reserved.
See http://bit.ly/12SHPLR for license details.

This script can be used to generate boxplot's of Mab's orbital element variations due to a 
perturber ring.
'''


################################################################################################### 
# Set up input deck
###################################################################################################

# Set absolute path to SQLite database with simulation data.
databasePath        = "/Users/kartikkumar/Documents/Education/PhD/Simulations/Tudat/" + \
                      "Workspace/tudatApplications/stochasticMigration/" + \
                      "stochasticMigrationResults.sqlite"

# Set random walk case names.
nominal             = "circular_equatorial_nominal"
light               = "circular_equatorial_light"
heavy               = "circular_equatorial_heavy"
sparse              = "circular_equatorial_sparse"
dense               = "circular_equatorial_dense"

# Set absolute path to output directory.
outputPath          = "/Users/kartikkumar/Desktop"

# Set figure dpi.
figureDPI           = 600

# Set flag whether outliers should be shown in boxplots.
showOutliers        = False

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

# Set query to retrieve output data for a given case, with the input data used as a look-up 
# table.
outputQuery = "SELECT * FROM random_walk_output \
               WHERE random_walk_output.monteCarloRunId \
               IN (SELECT \"monteCarloRunId\" FROM random_walk_input \
                   WHERE random_walk_input.randomWalkCaseId \
                   == (SELECT \"caseId\" FROM random_walk_case \
                       WHERE \"caseName\" == \"%s\"));"

# Connect to SQLite database.
database = sqlite3.connect(databasePath)

# Enter database and retrieve random walk case data.
with database:
    # Set cursor to scan through database.
    cursor = database.cursor()

    # Get output data for nominal case.
    cursor.execute(outputQuery % nominal)

    # Get column name list for all cases and store in a list.
    outputColumnNameList = [column_name[0] for column_name in cursor.description] 
    outputTypeList = [] 
    for name in outputColumnNameList:
        outputTypeList.append((name, '<f8'))

    rawOutputNominal = cursor.fetchall()
    outputNominal = numpy.array(rawOutputNominal, dtype=(outputTypeList))

    # Get output data for light case.
    cursor.execute(outputQuery % light)
    rawOutputLight = cursor.fetchall()
    outputLight = numpy.array(rawOutputLight, dtype=(outputTypeList))

    # Get output data for heavy case.
    cursor.execute(outputQuery % heavy)
    rawOutputHeavy = cursor.fetchall()
    outputHeavy = numpy.array(rawOutputHeavy, dtype=(outputTypeList))

    # Get output data for sparse case.
    cursor.execute(outputQuery % sparse)
    rawOutputSparse = cursor.fetchall()
    outputSparse = numpy.array(rawOutputSparse, dtype=(outputTypeList))

    # Get output data for dense case.
    cursor.execute(outputQuery % dense)
    rawOutputDense = cursor.fetchall()
    outputDense = numpy.array(rawOutputDense, dtype=(outputTypeList))

    # Get random walk case ID associated with output data.
    cursor.execute("SELECT testParticleCaseId FROM random_walk_case \
                    WHERE \"caseName\" == \"" + nominal + "\"")
    testParticleCaseId = cursor.fetchall()[0][0]

# Close database connection if it is open.    
if database:
    # Set cursor to scan through database and execute queries.
    database.close()             

###################################################################################################


###################################################################################################
# Plot boxplots of random walk output data as a function of perturber ring mass and density
###################################################################################################

# Set output path and case-prefix for files generated.
outputPathAndCasePrefix = outputPath + "/testParticleCase" + str(testParticleCaseId) + "_"

if showOutliers:
    outliers = '.'
else:
    outliers = ''

# Store plotting data.
longitudeResidualVsRingMass = []
longitudeResidualVsRingMass.append(numpy.rad2deg(outputLight['maximumLongitudeResidualChange']))
longitudeResidualVsRingMass.append(numpy.rad2deg(outputNominal['maximumLongitudeResidualChange']))
longitudeResidualVsRingMass.append(numpy.rad2deg(outputHeavy['maximumLongitudeResidualChange']))

eccentricityVsRingMass = []
eccentricityVsRingMass.append(outputLight['maximumEccentricityChange'])
eccentricityVsRingMass.append(outputNominal['maximumEccentricityChange'])
eccentricityVsRingMass.append(outputHeavy['maximumEccentricityChange'])

inclinationVsRingMass = []
inclinationVsRingMass.append(numpy.rad2deg(outputLight['maximumInclinationChange']))
inclinationVsRingMass.append(numpy.rad2deg(outputNominal['maximumInclinationChange']))
inclinationVsRingMass.append(numpy.rad2deg(outputHeavy['maximumInclinationChange']))

longitudeResidualVsRingDensity = []
longitudeResidualVsRingDensity.append(numpy.rad2deg(outputSparse['maximumLongitudeResidualChange']))
longitudeResidualVsRingDensity.append(numpy.rad2deg(outputNominal['maximumLongitudeResidualChange']))
longitudeResidualVsRingDensity.append(numpy.rad2deg(outputDense['maximumLongitudeResidualChange']))

eccentricityVsRingDensity = []
eccentricityVsRingDensity.append(outputSparse['maximumEccentricityChange'])
eccentricityVsRingDensity.append(outputNominal['maximumEccentricityChange'])
eccentricityVsRingDensity.append(outputDense['maximumEccentricityChange'])

inclinationVsRingDensity = []
inclinationVsRingDensity.append(numpy.rad2deg(outputSparse['maximumInclinationChange']))
inclinationVsRingDensity.append(numpy.rad2deg(outputNominal['maximumInclinationChange']))
inclinationVsRingDensity.append(numpy.rad2deg(outputDense['maximumInclinationChange']))

# Generate figure with subplots.
fig, ((axes1, axes2, axes3), (axes4, axes5, axes6)) = plt.subplots(nrows=2, ncols=3)

# Plot maximum longitude residual change [deg] vs. perturber ring mass [M_Mab].
plotOut = axes1.boxplot(longitudeResidualVsRingMass,sym=outliers)
plt.setp(plotOut['boxes'], color='black')
plt.setp(plotOut['whiskers'], color='black')
plt.setp(plotOut['medians'], color='black')
plt.setp(plotOut['fliers'], color='black')
axes1.xaxis.set_ticks([1, 2, 3])
axes1.set_xticklabels(['$\\frac{1}{3}$', '1', '3'])
axes1.set_xlabel('$M_{perturber}$ [$M_{Mab}$]')
axes1.set_ylabel('$\Delta L_{Mab,max}$ [deg]')

# Plot maximum eccentricity change [-] vs. perturber ring mass [M_Mab].
plotOut = axes2.boxplot(eccentricityVsRingMass,sym=outliers)
plt.setp(plotOut['boxes'], color='black')
plt.setp(plotOut['whiskers'], color='black')
plt.setp(plotOut['medians'], color='black')
plt.setp(plotOut['fliers'], color='black')
axes2.xaxis.set_ticks([1, 2, 3])
axes2.set_xticklabels(['$\\frac{1}{3}$', '1', '3'])
axes2.set_xlabel('$M_{perturber}$ [$M_{Mab}$]')
axes2.set_ylabel('$\Delta e_{Mab,max}$ [-]')
axes2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Plot maximum inclination change [deg] vs. perturber ring mass [M_Mab].
plotOut = axes3.boxplot(inclinationVsRingMass,sym=outliers)
plt.setp(plotOut['boxes'], color='black')
plt.setp(plotOut['whiskers'], color='black')
plt.setp(plotOut['medians'], color='black')
plt.setp(plotOut['fliers'], color='black')
axes3.xaxis.set_ticks([1, 2, 3])
axes3.set_xticklabels(['$\\frac{1}{3}$', '1', '3'])
axes3.set_xlabel('$M_{perturber}$ [$M_{Mab}$]')
axes3.set_ylabel('$\Delta i_{Mab,max}$ [deg]')
axes3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plotOut = axes4.boxplot(longitudeResidualVsRingDensity,sym=outliers)
plt.setp(plotOut['boxes'], color='black')
plt.setp(plotOut['whiskers'], color='black')
plt.setp(plotOut['medians'], color='black')
plt.setp(plotOut['fliers'], color='black')
axes4.xaxis.set_ticks([1, 2, 3])
axes4.set_xticklabels(['1', '$\\frac{10}{3}$', '10'])
axes4.set_xlabel(r'$\rho_{perturber}$ [N/$R_{Hill,Mab}$]')
axes4.set_ylabel('$\Delta L_{Mab,max}$ [deg]')

plotOut = axes5.boxplot(eccentricityVsRingDensity,sym=outliers)
plt.setp(plotOut['boxes'], color='black')
plt.setp(plotOut['whiskers'], color='black')
plt.setp(plotOut['medians'], color='black')
plt.setp(plotOut['fliers'], color='black')
axes5.xaxis.set_ticks([1, 2, 3])
axes5.set_xticklabels(['1', '$\\frac{10}{3}$', '10'])
axes5.set_xlabel(r'$\rho_{perturber}$ [N/$R_{Hill,Mab}$]')
axes5.set_ylabel('$\Delta e_{Mab,max}$ [-]')
axes5.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plotOut = axes6.boxplot(inclinationVsRingDensity,sym=outliers)
plt.setp(plotOut['boxes'], color='black')
plt.setp(plotOut['whiskers'], color='black')
plt.setp(plotOut['medians'], color='black')
plt.setp(plotOut['fliers'], color='black')
axes6.xaxis.set_ticks([1, 2, 3])
axes6.set_xticklabels(['1', '$\\frac{10}{3}$', '10'])
axes6.set_xlabel(r'$\rho_{perturber}$ [N/$R_{Hill,Mab}$]')
axes6.set_ylabel('$\Delta i_{Mab,max}$ [deg]')
axes6.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

fig.set_tight_layout(True)

plt.savefig(outputPathAndCasePrefix + "maximumChangesVsRingMass.pdf", \
            dpi = figureDPI, orientation='landscape')
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