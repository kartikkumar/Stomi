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
databasePath        = "/Users/kartikkumar/Documents/University/PhD/Simulations/Tudat/" + \
                      "Workspace/tudatApplications/stochasticMigration/" + \
                      "stochasticMigrationResults.sqlite"

# Set random walk case names.
nominal             = "circular_equatorial_nominal"
light               = "circular_equatorial_short"
heavy               = "circular_equatorial_long"
sparse              = "circular_equatorial_long"
dense               = "circular_equatorial_long"

# Set absolute path to output directory.
outputPath          = "/Users/kartikkumar/Desktop/results"

# Set figure dpi.
figureDPI           = 600

# Set flag whether outliers should be shown in boxplots.
showOutliers        = True

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

    # Get test particle case ID associated with output data.
    cursor.execute("SELECT testParticleCaseId FROM random_walk_case \
                    WHERE \"caseName\" == \"" + nominal + "\"")
    testParticleCaseId = cursor.fetchall()[0][0]

# Close database connection if it is open.    
if database:
    # Set cursor to scan through database and execute queries.
    database.close()             

###################################################################################################


###################################################################################################
# Plot boxplots of random walk output data as a function of perturber ring mass
###################################################################################################

# Set output path and case-prefix for files generated.
outputPathAndCasePrefix = outputPath + "/case" + str(testParticleCaseId) + "_"

if showOutliers:
    outliers = '.'
else:
    outliers = ''

# Plot maximum longitude residual change [deg] vs. perturber ring mass [M_Mab].
longitudeResidualVsRingMass = []
longitudeResidualVsRingMass.append(numpy.rad2deg(outputLight['maximumLongitudeResidualChange']))
longitudeResidualVsRingMass.append(numpy.rad2deg(outputNominal['maximumLongitudeResidualChange']))
longitudeResidualVsRingMass.append(numpy.rad2deg(outputHeavy['maximumLongitudeResidualChange']))

fig = plt.figure()
axes = fig.gca()
plotOut = plt.boxplot(longitudeResidualVsRingMass,sym=outliers)
plt.setp(plotOut['boxes'], color='black')
plt.setp(plotOut['whiskers'], color='black')
plt.setp(plotOut['medians'], color='black')
plt.setp(plotOut['fliers'], color='black')
plt.xticks([1, 2, 3], ['$\\frac{1}{3}$', '1', '3'])
plt.xlabel('Perturber ring mass [$M_{Mab}$]')
plt.ylabel('Maximum longitude residual change [deg]')

for line in plotOut['medians']:
    # Get position data for median line.
    medianX, medianY = line.get_xydata()[1]

    # Overlay median value.
    axes.annotate( '%1.1f' % medianY, xy=(medianX+0.11, medianY), xycoords='data',\
                  horizontalalignment='center', verticalalignment='center')

plt.savefig(outputPathAndCasePrefix + "maximumLongitudeResidualVsRingMass.pdf", \
            dpi = figureDPI)
plt.close()

# Plot maximum eccentricity change [-] vs. perturber ring mass [M_Mab].
eccentricityVsRingMass = []
eccentricityVsRingMass.append(outputLight['maximumEccentricityChange'])
eccentricityVsRingMass.append(outputNominal['maximumEccentricityChange'])
eccentricityVsRingMass.append(outputHeavy['maximumEccentricityChange'])

fig = plt.figure()
axes = fig.gca()
plotOut = plt.boxplot(eccentricityVsRingMass,sym=outliers)
plt.setp(plotOut['boxes'], color='black')
plt.setp(plotOut['whiskers'], color='black')
plt.setp(plotOut['medians'], color='black')
plt.setp(plotOut['fliers'], color='black')
plt.xticks([1, 2, 3], ['$\\frac{1}{3}$', '1', '3'])
plt.xlabel('Perturber ring mass [$M_{Mab}$]')
plt.ylabel('Maximum eccentricity change [-]')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

for line in plotOut['medians']:
    # Get position data for median line.
    medianX, medianY = line.get_xydata()[1]

    # Overlay median value.
    axes.annotate( '%1.1e' % medianY, xy=(medianX+0.18, medianY), xycoords='data',\
                  horizontalalignment='center', verticalalignment='center')

plt.savefig(outputPathAndCasePrefix + "maximumEccentricityVsRingMass.pdf", \
            dpi = figureDPI)
plt.close()

# Plot maximum inclination change [deg] vs. perturber ring mass [M_Mab].
inclinationVsRingMass = []
inclinationVsRingMass.append(numpy.rad2deg(outputLight['maximumInclinationChange']))
inclinationVsRingMass.append(numpy.rad2deg(outputNominal['maximumInclinationChange']))
inclinationVsRingMass.append(numpy.rad2deg(outputHeavy['maximumInclinationChange']))

fig = plt.figure()
axes = fig.gca()
plotOut = plt.boxplot(inclinationVsRingMass,sym=outliers)
plt.setp(plotOut['boxes'], color='black')
plt.setp(plotOut['whiskers'], color='black')
plt.setp(plotOut['medians'], color='black')
plt.setp(plotOut['fliers'], color='black')
plt.xticks([1, 2, 3], ['$\\frac{1}{3}$', '1', '3'])
plt.xlabel('Perturber ring mass [$M_{Mab}$]')
plt.ylabel('Maximum inclination change [deg]')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

for line in plotOut['medians']:
    # Get position data for median line.
    medianX, medianY = line.get_xydata()[1]

    # Overlay median value.
    axes.annotate( '%1.1e' % medianY, xy=(medianX+0.18, medianY), xycoords='data',\
                  horizontalalignment='center', verticalalignment='center')

plt.savefig(outputPathAndCasePrefix + "maximumInclinationVsRingMass.pdf", \
            dpi = figureDPI)
plt.close()

###################################################################################################


###################################################################################################
# Plot boxplots of random walk output data as a function of perturber density
###################################################################################################

# Plot maximum longitude residual change [deg] vs. perturber density [#/R_Hill,Mab].
longitudeResidualVsRingDensity = []
longitudeResidualVsRingDensity.append(numpy.rad2deg(outputSparse['maximumLongitudeResidualChange']))
longitudeResidualVsRingDensity.append(numpy.rad2deg(outputNominal['maximumLongitudeResidualChange']))
longitudeResidualVsRingDensity.append(numpy.rad2deg(outputDense['maximumLongitudeResidualChange']))

fig = plt.figure()
axes = fig.gca()
plotOut = plt.boxplot(longitudeResidualVsRingDensity,sym=outliers)
plt.setp(plotOut['boxes'], color='black')
plt.setp(plotOut['whiskers'], color='black')
plt.setp(plotOut['medians'], color='black')
plt.setp(plotOut['fliers'], color='black')
plt.xticks([1, 2, 3], ['1', '$\\frac{10}{3}$', '10'])
plt.xlabel('Perturber ring density [# per $R_{Hill,Mab}$]')
plt.ylabel('Maximum longitude residual change [deg]')

for line in plotOut['medians']:
    # Get position data for median line.
    medianX, medianY = line.get_xydata()[1]

    # Overlay median value.
    axes.annotate( '%1.1f' % medianY, xy=(medianX+0.11, medianY), xycoords='data',\
                  horizontalalignment='center', verticalalignment='center')

plt.savefig(outputPathAndCasePrefix + "maximumLongitudeResidualVsRingDensity.pdf", \
            dpi = figureDPI)
plt.close()

# Plot maximum eccentricity change [-] vs. perturber ring mass [M_Mab].
eccentricityVsRingDensity = []
eccentricityVsRingDensity.append(outputSparse['maximumEccentricityChange'])
eccentricityVsRingDensity.append(outputNominal['maximumEccentricityChange'])
eccentricityVsRingDensity.append(outputDense['maximumEccentricityChange'])

fig = plt.figure()
axes = fig.gca()
plotOut = plt.boxplot(eccentricityVsRingDensity,sym=outliers)
plt.setp(plotOut['boxes'], color='black')
plt.setp(plotOut['whiskers'], color='black')
plt.setp(plotOut['medians'], color='black')
plt.setp(plotOut['fliers'], color='black')
plt.xticks([1, 2, 3], ['1', '$\\frac{10}{3}$', '10'])
plt.xlabel('Perturber ring density [# per $R_{Hill,Mab}$]')
plt.ylabel('Maximum eccentricity change [-]')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

for line in plotOut['medians']:
    # Get position data for median line.
    medianX, medianY = line.get_xydata()[1]

    # Overlay median value.
    axes.annotate( '%1.1e' % medianY, xy=(medianX+0.18, medianY), xycoords='data',\
                  horizontalalignment='center', verticalalignment='center')

plt.savefig(outputPathAndCasePrefix + "maximumEccentricityVsRingDensity.pdf", \
            dpi = figureDPI)
plt.close()

# Plot maximum inclination change [deg] vs. perturber ring mass [M_Mab].
inclinationVsRingDensity = []
inclinationVsRingDensity.append(numpy.rad2deg(outputSparse['maximumInclinationChange']))
inclinationVsRingDensity.append(numpy.rad2deg(outputNominal['maximumInclinationChange']))
inclinationVsRingDensity.append(numpy.rad2deg(outputDense['maximumInclinationChange']))

fig = plt.figure()
axes = fig.gca()
plotOut = plt.boxplot(inclinationVsRingDensity,sym=outliers)
plt.setp(plotOut['boxes'], color='black')
plt.setp(plotOut['whiskers'], color='black')
plt.setp(plotOut['medians'], color='black')
plt.setp(plotOut['fliers'], color='black')
plt.xticks([1, 2, 3], ['1', '$\\frac{10}{3}$', '10'])
plt.xlabel('Perturber ring density [# per $R_{Hill,Mab}$]')
plt.ylabel('Maximum inclination change [deg]')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

for line in plotOut['medians']:
    # Get position data for median line.
    medianX, medianY = line.get_xydata()[1]

    # Overlay median value.
    axes.annotate( '%1.1e' % medianY, xy=(medianX+0.18, medianY), xycoords='data',\
                  horizontalalignment='center', verticalalignment='center')

plt.savefig(outputPathAndCasePrefix + "maximumInclinationVsRingDensity.pdf", \
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