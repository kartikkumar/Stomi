'''   
Copyright (c) 2010-2013, Delft University of Technology
Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
All rights reserved.
See http://bit.ly/12SHPLR for license details.

This script can be used to generate plots of the initial Keplerian elements used for 
test particle simulations.
'''

################################################################################################### 
# Set up input deck
###################################################################################################

# Set absolute path to SQLite database with simulation data.
databasePath    = "/Users/kartikkumar/Desktop/stochasticMigrationResults.sqlite.backup"

# Set case name.
caseName        = "circular_equatorial_large_sma"

# Set absolute path to output directory.
outputPath      = "/Users/kartikkumar/Desktop/plots"

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
if not os.path.exists( outputPath ):
    os.makedirs( outputPath )

###################################################################################################


################################################################################################### 
# Execute database operations
###################################################################################################

# Connect to SQLite database.
database = sqlite3.connect( databasePath )

# Enter database and retrieve desired data.
with database:
    # Set cursor to scan through database and execute queries.
    cursor = database.cursor( )
    
    # Select the case ID corresponding to the case name provided.
    cursor.execute( "SELECT caseId FROM test_particle_case \
                     WHERE caseName == \"" + caseName + "\";" )
    caseId = cursor.fetchall( )[ 0 ][ 0 ]  
    
    # Select all the case data associated with the case ID.
    cursor.execute( "SELECT * FROM test_particle_case \
                     WHERE caseId == " + str( caseId ) + ";" )
    rawCaseData = cursor.fetchall( )
    caseDataColumnNameList = [ column_name[ 0 ] for column_name in cursor.description ]  
    
    # Select all the input data associated with the case ID.
    cursor.execute( "SELECT * FROM test_particle_input \
                     WHERE caseId == " + str( caseId ) + ";" )
    inputDataColumnNameList = [ column_name[ 0 ] for column_name in cursor.description ]
    rawInputData = cursor.fetchall( )
    
# Store case data in dictionary, using column names from database.    
caseData = { } 
for i,name in enumerate( caseDataColumnNameList ):
    caseData[ name ] = rawCaseData[ 0 ][ i ]

# Store input data in structured array, using column names from database.
inputDataTypeList = [ ] 
for name in inputDataColumnNameList:
    inputDataTypeList.append( ( name, '<f8' ) )
    
inputData = numpy.array( rawInputData, dtype = numpy.dtype( inputDataTypeList ) )

###################################################################################################


###################################################################################################
# Plot histograms of input table data.
###################################################################################################

# Set output path and case-prefix for files generated.
outputPathAndCasePrefix = outputPath + "/case" + str( caseData[ 'caseId' ] ) + "_"

# Plot distribution of initial semi-major axes wrt to perturbed body.
plt.figure( )
plt.hist( ( inputData[ 'semiMajorAxis' ] - caseData[ 'perturbedBodySemiMajorAxisAtT0' ] ) \
            * constants.meterInKilometers, facecolor='w', edgecolor='k' )
plt.xlabel( "Initial semi-major axis with respect to perturbed body [km]" )
plt.ylabel( "Frequency [-]" )
plt.savefig( outputPathAndCasePrefix + "histogramInitialSemiMajorAxes.pdf", dpi = figureDPI )

# Plot distribution of initial eccentricities.
plt.figure( )
plt.hist( inputData[ 'eccentricity' ], facecolor='w', edgecolor='k' )
plt.xlabel( "Initial eccentricity [-]" )
plt.ylabel( "Frequency [-]" )
plt.savefig( outputPathAndCasePrefix + "histogramInitialEccentricities.pdf", dpi = figureDPI )

# Plot distribution of initial inclinations.
plt.figure( )
plt.hist( inputData[ 'inclination' ] * constants.radiansInDegrees, facecolor='w', edgecolor='k' )
plt.xlabel( "Initial inclination [deg]" )
plt.ylabel( "Frequency [-]" )
plt.savefig( outputPathAndCasePrefix + "histogramInitialInclinations.pdf", dpi = figureDPI )

# Plot distribution of initial argument of periapses.
plt.figure( )
plt.hist( inputData[ 'argumentOfPeriapsis' ] * constants.radiansInDegrees, \
          facecolor='w', edgecolor='k' )
plt.xlabel( "Initial argument of periapsis [deg]" )
plt.ylabel( "Frequency [-]" )
plt.savefig( outputPathAndCasePrefix + "histogramInitialArgumentOfPeriapsis.pdf", dpi = figureDPI )

# Plot distribution of initial longitude of ascending nodes.
plt.figure( )
plt.hist( inputData[ 'longitudeOfAscendingNode' ] * constants.radiansInDegrees, \
          facecolor='w', edgecolor='k' )
plt.xlabel( "Initial longitude of ascending node [deg]" )
plt.ylabel( "Frequency [-]" )
plt.savefig( outputPathAndCasePrefix + "histogramInitialLongitudesOfAscendingNode.pdf", \
             dpi = figureDPI )

# Plot distribution of initial true anomalies.
plt.figure( )
plt.hist( inputData[ 'trueAnomaly' ] * constants.radiansInDegrees, facecolor='w', edgecolor='k' )
plt.xlim( xmin = -180.0, xmax = 180.0 )
plt.xlabel( "Initial true anomaly [deg]" )
plt.ylabel( "Frequency [-]" )
plt.savefig( outputPathAndCasePrefix + "histogramInitialTrueAnomalies.pdf", dpi = figureDPI )

###################################################################################################


###################################################################################################
# Finalize timer and print elapsed time.
###################################################################################################

# Finalize timer.
endTime = time.time( )

# Print elapsed time for script [s].
print "This script took " + str(endTime-startTime) + "s"

###################################################################################################
    