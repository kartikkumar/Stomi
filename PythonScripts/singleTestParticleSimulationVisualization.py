'''   Copyright (c) 2010-2012, Delft University of Technology
      All rights reserved.
 
      Redistribution and use in source and binary forms, with or without modification, are
      permitted provided that the following conditions are met:
        - Redistributions of source code must retain the above copyright notice, this list of
          conditions and the following disclaimer.
        - Redistributions in binary form must reproduce the above copyright notice, this list of
          conditions and the following disclaimer in the documentation and/or other materials
          provided with the distribution.
        - Neither the name of the Delft University of Technology nor the names of its contributors
          may be used to endorse or promote products derived from this software without specific
          prior written permission.
  
      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
      OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
      MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
      COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
      EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
      GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
      AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
      NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
      OF THE POSSIBILITY OF SUCH DAMAGE.
  
      Changelog
        YYMMDD    Author            Comment
        120908    K. Kumar          File created based on MATLAB code.
  
      References
  
      Notes
  
'''

# Import necessary packages.
import numpy
import sqlite3
import matplotlib 
matplotlib.use( 'WXAgg' )   # generate pdf output 
import matplotlib.pyplot as pyplot
import time
import math

# Start timer.
startTime = time.time( )

''' 
Input deck
'''

# Set absolute path to SQLite database with simulation data.
databasePath = "/Users/kartikkumar/Documents/University/PhD/Simulations/Tudat/Workspace/tudatApplications/mabSimulations/PlanetaryRings/MabSimulations/data/MabSimulationFinalResults/testParticleData/case1/case1-hipparchos.db.testParticle"

# Set absolute path to directory with data files.
dataPath = "/Users/kartikkumar/Documents/University/PhD/Simulations/Tudat/Workspace/tudatApplications/mabSimulations/PlanetaryRings/MabSimulations/data/MabSimulationFinalResults/testParticleData/case1/individualTestParticleSimulations/5699"

# Set absolute path to output directory for plots.
outputPath = "/Users/kartikkumar/Documents/University/PhD/Simulations/Tudat/Workspace/tudatApplications/mabSimulations/PlanetaryRings/MabSimulations/data/MabSimulationFinalResults/testParticleData/case1/individualTestParticleSimulations/5699/plots"

# Set simulation number.
simulationNumber = 5699

# Show figures?
isShowFigures = 0

# Set figure dpi.
figureDPI = 600

''' 
Database operations
'''

# Connect to SQLite database.
database = sqlite3.connect( databasePath )

# Enter database and retrieve desired data.
with database:
    cursor = database.cursor( )
    cursor.execute( "SELECT * FROM cases" )
    rawCaseData = cursor.fetchall( )
    caseDataColumnNameList = [ column_name[ 0 ] for column_name in cursor.description ]
    cursor.execute("SELECT * FROM input_data WHERE \"simulation\" = " + str( simulationNumber ) )
    inputDataColumnNameList = [ column_name[ 0 ] for column_name in cursor.description ]
    inputData = cursor.fetchall( )
    cursor.execute("SELECT * FROM output_data WHERE \"simulation\" = " + str( simulationNumber ) )
    outputDataColumnNameList = [ column_name[ 0 ] for column_name in cursor.description ]
    outputData = cursor.fetchall( )
    
caseData = {} 
for i,name in enumerate( caseDataColumnNameList ):
    caseData[name] = rawCaseData[0][i]

inputDataTypeList = [ ] 
for name in inputDataColumnNameList:
    inputDataTypeList.append( ( name, '<f8' ) )
    
inputData = numpy.array( inputData, dtype = numpy.dtype( inputDataTypeList ) )

outputDataTypeList = [ ] 
for name in outputDataColumnNameList:
    outputDataTypeList.append( ( name, '<f8' ) )   
    
outputData = numpy.array( outputData, dtype = numpy.dtype( outputDataTypeList ) )

''' 
Input file operations
'''
    
# Read Mab propagation history from file and store in a structured array.    
mabPropagationHistoryFile = "/mabPropagationHistory_case" + str( caseData[ 'case' ] ) \
                             + "_simulation" + str( simulationNumber) + ".dat"  
                               
rawMabPropagationHistory = numpy.genfromtxt( dataPath + mabPropagationHistoryFile, delimiter=',' )

mabPropagationHistory = rawMabPropagationHistory.view( dtype = [ ( "epoch", '<f8' ),
                                                                 ( "xPosition", '<f8' ),
                                                                 ( "yPosition", '<f8' ),
                                                                 ( "zPosition", '<f8' ),
                                                                 ( "xVelocity", '<f8' ),
                                                                 ( "yVelocity", '<f8' ),
                                                                 ( "zVelocity", '<f8' ) ] )

# Read test particle propagation history from file and store in a structured array.    
testParticlePropagationHistoryFile = "/testParticlePropagationHistory_case" + str( caseData[ 'case' ] ) \
                                     + "_simulation" + str( simulationNumber) + ".dat"  
                               
rawTestParticlePropagationHistory = numpy.genfromtxt( dataPath + testParticlePropagationHistoryFile, \
                                                      delimiter=',' )
    
testParticlePropagationHistory = rawTestParticlePropagationHistory.view( dtype = [ ( "epoch", '<f8' ),
                                                                                   ( "xPosition", '<f8' ),
                                                                                   ( "yPosition", '<f8' ),
                                                                                   ( "zPosition", '<f8' ),
                                                                                   ( "xVelocity", '<f8' ),
                                                                                   ( "yVelocity", '<f8' ),
                                                                                   ( "zVelocity", '<f8' ) ] ) 

## Read Mab propagation history in Keplerian elements from file and store in a structured array.    
#mabPropagationHistoryInKeplerianElementsFile = "/mabPropagationHistoryInKeplerianElements_case" \
#                                                + str( caseData[ 'case' ] ) + "_simulation" \
#                                                + str( simulationNumber) + ".dat"  
#                               
#rawMabPropagationHistoryInKeplerianElements \
# = numpy.genfromtxt( dataPath + mabPropagationHistoryInKeplerianElementsFile, delimiter=',' )
#
#mabPropagationHistoryInKeplerianElements = rawMabPropagationHistoryInKeplerianElements.view( \
#    dtype = [ ( "epoch", '<f8' ), ( "semiMajorAxis", '<f8' ), ( "eccentricity", '<f8' ), \
#              ( "inclination", '<f8' ), ( "argumentOfPeriapsis", '<f8' ), \
#              ( "longitudeOfAscendingNode", '<f8' ), ( "trueAnomaly", '<f8' ) ] )

# Read test particle propagation history in Keplerian elements from file and store in a structured array.    
testParticlePropagationHistoryInKeplerianElementsFile = "/testParticlePropagationHistoryInKeplerianElements_case" \
                                                        + str( caseData[ 'case' ] ) + "_simulation" \
                                                        + str( simulationNumber) + ".dat"  
                               
rawTestParticlePropagationHistoryInKeplerianElements \
 = numpy.genfromtxt( dataPath + testParticlePropagationHistoryInKeplerianElementsFile, delimiter=',' )

testParticlePropagationHistoryInKeplerianElements \
= rawTestParticlePropagationHistoryInKeplerianElements.view( \
    dtype = [ ( "epoch", '<f8' ), ( "semiMajorAxis", '<f8' ), ( "eccentricity", '<f8' ), \
              ( "inclination", '<f8' ), ( "argumentOfPeriapsis", '<f8' ), \
              ( "longitudeOfAscendingNode", '<f8' ), ( "trueAnomaly", '<f8' ) ] )

# Read distancehistory from file and store in a structured array.    
distanceHistoryFile = "/distanceHistory_case" + str( caseData[ 'case' ] ) \
 + "_simulation" + str( simulationNumber ) + ".dat"  
                               
rawDistanceHistory = numpy.genfromtxt( dataPath + distanceHistoryFile, delimiter=',' )
    
distanceHistory = rawDistanceHistory.view( dtype = [ ( "epoch", '<f8' ), ( "distance", '<f8' ) ] ) 

''' 
Data processing
'''

# Set length of Julian year [s].
julianYear = 365.25 * 86400.0
   
# Set minimum crossing distance [km].
minimumCrossingDistance = 2.0e4

# Set maximum crossing distance [km].
maximumCrossingDistance = 1.8e5

# Compute Jacobi energy.
#jacobiIntegral = 2.0 * ( caseData[ 'uranusGravitationalParameter' ] /  )

# Compute energy history.
#plottingData = numpy.array( numpy.zeros( len( distanceHistory ) ),
#                            dtype = numpy.dtype( [ ( 'semiMajorAxis', '<f8' ) ] ) )
#plottingData[ 'semiMajorAxis' ] = testParticlePropagationHistoryInKeplerianElements[ 'semiMajorAxis' ] \
#                                    - testParticlePropagationHistoryInKeplerianElements[ 'semiMajorAxis' ][ 0 ]

''' 
Visualization
'''

# Set output path and case-prefix for files generated.
outputPathAndCasePrefix = outputPath + "/case" + str( caseData[ 'case' ] ) + "_"

###################################################################################################  
# Plot distance history.
###################################################################################################  
    
figure = pyplot.figure( )
figureAxes = pyplot.gca( )   
pyplot.plot( distanceHistory[ 'epoch' ] / julianYear, distanceHistory[ 'distance' ] * 1.0e-8, 'k' ) 
pyplot.plot( outputData[ 'encounterEpoch' ] / julianYear, outputData[ 'encounterDistance' ] \
             * 1.0e-8, 'o', markerfacecolor = '0.75' ) 
pyplot.plot( outputData[ 'preEncounterEpoch' ] / julianYear, outputData[ 'preEncounterDistance' ] \
             * 1.0e-8, 'wo' ) 
pyplot.plot( outputData[ 'postEncounterEpoch' ] / julianYear, \
             outputData[ 'postEncounterDistance' ] * 1.0e-8, 'wo' ) 
pyplot.plot( [ distanceHistory[ 'epoch' ][ 0 ] / julianYear, \
               distanceHistory[ 'epoch' ][ -1 ] / julianYear ] , \
             [ minimumCrossingDistance * 1.0e-5, minimumCrossingDistance * 1.0e-5 ], 'k', \
             linestyle = 'dashed' ) 
pyplot.plot( [ distanceHistory[ 'epoch' ][ 0 ] / julianYear, \
               distanceHistory[ 'epoch' ][ -1 ] / julianYear ] , \
             [ maximumCrossingDistance * 1.0e-5, maximumCrossingDistance * 1.0e-5 ], 'k', \
             linestyle = 'dashed' ) 
pyplot.xlabel( "Simulation epoch [Julian year]" )
pyplot.ylabel( "Distance between test particle and Mab [$10^{5}$ km]" )
pyplot.xlim( xmin = -10.0, xmax = 60.0 )
pyplot.ylim( ymin = -0.05 )
pyplot.savefig( outputPathAndCasePrefix + "distanceHistory_case" + str( caseData[ 'case' ] ) \
                + "_simulation" + str( simulationNumber ) + ".pdf", dpi = figureDPI )
if ( isShowFigures == 1 ):
    pyplot.show( )
pyplot.close( )

###################################################################################################  
# Plot semi-major axis, eccentricity and inclination histories of test particle.
###################################################################################################  

figure = pyplot.figure( )

# Three subplots sharing x-axes.
figure, (subPlot1Axes, subPlot2Axes, subPlot3Axes) = pyplot.subplots( 3, sharex=True )
pyplot.xlabel( "Simulation epoch [Julian year]" )
pyplot.xlim( xmin = 0.0, xmax = 50.0 )
subPlot1Axes.plot( distanceHistory[ 'epoch' ] / julianYear, \
             ( testParticlePropagationHistoryInKeplerianElements[ 'semiMajorAxis' ] \
             - caseData[ 'mabSemiMajorAxisAtT0' ] ) * 0.001, 'k' )
subPlot1Axes.set_ylabel( "$\Delta a_{Mab}$ [km]" )
subPlot1Axes.text( -6.0, 600.0, "(a)", ha="center", va="center", size=14 )
subPlot2Axes.plot( distanceHistory[ 'epoch' ] / julianYear, \
                   testParticlePropagationHistoryInKeplerianElements[ 'eccentricity' ] * 1.0e3, \
                   'k' ) 
subPlot2Axes.set_ylabel( "e [$10^{-3}$]" )
subPlot2Axes.text( -6.0, 10.0, "(b)", ha="center", va="center", size=14 )
subPlot2Axes.set_ylim( ymin = 0.0, ymax = 10.0 )
subPlot3Axes.plot( distanceHistory[ 'epoch' ] / julianYear, \
                   testParticlePropagationHistoryInKeplerianElements[ 'inclination' ] \
                   / math.pi * 180.0, 'k' ) 
subPlot3Axes.set_ylabel( "i [deg]" )
subPlot3Axes.text( -6.0, 0.20, "(c)", ha="center", va="center", size=14 )
subPlot3Axes.set_ylim( ymin = 0.0, ymax = 0.2 )
# Fine-tune figure and hide x ticks for all but bottom plot.
figure.subplots_adjust( hspace = 0.35 )
pyplot.setp( [ a.get_xticklabels( ) for a in figure.axes[ :-1 ] ], visible = False )
pyplot.savefig( outputPathAndCasePrefix + "orbitalElementHistory_case" \
                + str( caseData[ 'case' ] ) + "_simulation" + str( simulationNumber ) \
                + ".pdf", dpi = figureDPI )
if ( isShowFigures == 1 ):
    pyplot.show( )
pyplot.close( )
    
''' 
End of file
'''

endTime = time.time( )
print endTime - startTime
