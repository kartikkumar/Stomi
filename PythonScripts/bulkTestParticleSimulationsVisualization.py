'''   
Copyright (c) 2010-2013, Delft University of Technology      
Distributed under the BSD 3-Clause License 
(see COPYING or visit http://opensource.org/licenses/BSD-3-Clause )
'''

# Import necessary external modules.
import math
import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import pylab
import sqlite3
import time

# Import user modules.
import constants

# Start timer.
startTime = time.time( )

################################################################################################### 
# Input deck
###################################################################################################

# Set root path.
rootPath = "/Users/kartikkumar/Desktop/data/case3/"

# Set absolute path to SQLite database with simulation data.
databasePath = rootPath + "/case3-hipparchos.db.testParticle"

# Set absolute path to output directory.
outputPath = rootPath + "/plots"

# Set cut-off for big kicks included in plots.
bigKicksCutOff = 5000

# Show figures?
isShowFigures = 0

# Set figure dpi.
figureDPI = 600

# Show input histograms?
isShowInputHistograms = 0

################################################################################################### 
# Database operations
###################################################################################################

# Connect to SQLite database.
database = sqlite3.connect( databasePath )

# Enter database and retrieve desired data.
with database:
    cursor = database.cursor( )
    cursor.execute( "SELECT * FROM cases" )
    rawCaseData = cursor.fetchall( )
    caseDataColumnNameList = [ column_name[ 0 ] for column_name in cursor.description ]
    cursor.execute( "SELECT * FROM input_data" )
    inputDataColumnNameList = [ column_name[ 0 ] for column_name in cursor.description ]
    inputData = cursor.fetchall( )
    cursor.execute( "SELECT * FROM output_data" )
    outputDataColumnNameList = [ column_name[ 0 ] for column_name in cursor.description ]
    outputData = cursor.fetchall( )
    
# Store case data in dictionary, using column names from database.    
caseData = { } 
for i,name in enumerate( caseDataColumnNameList ):
    caseData[ name ] = rawCaseData[ 0 ][ i ]

# Store input data in structured array, using column names from database.
inputDataTypeList = [ ] 
for name in inputDataColumnNameList:
    inputDataTypeList.append( ( name, '<f8' ) )
    
inputData = numpy.array( inputData, dtype = numpy.dtype( inputDataTypeList ) )

# Store output data in structured array, using column names from database.
outputDataTypeList = [ ] 
for name in outputDataColumnNameList:
    outputDataTypeList.append( ( name, '<f8' ) )
    
outputData = numpy.array( outputData, dtype = numpy.dtype( outputDataTypeList ) )    
    
##################################################################################################
# Data processing
###################################################################################################
   
# Declare structured array containing data using to generate plots.   
plottingData = numpy.array( numpy.zeros( len( outputData ) ),
                            dtype = numpy.dtype( [ ( 'semiMajorAxisKick', '<f8' ),
                                                   ( 'eccentricityKick', '<f8' ),
                                                   ( 'inclinationKick', '<f8' ),
                                                   ( 'semiMajorAxisKickMagnitude', '<f8' ),
                                                   ( 'eccentricityKickMagnitude', '<f8' ),
                                                   ( 'inclinationKickMagnitude', '<f8' ),
                                                   ( 'preEncounterSemiMajorAxis', '<f8' ),
                                                   ( 'preEncounterEccentricity', '<f8' ),
                                                   ( 'preEncounterInclination', '<f8' ),
                                                   ( 'encounterDuration', '<f8' ),
                                                   ( 'encounterDistance', '<f8' ),
                                                   ( 'simulation', '<i4' ),
                                                   ( 'initialSemiMajorAxis', '<f8' ) ] ) )

# Store semi-major axis kicks [km].
plottingData[ 'semiMajorAxisKick' ] = ( outputData['postEncounterSemiMajorAxis'] \
                                      - outputData['preEncounterSemiMajorAxis'] ) \
                                      * constants.meterInKilometers
     
# Store eccentricity kicks [-].                                    
plottingData[ 'eccentricityKick' ] = ( outputData[ 'postEncounterEccentricity' ] - \
                                       outputData[ 'preEncounterEccentricity' ] )

# Store inclination kicks [deg].
plottingData[ 'inclinationKick' ] = ( outputData[ 'postEncounterInclination' ] - \
                                      outputData[ 'preEncounterInclination' ] ) \
                                      * constants.radiansInDegrees

# Store (magnitude) semi-major axis kicks [km].
plottingData[ 'semiMajorAxisKickMagnitude' ] = numpy.abs( plottingData[ 'semiMajorAxisKick' ] )
     
# Store (magnitude) eccentricity kicks [-].                                    
plottingData[ 'eccentricityKickMagnitude' ] = numpy.abs( plottingData[ 'eccentricityKick' ] ) 

# Store (magnitude) inclination kicks [deg].
plottingData[ 'inclinationKickMagnitude' ] = numpy.abs( plottingData[ 'inclinationKick' ] )                                    

# Store pre-encounter semi-major axes [km].                                               
plottingData[ 'preEncounterSemiMajorAxis' ] = ( outputData[ 'preEncounterSemiMajorAxis' ] - \
                                                caseData[ 'mabSemiMajorAxisAtT0' ] ) \
                                              * constants.meterInKilometers

# Store pre-encounter eccentricity [-].
plottingData[ 'preEncounterEccentricity' ] = outputData[ 'preEncounterEccentricity' ] 

# Store pre-encounter eccentricity [km].
plottingData[ 'preEncounterInclination' ] = outputData[ 'preEncounterInclination' ] \
                                            * constants.radiansInDegrees                                                                                           

# Store encounter duration [Julian days].
plottingData[ 'encounterDuration' ] = numpy.abs( outputData[ 'encounterDuration' ] \
                                                 * constants.secondsInJulianDays )

# Store encounter distance [km].
plottingData[ 'encounterDistance' ] = outputData[ 'encounterDistance' ] \
                                        * constants.meterInKilometers

# Store simulation numbers.
plottingData[ 'simulation' ] = outputData[ 'simulation' ]

# Store initial semi-major axes with respect to Mab [km].
plottingData[ 'initialSemiMajorAxis' ] \
= ( inputData[ plottingData[ 'simulation' ] - 1 ][ 'semiMajorAxis' ] \
    - caseData[ 'mabSemiMajorAxisAtT0' ] ) * constants.meterInKilometers

# Sort plotting data in descending order based on semi-major axis kick (magnitude).
semiMajorAxisKickSortedData \
= plottingData[ numpy.argsort( plottingData, order='semiMajorAxisKickMagnitude' ) ][ : :-1 ]

# Sort plotting data in descending order based on eccentricity kick (magnitude).
eccentricitySortedData \
= plottingData[ numpy.argsort( plottingData, order='eccentricityKickMagnitude' ) ][ : :-1 ]

# Sort plotting data in descending order based on inclination kick (magnitude).
inclinationSortedData \
= plottingData[ numpy.argsort( plottingData, order='inclinationKickMagnitude' ) ][ : :-1 ]

# Compute Hill sphere radius.
hillRadius = caseData[ 'mabSemiMajorAxisAtT0' ] \
             * ( caseData[ 'mabGravitationalParameter' ] \
                 / ( 3.0 * caseData[ 'uranusGravitationalParameter' ]  ) )**( 1.0 / 3.0 )
                 
# Compute Hill velocity using Hill radius
hillVelocity = hillRadius * math.sqrt( caseData[ 'uranusGravitationalParameter' ] \
                 / caseData[ 'mabSemiMajorAxisAtT0' ]**3.0 )
                 
# Extract indices for horseshoe orbits.
horseshoeOrbitIndices = ( numpy.abs( plottingData[ 'initialSemiMajorAxis' ] ) < 40.0 ) \
                        & ( plottingData[ 'semiMajorAxisKickMagnitude' ] > 1.0 )
                        
# Extract unique simulation numbers for horseshoe orbits.
horseshoeOrbitSimulationNumbers \
= numpy.unique( plottingData[ 'simulation' ][ horseshoeOrbitIndices ] ) 

###################################################################################################
# Visualization
###################################################################################################

# Set output path and case-prefix for files generated.
outputPathAndCasePrefix = outputPath + "/case" + str( caseData[ 'case' ] ) + "_"

'''
Plot histograms of initial conditions (if isShowInputHistograms = 1).
'''

if ( isShowInputHistograms == 1 ):
  pyplot.figure( )
  pyplot.hist( ( inputData[ 'semiMajorAxis' ] - caseData[ 'mabSemiMajorAxisAtT0' ] ) \
              * constants.meterInKilometers, facecolor='w', edgecolor='k' )
  pyplot.xlabel( "Initial semi-major axis with respect to Mab [km]" )
  pyplot.ylabel( "Frequency [-]" )
  pyplot.savefig( outputPathAndCasePrefix +  "histogramInitialSemiMajorAxes.pdf", dpi = figureDPI )
  pyplot.close( )

  pyplot.figure( )
  pyplot.hist( inputData[ 'eccentricity' ], facecolor='w', edgecolor='k' )
  pyplot.xlabel( "Initial eccentricity [-]" )
  pyplot.ylabel( "Frequency [-]" )
  pyplot.savefig( outputPathAndCasePrefix + "histogramInitialEccentricities.pdf", dpi = figureDPI )
  pyplot.close( )

  pyplot.figure( )
  pyplot.hist( inputData[ 'inclination' ] * constants.radiansInDegrees, \
               facecolor='w', edgecolor='k' )
  pyplot.xlabel( "Initial inclination [deg]" )
  pyplot.ylabel( "Frequency [-]" )
  pyplot.savefig( outputPathAndCasePrefix + "histogramInitialInclinations.pdf", dpi = figureDPI )
  pyplot.close( )

  pyplot.figure( )
  pyplot.hist( inputData[ 'argumentOfPeriapsis' ] * constants.radiansInDegrees, 
               facecolor='w', edgecolor='k' )
  pyplot.xlabel( "Initial argument of periapsis [deg]" )
  pyplot.ylabel( "Frequency [-]" )
  pyplot.savefig( outputPathAndCasePrefix + "histogramInitialArgumentOfPeriapsis.pdf", \
                 dpi = figureDPI )
  pyplot.close( )

  pyplot.figure( )
  pyplot.hist( inputData[ 'longitudeOfAscendingNode' ] * constants.radiansInDegrees, \
              facecolor='w', edgecolor='k' )
  pyplot.xlabel( "Initial longitude of ascending node [deg]" )
  pyplot.ylabel( "Frequency [-]" )
  pyplot.savefig( outputPathAndCasePrefix + "histogramInitialLongitudesOfAscendingNode.pdf", \
                 dpi = figureDPI )
  pyplot.close( )

  pyplot.figure( )
  pyplot.hist( inputData[ 'trueAnomaly' ] * constants.radiansInDegrees, 
               facecolor='w', edgecolor='k' )
  pyplot.xlim( xmin = -180.0, xmax = 180.0 )
  pyplot.xlabel( "Initial true anomaly [deg]" )
  pyplot.ylabel( "Frequency [-]" )
  pyplot.savefig( outputPathAndCasePrefix + "histogramInitialTrueAnomalies.pdf", dpi = figureDPI )
  pyplot.close( )

'''
Plot eccentricity kicks vs. semi-major axis kicks.
'''

# Plot eccentricity vs. semi-major axis kicks (magnitude).
pyplot.figure( )   
pyplot.plot( plottingData[ 'semiMajorAxisKickMagnitude' ], \
             plottingData[ 'eccentricityKickMagnitude' ], '.k', rasterized = True )
pyplot.xscale( 'log' )
pyplot.yscale( 'log' )
pyplot.xlabel( "Semi-major axis kick [km]" )
pyplot.ylabel( "Eccentricity kick [-]" )
pyplot.savefig( outputPathAndCasePrefix + "eccentricityVsSemiMajorAxisKicks.pdf", \
               dpi = figureDPI )
pyplot.close( )

# Plot eccentricity vs. semi-major axis kicks (big kicks plot). 
pyplot.figure( ) 
figureAxes = pyplot.gca( )        
pyplot.plot( semiMajorAxisKickSortedData[ 'semiMajorAxisKick' ][ :bigKicksCutOff ], \
            semiMajorAxisKickSortedData[ 'eccentricityKick' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( eccentricitySortedData[ 'semiMajorAxisKick' ][ :bigKicksCutOff ], \
            eccentricitySortedData[ 'eccentricityKick' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( inclinationSortedData[ 'semiMajorAxisKick' ][ :bigKicksCutOff ], \
            inclinationSortedData[ 'eccentricityKick' ][ :bigKicksCutOff ], '.k' )
pyplot.xlabel( "Semi-major axis kick [km]" )
pyplot.ylabel( "Eccentricity kick [-]" )
pyplot.savefig( outputPathAndCasePrefix + "eccentricityVsSemiMajorAxisBigKicks.pdf", \
               dpi = figureDPI )    
pyplot.close( )

# Plot eccentricity vs. semi-major axis kicks (magnitudes) (big kicks plot). 
pyplot.figure( ) 
figureAxes = pyplot.gca( )        
pyplot.plot( semiMajorAxisKickSortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], \
            semiMajorAxisKickSortedData[ 'eccentricityKickMagnitude' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( eccentricitySortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], \
            eccentricitySortedData[ 'eccentricityKickMagnitude' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( inclinationSortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], \
            inclinationSortedData[ 'eccentricityKickMagnitude' ][ :bigKicksCutOff ], '.k' )
figureAxes.yaxis.major.formatter.set_powerlimits( ( 0, 0 ) ) 
pyplot.xlabel( "Magnitude of semi-major axis kick [km]" )
pyplot.ylabel( "Magnitude of eccentricity kick [-]" )
pyplot.xlim( xmin = 0.0 )
pyplot.ylim( ymin = 0.0 )
pyplot.savefig( outputPathAndCasePrefix + "eccentricityVsSemiMajorAxisBigKickMagnitudes.pdf", \
               dpi = figureDPI )    
pyplot.close( )

# Plot eccentricity vs. semi-major axis kicks (magnitudes) (big kicks zoom plot).      
pyplot.figure( ) 
figureAxes = pyplot.gca( )   
pyplot.plot( semiMajorAxisKickSortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], \
            semiMajorAxisKickSortedData[ 'eccentricityKickMagnitude' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( eccentricitySortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], \
            eccentricitySortedData[ 'eccentricityKickMagnitude' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( inclinationSortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], \
            inclinationSortedData[ 'eccentricityKickMagnitude' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( plottingData[ 'semiMajorAxisKickMagnitude' ][ horseshoeOrbitIndices ], \
            plottingData[ 'eccentricityKickMagnitude' ][ horseshoeOrbitIndices ], \
            marker = '.', color = '0.5', linestyle = 'None' )
figureAxes.yaxis.major.formatter.set_powerlimits( ( 0, 0 ) ) 
pyplot.xlabel( "Magnitude of semi-major axis kick [km]" )
pyplot.ylabel( "Magnitude of eccentricity kick [-]" )
pyplot.xlim( xmin = 20.0, xmax = 60.0 )
pyplot.ylim( ymin = 0.0, ymax = 3.0e-4 )    
pyplot.savefig( outputPathAndCasePrefix + "eccentricityVsSemiMajorAxisKickMagnitudesZoom.pdf", \
               dpi = figureDPI )
pyplot.close( )

# Plot eccentricity vs. semi-major axis kicks (magnitudes) (big kicks interactive plot). 
labels = ['{0}'.format( simulation ) for simulation in \
          semiMajorAxisKickSortedData[ 'simulation' ][ :bigKicksCutOff ] ]  
pyplot.figure( )  
figureAxes = pyplot.gca( )        
pyplot.plot( semiMajorAxisKickSortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], \
            semiMajorAxisKickSortedData[ 'eccentricityKickMagnitude' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( eccentricitySortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], \
            eccentricitySortedData[ 'eccentricityKickMagnitude' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( inclinationSortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], \
            inclinationSortedData[ 'eccentricityKickMagnitude' ][ :bigKicksCutOff ], '.k' )
figureAxes.yaxis.major.formatter.set_powerlimits( ( 0, 0 ) ) 
for label, semiMajorAxisKick, eccentricityKick in \
   zip( labels, semiMajorAxisKickSortedData[ 'semiMajorAxisKickMagnitude' ], \
        semiMajorAxisKickSortedData[ 'eccentricityKickMagnitude' ] ):
   pyplot.annotate( label, xy = ( semiMajorAxisKick, eccentricityKick ), xytext = ( -5, 5 ), \
                    textcoords = 'offset points', ha = 'right', va = 'bottom', size = 'small' )
pyplot.xlabel( "Magnitude of semi-major axis kick [km]" )
pyplot.ylabel( "Magnitude of eccentricity kick [-]" )
pyplot.xlim( xmin = 0.0 )
pyplot.ylim( ymin = 0.0 )
if ( isShowFigures == 1 ):
   pyplot.show( )
pyplot.close( )

'''
Plot semi-major axis kicks vs. initial semi-major axes 
'''

# Plot semi-major axis kicks (magnitude) vs. initial semi-major axes with respect to Mab.
figure = pyplot.figure( )  
axis1 = figure.add_subplot( 111 )
axis1.plot( plottingData[ 'initialSemiMajorAxis' ], plottingData[ 'semiMajorAxisKickMagnitude' ], \
            '.k', rasterized = True )
axis1.plot( plottingData[ 'initialSemiMajorAxis' ][ horseshoeOrbitIndices ], \
            plottingData[ 'semiMajorAxisKickMagnitude' ][ horseshoeOrbitIndices ], \
            marker = '.', color = '0.5', linestyle = 'None' )
axis1.set_yscale('log')
axis1.set_xlabel( "Initial semi-major axis relative to Mab [km]" )
axis1.set_ylabel( "Semi-major axis kick magnitude [km]" )
axis1.set_xlim( -1000.0, 1000.0 )
axis1.set_ylim( ymin = 0.0 )
axis2 = axis1.twiny( )
axis2.set_xlim( -1000.0, 1000.0 )
axis2.set_xticks( hillRadius * constants.meterInKilometers * numpy.arange( -25.0, 30.0, 5.0 ) )
axis2.set_xticklabels( numpy.arange( -25.0, 30.0, 5.0 ) )
axis2.set_xlabel( "Initial semi-major axis relative to Mab [Hill radii]" )
pyplot.savefig( outputPathAndCasePrefix + "semiMajorAxisKickMagnitudeVsInitialSemiMajorAxes.pdf", \
               dpi = figureDPI )
pyplot.close( )

# Plot semi-major axis kicks (magnitude) vs. initial semi-major axes with respect to Mab 
# (interactive).
figure = pyplot.figure( )  
axis1 = figure.add_subplot( 111 )
axis1.plot( plottingData[ 'initialSemiMajorAxis' ], plottingData[ 'semiMajorAxisKickMagnitude' ], \
           '.k', rasterized = True )
axis1.set_yscale('log')
axis1.set_xlabel( "Initial semi-major axis relative to Mab [km]" )
axis1.set_ylabel( "Semi-major axis kick magnitude [km]" )
axis1.set_xlim( -1000.0, 1000.0 )
axis1.set_ylim( ymin = 0.0 )
axis2 = axis1.twiny( )
axis2.set_xlim( -1000.0, 1000.0 )
axis2.set_xticks( hillRadius * constants.meterInKilometers * numpy.arange( -25.0, 30.0, 5.0 ) )
axis2.set_xticklabels( numpy.arange( -25.0, 30.0, 5.0 ) )
axis2.set_xlabel( "Initial semi-major axis relative to Mab [Hill radii]" )
if ( isShowFigures == 1 ):
   pyplot.show( )
pyplot.close( )

# Plot semi-major axis kicks (magnitude) vs. initial semi-major axes with respect to Mab 
# (big kicks interactive plot). 
labels = ['{0}'.format( simulation ) for simulation in \
         semiMajorAxisKickSortedData[ 'simulation' ][ :bigKicksCutOff ] ]  
pyplot.figure( )  
pyplot.plot( semiMajorAxisKickSortedData[ 'initialSemiMajorAxis' ][ :bigKicksCutOff ], \
            semiMajorAxisKickSortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( eccentricitySortedData[ 'initialSemiMajorAxis' ][ :bigKicksCutOff ], \
            eccentricitySortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( inclinationSortedData[ 'initialSemiMajorAxis' ][ :bigKicksCutOff ], \
            inclinationSortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], '.k' )
for label, initialSemiMajorAxis, semiMajorAxisKickMagnitude in \
   zip( labels, semiMajorAxisKickSortedData[ 'initialSemiMajorAxis' ], \
        semiMajorAxisKickSortedData[ 'semiMajorAxisKickMagnitude' ] ):
   pyplot.annotate( label, xy = ( initialSemiMajorAxis, semiMajorAxisKickMagnitude ), \
                    xytext = ( -5, 5 ), textcoords = 'offset points', ha = 'right', \
                    va = 'bottom', size = 'small' )
pyplot.yscale('log')
pyplot.xlabel( "Initial semi-major axis with respect to Mab [km]" )
pyplot.ylabel( "Magnitude of semi-major axis kick [km]" )
pyplot.xlim( -1000.0, 1000.0 )
pyplot.ylim( ymin = 0.0 )
if ( isShowFigures == 1 ):
   pyplot.show( )
pyplot.close( )

'''
Plot eccentricity kicks vs. initial semi-major axes 
'''

# Plot eccentricity kicks (magnitude) vs. initial semi-major axes with respect to Mab.
figure = pyplot.figure( )  
axis1 = figure.add_subplot( 111 )
axis1.plot( plottingData[ 'initialSemiMajorAxis' ], plottingData[ 'eccentricityKickMagnitude' ], \
            '.k', rasterized = True )
axis1.plot( plottingData[ 'initialSemiMajorAxis' ][ horseshoeOrbitIndices ], \
            plottingData[ 'eccentricityKickMagnitude' ][ horseshoeOrbitIndices ], \
            marker = '.', color = '0.5', linestyle = 'None' )
axis1.set_yscale('log')
axis1.set_xlabel( "Initial semi-major axis relative to Mab [km]" )
axis1.set_ylabel( "Eccentricity kick magnitude [-]" )
axis1.set_xlim( -1000.0, 1000.0 )
axis1.set_ylim( ymin = 0.0 )
axis2 = axis1.twiny( )
axis2.set_xlim( -1000.0, 1000.0 )
axis2.set_xticks( hillRadius * constants.meterInKilometers * numpy.arange( -25.0, 30.0, 5.0 ) )
axis2.set_xticklabels( numpy.arange( -25.0, 30.0, 5.0 ) )
axis2.set_xlabel( "Initial semi-major axis relative to Mab [Hill radii]" )
pyplot.savefig( outputPathAndCasePrefix + "eccentricityKickMagnitudeVsInitialSemiMajorAxes.pdf", \
               dpi = figureDPI )
pyplot.close( )

''' 
Plot semi-major axis kicks vs. pre-encounter semi-major axes.
'''

# Plot semi-major axis kicks (magnitude) vs. pre-encounter semi-major axes with respect to Mab.   
figure = pyplot.figure( ) 
axis1 = figure.add_subplot( 111 )
pyplot.plot( plottingData[ 'preEncounterSemiMajorAxis' ], \
            plottingData[ 'semiMajorAxisKickMagnitude' ], '.k', rasterized = True )
axis1.plot( plottingData[ 'preEncounterSemiMajorAxis' ][ horseshoeOrbitIndices ], \
           plottingData[ 'semiMajorAxisKickMagnitude' ][ horseshoeOrbitIndices ], \
           marker = '.', color = '0.5', linestyle = 'None' )
axis1.set_yscale('log')
axis1.set_xlabel( "Pre-encounter semi-major axis relative to Mab [km]" )
axis1.set_ylabel( "Semi-major axis kick magnitude [km]" )
axis1.set_xlim( -1000.0, 1000.0 )
axis1.set_ylim( ymin = 0.0 )
axis2 = axis1.twiny( )
axis2.set_xlim( -1000.0, 1000.0 )
axis2.set_xticks( hillRadius * constants.meterInKilometers * numpy.arange( -25.0, 30.0, 5.0 ) )
axis2.set_xticklabels( numpy.arange( -25.0, 30.0, 5.0 ) )
axis2.set_xlabel( "Pre-encounter semi-major axis relative to Mab [Hill radii]" )
pyplot.savefig( outputPathAndCasePrefix + "semiMajorAxisKickVsPreEncounterSemiMajorAxes.pdf", \
               dpi = figureDPI )
pyplot.close( )

'''
Plot kick duration vs. semi-major axis kicks (magnitude) 
'''

# Plot kick duration vs. semi-major axis kicks (magnitude) (log-log plot).
pyplot.figure( )
pyplot.plot( plottingData[ 'semiMajorAxisKickMagnitude' ], plottingData[ 'encounterDuration' ], \
            '.k', rasterized = True )
pyplot.xscale( 'log' )
pyplot.yscale( 'log' )
pyplot.xlabel( "Magnitude of semi-major axis kick [km]" )
pyplot.ylabel( "Relative kick duration [-]" )
pyplot.xlim( xmin = 0.0 )
pyplot.ylim( ymin = 0.0 )
pyplot.savefig( outputPathAndCasePrefix + "kickDurationVsSemiMajorAxisKickMagnitudesLogLog.pdf", \
               dpi = figureDPI )    
pyplot.close( )

# Plot kick duration vs. semi-major axis kicks (magnitude) (big kicks plot).
pyplot.figure( )
figureAxes = pyplot.gca( )   
pyplot.plot( semiMajorAxisKickSortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], \
            semiMajorAxisKickSortedData[ 'encounterDuration' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( eccentricitySortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], \
            eccentricitySortedData[ 'encounterDuration' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( inclinationSortedData[ 'semiMajorAxisKickMagnitude' ][ :bigKicksCutOff ], \
            inclinationSortedData[ 'encounterDuration' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( plottingData[ 'semiMajorAxisKickMagnitude' ][ horseshoeOrbitIndices ], \
            plottingData[ 'encounterDuration' ][ horseshoeOrbitIndices ], \
            marker = '.', color = '0.5', linestyle = 'None' )
pyplot.xlabel( "Magnitude of semi-major axis kick [km]" )
pyplot.ylabel( "Relative kick duration [-]" )
pyplot.xlim( xmin = 0.0 )
pyplot.ylim( ymin = 0.0 )
pyplot.savefig( outputPathAndCasePrefix + "kickDurationVsSemiMajorAxisKickMagnitudeLogLogBigKicks.pdf", \
               dpi = figureDPI )
pyplot.close( )

'''
Plot kick duration vs. initial semi-major axis.
'''

# Plot kick duration vs. initial semi-major axis.
pyplot.figure( )
pyplot.plot( plottingData[ 'initialSemiMajorAxis' ], plottingData[ 'encounterDuration' ], \
            '.k', rasterized = True )
pyplot.yscale( 'log' )
pyplot.xlabel( "Initial semi-major axis relative to Mab [km]" )
pyplot.ylabel( "Relative kick duration [-]" )
pyplot.ylim( ymin = 0.0 )
pyplot.savefig( outputPathAndCasePrefix + "kickDurationVsInitialSemiMajorAxis.pdf", \
               dpi = figureDPI )    
pyplot.close( )

# Plot kick duration vs. semi-major axis kicks (magnitude) (big kicks plot).
pyplot.figure( )
figureAxes = pyplot.gca( )   
pyplot.plot( semiMajorAxisKickSortedData[ 'initialSemiMajorAxis' ][ :bigKicksCutOff ], \
            semiMajorAxisKickSortedData[ 'encounterDuration' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( eccentricitySortedData[ 'initialSemiMajorAxis' ][ :bigKicksCutOff ], \
            eccentricitySortedData[ 'encounterDuration' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( inclinationSortedData[ 'initialSemiMajorAxis' ][ :bigKicksCutOff ], \
            inclinationSortedData[ 'encounterDuration' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( plottingData[ 'initialSemiMajorAxis' ][ horseshoeOrbitIndices ], \
            plottingData[ 'encounterDuration' ][ horseshoeOrbitIndices ], \
            marker = '.', color = '0.5', linestyle = 'None' )
pyplot.xlabel( "Initial semi-major axis relative to Mab [km]" )
pyplot.ylabel( "Relative kick duration [-]" )
pyplot.ylim( ymin = 0.0 )
pyplot.savefig( outputPathAndCasePrefix + "kickDurationVsInitialSemiMajorAxisBigKicks.pdf", \
               dpi = figureDPI )
pyplot.close( )

'''
Plot kick distance vs. initial semi-major axis.
'''

# Plot kick distance vs. initial semi-major axis.
pyplot.figure( )
pyplot.plot( plottingData[ 'initialSemiMajorAxis' ], plottingData[ 'encounterDistance' ], \
            '.k', rasterized = True )
pyplot.yscale( 'log' )
pyplot.xlabel( "Initial semi-major axis relative to Mab [km]" )
pyplot.ylabel( "Kick distance [km]" )
pyplot.ylim( ymin = 0.0 )
pyplot.savefig( outputPathAndCasePrefix + "kickDistanceVsInitialSemiMajorAxis.pdf", \
               dpi = figureDPI )    
pyplot.close( )

# Plot kick distance vs. semi-major axis kicks (magnitude) (big kicks plot).
pyplot.figure( )
figureAxes = pyplot.gca( )   
pyplot.plot( semiMajorAxisKickSortedData[ 'initialSemiMajorAxis' ][ :bigKicksCutOff ], \
            semiMajorAxisKickSortedData[ 'encounterDistance' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( eccentricitySortedData[ 'initialSemiMajorAxis' ][ :bigKicksCutOff ], \
            eccentricitySortedData[ 'encounterDistance' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( inclinationSortedData[ 'initialSemiMajorAxis' ][ :bigKicksCutOff ], \
            inclinationSortedData[ 'encounterDistance' ][ :bigKicksCutOff ], '.k' )
pyplot.plot( plottingData[ 'initialSemiMajorAxis' ][ horseshoeOrbitIndices ], \
            plottingData[ 'encounterDistance' ][ horseshoeOrbitIndices ], \
            marker = '.', color = '0.5', linestyle = 'None' )
pyplot.xlabel( "Initial semi-major axis relative to Mab [km]" )
pyplot.ylabel( "Kick distance [km]" )
pyplot.ylim( ymin = 0.0 )
pyplot.savefig( outputPathAndCasePrefix + "kickDistanceVsInitialSemiMajorAxisBigKicks.pdf", \
               dpi = figureDPI )
pyplot.close( )

'''
Plot pre-encounter orbital velocity vs. pre-encounter semi-major axis
'''

#pyplot.figure( )
#figureAxes = pyplot.gca( )  
#pyplot.xlabel( "Pre-encounter semi-major axis relative to Mab [km]" )
#pyplot.ylabel( "Orbital velocity [Hill velocity]" )
#pyplot.ylim( ymin = 0.0 )
#pyplot.savefig( outputPathAndCasePrefix + "kickDistanceVsInitialSemiMajorAxisBigKicks.pdf", \
#                dpi = figureDPI )
#pyplot.close( )

''' 
End of file
'''

endTime = time.time( )
print endTime - startTime
    