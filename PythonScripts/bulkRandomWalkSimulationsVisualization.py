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
        121011    K. Kumar          File created.
  
      References
  
      Notes
  
'''

# Import necessary packages.
import math
import matplotlib
from enable.enable_traits import LineStyle
# Set backend to use to generate PDF output.
matplotlib.use( 'WXAgg' )
import matplotlib.pyplot as pyplot
import numpy
import sqlite3
import time

# Start timer.
startTime = time.time( )

################################################################################################### 
# Input deck
###################################################################################################

# Set absolute path to SQLite database with simulation data.
databasePath = "/Users/kartikkumar/Desktop/MabSimulationFinalResults/randomWalkData/case2_horseshoe/case2_horseshoe-hipparchos.db.randomWalk"

# Set absolute path to output directory.
outputPath = "/Users/kartikkumar/Desktop/MabSimulationFinalResults/randomWalkData/case2_horseshoe/plots/p=1000,m=equal,0d001/"

# Set perturber population.
perturberPopulation = 1000

# Set mass distribution type.
massDistributionType = "EQUAL"

# Set mass distribution parameter 1.
massDistributionParameter1 = 0.001

# Set mass distribution parameter 2.
massDistributionParameter2 = 0.0

# Set number of histogram bins.
histogramBins = 50

# Show figures?
isShowFigures = 0

# Set figure dpi.
figureDPI = 600

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
    cursor.execute( "SELECT * FROM random_walk WHERE \"run\" IN " \
                    "(SELECT \"run\" FROM random_walk_runs WHERE \"perturberPopulation\" = " \
                   + str( perturberPopulation ) + ")" )
    rawRandomWalkData = cursor.fetchall( )
    randomWalkDataColumnNameList = [ column_name[ 0 ] for column_name in cursor.description ]

# Store case data in dictionary, using column names from database.    
caseData = { } 
for i,name in enumerate( caseDataColumnNameList ):
    caseData[ name ] = rawCaseData[ 0 ][ i ]
    
# Store random walk data in structured array, using column names from database.
randomWalkDataTypeList = [ ] 
for name in randomWalkDataColumnNameList:
    randomWalkDataTypeList.append( ( name, '<f8' ) )

randomWalkData = numpy.array( rawRandomWalkData, dtype = numpy.dtype( randomWalkDataTypeList ) )

###################################################################################################
# Visualization
###################################################################################################

# Set output path and case-prefix for files generated.
outputPathAndCasePrefix = outputPath + "/case" + str( caseData[ 'case' ] ) + "_"
   
# Plot histograms.   
pyplot.figure( ) 
figureAxes = pyplot.gca( ) 
pyplot.hist( randomWalkData[ 'maximumLongitudeResidualChange' ] / math.pi * 180.0, \
             bins = histogramBins, facecolor = 'w', edgecolor = 'k' )
pyplot.xlabel( "Maximum longitude residual change [deg]" )
pyplot.ylabel( "Frequency [-]" )
pyplot.savefig( outputPathAndCasePrefix + "maximumLongitudeResidualChanges.png", dpi = figureDPI )
pyplot.close( )      
 
pyplot.figure( ) 
figureAxes = pyplot.gca( ) 
pyplot.hist( randomWalkData[ 'maximumEccentricityChange' ], \
             bins = histogramBins, facecolor = 'w', edgecolor = 'k' )
figureAxes.xaxis.major.formatter.set_powerlimits( ( 0, 0 ) ) 
pyplot.xlabel( "Maximum eccentricity change [-]" )
pyplot.ylabel( "Frequency [-]" )
pyplot.savefig( outputPathAndCasePrefix + "maximumEccentricityChanges.png", dpi = figureDPI )
pyplot.close( )

pyplot.figure( ) 
figureAxes = pyplot.gca( ) 
pyplot.hist( randomWalkData[ 'maximumInclinationChange' ] / math.pi * 180.0, \
             bins = histogramBins, facecolor = 'w', edgecolor = 'k' )
pyplot.xlabel( "Maximum inclination change [deg]" )
pyplot.ylabel( "Frequency [-]" )
pyplot.savefig( outputPathAndCasePrefix + "maximumInclinationChanges.png", dpi = figureDPI )
pyplot.close( )      

''' 
End of file
'''

endTime = time.time( )
print endTime - startTime