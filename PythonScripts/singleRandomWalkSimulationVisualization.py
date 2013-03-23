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
        121018    K. Kumar          File created.
  
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

# Set absolute path to directory with data files.
outputDataPath = "/Users/kartikkumar/Desktop/"

# Set absolute path to output directory for plots.
outputPlotPath = "/Users/kartikkumar/Desktop/"

# Set case number.
caseNumber = 2

# Set perturber population.
perturberPopulation = 1000

# Show figures?
isShowFigures = 1

# Set figure dpi.
figureDPI = 600

''' 
Input file operations
'''

# Read Mab propagation history from file and store in a structured array.    
mabPropagationHistoryFile = "/keplerianActionElementsHistory_case" + str( caseNumber ) \
                             + "_perturbers" + str( perturberPopulation ) + ".dat"  
                               
rawMabPropagationHistory = numpy.genfromtxt( outputDataPath + mabPropagationHistoryFile, delimiter=',' )

mabPropagationHistory = rawMabPropagationHistory.view( dtype = [ ( "epoch", '<f8' ),
                                                                 ( "semiMajorAxis", '<f8' ),
                                                                 ( "eccentricity", '<f8' ),
                                                                 ( "inclination", '<f8' ) ] )

''' 
Data processing
'''

# Set length of Julian year [s].
julianYear = 365.25 * 86400.0

''' 
Visualization
'''

###################################################################################################  
# Plot Keplerian elements history.
###################################################################################################  

# Three subplots sharing x-axes.
figure, (subPlot1Axes, subPlot2Axes, subPlot3Axes) = pyplot.subplots( 3, sharex=True )
pyplot.xlabel( "Simulation epoch [Julian year]" )
pyplot.xlim( xmin = 0.0, xmax = 50.0 )
subPlot1Axes.plot( mabPropagationHistory[ 'epoch' ] / julianYear, \
                   ( mabPropagationHistory[ 'semiMajorAxis' ] - 97.736e6 ) * 0.001, 'k' )
subPlot1Axes.set_ylabel( "$\Delta a$ [km]" )
subPlot2Axes.plot( mabPropagationHistory[ 'epoch' ] / julianYear, \
                   mabPropagationHistory[ 'eccentricity' ], 'k' ) 
subPlot2Axes.yaxis.major.formatter.set_powerlimits( ( 0, 0 ) )
subPlot2Axes.set_ylabel( "e [-]" )
#subPlot2Axes.text( -200.0, 6.8e-3, "(b)", ha="center", va="center", size=14 )
#subPlot2Axes.set_ylim( ymin = 0.0 )
subPlot3Axes.plot( mabPropagationHistory[ 'epoch' ] / julianYear, \
                   mabPropagationHistory[ 'inclination' ] / math.pi * 180.0, 'k' ) 
subPlot3Axes.set_ylabel( "i [deg]" )
#subPlot3Axes.text( -200.0, 0.35, "(c)", ha="center", va="center", size=14 )
#subPlot3Axes.set_ylim( ymin = 0.0 )
# Fine-tune figure and hide x ticks for all but bottom plot.
figure.subplots_adjust( hspace = 0.35 )
pyplot.setp( [ a.get_xticklabels( ) for a in figure.axes[ :-1 ] ], visible = False )
pyplot.savefig( outputPlotPath + "/orbitalElementHistory_case" + str( caseNumber ) \
                + "_perturbers" + str( perturberPopulation ) + ".pdf", dpi = figureDPI )
if ( isShowFigures == 1 ):
    pyplot.show( )
pyplot.close( )

''' 
End of file
'''

endTime = time.time( )
print endTime - startTime
