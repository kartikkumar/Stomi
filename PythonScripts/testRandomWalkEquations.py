'''   
Copyright (c) 2010-2013, Delft University of Technology
Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
All rights reserved.
See http://bit.ly/12SHPLR for license details.

This file can be used to compute values for unit-testing the executeKick() random walk function
in the Stochastic Migration project. 
'''


################################################################################################### 
# Import Python packages and modules
###################################################################################################

# Import necessary external packages.
import math
import numpy
import time

###################################################################################################


################################################################################################### 
# Set up input deck
###################################################################################################

# Perturbed body initial semi-major axis [m].
semiMajorAxisAtT0 = 97.736e6

# Perturbed body initial eccentricity [-].
eccentricityAtT0 = 2.54e-3

# Perturbed body initial inclination [rad].
inclinationAtT0 = numpy.deg2rad(0.14)

# Perturber semi-major axis before kick [m].
perturberSemiMajorAxisBeforeKick = 98.777e6

# Perturber eccentricity before kick [-].
perturberEccentricityBeforeKick = 2.0e-4

# Perturber inclination before kick [deg].
perturberInclinationBeforeKick = numpy.deg2rad(0.11)

# Perturber semi-major axis after kick [m].
perturberSemiMajorAxisAfterKick = 99.123e6

# Perturber eccentricity after kick [-].
perturberEccentricityAfterKick = 8.8e-4

# Perturber inclination after kick [rad].
perturberInclinationAfterKick = numpy.deg2rad(0.09)

# Mass ratio between perturber and perturbed body [-].
massRatio = 0.01

###################################################################################################


'''
    
                            DO NOT EDIT PARAMETERS BEYOND THIS POINT!!!

'''


################################################################################################### 
# Start timer
###################################################################################################

startTime = time.time()

###################################################################################################


################################################################################################### 
# Print input data for computation
###################################################################################################

# Set print format to scientific with 16 decimal places.
printFormat = '%.15e'

# Print statements.
print "Perturbed body semi-major axis before kick: " + str(printFormat % semiMajorAxisAtT0) \
      + " [m]"
print "Perturbed body eccentricity before kick: " + str(printFormat % eccentricityAtT0) + " [-]"
print "Perturbed body inclination before kick: " + str(printFormat % inclinationAtT0) + " [rad]"

print "Perturber semi-major axis before kick: " \
      + str(printFormat % perturberSemiMajorAxisBeforeKick) + " [m]"
print "Perturber eccentricity before kick: " \
      + str(printFormat % perturberEccentricityBeforeKick) + " [-]"
print "Perturber inclination before kick: " + str(printFormat % perturberInclinationBeforeKick) \
      + " [rad]"

print "Perturber semi-major axis after kick: " \
      + str(printFormat % perturberSemiMajorAxisAfterKick) + " [m]"
print "Perturber eccentricity after kick: " \
      + str(printFormat % perturberEccentricityAfterKick) + " [-]"
print "Perturber inclination after kick: " \
      + str(printFormat % perturberInclinationAfterKick) + " [rad]"

###################################################################################################


################################################################################################### 
# Compute perturbed body Keplerian elements after kick is executed
###################################################################################################

# Compute perturbed body semi-major axis after kick [m].
perturbedBodySemiMajorAxisAfterKick \
    = 1.0 / (1.0 / semiMajorAxisAtT0 \
             + massRatio * (1.0 / perturberSemiMajorAxisBeforeKick \
                            - 1.0 / perturberSemiMajorAxisAfterKick)) 

print "Perturbed body semi-major axis after kick is: " \
      + str(printFormat % perturbedBodySemiMajorAxisAfterKick) + " [m]"

# Compute perturbed body eccentricity after kick [-].
perturbedBodyEccentricityAfterKick = math.sqrt( \
    1.0 \
    - (math.sqrt(semiMajorAxisAtT0 * (1.0 - eccentricityAtT0**2.0)) \
       + massRatio * (math.sqrt(perturberSemiMajorAxisBeforeKick \
                                * (1.0 - perturberEccentricityBeforeKick**2.0)) \
                      - math.sqrt(perturberSemiMajorAxisAfterKick \
                                  * (1.0 - perturberEccentricityAfterKick**2.0))))**2.0 \
         / perturbedBodySemiMajorAxisAfterKick)      

print "Perturbed body eccentricity after kick is: " \
      + str(printFormat % perturbedBodyEccentricityAfterKick) + " [-]"

# Compute perturbed body inclination after kick [rad].
perturbedBodyInclinationAfterKick = math.acos( \
    (math.sqrt(semiMajorAxisAtT0 * (1.0 - eccentricityAtT0**2.0)) * math.cos(inclinationAtT0) \
     + massRatio * (math.sqrt(perturberSemiMajorAxisBeforeKick \
                              * (1.0 - perturberEccentricityBeforeKick**2.0)) \
                    * math.cos(perturberInclinationBeforeKick) \
                    - math.sqrt(perturberSemiMajorAxisAfterKick \
                                * (1.0 - perturberEccentricityAfterKick**2.0)) \
                    * math.cos(perturberInclinationAfterKick))) \
    / math.sqrt(perturbedBodySemiMajorAxisAfterKick \
                * (1.0 - perturbedBodyEccentricityAfterKick**2.0)))

print "Perturbed body inclination after kick is: " \
      + str(printFormat % perturbedBodyInclinationAfterKick) + " [rad]"                      

###################################################################################################


###################################################################################################
# Finalize timer and print elapsed time.
###################################################################################################

# Finalize timer.
endTime = time.time()

# Print elapsed time for script [s].
print "This script took " + str(endTime - startTime) + "s"

###################################################################################################