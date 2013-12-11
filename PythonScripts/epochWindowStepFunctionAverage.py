'''   
Copyright (c) 2010-2013, Delft University of Technology
Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
All rights reserved.
See http://bit.ly/12SHPLR for license details.

This script can be used to compute the weighted average of a step function within an epoch window.
'''

def computeStepFunctionAverage(independentVariable, dependentVariable, lowerBound, upperBound):
    for i in xrange(0,len(independentVariable)):
        if independentVariable[i] > lowerBound:
            lowerBoundIndex = i
            break

    upperBoundIndex = len(independentVariable) - 1
    for i in xrange(0,len(independentVariable)):
        if independentVariable[i] > upperBound:
            upperBoundIndex = i-1
            break

    firstMoment = 0.0
    for i in xrange(lowerBoundIndex, upperBoundIndex):
        firstMoment += dependentVariable[i] * (independentVariable[i+1] - independentVariable[i])

    firstMoment += dependentVariable[lowerBoundIndex-1] \
                   * (independentVariable[lowerBoundIndex]-lowerBound)
    firstMoment += dependentVariable[upperBoundIndex] \
                   * (upperBound-independentVariable[upperBoundIndex])

    return firstMoment / (upperBound-lowerBound)