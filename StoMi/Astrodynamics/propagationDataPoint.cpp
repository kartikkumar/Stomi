/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */


#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/comparisonFunctions.h>

#include "StoMi/Astrodynamics/propagationDataPoint.h"

namespace stomi
{
namespace astrodynamics
{

//! Overload assignment-operator.
PropagationDataPoint& PropagationDataPoint::operator=( 
    const PropagationDataPoint& sourceDataPoint )
{
    // Check for self-assignment by comparing the address of the
    // implicit object and the parameter.
    if ( this == &sourceDataPoint )
    {
        return *this;
    }
 
    // Else, do the copy.
    epoch = sourceDataPoint.epoch;
    mutualDistance = sourceDataPoint.mutualDistance;
    testParticleStateInKeplerianElements 
        = sourceDataPoint.testParticleStateInKeplerianElements;
    perturbedBodyStateInKeplerianElements 
        = sourceDataPoint.perturbedBodyStateInKeplerianElements;

    // Return the existing object.
    return *this;
}

//! Overload == operator.
bool operator==( const PropagationDataPoint& propagationDataPoint1,
                 const PropagationDataPoint& propagationDataPoint2 )
{
    return ( propagationDataPoint1.epoch == propagationDataPoint2.epoch
             && propagationDataPoint1.mutualDistance == propagationDataPoint2.mutualDistance
             && propagationDataPoint1.testParticleStateInKeplerianElements 
             == propagationDataPoint2.testParticleStateInKeplerianElements 
             && propagationDataPoint1.perturbedBodyStateInKeplerianElements 
             == propagationDataPoint2.perturbedBodyStateInKeplerianElements );
}

//! Overload < operator.
bool operator<( const PropagationDataPoint& propagationDataPoint1,
                const PropagationDataPoint& propagationDataPoint2 )
{
    return propagationDataPoint1.epoch < propagationDataPoint2.epoch;
}

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, const PropagationDataPoint& PropagationDataPoint )
{
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_astrodynamics::unit_conversions;
    using namespace assist::astrodynamics;

    // Write contents of PropagationDataPointPointer object to output stream.
    outputStream << "Event epoch (since start): "
                 << convertSecondsToJulianYears( PropagationDataPoint.epoch ) << std::endl;
    outputStream << "Event distance [km]: "
                 << convertMetersToKilometers( PropagationDataPoint.mutualDistance ) << std::endl;
    outputStream << "Test particle semi-major axis at event [km]: "
                 << convertMetersToKilometers(
                        PropagationDataPoint.testParticleStateInKeplerianElements( 
                            semiMajorAxisIndex ) )
                 << std::endl;
    outputStream << "Test particle eccentricity at event [-]: "
                 << PropagationDataPoint.testParticleStateInKeplerianElements( 
                        eccentricityIndex ) << std::endl;
    outputStream << "Test particle inclination at event [deg]: "
                 << convertRadiansToDegrees(
                        PropagationDataPoint.testParticleStateInKeplerianElements(
                            inclinationIndex ) )
                 << std::endl;                                  
    outputStream << "Test particle argument of periapsis at event [deg]: "
                 << convertRadiansToDegrees(
                        PropagationDataPoint.testParticleStateInKeplerianElements( 
                            argumentOfPeriapsisIndex ) )
                 << std::endl;
    outputStream << "Test particle longitude of ascending node at event [deg]: "
                 << convertRadiansToDegrees(
                        PropagationDataPoint.testParticleStateInKeplerianElements( 
                            longitudeOfAscendingNodeIndex ) )
                 << std::endl;
    outputStream << "Test particle true anomaly at event [deg]: "
                 << convertRadiansToDegrees(
                        PropagationDataPoint.testParticleStateInKeplerianElements( 
                            trueAnomalyIndex ) )
                 << std::endl;
    outputStream << "Perturbed body semi-major axis at event [km]: "
                 << convertMetersToKilometers(
                        PropagationDataPoint.perturbedBodyStateInKeplerianElements( 
                            semiMajorAxisIndex ) )
                 << std::endl;
    outputStream << "Perturbed body eccentricity at event [-]: "
                 << PropagationDataPoint.perturbedBodyStateInKeplerianElements( 
                        eccentricityIndex ) << std::endl;
    outputStream << "Perturbed body inclination at event [deg]: "
                 << convertRadiansToDegrees(
                        PropagationDataPoint.perturbedBodyStateInKeplerianElements( 
                            inclinationIndex ) )
                 << std::endl;                                  
    outputStream << "Perturbed body argument of periapsis at event [deg]: "
                 << convertRadiansToDegrees(
                        PropagationDataPoint.perturbedBodyStateInKeplerianElements( 
                            argumentOfPeriapsisIndex ) )
                 << std::endl;
    outputStream << "Perturbed body longitude of ascending node at event [deg]: "
                 << convertRadiansToDegrees(
                        PropagationDataPoint.perturbedBodyStateInKeplerianElements( 
                            longitudeOfAscendingNodeIndex ) )
                 << std::endl;
    outputStream << "Perturbed body true anomaly at event [deg]: "
                 << convertRadiansToDegrees(
                        PropagationDataPoint.perturbedBodyStateInKeplerianElements( 
                            trueAnomalyIndex ) )
                 << std::endl;                 

    // Return output stream.
    return outputStream;
}    

} // namespace astrodynamics
} // namespace stomi