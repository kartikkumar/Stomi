/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/comparisonFunctions.h>

#include "Stomi/Database/testParticleKick.h"

namespace stomi
{
namespace database
{

using namespace assist::basics;
using namespace tudat::basic_astrodynamics::orbital_element_conversions;

//! Default constructor, taking input values for all elements of kick.
TestParticleKick::TestParticleKick( const int aKickId,
                                    const int aTestParticleSimulationId,
                                    const double aConjunctionEpoch,
                                    const double aConjunctionDistance,
                                    const double aPreConjunctionEpoch,
                                    const double aPreConjunctionDistance,
                                    const tudat::basic_mathematics::Vector6d 
                                            aPreConjunctionStateInKeplerianElements,
                                    const double aPostConjunctionEpoch,
                                    const double aPostConjunctionDistance,
                                    const tudat::basic_mathematics::Vector6d 
                                            aPostConjunctionStateInKeplerianElements )
    : kickId( checkPositive( aKickId, "Kick ID" ) ),
      testParticleSimulationId( checkPositive( aTestParticleSimulationId, 
                                               "Test particle simulation ID" ) ),
      conjunctionEpoch( checkPositive( aConjunctionEpoch, "Conjunction epoch" ) ),
      conjunctionDistance( checkPositive( aConjunctionDistance, "Conjunction distance [m]" ) ),
      preConjunctionEpoch( checkLessThan( aPreConjunctionEpoch, "Pre-conjunction epoch [s]",
                                          conjunctionEpoch ) ),
      preConjunctionDistance( 
        checkGreaterThan( checkPositive( aPreConjunctionDistance, "Pre-conjunction distance [m]" ),
                          "Pre-conjunction distance [m]", conjunctionDistance ) ),
      preConjunctionStateInKeplerianElements( 
        ( Eigen::VectorXd( 6 ) 
            << checkPositive( aPreConjunctionStateInKeplerianElements( semiMajorAxisIndex ),
                              "Pre-conjunction semi-major axis [m]" ),
               checkInRange( aPreConjunctionStateInKeplerianElements( eccentricityIndex ),
                              "Pre-conjunction eccentricity [-]", 0.0, 1.0 ),
               checkPositive( aPreConjunctionStateInKeplerianElements( inclinationIndex ),
                              "Pre-conjunction inclination [rad]" ),
               aPreConjunctionStateInKeplerianElements.segment( 3, 3 ) ).finished( ) ),
      postConjunctionEpoch( checkGreaterThan( aPostConjunctionEpoch, "Post-conjunction epoch [s]",
                                              conjunctionEpoch ) ),
      postConjunctionDistance( checkGreaterThan(
                                 checkPositive( aPostConjunctionDistance, 
                                                "Post-conjunction distance [m]" ),
                                 "Post-conjunction distance [m]", conjunctionDistance ) ),
      postConjunctionStateInKeplerianElements( 
        ( Eigen::VectorXd( 6 ) 
            << checkPositive( aPostConjunctionStateInKeplerianElements( semiMajorAxisIndex ),
                              "Post-conjunction semi-major axis [m]" ),
               checkInRange( aPostConjunctionStateInKeplerianElements( eccentricityIndex ),
                              "Post-conjunction eccentricity [-]", 0.0, 1.0 ),
               checkPositive( aPostConjunctionStateInKeplerianElements( inclinationIndex ),
                              "Post-conjunction inclination [rad]" ),
               aPostConjunctionStateInKeplerianElements.segment( 3, 3 ) ).finished( ) )
{ }

//! Overload == operator.
bool operator==( const TestParticleKick& testParticleKick1,
                 const TestParticleKick& testParticleKick2 )
{
    return ( testParticleKick1.kickId == testParticleKick2.kickId
             && testParticleKick1.testParticleSimulationId == testParticleKick2.testParticleSimulationId
             && testParticleKick1.conjunctionEpoch == testParticleKick2.conjunctionEpoch
             && testParticleKick1.conjunctionDistance == testParticleKick2.conjunctionDistance
             && testParticleKick1.preConjunctionEpoch == testParticleKick2.preConjunctionEpoch
             && testParticleKick1.preConjunctionDistance
             == testParticleKick2.preConjunctionDistance
             && testParticleKick1.preConjunctionStateInKeplerianElements
             == testParticleKick2.preConjunctionStateInKeplerianElements
             && testParticleKick1.postConjunctionEpoch == testParticleKick2.postConjunctionEpoch
             && testParticleKick1.postConjunctionDistance
             == testParticleKick2.postConjunctionDistance
             && testParticleKick1.postConjunctionStateInKeplerianElements
             == testParticleKick2.postConjunctionStateInKeplerianElements );
}

//! Overload < operator.
bool operator<( const TestParticleKick& testParticleKick1,
                const TestParticleKick& testParticleKick2 )
{
    return testParticleKick1.conjunctionEpoch < testParticleKick2.conjunctionEpoch;
}

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, const TestParticleKick& testParticleKick )
{
    using std::endl;
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_astrodynamics::unit_conversions;
    using namespace assist::astrodynamics;

    // Write contents of TestParticleKickPointer object to output stream.
    outputStream << "Test particle kick ID: "
                 << testParticleKick.kickId << endl;
    outputStream << "Test particle simulation ID: "
                 << testParticleKick.testParticleSimulationId << endl;
    outputStream << "Conjunction epoch (since start) [Jyr]: "
                 << convertSecondsToJulianYears( testParticleKick.conjunctionEpoch ) << endl;
    outputStream << "Conjunction distance [m]: "
                 << testParticleKick.conjunctionDistance << endl;
    outputStream << "Pre-conjunction epoch (since start) [Jyr]: "
                 << convertSecondsToJulianYears( testParticleKick.preConjunctionEpoch )
                 << endl;
    outputStream << "Pre-conjunction distance [m]: "
                 << testParticleKick.preConjunctionDistance << endl;
    outputStream << "Pre-conjunction semi-major axis at T0 [m]: "
                 << testParticleKick.preConjunctionStateInKeplerianElements( 
                        semiMajorAxisIndex ) << endl;
    outputStream << "Pre-conjunction eccentricity at T0 [-]: "
                 << testParticleKick.preConjunctionStateInKeplerianElements( eccentricityIndex )
                 << endl;
    outputStream << "Pre-conjunction inclination at T0 [deg]: "
                 << convertRadiansToDegrees(
                        testParticleKick.preConjunctionStateInKeplerianElements(
                            inclinationIndex ) ) << endl;
    outputStream << "Pre-conjunction argument of periapsis at T0 [deg]: "
                 << convertRadiansToDegrees(
                        testParticleKick.preConjunctionStateInKeplerianElements(
                            argumentOfPeriapsisIndex ) ) << endl;
    outputStream << "Pre-conjunction longitude of ascending node at T0 [deg]: "
                 << convertRadiansToDegrees(
                        testParticleKick.preConjunctionStateInKeplerianElements(
                            longitudeOfAscendingNodeIndex ) ) << endl;
    outputStream << "Pre-conjunction true anomaly at T0 [deg]: "
                 << convertRadiansToDegrees(
                        testParticleKick.preConjunctionStateInKeplerianElements(
                            trueAnomalyIndex ) ) << endl;  
    outputStream << "Post-conjunction epoch (since start) [Jyr]: "
                 << convertSecondsToJulianYears( testParticleKick.postConjunctionEpoch )
                 << endl;
    outputStream << "Post-conjunction distance [m]: "
                 << testParticleKick.postConjunctionDistance << endl;
    outputStream << "Post-conjunction semi-major axis at T0 [m]: "
                 << testParticleKick.postConjunctionStateInKeplerianElements( 
                        semiMajorAxisIndex ) << endl;
    outputStream << "Post-conjunction eccentricity at T0 [-]: "
                 << testParticleKick.postConjunctionStateInKeplerianElements( eccentricityIndex )
                 << endl;
    outputStream << "Post-conjunction inclination at T0 [deg]: "
                 << convertRadiansToDegrees(
                        testParticleKick.postConjunctionStateInKeplerianElements(
                            inclinationIndex ) ) << endl;
    outputStream << "Post-conjunction argument of periapsis at T0 [deg]: "
                 << convertRadiansToDegrees(
                        testParticleKick.postConjunctionStateInKeplerianElements(
                            argumentOfPeriapsisIndex ) ) << endl;
    outputStream << "Post-conjunction longitude of ascending node at T0 [deg]: "
                 << convertRadiansToDegrees(
                        testParticleKick.postConjunctionStateInKeplerianElements(
                            longitudeOfAscendingNodeIndex ) ) << endl;
    outputStream << "Post-conjunction true anomaly at T0 [deg]: "
                 << convertRadiansToDegrees(
                        testParticleKick.postConjunctionStateInKeplerianElements(
                            trueAnomalyIndex ) ) << endl;                                          

    // Return output stream.
    return outputStream;
}

} // namespace database
} // namespace stomi
