/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130218    K. Kumar          File created.
 *      130308    K. Kumar          Moved constructor implementation to source file and added
 *                                  sanity checks to ensure that test particle kick is valid.
 *      130309    K. Kumar          Moved all existing operator overloads to non-member functions
 *                                  and added new ones, including for pointer comparisons.
 *      130328    K. Kumar          Moved standard operator overload functions to Assist; used  
 *                                  comparison functions in Assist for checks in constructor.
 *      130708    K. Kumar          Added kick ID, Tisserance parameter, and energy & angular 
 *                                  momentum relative errors; removed mass ratio.
 *
 *    References
 *
 *    Notes
 *
 */

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/comparisonFunctions.h>

#include "StochasticMigration/Database/testParticleKick.h"

namespace stochastic_migration
{
namespace database
{

using namespace assist::basics;

//! Default constructor, taking input values for all elements of kick.
TestParticleKick::TestParticleKick( const int aKickId,
                                    const int aSimulationNumber,
                                    const double aConjunctionEpoch,
                                    const double aConjunctionDistance,
                                    const double aPreConjunctionEpoch,
                                    const double aPreConjunctionDistance,
                                    const double aPreConjunctionSemiMajorAxis,
                                    const double aPreConjunctionEccentricity,
                                    const double aPreConjunctionInclination,
                                    const double aPostConjunctionEpoch,
                                    const double aPostConjunctionDistance,
                                    const double aPostConjunctionSemiMajorAxis,
                                    const double aPostConjunctionEccentricity,
                                    const double aPostConjunctionInclination )
    : kickId( checkPositive( aKickId, "Kick ID" ) ),
      simulationNumber( checkPositive( aSimulationNumber, "Simulation number" ) ),
      conjunctionEpoch( checkPositive( aConjunctionEpoch, "Conjunction epoch" ) ),
      conjunctionDistance( checkPositive( aConjunctionDistance, "Conjunction distance [m]" ) ),
      preConjunctionEpoch( checkLessThan( aPreConjunctionEpoch, "Pre-conjunction epoch [s]",
                                          conjunctionEpoch ) ),
      preConjunctionDistance( 
        checkGreaterThan( checkPositive( aPreConjunctionDistance, "Pre-conjunction distance [m]" ),
                          "Pre-conjunction distance [m]", conjunctionDistance ) ),
      preConjunctionSemiMajorAxis( 
        checkPositive( aPreConjunctionSemiMajorAxis, "Pre-conjunction semi-major axis [m]" ) ),
      preConjunctionEccentricity( checkInRange( aPreConjunctionEccentricity,
                                                "Pre-conjunction eccentricity [-]", 0.0, 1.0 ) ),
      preConjunctionInclination( checkPositive( aPreConjunctionInclination,
                                                "Pre-conjunction inclination [rad]" ) ),
      postConjunctionEpoch( checkGreaterThan( aPostConjunctionEpoch, "Post-conjunction epoch [s]",
                                              conjunctionEpoch ) ),
      postConjunctionDistance( checkGreaterThan(
                                 checkPositive( aPostConjunctionDistance, 
                                                "Post-conjunction distance [m]" ),
                                 "Post-conjunction distance [m]", conjunctionDistance ) ),
      postConjunctionSemiMajorAxis( checkPositive( aPostConjunctionSemiMajorAxis,
                                                   "Post-conjunction semi-major axis [m]" ) ),
      postConjunctionEccentricity( checkInRange( aPostConjunctionEccentricity,
                                                 "Post-conjunction eccentricity [-]", 0.0, 1.0 ) ),
      postConjunctionInclination( checkPositive( aPostConjunctionInclination,
                                                 "Post-conjunction inclination [rad]" ) )

{ }

//! Overload == operator.
bool operator==( const TestParticleKick& testParticleKick1,
                 const TestParticleKick& testParticleKick2 )
{
    return ( testParticleKick1.kickId == testParticleKick2.kickId
             && testParticleKick1.simulationNumber == testParticleKick2.simulationNumber
             && testParticleKick1.conjunctionEpoch == testParticleKick2.conjunctionEpoch
             && testParticleKick1.conjunctionDistance == testParticleKick2.conjunctionDistance
             && testParticleKick1.preConjunctionEpoch == testParticleKick2.preConjunctionEpoch
             && testParticleKick1.preConjunctionDistance
             == testParticleKick2.preConjunctionDistance
             && testParticleKick1.preConjunctionSemiMajorAxis
             == testParticleKick2.preConjunctionSemiMajorAxis
             && testParticleKick1.preConjunctionEccentricity
             == testParticleKick2.preConjunctionEccentricity
             && testParticleKick1.preConjunctionInclination
             == testParticleKick2.preConjunctionInclination
             && testParticleKick1.postConjunctionEpoch == testParticleKick2.postConjunctionEpoch
             && testParticleKick1.postConjunctionDistance
             == testParticleKick2.postConjunctionDistance
             && testParticleKick1.postConjunctionSemiMajorAxis
             == testParticleKick2.postConjunctionSemiMajorAxis
             && testParticleKick1.postConjunctionEccentricity
             == testParticleKick2.postConjunctionEccentricity
             && testParticleKick1.postConjunctionInclination
             == testParticleKick2.postConjunctionInclination );
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
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_astrodynamics::unit_conversions;
    using namespace assist::astrodynamics;

    // Write contents of TestParticleKickPointer object to output stream.
    outputStream << "Test particle kick ID: "
                 << testParticleKick.kickId << std::endl;
    outputStream << "Test particle simulation number: "
                 << testParticleKick.simulationNumber << std::endl;
    outputStream << "Conjunction epoch (since start) [Jyr]: "
                 << convertSecondsToJulianYears( testParticleKick.conjunctionEpoch ) << std::endl;
    outputStream << "Conjunction distance [m]: "
                 << testParticleKick.conjunctionDistance << std::endl;
    outputStream << "Pre-conjunction epoch (since start) [Jyr]: "
                 << convertSecondsToJulianYears( testParticleKick.preConjunctionEpoch )
                 << std::endl;
    outputStream << "Pre-conjunction distance [m]: "
                 << testParticleKick.preConjunctionDistance << std::endl;
    outputStream << "Pre-conjunction semi-major axis [m]: "
                 << testParticleKick.preConjunctionSemiMajorAxis << std::endl;
    outputStream << "Pre-conjunction eccentricity [-]: "
                 << testParticleKick.preConjunctionEccentricity << std::endl;
    outputStream << "Pre-conjunction inclination [deg]: "
                 << convertRadiansToDegrees( testParticleKick.preConjunctionInclination )
                 << std::endl;
    outputStream << "Post-conjunction epoch (since start) [Jyr]: "
                 << convertSecondsToJulianYears( testParticleKick.postConjunctionEpoch )
                 << std::endl;
    outputStream << "Post-conjunction distance [m]: "
                 << testParticleKick.postConjunctionDistance << std::endl;
    outputStream << "Post-conjunction semi-major axis [m]: "
                 << testParticleKick.postConjunctionSemiMajorAxis << std::endl;
    outputStream << "Post-conjunction eccentricity [-]: "
                 << testParticleKick.postConjunctionEccentricity << std::endl;
    outputStream << "Post-conjunction inclination [deg]: "
                 << convertRadiansToDegrees( testParticleKick.postConjunctionInclination )
                 << std::endl;                                              

    // Return output stream.
    return outputStream;
}

} // namespace database
} // namespace stochastic_migration
