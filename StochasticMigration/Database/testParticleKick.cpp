/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130218    K. Kumar          File created.
 *      130308    K. Kumar          Moved constructor implementation to source file and added
 *                                  sanity checks to ensure that test particle kick is valid.
 *      130309    K. Kumar          Moved all existing operator overloads to non-member functions
 *                                  and added new ones, including for pointer comparisons.
 *
 *    References
 *
 *    Notes
 *
 */

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/basics.h>
#include <Assist/Basics/comparisonFunctions.h>

#include "StochasticMigration/Database/testParticleKick.h"

namespace stochastic_migration
{
namespace database
{

using namespace assist::basics;

//! Default constructor, taking input values for all elements of kick.
TestParticleKick::TestParticleKick( const int aSimulationNumber,
                                    const double aConjunctionEpoch,
                                    const double aConjunctionDistance,
                                    const double aConjunctionDuration,
                                    const double aPreConjunctionEpoch,
                                    const double aPreConjunctionDistance,
                                    const double aPreConjunctionSemiMajorAxis,
                                    const double aPreConjunctionEccentricity,
                                    const double aPreConjunctionInclination,
                                    const double aPostConjunctionEpoch,
                                    const double aPostConjunctionDistance,
                                    const double aPostConjunctionSemiMajorAxis,
                                    const double aPostConjunctionEccentricity,
                                    const double aPostConjunctionInclination,
                                    const double aMassRatio )
    : simulationNumber( checkPositive( aSimulationNumber, "Simulation number" ) ),
      conjunctionEpoch( aConjunctionEpoch ),
      conjunctionDistance( checkPositive( aConjunctionDistance, "Conjunction distance [m]" ) ),
      conjunctionDuration( checkLessThan( 
                            checkPositive( aConjunctionDuration, "Conjunction duration [s]" ),
                            "Conjunction duration [s]", 
                            aPostConjunctionEpoch - aPreConjunctionEpoch ) ),
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
                                                 "Post-conjunction inclination [rad]" ) ),
      massRatio( checkInRange( aMassRatio, "Mass ratio [-]", 0.0, 1.0 ) )
{ }

//! Overload == operator.
bool operator==( const TestParticleKick& testParticleKick1,
                 const TestParticleKick& testParticleKick2 )
{
    return ( testParticleKick1.simulationNumber == testParticleKick2.simulationNumber
             && testParticleKick1.conjunctionEpoch == testParticleKick2.conjunctionEpoch
             && testParticleKick1.conjunctionDistance == testParticleKick2.conjunctionDistance
             && testParticleKick1.conjunctionDuration == testParticleKick2.conjunctionDuration
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
             == testParticleKick2.postConjunctionInclination
             && testParticleKick1.massRatio == testParticleKick2.massRatio );
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
    using namespace assist::basics;
    using namespace assist::astrodynamics;

    // Write contents of TestParticleKickPointer object to output stream.
    outputStream << "Test particle simulation number: "
                 << testParticleKick.simulationNumber << std::endl;
    outputStream << "Conjunction epoch (since start) [Jyr]: "
                 << convertSecondsToJulianYears( testParticleKick.conjunctionEpoch ) << std::endl;
    outputStream << "Conjunction distance [m]: "
                 << testParticleKick.conjunctionDistance << std::endl;
    outputStream << "Conjunction duration [Jdays]: "
                 << convertSecondsToJulianYears( testParticleKick.conjunctionDuration )
                 << std::endl;
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
    outputStream << "Mass ratio [-]: " << testParticleKick.massRatio << std::endl;

    // Return output stream.
    return outputStream;
}

} // namespace database
} // namespace stochastic_migration
