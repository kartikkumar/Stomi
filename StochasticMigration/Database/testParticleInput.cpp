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
 *      130309    K. Kumar          Moved constructor implementation to source file and added
 *                                  sanity checks to ensure that test particle kick is valid; Moved
 *                                  all existing operator overloads to non-member functions and
 *                                  added new ones, including for pointer comparisons.
 *      130328    K. Kumar          Moved standard operator overload functions to Assist; used  
 *                                  comparison functions in Assist for checks in constructor.
 *
 *    References
 *
 *    Notes
 *
 */

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include <Assist/Basics/comparisonFunctions.h>

#include "StochasticMigration/Database/testParticleInput.h"

namespace stochastic_migration
{
namespace database
{

using namespace assist::basics;
using namespace tudat::basic_astrodynamics::orbital_element_conversions;

//! Default constructor, initializing class members.
TestParticleInput::TestParticleInput(
        const int aSimulationNumber,
        const bool aCompletedFlag,
        const tudat::basic_mathematics::Vector6d& anInitialStateInKeplerianElements )
    : simulationNumber( checkPositive( aSimulationNumber, "Simulation number" ) ),
      isCompleted( aCompletedFlag ),
      initialStateInKeplerianElements( 
        ( Eigen::VectorXd( 6 ) 
            << checkPositive( anInitialStateInKeplerianElements( semiMajorAxisIndex ),
                              "Initial semi-major axis [m]" ),
               checkInRange( anInitialStateInKeplerianElements( eccentricityIndex ), 
                             "Initial eccentricity [-]", 0.0, 1.0 ),
               checkPositive( anInitialStateInKeplerianElements( inclinationIndex ),
                              "Initial inclination [rad]" ),
               anInitialStateInKeplerianElements.segment( 3, 3 ) ).finished( ) )
{ }

//! Overload == operator.
bool operator==( const TestParticleInput& testParticleInput1,
                 const TestParticleInput& testParticleInput2 )
{
    return ( testParticleInput1.simulationNumber == testParticleInput2.simulationNumber
             && testParticleInput1.isCompleted == testParticleInput2.isCompleted
             && testParticleInput1.initialStateInKeplerianElements
             == testParticleInput2.initialStateInKeplerianElements );
}

//! Overload < operator.
bool operator<( const TestParticleInput& testParticleInput1,
                const TestParticleInput& testParticleInput2 )
{
    return testParticleInput1.simulationNumber < testParticleInput2.simulationNumber;
}

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, const TestParticleInput& testParticleInput )
{
    using namespace tudat::basic_astrodynamics::unit_conversions;

    // Set string indicating whether test particle simulation has been completed or not.
    std::string completedStatus;

    if ( testParticleInput.isCompleted == true )
    {
        completedStatus = "TRUE";
    }

    else
    {
        completedStatus = "FALSE";
    }

    // Write contents of TestParticleInput object to output stream.
    outputStream << "Test particle simulation number: "
                 << testParticleInput.simulationNumber << std::endl;
    outputStream << "Test particle simulation completed?: " << completedStatus << std::endl;
    outputStream << "Test particle's initial semi-major axis [m]: "
                 << testParticleInput.initialStateInKeplerianElements( semiMajorAxisIndex )
                 << std::endl;
    outputStream << "Test particle's initial eccentricity [-]: "
                 << testParticleInput.initialStateInKeplerianElements( eccentricityIndex )
                 << std::endl;
    outputStream << "Test particle's initial inclination [deg]: "
                 << convertRadiansToDegrees(
                        testParticleInput.initialStateInKeplerianElements( inclinationIndex ) )
                 << std::endl;
    outputStream << "Test particle's initial argument of periapsis [deg]: "
                 << convertRadiansToDegrees( testParticleInput.initialStateInKeplerianElements(
                                                 argumentOfPeriapsisIndex ) ) << std::endl;
    outputStream << "Test particle's initial longitude of ascending node [deg]: "
                 << convertRadiansToDegrees( testParticleInput.initialStateInKeplerianElements(
                                                 longitudeOfAscendingNodeIndex ) ) << std::endl;
    outputStream << "Test particle's initial true anomaly [deg]: "
                 << convertRadiansToDegrees(
                        testParticleInput.initialStateInKeplerianElements( trueAnomalyIndex ) )
                 << std::endl;

    // Return output stream.
    return outputStream;
}

} // namespace database
} // namespace stochastic_migration
