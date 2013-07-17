/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
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
 *      130705    K. Kumar          Added caseId variable.
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
        const int aSimulationId,
        const int aCaseId,
        const bool aCompletedFlag,
        const tudat::basic_mathematics::Vector6d& anInitialStateInKeplerianElements )
    : simulationId( checkPositive( aSimulationId, "Simulation ID" ) ),
      caseId( checkPositive( aCaseId, "Case ID" ) ),
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
    return ( testParticleInput1.simulationId == testParticleInput2.simulationId
             && testParticleInput1.caseId == testParticleInput2.caseId
             && testParticleInput1.isCompleted == testParticleInput2.isCompleted
             && testParticleInput1.initialStateInKeplerianElements
             == testParticleInput2.initialStateInKeplerianElements );
}

//! Overload < operator.
bool operator<( const TestParticleInput& testParticleInput1,
                const TestParticleInput& testParticleInput2 )
{
    return testParticleInput1.simulationId < testParticleInput2.simulationId;
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
    outputStream << "Test particle simulation ID: "
                 << testParticleInput.simulationId << std::endl;
    outputStream << "Case ID: " << testParticleInput.caseId << std::endl;
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
