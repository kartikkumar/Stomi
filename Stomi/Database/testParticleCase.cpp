/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#define NEGATIVE_ZERO -std::numeric_limits< double >::min( )

#include <limits>
#include <sstream>
#include <stdexcept>

#include <boost/algorithm/string/predicate.hpp>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/comparisonFunctions.h>

#include "Stomi/Database/testParticleCase.h"

namespace stomi
{
namespace database
{

using namespace assist::basics;
using namespace tudat::basic_astrodynamics::orbital_element_conversions;

//! Constructor taking all test particle case data as input.
TestParticleCase::TestParticleCase(
    // Required parameters.
    const int aTestParticleCaseId,
    const std::string& aTestParticleCaseName,
    const double aRandomWalkSimulationPeriod,
    const double aCentralBodyGravitationalParameter,
    const double aPerturbedBodyRadius,
    const double aPerturbedBodyBulkDensity,
    const tudat::basic_mathematics::Vector6d& aPerturbedBodyStateInKeplerianElementsAtT0,
    const double aSemiMajorAxisDistributionLimit,
    // Optional parameters.
    const double aSynodicPeriodMaximum,
    const double aStartUpIntegrationPeriod,
    const double aCentralBodyJ2GravityCoefficient,
    const double aCentralBodyEquatorialRadius,            
    const double aConjunctionEventDetectionDistance,
    const double anOppositionEventDetectionDistance,
    const double anEccentricityDistributionMean,
    const double anEccentricityDistributionStandardDeviation,
    const double anInclinationDistributionMean,
    const double anInclinationDistributionStandardDeviation,
    const std::string& aNumericalIntegratorType,
    const double aNumericalIntegratorInitialStepSize,
    const double aNumericalIntegratorRelativeTolerance,
    const double aNumericalIntegratorAbsoluteTolerance )
        : testParticleCaseId( checkPositive( aTestParticleCaseId, "Test particle case ID" ) ),
          testParticleCaseName( aTestParticleCaseName ),
          randomWalkSimulationPeriod( 
            checkPositive( aRandomWalkSimulationPeriod, "Random walk period [s]" ) ),
          centralBodyGravitationalParameter( 
            checkPositive( aCentralBodyGravitationalParameter, 
                           "Central body gravitational parameter [m^3 s^-2]" ) ),
          perturbedBodyRadius( 
            checkPositive( aPerturbedBodyRadius, "Perturbed body radius [m]" ) ),
          perturbedBodyBulkDensity( 
            checkPositive( aPerturbedBodyBulkDensity, "Perturbed body bulk density [kg m^-3]" ) ),
          perturbedBodyStateInKeplerianElementsAtT0( 
            ( Eigen::VectorXd( 6 ) 
                << checkPositive( aPerturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ),
                                  "Perturbed body initial semi-major axis [m]" ),
                   checkPositive( aPerturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex ),
                                  "Perturbed body initial eccentricity [-]" ),
                   checkPositive( aPerturbedBodyStateInKeplerianElementsAtT0( inclinationIndex ),
                                  "Perturbed body initial inclination [rad]" ),
                   aPerturbedBodyStateInKeplerianElementsAtT0.segment( 3, 3 ) ).finished( ) ),
          semiMajorAxisDistributionLimit( 
            checkPositive( aSemiMajorAxisDistributionLimit,
                           "Semi-major axis distribution limit [m]" ) ),
          synodicPeriodMaximum( 
            checkPositive( aSynodicPeriodMaximum, "Synodic period limit [s]" ) ),
          startUpIntegrationPeriod( 
            checkPositive( aStartUpIntegrationPeriod, "Startup integration period [s]" ) ),
          centralBodyJ2GravityCoefficient( 
            checkGreaterThan( aCentralBodyJ2GravityCoefficient, 
                              "Central body J2 gravity coefficient [-]", NEGATIVE_ZERO ) ),
          centralBodyEquatorialRadius( 
            checkGreaterThan( aCentralBodyEquatorialRadius, 
                              "Central body equatorial radius [m]", NEGATIVE_ZERO ) ),
          conjunctionEventDetectionDistance( 
            checkPositive( aConjunctionEventDetectionDistance, 
                           "Conjunction event detection distance [m]" ) ),
          oppositionEventDetectionDistance( 
            checkPositive( anOppositionEventDetectionDistance,
                           "Opposition event detection distance [m]" ) ),
          eccentricityDistributionMean( 
            checkGreaterThan( anEccentricityDistributionMean, 
                              "Eccentricity distribution mean [-]", NEGATIVE_ZERO ) ),
          eccentricityDistributionStandardDeviation( 
            checkGreaterThan( anEccentricityDistributionStandardDeviation, 
                              "Eccentricity distribution standard deviation [-]",
                              NEGATIVE_ZERO ) ),
          inclinationDistributionMean( 
            checkGreaterThan( anInclinationDistributionMean, 
                              "Inclination distribution mean [rad]", NEGATIVE_ZERO ) ),
          inclinationDistributionStandardDeviation( 
            checkGreaterThan( anInclinationDistributionStandardDeviation, 
                              "Inclination distribution standard deviation [rad]",
                              NEGATIVE_ZERO ) ),
          numericalIntegratorInitialStepSize( checkPositive( 
            aNumericalIntegratorInitialStepSize, "Numerical integrator initial step size [s]" ) ),
          numericalIntegratorRelativeTolerance( 
            checkPositive( aNumericalIntegratorRelativeTolerance, 
                           "Numerical integrator relative tolerance [-]" ) ),
          numericalIntegratorAbsoluteTolerance( 
            checkPositive( aNumericalIntegratorAbsoluteTolerance, 
                           "Numerical integrator absolute tolerance [-]" ) )
{
    if ( testParticleCaseName.empty( ) )
    {
        throw std::runtime_error( "Test particle case name is empty!" );
    } 

    if ( aNumericalIntegratorType.empty( ) )
    {
        throw std::runtime_error( "Numerical integrator type is empty!" );
    }

    else if ( boost::iequals( aNumericalIntegratorType, "DOPRI853" ) )
    {
        numericalIntegratorType = DOPRI853;   
    } 

    else if ( boost::iequals( aNumericalIntegratorType, "RKF78" ) )
    {
        numericalIntegratorType = RKF78;   
    } 

    else
    {
        throw std::runtime_error( "Invalid numerical integrator type requested!" );
    }
}

//! Overload == operator.
bool operator==( const TestParticleCase& testParticleCase1,
                 const TestParticleCase& testParticleCase2 )
{
    return ( testParticleCase1.testParticleCaseId == testParticleCase2.testParticleCaseId
             && testParticleCase1.testParticleCaseName == testParticleCase2.testParticleCaseName
             && testParticleCase1.randomWalkSimulationPeriod
             == testParticleCase2.randomWalkSimulationPeriod
             && testParticleCase1.centralBodyGravitationalParameter
             == testParticleCase2.centralBodyGravitationalParameter
             && testParticleCase1.perturbedBodyRadius == testParticleCase2.perturbedBodyRadius
             && testParticleCase1.perturbedBodyBulkDensity 
             == testParticleCase2.perturbedBodyBulkDensity
             && testParticleCase1.perturbedBodyStateInKeplerianElementsAtT0
             == testParticleCase2.perturbedBodyStateInKeplerianElementsAtT0
             && testParticleCase1.semiMajorAxisDistributionLimit 
             == testParticleCase2.semiMajorAxisDistributionLimit
             && testParticleCase1.synodicPeriodMaximum == testParticleCase2.synodicPeriodMaximum
             && testParticleCase1.startUpIntegrationPeriod
             == testParticleCase2.startUpIntegrationPeriod
             && testParticleCase1.centralBodyJ2GravityCoefficient
             == testParticleCase2.centralBodyJ2GravityCoefficient
             && testParticleCase1.centralBodyEquatorialRadius
             == testParticleCase2.centralBodyEquatorialRadius             
             && testParticleCase1.conjunctionEventDetectionDistance
             == testParticleCase2.conjunctionEventDetectionDistance
             && testParticleCase1.oppositionEventDetectionDistance
             == testParticleCase2.oppositionEventDetectionDistance
             && testParticleCase1.eccentricityDistributionMean 
             == testParticleCase2.eccentricityDistributionMean
             && testParticleCase1.eccentricityDistributionStandardDeviation
             == testParticleCase2.eccentricityDistributionStandardDeviation
             && testParticleCase1.inclinationDistributionMean 
             == testParticleCase2.inclinationDistributionMean
             && testParticleCase1.inclinationDistributionStandardDeviation
             == testParticleCase2.inclinationDistributionStandardDeviation
             && testParticleCase1.numericalIntegratorType
             == testParticleCase2.numericalIntegratorType
             && testParticleCase1.numericalIntegratorInitialStepSize 
             == testParticleCase2.numericalIntegratorInitialStepSize
             && testParticleCase1.numericalIntegratorRelativeTolerance
             == testParticleCase2.numericalIntegratorRelativeTolerance
             && testParticleCase1.numericalIntegratorAbsoluteTolerance
             == testParticleCase2.numericalIntegratorAbsoluteTolerance );
}

//! Overload < operator.
bool operator<( const TestParticleCase& testParticleCase1,
                const TestParticleCase& testParticleCase2 )
{
    return testParticleCase1.testParticleCaseId < testParticleCase2.testParticleCaseId;
}

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, const TestParticleCase& testParticleCase )
{
    using std::endl;
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_astrodynamics::unit_conversions;
    using namespace assist::astrodynamics;
    
    // Write contents of TestParticleCase object to output stream.

    // Write required parameters to stream.
    outputStream << endl;
    outputStream << "********************************************************************" << endl;
    outputStream << "Required parameters" << endl;
    outputStream << "********************************************************************" << endl;
    outputStream << endl;

    outputStream << "Test particle case ID: " << testParticleCase.testParticleCaseId << endl;
    outputStream << "Test particle case: " << testParticleCase.testParticleCaseName << endl;    
    outputStream << "Random walk simulation period [Jyr]: "
                 << convertSecondsToJulianYears( 
                        testParticleCase.randomWalkSimulationPeriod ) << endl;
    outputStream << "Central body gravitational parameter [m^3 s^-2]: "
                 << testParticleCase.centralBodyGravitationalParameter << endl;
    outputStream << "Semi-major axis distribution limit [m]: "
                 << testParticleCase.semiMajorAxisDistributionLimit << endl;                 
    outputStream << "Perturbed body's radius [km]: "
                 << convertMetersToKilometers( testParticleCase.perturbedBodyRadius ) << endl;
    outputStream << "Perturbed body's bulk density [kg m^-3]: "
                 << convertMetersToKilometers( testParticleCase.perturbedBodyBulkDensity ) 
                 << endl;                 
    outputStream << "Perturbed body's semi-major axis at T0 [m]: "
                 << testParticleCase.perturbedBodyStateInKeplerianElementsAtT0( 
                        semiMajorAxisIndex ) << endl;
    outputStream << "Perturbed body's eccentricity at T0 [-]: "
                 << testParticleCase.perturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex )
                 << endl;
    outputStream << "Perturbed body's inclination at T0 [deg]: "
                 << convertRadiansToDegrees(
                        testParticleCase.perturbedBodyStateInKeplerianElementsAtT0(
                            inclinationIndex ) ) << endl;
    outputStream << "Perturbed body's argument of periapsis at T0 [deg]: "
                 << convertRadiansToDegrees(
                        testParticleCase.perturbedBodyStateInKeplerianElementsAtT0(
                            argumentOfPeriapsisIndex ) ) << endl;
    outputStream << "Perturbed body's longitude of ascending node at T0 [deg]: "
                 << convertRadiansToDegrees(
                        testParticleCase.perturbedBodyStateInKeplerianElementsAtT0(
                            longitudeOfAscendingNodeIndex ) ) << endl;
    outputStream << "Perturbed body's true anomaly at T0 [deg]: "
                 << convertRadiansToDegrees(
                        testParticleCase.perturbedBodyStateInKeplerianElementsAtT0(
                            trueAnomalyIndex ) ) << endl;                 

    // Write optional parameters to stream.
    outputStream << endl;
    outputStream << "********************************************************************" << endl;
    outputStream << "Optional parameters" << endl;
    outputStream << "********************************************************************" << endl;
    outputStream << endl;

    outputStream << "Synodic period maximum [Jyr]: "
                 << convertSecondsToJulianYears( testParticleCase.synodicPeriodMaximum )
                 << endl;
    outputStream << "Startup integration period [Jyr]: "
                 << convertSecondsToJulianYears( testParticleCase.startUpIntegrationPeriod )
                 << endl;
    outputStream << "Central body J2 gravity coefficient: "
                 << testParticleCase.centralBodyJ2GravityCoefficient << endl;
    outputStream << "Central body equatorial radius [m]: "
                 << testParticleCase.centralBodyEquatorialRadius << endl;                 
    outputStream << "Conjunction event detection distance [m]: "
                 << testParticleCase.conjunctionEventDetectionDistance << endl;
    outputStream << "Opposition event detection distance [m]: "
                 << testParticleCase.oppositionEventDetectionDistance << endl;
    outputStream << "Eccentricity distribution mean [-]: "
                 << testParticleCase.eccentricityDistributionMean << endl;
    outputStream << "Eccentricity distribution standard deviation [-]: "
                 << testParticleCase.eccentricityDistributionStandardDeviation << endl;
    outputStream << "Inclination distribution mean [deg]: "
                 << convertRadiansToDegrees( testParticleCase.inclinationDistributionMean ) 
                 << endl;           
    outputStream << "Inclination distribution standard deviation [deg]: "
                 << convertRadiansToDegrees(
                        testParticleCase.inclinationDistributionStandardDeviation )
                 << endl;

    outputStream << "Numerical integration type: ";
    if ( testParticleCase.numericalIntegratorType == DOPRI853 )
    {
        outputStream << "DOPRI853" << endl;
    }                 

    else if ( testParticleCase.numericalIntegratorType == RKF78 )
    {
        outputStream << "RKF78" << endl;
    }                 
                 
    outputStream << "Numerical integrator initial step size [s]: "
                 << testParticleCase.numericalIntegratorInitialStepSize << endl;                 
    outputStream << "Numerical integrator relative tolerance: "
                 << testParticleCase.numericalIntegratorRelativeTolerance << endl;
    outputStream << "Numerical integrator absolute tolerance: "
                 << testParticleCase.numericalIntegratorAbsoluteTolerance << endl;

    // Return output stream.
    return outputStream;
}

} // namespace database
} // namespace stomi
