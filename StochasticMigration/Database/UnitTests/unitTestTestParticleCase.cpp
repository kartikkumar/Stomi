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
 *      130329    K. Kumar          File created.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <iostream>
#include <stdexcept>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Basics/testMacros.h>

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include <Assist/Basics/operatorOverloadFunctions.h>

#include "StochasticMigration/Database/testParticleCase.h"

namespace stochastic_migration
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics::orbital_element_conversions;

//! Test fixture used to test the TestParticleCase struct.
struct TestParticleCaseFixture
{
public:

    // Set Eigen macro to correctly align class with fixed-size vectorizable types.
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! ructor initializing valid parameters.
    TestParticleCaseFixture( )
        : caseNumber( 1 ),
          randomWalkSimulationDuration( 1.2345e6 ),
          synodicPeriodLimit( 1.23e4 ),
          outputInterval( 3600.0 ),
          startUpIntegrationDuration( 0.0 ),
          conjunctionEventDetectionDistance( 1.234e5 ),
          oppositionEventDetectionDistance( 1.234e5 ),
          centralBodyGravitationalParameter( 2.345e10 ),
          centralBodyJ2GravityCoefficient( 0.0 ),
          centralBodyEquatorialRadius( 0.0 ),
          semiMajorAxisLimit( 1000.0e3 ),
          meanEccentricity( 1.0e-3 ),
          fullWidthHalfMaxmimumEccentricityDistribution( 1.0e-4 ),
          meanInclination( 1.0e-3 ),
          fullWidthHalfMaxmimumInclinationDistribution( 1.0e-4 ),          
          perturbedBodyGravitationalParameter( 2.345e6 ),
          perturbedBodyStateInKeplerianElementsAtT0(
              ( Eigen::VectorXd( 6 ) << 1234.5, 0.55, 0.975, 1.234, 4.567, 3.456 ).finished( ) ),
          numericalIntegratorType( "RKF78" ),
          numericalIntegratorRelativeTolerance( 1.0e-12 ),
          numericalIntegratorAbsoluteTolerance( 1.0e-15 ),
          initialStepSize( 60.0 )
    { }

    //! Declare parameters of test particle case.
    //! Case number.
    int caseNumber;

    //! Random walk simulation duration [s].
    double randomWalkSimulationDuration;

    //! Synodic period limit [s].
    double synodicPeriodLimit;

    //! Output interval, determining output frequency for data files [s].
    double outputInterval;

    //! Startup integration duration [s].
    double startUpIntegrationDuration;

    //! Mutual distance used to detect start and end of conjunction events [m].
    double conjunctionEventDetectionDistance;

    //! Distance used to detect start and end of opposition events [m].
    double oppositionEventDetectionDistance;

    //! Central body gravitational parameter [m^3 s^-2].
    double centralBodyGravitationalParameter;

    //! Central body J2 gravity field coefficient.
    double centralBodyJ2GravityCoefficient;

    //! Central body equatorial radius [m].
    double centralBodyEquatorialRadius;

    //! Limits on maximum semi-major axis values wrt perturbed body [m].
    double semiMajorAxisLimit;

    //! Mean eccentricity value for distribution.
    double meanEccentricity;

    //! FWHM eccentricity value for distribution.
    double fullWidthHalfMaxmimumEccentricityDistribution;

    //! Mean inclination value for distribution.
    double meanInclination;

    //! FWHM inclination value for distribution.
    double fullWidthHalfMaxmimumInclinationDistribution;

    //! Perturbed gravitational parameter [m^3 s^-2].
    double perturbedBodyGravitationalParameter;

    //! Perturbed body state in Keplerian elements at T0.
    tudat::basic_mathematics::Vector6d perturbedBodyStateInKeplerianElementsAtT0;

    //! Numerical integrator type.
    std::string numericalIntegratorType;

    //! Relative tolerance for numerical integrator.
    double numericalIntegratorRelativeTolerance;

    //! Absolute tolerance for numerical integrator.
    double numericalIntegratorAbsoluteTolerance;

    //! Initial step size for numerical integrator.
    double initialStepSize;

    //! Get test particle case created from specified parameters.
    database::TestParticleCasePointer getTestParticleCase( )
    {
        return boost::make_shared< database::TestParticleCase >(
            database::TestParticleCase( 
                caseNumber, randomWalkSimulationDuration, synodicPeriodLimit,
                outputInterval, startUpIntegrationDuration, conjunctionEventDetectionDistance,
                oppositionEventDetectionDistance, centralBodyGravitationalParameter,
                centralBodyJ2GravityCoefficient, centralBodyEquatorialRadius, semiMajorAxisLimit,
                meanEccentricity, fullWidthHalfMaxmimumEccentricityDistribution, meanInclination,
                fullWidthHalfMaxmimumInclinationDistribution, perturbedBodyGravitationalParameter,
                perturbedBodyStateInKeplerianElementsAtT0, numericalIntegratorType,
                numericalIntegratorRelativeTolerance, numericalIntegratorAbsoluteTolerance,
                initialStepSize ) );
    }

protected:

private:
};    

BOOST_FIXTURE_TEST_SUITE( test_test_particle_case, TestParticleCaseFixture )

//! Test correct construction of test particle case.
BOOST_AUTO_TEST_CASE( testTestParticleCaseStructContruction )
{
    // Create test particle case.
    database::TestParticleCasePointer testParticleCase = getTestParticleCase( );

    // Check that the test particle case created contains all the data as required.
    BOOST_CHECK_EQUAL( testParticleCase->caseNumber, caseNumber );     
    BOOST_CHECK_EQUAL( testParticleCase->randomWalkSimulationDuration, 
                       randomWalkSimulationDuration );
    BOOST_CHECK_EQUAL( testParticleCase->synodicPeriodLimit, synodicPeriodLimit );
    BOOST_CHECK_EQUAL( testParticleCase->outputInterval, outputInterval );
    BOOST_CHECK_EQUAL( testParticleCase->startUpIntegrationDuration, startUpIntegrationDuration );
    BOOST_CHECK_EQUAL( testParticleCase->conjunctionEventDetectionDistance, 
                       conjunctionEventDetectionDistance );
    BOOST_CHECK_EQUAL( testParticleCase->oppositionEventDetectionDistance, 
                       oppositionEventDetectionDistance );
    BOOST_CHECK_EQUAL( testParticleCase->centralBodyGravitationalParameter, 
                       centralBodyGravitationalParameter );
    BOOST_CHECK_EQUAL( testParticleCase->centralBodyJ2GravityCoefficient, 
                       centralBodyJ2GravityCoefficient );
    BOOST_CHECK_EQUAL( testParticleCase->centralBodyEquatorialRadius, 
                       centralBodyEquatorialRadius );
    BOOST_CHECK_EQUAL( testParticleCase->semiMajorAxisLimit, semiMajorAxisLimit );
    BOOST_CHECK_EQUAL( testParticleCase->meanEccentricity, meanEccentricity );
    BOOST_CHECK_EQUAL( testParticleCase->fullWidthHalfMaxmimumEccentricityDistribution,
                       fullWidthHalfMaxmimumEccentricityDistribution );
    BOOST_CHECK_EQUAL( testParticleCase->fullWidthHalfMaxmimumInclinationDistribution, 
                       fullWidthHalfMaxmimumInclinationDistribution );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyGravitationalParameter, 
                       perturbedBodyGravitationalParameter );

    {
        TUDAT_CHECK_MATRIX_BASE( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0,
                                 perturbedBodyStateInKeplerianElementsAtT0 )
                BOOST_CHECK_EQUAL(
                    testParticleCase->perturbedBodyStateInKeplerianElementsAtT0.coeff( row, col ),
                    perturbedBodyStateInKeplerianElementsAtT0.coeff( row, col ) );
    }

    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorType, numericalIntegratorType );
    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorRelativeTolerance, 
                       numericalIntegratorRelativeTolerance );
    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorAbsoluteTolerance, 
                       numericalIntegratorAbsoluteTolerance );
    BOOST_CHECK_EQUAL( testParticleCase->initialStepSize, initialStepSize );
}

//! Test initialization of test particle case with non-positive case number.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveCaseNumberError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set case number to invalid (non-positive) number.
    caseNumber = -1;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive random walk simulation duration.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveRandomWalkSimulationDurationError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set random walk duration to invalid (non-positive) value.
    randomWalkSimulationDuration = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive synodic period limit.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveSynodicPeriodLimitError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set synodic period limit to invalid (non-positive) value.
    synodicPeriodLimit = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive output interval.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveOutputIntervalError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set output interval to invalid (non-positive) value.
    outputInterval = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive start-up integration duration.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveStartUpIntegrationDurationError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set start-up integration duration to invalid (non-positive) value.
    startUpIntegrationDuration = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive conjunction event detection. 
//! distance.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveConjunctionEventDetectionDistanceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set conjunction event detection distance to invalid (non-positive) value.
    conjunctionEventDetectionDistance = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive opposition event detection.
//! distance.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveOppositionEventDetectionDistanceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set opposition event detection distance to invalid (non-positive) value.
    oppositionEventDetectionDistance = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive central body gravitational 
//! parameter.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveCentralBodyGravitationalParameterError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set central body gravitational parameter to invalid (non-positive) value.
    centralBodyGravitationalParameter = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive central body J2 gravity 
//! coefficient. 
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveCentralBodyJ2GravityCoefficientError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set central body J2 gravity coefficient to invalid (non-positive) value.
    centralBodyJ2GravityCoefficient = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive central body equatorial radius.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveCentralBodyEquatorialRadiusError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set central body equatorial radius to invalid (non-positive) value.
    centralBodyEquatorialRadius = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive semi-major axis limit.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveSemiMajorAxisLimitError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set semi-major axis limit to invalid (non-positive) value.
    semiMajorAxisLimit = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive mean eccentricity.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveMeanEccentricityError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set mean eccentricity to invalid (non-positive) value.
    meanEccentricity = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive FWHM eccentricity distribution.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveFWHMEccentricityDistributionError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set FWHM eccentricity distribution to invalid (non-positive) value.
    fullWidthHalfMaxmimumEccentricityDistribution = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive mean inclination.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveMeanInclinationError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set mean inclination to invalid (non-positive) value.
    meanInclination = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive FWHM inclination distribution.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveFWHMInclinationDistributionError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set FWHM inclination distribution to invalid (non-positive) value.
    fullWidthHalfMaxmimumInclinationDistribution = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive perturbed body gravitational 
//! parameter.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositivePerturbedBodyGravitationalParameterError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set perturbed body gravitational parameter to invalid (non-positive) value.
    perturbedBodyGravitationalParameter = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive perturbed body initial semi-major 
//! axis.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositivePerturbedBodyInitialSemiMajorAxisError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set perturbed body initial semi-major axis to invalid (non-positive) value.
    perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-elliptical perturbed body initial
//! eccentricity axis.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositivePerturbedBodyInitialEccentricityError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set perturbed body initial eccentricity to invalid (non-elliptical) value.
    perturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex ) = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive perturbed body initial inclination
//! axis.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositivePerturbedBodyInitialInclinationError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set perturbed body initial inclination to invalid (non-positive) value.
    perturbedBodyStateInKeplerianElementsAtT0( inclinationIndex ) = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with empty numerical integrator type.
BOOST_AUTO_TEST_CASE( testTestParticleCaseEmptyNumericalIntegratorTypeError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set numerical integrator type to invalid (non-positive) value.
    numericalIntegratorType = "";

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive numerical integrator relative 
//! tolerance.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveNumericalIntegratorRelativeToleranceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set numerical integrator relative tolerance to invalid (non-positive) value.
    numericalIntegratorRelativeTolerance = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case with non-positive numerical integrator absolute 
//! tolerance.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonPositiveNumericalIntegratorAbsoluteToleranceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set numerical integrator absolute tolerance to invalid (non-positive) value.
    numericalIntegratorAbsoluteTolerance = -1.0;

    // Try to create test particle case.
    try { database::TestParticleCasePointer testParticleCase = getTestParticleCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test comparison of TestParticleCase pointers using overloaded == operator.
BOOST_AUTO_TEST_CASE( testTestParticleCaseEqualComparison )
{
    // Create TestParticleCase pointers.
    database::TestParticleCasePointer testParticleCase = getTestParticleCase( );
    database::TestParticleCasePointer testParticleCase2 = getTestParticleCase( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleCase == *testParticleCase2 );
}

//! Test comparison of TestParticleCase pointers using overloaded != operator.
BOOST_AUTO_TEST_CASE( testTestParticleCaseNonEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleCase pointers.
    database::TestParticleCasePointer testParticleCase = getTestParticleCase( );

    caseNumber = 2;
    database::TestParticleCasePointer testParticleCase2 = getTestParticleCase( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleCase != *testParticleCase2 );
}

//! Test comparison of TestParticleCase pointers using overloaded < operator.
BOOST_AUTO_TEST_CASE( testTestParticleCaseLessThanComparison )
{
    // Create TestParticleCase pointers.
    database::TestParticleCasePointer testParticleCase = getTestParticleCase( );

    caseNumber = 2;
    database::TestParticleCasePointer testParticleCase2 = getTestParticleCase( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleCase < *testParticleCase2 );
}

//! Test comparison of TestParticleCase pointers using overloaded > operator.
BOOST_AUTO_TEST_CASE( testTestParticleCaseGreaterThanComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleCase pointers.
    database::TestParticleCasePointer testParticleCase = getTestParticleCase( );

    caseNumber = 2;
    database::TestParticleCasePointer testParticleCase2 = getTestParticleCase( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleCase2 > *testParticleCase );
}

//! Test comparison of TestParticleCase pointers using overloaded <= operator.
BOOST_AUTO_TEST_CASE( testTestParticleCaseLessThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleCase pointers.
    database::TestParticleCasePointer testParticleCase = getTestParticleCase( );

    caseNumber = 2;
    database::TestParticleCasePointer testParticleCase2 = getTestParticleCase( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleCase <= *testParticleCase2 );
}

//! Test comparison of TestParticleCase pointers using overloaded >= operator.
BOOST_AUTO_TEST_CASE( testTestParticleCaseGreaterThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleCase pointers.
    database::TestParticleCasePointer testParticleCase = getTestParticleCase( );

    caseNumber = 2;
    database::TestParticleCasePointer testParticleCase2 = getTestParticleCase( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleCase2 >= *testParticleCase );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stochastic_migration    
