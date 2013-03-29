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
 *      130212    K. Kumar          File created.
 *      130217    K. Kumar          Updated "mab simulations" references to "stochastic migration";
 *                                  renamed file.
 *      130309    K. Kumar          Expanded unit tests to test all functionality of
 *                                  TestParticleKick objects.
 *      130328    K. Kumar          Updated unit tests to use boiler plate operator overload 
 *                                  functions in Assist.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <set>
#include <typeinfo>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include <Assist/Basics/operatorOverloadFunctions.h>

#include "StochasticMigration/Database/testParticleKick.h"

namespace stochastic_migration
{
namespace unit_tests
{

//! Test fixture used to test the TestParticleKick struct.
struct TestParticleKickFixture
{
public:

    //! Constructor initializing valid parameters.
    TestParticleKickFixture( )
        : simulationNumber( 1 ),
          conjunctionEpoch( 1.23e4 ),
          conjunctionDistance( 2.0e7 ),
          preConjunctionEpoch( 5.67e2 ),
          preConjunctionDistance( 789.54e8 ),
          preConjunctionSemiMajorAxis( 112233.44 ),
          preConjunctionEccentricity( 0.456 ),
          preConjunctionInclination( 1.234 ),
          postConjunctionEpoch( 9.87e5 ),
          postConjunctionDistance( 987.65e8 ),
          postConjunctionSemiMajorAxis( 112244.55 ),
          postConjunctionEccentricity( 0.567 ),
          postConjunctionInclination( 1.345 ),
          massRatio( 0.222 )
    { }

    //! Declare parameters of test particle kick.
    int simulationNumber;
    double conjunctionEpoch;
    double conjunctionDistance;
    double conjunctionDuration;
    double preConjunctionEpoch;
    double preConjunctionDistance;
    double preConjunctionSemiMajorAxis;
    double preConjunctionEccentricity;
    double preConjunctionInclination;
    double postConjunctionEpoch;
    double postConjunctionDistance;
    double postConjunctionSemiMajorAxis;
    double postConjunctionEccentricity;
    double postConjunctionInclination;
    double massRatio;

    //! Get test particle kick created from specified parameters.
    database::TestParticleKickPointer getTestParticleKick( )
    {
        return boost::make_shared< database::TestParticleKick >(
                    database::TestParticleKick(
                        simulationNumber, conjunctionEpoch, conjunctionDistance,
                        conjunctionDuration, preConjunctionEpoch, preConjunctionDistance,
                        preConjunctionSemiMajorAxis, preConjunctionEccentricity,
                        preConjunctionInclination, postConjunctionEpoch, postConjunctionDistance,
                        postConjunctionSemiMajorAxis, postConjunctionEccentricity,
                        postConjunctionInclination, massRatio ) );
    }

protected:

private:
};

BOOST_FIXTURE_TEST_SUITE( test_test_particle_kick, TestParticleKickFixture )

//! Test definition of typedef for test particle kick database struct.
BOOST_AUTO_TEST_CASE( testTestParticleKickTypedef )
{
    BOOST_CHECK( typeid( database::TestParticleKickPointerTable )
                 == typeid( std::set< database::TestParticleKickPointer > ) );
}

//! Test correct construction of test particle kick.
BOOST_AUTO_TEST_CASE( testTestParticleKickStructContruction )
{
    // Create test particle kick.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );

    // Check that the test particle kick created contains all the data as required.
    BOOST_CHECK_EQUAL( testParticleKick->simulationNumber, simulationNumber );
    BOOST_CHECK_EQUAL( testParticleKick->conjunctionEpoch, conjunctionEpoch );
    BOOST_CHECK_EQUAL( testParticleKick->conjunctionDistance, conjunctionDistance );
    BOOST_CHECK_EQUAL( testParticleKick->conjunctionDuration, conjunctionDuration );
    BOOST_CHECK_EQUAL( testParticleKick->preConjunctionEpoch, preConjunctionEpoch );
    BOOST_CHECK_EQUAL( testParticleKick->preConjunctionDistance, preConjunctionDistance );
    BOOST_CHECK_EQUAL( testParticleKick->preConjunctionSemiMajorAxis,
                       preConjunctionSemiMajorAxis );
    BOOST_CHECK_EQUAL( testParticleKick->preConjunctionEccentricity, preConjunctionEccentricity );
    BOOST_CHECK_EQUAL( testParticleKick->preConjunctionInclination, preConjunctionInclination );
    BOOST_CHECK_EQUAL( testParticleKick->postConjunctionEpoch, postConjunctionEpoch );
    BOOST_CHECK_EQUAL( testParticleKick->postConjunctionDistance, postConjunctionDistance );
    BOOST_CHECK_EQUAL( testParticleKick->postConjunctionSemiMajorAxis,
                       postConjunctionSemiMajorAxis );
    BOOST_CHECK_EQUAL( testParticleKick->postConjunctionEccentricity,
                       postConjunctionEccentricity );
    BOOST_CHECK_EQUAL( testParticleKick->postConjunctionInclination, postConjunctionInclination );
    BOOST_CHECK_EQUAL( testParticleKick->massRatio, massRatio );
}

//! Test initialization of test particle kick with non-positive simulation number.
BOOST_AUTO_TEST_CASE( testTestParticleKickNonPositiveSimulationNumberError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set simulation number to invalid (non-positive) number.
    simulationNumber = -1;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative conjunction distance.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativeConjunctionDistanceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set conjunction distance to invalid (negative) number.
    conjunctionDistance = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative conjunction duration.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativeConjunctionDurationError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set conjunction duration to invalid (negative) number.
    conjunctionDistance = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative pre-conjunction distance.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePreConjunctionDistanceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction distance to invalid (negative) number.
    preConjunctionDistance = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative pre-conjunction semi-major axis.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePreConjunctionSemiMajorAxisError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction semi-major axis to invalid (negative) number.
    preConjunctionSemiMajorAxis = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative pre-conjunction eccentricity.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePreConjunctionEccentricityError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction eccentricity to invalid (negative) number.
    preConjunctionEccentricity = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with >1.0 pre-conjunction eccentricity.
BOOST_AUTO_TEST_CASE( testTestParticleKickGreaterThanOnePreConjunctionEccentricityError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction eccentricity to invalid (greater than 1.0) number.
    preConjunctionEccentricity = 1.1;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative pre-conjunction inclination.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePreConjunctionInclinationError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction eccentricity to invalid (negative) number.
    preConjunctionInclination = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative post-conjunction distance.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePostConjunctionDistanceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction distance to invalid (negative) number.
    postConjunctionDistance = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative post-conjunction semi-major axis.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePostConjunctionSemiMajorAxisError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction semi-major axis to invalid (negative) number.
    postConjunctionSemiMajorAxis = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative post-conjunction eccentricity.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePostConjunctionEccentricityError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction eccentricity to invalid (negative) number.
    postConjunctionEccentricity = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with >1.0 post-conjunction eccentricity.
BOOST_AUTO_TEST_CASE( testTestParticleKickGreaterThanOnePostConjunctionEccentricityError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction eccentricity to invalid (greater than 1.0) number.
    postConjunctionEccentricity = 1.1;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative post-conjunction inclination.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePostConjunctionInclinationError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction eccentricity to invalid (negative) number.
    postConjunctionInclination = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative mass factor.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativemassRatioError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set mass factor to invalid (negative) number.
    massRatio = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with >1.0 mass factor.
BOOST_AUTO_TEST_CASE( testTestParticleKickGreaterThanOnemassRatioError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set mass factor to invalid (greater than 1.0) number.
    massRatio = 1.1;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with pre-conjunction epoch after conjunction epoch.
BOOST_AUTO_TEST_CASE( testTestParticleKickPreConjunctionAfterConjunctionEpochError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction epoch to value after conjunction epoch.
    preConjunctionEpoch = 2.34e5;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick woth conjunction epoch after post-conjunction epoch.
BOOST_AUTO_TEST_CASE( testTestParticleKickConjunctionAfterPostConjunctionEpochError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction epoch to value before conjunction epoch.
    postConjunctionEpoch = 0.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with pre-conjunction distance less than conjunction
//! distance.
BOOST_AUTO_TEST_CASE( testTestParticleKickPreConjunctionLessThanConjunctionDistanceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction distance to value less than conjunction distance.
    preConjunctionDistance = 1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with post-conjunction distance less than conjunction
//! distance.
BOOST_AUTO_TEST_CASE( testTestParticleKickPostConjunctionLessThanConjunctionDistanceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction distance to value less than conjunction distance.
    postConjunctionDistance = 1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with conjunction duration greater than
//! post-conjunction - pre-conjunction epochs.
BOOST_AUTO_TEST_CASE( testTestParticleKickConjunctionDurationGreaterThanPostMinusPreError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set conjunction duration to value less than post-conjunction duration - pre-conjunction
    // duration.
    conjunctionDistance = 999.e9;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test comparison of TestParticleKick pointers using overloaded == operator.
BOOST_AUTO_TEST_CASE( testTestParticleKickEqualComparison )
{
    // Create TestParticleKick pointers.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );
    database::TestParticleKickPointer testParticleKick2 = getTestParticleKick( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleKick == *testParticleKick2 );
}

//! Test comparison of TestParticleKick pointers using overloaded != operator.
BOOST_AUTO_TEST_CASE( testTestParticleKickNonEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleKick pointers.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );

    simulationNumber = 2;
    database::TestParticleKickPointer testParticleKick2 = getTestParticleKick( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleKick != *testParticleKick2 );
}

//! Test comparison of TestParticleKick pointers using overloaded < operator.
BOOST_AUTO_TEST_CASE( testTestParticleKickLessThanComparison )
{
    // Create TestParticleKick pointers.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );

    conjunctionEpoch = 1.24e4;
    database::TestParticleKickPointer testParticleKick2 = getTestParticleKick( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleKick < *testParticleKick2 );
}

//! Test comparison of TestParticleKick pointers using overloaded > operator.
BOOST_AUTO_TEST_CASE( testTestParticleKickGreaterThanComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleKick pointers.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );

    conjunctionEpoch = 1.24e4;
    database::TestParticleKickPointer testParticleKick2 = getTestParticleKick( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleKick2 > *testParticleKick );
}

//! Test comparison of TestParticleKick pointers using overloaded <= operator.
BOOST_AUTO_TEST_CASE( testTestParticleKickLessThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleKick pointers.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );

    conjunctionEpoch = 1.24e4;
    database::TestParticleKickPointer testParticleKick2 = getTestParticleKick( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleKick <= *testParticleKick2 );
}

//! Test comparison of TestParticleKick pointers using overloaded >= operator.
BOOST_AUTO_TEST_CASE( testTestParticleKickGreaterThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleKick pointers.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );

    conjunctionEpoch = 1.24e4;
    database::TestParticleKickPointer testParticleKick2 = getTestParticleKick( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleKick2 >= *testParticleKick );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stochastic_migration
