/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <string>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include <Assist/Basics/operatorOverloadFunctions.h>

#include "Stomi/Database/randomWalkRun.h"
 
namespace stomi
{
namespace unit_tests
{

//! Test fixture used to test the RandomWalkRun struct.
struct RandomWalkRunFixture
{
public:

    //! Constructor initializing valid parameters.
    RandomWalkRunFixture( )
        : randomWalkRunId( 1 ),
          randomWalkRunName( "test_run" ),
          testParticleCaseId( 1 ),
          perturberRingNumberDensity( 1.23 ),
          perturberRingMass( 4.56 ),
          observationPeriod( 1234567.89 ),
          numberOfEpochWindows( 4 ),
          epochWindowSize( 987.65 )
    { }

    //! Declare parameters of random walk run.

    // Required parameters.

    //! Random walk run ID.
    int randomWalkRunId;

    //! Random walk run name.
    std::string randomWalkRunName;

    //! Test particle run ID.
    int testParticleCaseId;

    //! Perturber ring number density [Number of perturber per Hill radius of perturbed body].
    double perturberRingNumberDensity;

    //! Perturber ring mass [M_perturbedBody].
    double perturberRingMass;

    //! Observation period [s].
    double observationPeriod;

    //! Number of epoch windows within observation period.
    int numberOfEpochWindows;

    //! Epoch window size [s].
    double epochWindowSize;

    //! Get random walk run created from specified parameters.
    database::RandomWalkRunPointer getRandomWalkRun( )
    {
        return boost::make_shared< database::RandomWalkRun >(
            randomWalkRunId, randomWalkRunName, testParticleCaseId, perturberRingNumberDensity, 
            perturberRingMass, observationPeriod, numberOfEpochWindows, epochWindowSize );
    }

protected:

private:
};    

BOOST_FIXTURE_TEST_SUITE( test_random_walk_run, RandomWalkRunFixture )

//! Test construction of random walk run.
BOOST_AUTO_TEST_CASE( testRandomWalkRunStructContruction )
{
    using namespace database;

    // Create RandomWalkRun object.
    database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( );

    // Check that the RandomWalkRun object created contains all the data as required.
    BOOST_CHECK_EQUAL( RandomWalkRun->randomWalkRunId, randomWalkRunId );    
    BOOST_CHECK_EQUAL( RandomWalkRun->randomWalkRunName, randomWalkRunName );      
    BOOST_CHECK_EQUAL( RandomWalkRun->testParticleCaseId, testParticleCaseId );      
    BOOST_CHECK_EQUAL( RandomWalkRun->perturberRingNumberDensity, perturberRingNumberDensity );      
    BOOST_CHECK_EQUAL( RandomWalkRun->perturberRingMass, perturberRingMass );      
    BOOST_CHECK_EQUAL( RandomWalkRun->observationPeriod, observationPeriod );      
    BOOST_CHECK_EQUAL( RandomWalkRun->numberOfEpochWindows, numberOfEpochWindows );      
    BOOST_CHECK_EQUAL( RandomWalkRun->epochWindowSize, epochWindowSize );      
}

//! Test initialization of random walk run with empty name.
BOOST_AUTO_TEST_CASE( testRandomWalkRunNameEmpty )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set run name to empty string.
    randomWalkRunName = "";

    // Try to create random walk run
    try { database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk run ID with non-positive number.
BOOST_AUTO_TEST_CASE( testRandomWalkRunIdNonPositiveNumberError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set random walk run ID to invalid (non-positive) number.
    randomWalkRunId = -1;

    // Try to create random walk run.
    try { database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle run ID with non-positive number.
BOOST_AUTO_TEST_CASE( testTestParticleCaseIdNonPositiveNumberError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set test particle case ID to invalid (non-positive) number.
    testParticleCaseId = -1;

    // Try to create random walk run.
    try { database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of perturber density with non-positive number.
BOOST_AUTO_TEST_CASE( testPerturberDensityNonPositiveError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set perturber density to invalid (non-positive) number.
    perturberRingNumberDensity = -1.0;

    // Try to create random walk run.
    try { database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of perturber ring mass with non-positive number.
BOOST_AUTO_TEST_CASE( testPerturberRingMassNonPositiveError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set perturber ring mass to invalid (non-positive) number.
    perturberRingMass = -1.0;

    // Try to create random walk run.
    try { database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of observation period with non-positive number.
BOOST_AUTO_TEST_CASE( testObservationPeriodNonPositiveError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set observation period to invalid (non-positive) number.
    observationPeriod = -1.2345;

    // Try to create random walk run.
    try { database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of number of epoch windows with non-positive number.
BOOST_AUTO_TEST_CASE( testNumberOfEpochWindowsNonPositiveError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set number of epoch windows to invalid (non-positive) number.
    numberOfEpochWindows = -2;

    // Try to create random walk run.
    try { database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of epoch window size with non-positive number.
BOOST_AUTO_TEST_CASE( testEpochWindowSizeNonPositiveError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set epoch window size to invalid (non-positive) number.
    epochWindowSize = -2.345;

    // Try to create random walk run.
    try { database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of total epoch window size greater than observation period.
BOOST_AUTO_TEST_CASE( testTotalEpochWindowSizeError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set observation period.
    observationPeriod = 10.3456;

    // Set number of epoch windows.
    numberOfEpochWindows = 5;

    // Set epoch window size.
    epochWindowSize = 5.67;

    // Try to create random walk run.
    try { database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk run failed.
    BOOST_CHECK( isError );
}

//! Test comparison of RandomWalkRun pointers using overloaded == operator.
BOOST_AUTO_TEST_CASE( testRandomWalkRunEqualComparison )
{
    // Create RandomWalkRun pointers.
    database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( );
    database::RandomWalkRunPointer RandomWalkRun2 = getRandomWalkRun( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *RandomWalkRun == *RandomWalkRun2 );
}

//! Test comparison of RandomWalkRun pointers using overloaded != operator.
BOOST_AUTO_TEST_CASE( testRandomWalkRunNonEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkRun pointers.
    database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( );

    randomWalkRunId = 2;
    database::RandomWalkRunPointer RandomWalkRun2 = getRandomWalkRun( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *RandomWalkRun != *RandomWalkRun2 );
}

//! Test comparison of RandomWalkRun pointers using overloaded < operator.
BOOST_AUTO_TEST_CASE( testRandomWalkRunLessThanComparison )
{
    // Create RandomWalkRun pointers.
    database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( );

    randomWalkRunId = 2;
    database::RandomWalkRunPointer RandomWalkRun2 = getRandomWalkRun( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *RandomWalkRun < *RandomWalkRun2 );
}

//! Test comparison of RandomWalkRun pointers using overloaded > operator.
BOOST_AUTO_TEST_CASE( testRandomWalkRunGreaterThanComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkRun pointers.
    database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( );

    randomWalkRunId = 2;
    database::RandomWalkRunPointer RandomWalkRun2 = getRandomWalkRun( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *RandomWalkRun2 > *RandomWalkRun );
}

//! Test comparison of RandomWalkRun pointers using overloaded <= operator.
BOOST_AUTO_TEST_CASE( testRandomWalkRunLessThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkRun pointers.
    database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( );

    randomWalkRunId = 2;
    database::RandomWalkRunPointer RandomWalkRun2 = getRandomWalkRun( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *RandomWalkRun <= *RandomWalkRun2 );
}

//! Test comparison of RandomWalkRun pointers using overloaded >= operator.
BOOST_AUTO_TEST_CASE( testRandomWalkRunGreaterThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkRun pointers.
    database::RandomWalkRunPointer RandomWalkRun = getRandomWalkRun( );

    randomWalkRunId = 2;
    database::RandomWalkRunPointer RandomWalkRun2 = getRandomWalkRun( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *RandomWalkRun2 >= *RandomWalkRun );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stomi  
