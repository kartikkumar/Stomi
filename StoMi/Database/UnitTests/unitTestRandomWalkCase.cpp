/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#define BOOST_TEST_MAIN

#include <string>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include <Assist/Basics/operatorOverloadFunctions.h>

#include "StoMi/Database/randomWalkCase.h"
 
namespace stochastic_migration
{
namespace unit_tests
{

//! Test fixture used to test the RandomWalkCase struct.
struct RandomWalkCaseFixture
{
public:

    //!ructor initializing valid parameters.
    RandomWalkCaseFixture( )
        : caseId( 1 ),
          caseName( "test_case" ),
          testParticleCaseId( 1 ),
          perturberDensity( 1.23 ),
          perturberRingMass( 4.56 ),
          observationPeriod( 1234567.89 ),
          numberOfEpochWindows( 4 ),
          epochWindowSize( 987.65 )
    { }

    //! Declare parameters of test particle case.

    // Required parameters.

    //! Case ID.
    int caseId;

    //! Case name.
    std::string caseName;

    //! Test particle case ID.
    int testParticleCaseId;

    //! Perturber density [Number of perturber per Hill radius of perturbed body].
    double perturberDensity;

    //! Perturber ring mass [M_perturbedBody].
    double perturberRingMass;

    //! Observation period [s].
    double observationPeriod;

    //! Number of epoch windows within observation period.
    int numberOfEpochWindows;

    //! Epoch window size [s].
    double epochWindowSize;

    //! Get random walk case created from specified parameters.
    database::RandomWalkCasePointer getRandomWalkCase( )
    {
        return boost::make_shared< database::RandomWalkCase >(
            caseId, caseName, testParticleCaseId, perturberDensity, perturberRingMass,
            observationPeriod, numberOfEpochWindows, epochWindowSize );
    }

protected:

private:
};    

BOOST_FIXTURE_TEST_SUITE( test_random_walk_case, RandomWalkCaseFixture )

//! Test correctruction of random walk case.
BOOST_AUTO_TEST_CASE( testRandomWalkCaseStructContruction )
{
    using namespace database;

    // Create test particle case.
    database::RandomWalkCasePointer RandomWalkCase = getRandomWalkCase( );

    // Check that the test particle case created contains all the data as required.
    BOOST_CHECK_EQUAL( RandomWalkCase->caseId, caseId );    
    BOOST_CHECK_EQUAL( RandomWalkCase->caseName, caseName );      
    BOOST_CHECK_EQUAL( RandomWalkCase->testParticleCaseId, testParticleCaseId );      
    BOOST_CHECK_EQUAL( RandomWalkCase->perturberDensity, perturberDensity );      
    BOOST_CHECK_EQUAL( RandomWalkCase->perturberRingMass, perturberRingMass );      
    BOOST_CHECK_EQUAL( RandomWalkCase->observationPeriod, observationPeriod );      
    BOOST_CHECK_EQUAL( RandomWalkCase->numberOfEpochWindows, numberOfEpochWindows );      
    BOOST_CHECK_EQUAL( RandomWalkCase->epochWindowSize, epochWindowSize );      
}

//! Test initialization of random walk case with empty name.
BOOST_AUTO_TEST_CASE( testRandomWalkCaseNameEmpty )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set case name to empty string.
    caseName = "";

    // Try to create test particle case.
    try { database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk case ID with non-positive number.
BOOST_AUTO_TEST_CASE( testRandomWalkCaseIdNonPositiveNumberError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set case ID to invalid (non-positive) number.
    caseId = -1;

    // Try to create random walk case.
    try { database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle case ID with non-positive number.
BOOST_AUTO_TEST_CASE( testrandomWalkCaseIdNonPositiveNumberError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set test particle case ID to invalid (non-positive) number.
    testParticleCaseId = -1;

    // Try to create random walk case.
    try { database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of perturber density with non-positive number.
BOOST_AUTO_TEST_CASE( testPerturberDensityNonPositiveError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set perturber density to invalid (non-positive) number.
    perturberDensity = -1.0;

    // Try to create random walk case.
    try { database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of perturber ring mass with non-positive number.
BOOST_AUTO_TEST_CASE( testPerturberRingMassNonPositiveError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set perturber ring mass to invalid (non-positive) number.
    perturberRingMass = -1.0;

    // Try to create random walk case.
    try { database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of observation period with non-positive number.
BOOST_AUTO_TEST_CASE( testObservationPeriodNonPositiveError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set observation period to invalid (non-positive) number.
    observationPeriod = -1.2345;

    // Try to create random walk case.
    try { database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of number of epoch windows with non-positive number.
BOOST_AUTO_TEST_CASE( testNumberOfEpochWindowsNonPositiveError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set number of epoch windows to invalid (non-positive) number.
    numberOfEpochWindows = -2;

    // Try to create random walk case.
    try { database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk case failed.
    BOOST_CHECK( isError );
}

//! Test initialization of epoch window size with non-positive number.
BOOST_AUTO_TEST_CASE( testEpochWindowSizeNonPositiveError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set epoch window size to invalid (non-positive) number.
    epochWindowSize = -2.345;

    // Try to create random walk case.
    try { database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk case failed.
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

    // Try to create random walk case.
    try { database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk case failed.
    BOOST_CHECK( isError );
}

//! Test comparison of RandomWalkCase pointers using overloaded == operator.
BOOST_AUTO_TEST_CASE( testRandomWalkCaseEqualComparison )
{
    // Create RandomWalkCase pointers.
    database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( );
    database::RandomWalkCasePointer randomWalkCase2 = getRandomWalkCase( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkCase == *randomWalkCase2 );
}

//! Test comparison of RandomWalkCase pointers using overloaded != operator.
BOOST_AUTO_TEST_CASE( testRandomWalkCaseNonEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkCase pointers.
    database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( );

    caseId = 2;
    database::RandomWalkCasePointer randomWalkCase2 = getRandomWalkCase( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkCase != *randomWalkCase2 );
}

//! Test comparison of RandomWalkCase pointers using overloaded < operator.
BOOST_AUTO_TEST_CASE( testRandomWalkCaseLessThanComparison )
{
    // Create RandomWalkCase pointers.
    database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( );

    caseId = 2;
    database::RandomWalkCasePointer randomWalkCase2 = getRandomWalkCase( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkCase < *randomWalkCase2 );
}

//! Test comparison of RandomWalkCase pointers using overloaded > operator.
BOOST_AUTO_TEST_CASE( testRandomWalkCaseGreaterThanComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkCase pointers.
    database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( );

    caseId = 2;
    database::RandomWalkCasePointer randomWalkCase2 = getRandomWalkCase( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkCase2 > *randomWalkCase );
}

//! Test comparison of RandomWalkCase pointers using overloaded <= operator.
BOOST_AUTO_TEST_CASE( testRandomWalkCaseLessThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkCase pointers.
    database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( );

    caseId = 2;
    database::RandomWalkCasePointer randomWalkCase2 = getRandomWalkCase( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkCase <= *randomWalkCase2 );
}

//! Test comparison of RandomWalkCase pointers using overloaded >= operator.
BOOST_AUTO_TEST_CASE( testRandomWalkCaseGreaterThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkCase pointers.
    database::RandomWalkCasePointer randomWalkCase = getRandomWalkCase( );

    caseId = 2;
    database::RandomWalkCasePointer randomWalkCase2 = getRandomWalkCase( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkCase2 >= *randomWalkCase );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stochastic_migration  