/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <vector>

#include <boost/assign/list_of.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include <Assist/Basics/operatorOverloadFunctions.h>

#include "Stomi/Database/randomWalkInput.h"

namespace stomi
{
namespace unit_tests
{

//! Test fixture used to test the RandomWalkInput struct.
struct RandomWalkInputFixture
{
public:

    //! Constructor initializing valid parameters.
    RandomWalkInputFixture( )
        : randomWalkSimulationId( 1 ),
          randomWalkRunId( 1 ),
          isCompleted( false ),
          observationPeriodStartEpoch( 1.23 ),
          testParticleSimulationIds( boost::assign::list_of( 1 )( 2 )
                                        .convert_to_container<std::vector< int > >( ) )
    { }

    //! Declare parameters of random walk input.

    // Required parameters.

    //! Random walk simulation ID.
    int randomWalkSimulationId;

    //! Random walk run ID.
    int randomWalkRunId;

    //! Flag indicating if simulation has been completed/executed and stored in database.
    bool isCompleted;    

    //! Epoch at start of observation period [s].
    double observationPeriodStartEpoch;

    //! List of test particle simulation IDs, used to generate perturbers.
    std::vector< int > testParticleSimulationIds;

    //! Get random walk input created from specified parameters.
    database::RandomWalkInputPointer getRandomWalkInput( )
    {
        return boost::make_shared< database::RandomWalkInput >(
            randomWalkSimulationId, randomWalkRunId, isCompleted, observationPeriodStartEpoch, 
            testParticleSimulationIds );
    }

protected:

private:
};       

BOOST_FIXTURE_TEST_SUITE( test_random_walk_input, RandomWalkInputFixture )

//! Test correct construction of random walk input.
BOOST_AUTO_TEST_CASE( testRandomWalkInputStructContruction )
{
    using namespace database;

    // Create random walk input.
    RandomWalkInputPointer randomWalkInput = getRandomWalkInput( );

    // Check that the random walk input created contains all the data as required.
    BOOST_CHECK_EQUAL( randomWalkInput->randomWalkSimulationId, randomWalkSimulationId );    
    BOOST_CHECK_EQUAL( randomWalkInput->randomWalkRunId, randomWalkRunId );      
    BOOST_CHECK_EQUAL( randomWalkInput->isCompleted, isCompleted );      
    BOOST_CHECK_EQUAL( randomWalkInput->observationPeriodStartEpoch, 
                       observationPeriodStartEpoch ); 

    for ( unsigned int i = 0; i < testParticleSimulationIds.size( ); i++ )
    {
        BOOST_CHECK_EQUAL( randomWalkInput->testParticleSimulationIds.at( i ), 
                           testParticleSimulationIds.at( i ) );         
    }     
}

//! Test initialization of random walk input with non-positive number random walk simulation ID.
BOOST_AUTO_TEST_CASE( testRandomWalkInputRandomWalkRunIdNonPositiveNumberError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set random walk simulation ID to invalid (non-positive) number.
    randomWalkSimulationId = -1;

    // Try to create random walk input.
    try { database::RandomWalkInputPointer randomWalkInput = getRandomWalkInput( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk input failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk input with non-positive number random walk run ID.
BOOST_AUTO_TEST_CASE( testRandomWalkInputCaseIdNonPositiveNumberError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set random walk run ID to invalid (non-positive) number.
    randomWalkRunId = -1;

    // Try to create random walk input.
    try { database::RandomWalkInputPointer randomWalkInput = getRandomWalkInput( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk input failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk input observation period epoch with non-positive number.
BOOST_AUTO_TEST_CASE( testRandomWalkInputObservationPeriodEpochNonPositiveNumberError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set observation period epoch to invalid (non-positive) number.
    observationPeriodStartEpoch = -1;

    // Try to create random walk input.
    try { database::RandomWalkInputPointer randomWalkInput = getRandomWalkInput( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk input failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk input empty test particle simulation IDs list.
BOOST_AUTO_TEST_CASE( testRandomWalkInputTestParticleSimulationIdsEmptyError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set test particle simulations ID list to empty.
    testParticleSimulationIds = std::vector< int >( );

    // Try to create random walk input.
    try { database::RandomWalkInputPointer randomWalkInput = getRandomWalkInput( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk input failed.
    BOOST_CHECK( isError );
}

//! Test comparison of RandomWalkInput pointers using overloaded == operator.
BOOST_AUTO_TEST_CASE( testRandomWalkInputEqualComparison )
{
    // Create RandomWalkInput pointers.
    database::RandomWalkInputPointer randomWalkInput = getRandomWalkInput( );
    database::RandomWalkInputPointer randomWalkInput2 = getRandomWalkInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkInput == *randomWalkInput2 );
}

//! Test comparison of RandomWalkInput pointers using overloaded != operator.
BOOST_AUTO_TEST_CASE( testRandomWalkInputNonEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkInput pointers.
    database::RandomWalkInputPointer randomWalkInput = getRandomWalkInput( );

    randomWalkSimulationId = 2;
    database::RandomWalkInputPointer randomWalkInput2 = getRandomWalkInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkInput != *randomWalkInput2 );
}

//! Test comparison of RandomWalkInput pointers using overloaded < operator.
BOOST_AUTO_TEST_CASE( testRandomWalkInputLessThanComparison )
{
    // Create RandomWalkInput pointers.
    database::RandomWalkInputPointer randomWalkInput = getRandomWalkInput( );

    randomWalkSimulationId = 2;
    database::RandomWalkInputPointer randomWalkInput2 = getRandomWalkInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkInput < *randomWalkInput2 );
}

//! Test comparison of RandomWalkInput pointers using overloaded > operator.
BOOST_AUTO_TEST_CASE( testRandomWalkInputGreaterThanComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkInput pointers.
    database::RandomWalkInputPointer randomWalkInput = getRandomWalkInput( );

    randomWalkSimulationId = 2;
    database::RandomWalkInputPointer randomWalkInput2 = getRandomWalkInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkInput2 > *randomWalkInput );
}

//! Test comparison of RandomWalkInput pointers using overloaded <= operator.
BOOST_AUTO_TEST_CASE( testRandomWalkInputLessThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkInput pointers.
    database::RandomWalkInputPointer randomWalkInput = getRandomWalkInput( );

    randomWalkSimulationId = 2;
    database::RandomWalkInputPointer randomWalkInput2 = getRandomWalkInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkInput <= *randomWalkInput2 );
}

//! Test comparison of RandomWalkInput pointers using overloaded >= operator.
BOOST_AUTO_TEST_CASE( testRandomWalkInputGreaterThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkInput pointers.
    database::RandomWalkInputPointer randomWalkInput = getRandomWalkInput( );

    randomWalkSimulationId = 2;
    database::RandomWalkInputPointer randomWalkInput2 = getRandomWalkInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkInput2 >= *randomWalkInput );
}

BOOST_AUTO_TEST_SUITE_END( )
 
} // namespace unit_tests
} // namespace stomi   