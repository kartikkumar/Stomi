/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <vector>

#include <boost/assign/list_of.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "StoMi/Database/randomWalkInput.h"

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
        : monteCarloRunId( 1 ),
          randomWalkCaseId( 1 ),
          isCompleted( false ),
          observationPeriodStartEpoch( 1.23 ),
          testParticleSimulationIds( boost::assign::list_of( 1 )( 2 )
                                        .convert_to_container<std::vector< int > >( ) )
    { }

    //! Declare parameters of random walk input.

    // Required parameters.

    //! Monte Carlo run ID.
    int monteCarloRunId;

    //! Random walk case ID.
    int randomWalkCaseId;

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
            monteCarloRunId, randomWalkCaseId, isCompleted, observationPeriodStartEpoch, 
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
    BOOST_CHECK_EQUAL( randomWalkInput->monteCarloRunId, monteCarloRunId );    
    BOOST_CHECK_EQUAL( randomWalkInput->randomWalkCaseId, randomWalkCaseId );      
    BOOST_CHECK_EQUAL( randomWalkInput->isCompleted, isCompleted );      
    BOOST_CHECK_EQUAL( randomWalkInput->observationPeriodStartEpoch, 
                       observationPeriodStartEpoch ); 

    for ( unsigned int i = 0; i < testParticleSimulationIds.size( ); i++ )
    {
        BOOST_CHECK_EQUAL( randomWalkInput->testParticleSimulationIds.at( i ), 
                           testParticleSimulationIds.at( i ) );         
    }     
}

BOOST_AUTO_TEST_SUITE_END( )
 
} // namespace unit_tests
} // namespace stomi   