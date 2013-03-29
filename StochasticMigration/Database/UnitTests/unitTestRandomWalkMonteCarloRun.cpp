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

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include <Assist/Basics/operatorOverloadFunctions.h>

#include "StochasticMigration/Database/randomWalkMonteCarloRun.h" 

namespace stochastic_migration
{
namespace unit_tests
{

//! Test fixture used to test the RandomWalkMonteCarloRun struct.
struct RandomWalkMonteCarloRunFixture
{
public:

    //! Constructor initializing valid parameters.
    RandomWalkMonteCarloRunFixture( )
        : monteCarloRun( 1 ),
          perturberPopulation( 10 ),
          massDistributionType( "EQUAL" ),
          massDistributionParameter1( 0.1 ),
          massDistributionParameter2( 0.0 ),
          observationPeriod( 4.0 * 365.25 * 86400.0 ),
          epochWindowSize( 30.0 * 86400.0 ),
          numberOfEpochWindows( 4 )
    { }

    //! Declare parameters of random walk Monte Carlo run.
    //! Monte Carlo run.
    int monteCarloRun;

    //! Perturber population.
    int perturberPopulation;

    //! Mass distribution type used for random walk simulation (can be equal, uniform, linear, 
    //! or power-law).
    std::string massDistributionType;

    //! Mass distribution parameter 1.
    double massDistributionParameter1;

    //! Mass distribution parameter 2.
    double massDistributionParameter2;

    //! Observation period [s].
    double observationPeriod;

    //! Epoch window size used to average observation data [s].
    double epochWindowSize;

    //! Number of epoch windows within prescribed observation period.
    int numberOfEpochWindows;

    //! Get random walk Monte Carlo run created from specified parameters.
    database::RandomWalkMonteCarloRunPointer getRandomWalkMonteCarloRun( )
    {
        return boost::make_shared< database::RandomWalkMonteCarloRun >( 
            monteCarloRun, perturberPopulation, massDistributionType, massDistributionParameter1,
            massDistributionParameter2, observationPeriod, epochWindowSize, numberOfEpochWindows );
    }

protected:

private:
};        

BOOST_FIXTURE_TEST_SUITE( test_random_walk_monte_carlo_run, RandomWalkMonteCarloRunFixture )    

//! Test correct construction of random walk Monte Carlo run.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunStructContruction )
{
    // Create random walk Monte Carlo run.
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarlorun 
        = getRandomWalkMonteCarloRun( );

    // Check that the random walk Monte Carlo run created contains all the data as required.
    BOOST_CHECK_EQUAL( randomWalkMonteCarlorun->monteCarloRun, monteCarloRun );
    BOOST_CHECK_EQUAL( randomWalkMonteCarlorun->perturberPopulation, perturberPopulation );
    BOOST_CHECK_EQUAL( randomWalkMonteCarlorun->massDistributionType, massDistributionType );
    BOOST_CHECK_EQUAL( randomWalkMonteCarlorun->massDistributionParameter1,
                       massDistributionParameter1 );
    BOOST_CHECK_EQUAL( randomWalkMonteCarlorun->massDistributionParameter2,
                       massDistributionParameter2 );
    BOOST_CHECK_EQUAL( randomWalkMonteCarlorun->observationPeriod, observationPeriod );
    BOOST_CHECK_EQUAL( randomWalkMonteCarlorun->epochWindowSize, epochWindowSize );
    BOOST_CHECK_EQUAL( randomWalkMonteCarlorun->numberOfEpochWindows, numberOfEpochWindows );
}

//! Test initialization of random walk Monte Carlo run with non-positive Monte Carlo run.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunNonPositiveMonteCarloRunError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set Monte Carlo run to invalid (non-positive) number.
    monteCarloRun = -1;

    // Try to create random walk Monte Carlo run.
    try 
    { 
        database::RandomWalkMonteCarloRunPointer randomWalkMonteCarlorun 
            = getRandomWalkMonteCarloRun( );
    }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk Monte Carlo run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk Monte Carlo run with non-positive perturber population.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunNonPositivePerturberPopulationError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pertuber population to invalid (non-positive) number.
    perturberPopulation = -1;

    // Try to create random walk Monte Carlo run.
    try 
    { 
        database::RandomWalkMonteCarloRunPointer randomWalkMonteCarlorun 
            = getRandomWalkMonteCarloRun( );
    }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk Monte Carlo run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk Monte Carlo run with empty mass-distribution type.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunEmptyMassDistributionTypeError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set mass-distribution type to empty.
    massDistributionType = "";

    // Try to create random walk Monte Carlo run.
    try 
    { 
        database::RandomWalkMonteCarloRunPointer randomWalkMonteCarlorun 
            = getRandomWalkMonteCarloRun( );
    }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk Monte Carlo run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk Monte Carlo run with non-positive observation period.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunNonPositiveObservationPeriodError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set observation period to invalid (non-positive) number.
    perturberPopulation = -1.234e5;

    // Try to create random walk Monte Carlo run.
    try 
    { 
        database::RandomWalkMonteCarloRunPointer randomWalkMonteCarlorun 
            = getRandomWalkMonteCarloRun( );
    }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk Monte Carlo run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk Monte Carlo run with non-positive epoch window size.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunNonPositiveEpochWindowSizeError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set epoch window size to invalid (non-positive) number.
    epochWindowSize = -1.23e4;

    // Try to create random walk Monte Carlo run.
    try 
    { 
        database::RandomWalkMonteCarloRunPointer randomWalkMonteCarlorun 
            = getRandomWalkMonteCarloRun( );
    }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk Monte Carlo run failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk Monte Carlo run with non-positive number of epoch windows.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunNonPositiveNumberOfEpochWindowsError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set number of epoch windows to invalid (non-positive) number.
    numberOfEpochWindows = -1;

    // Try to create random walk Monte Carlo run.
    try 
    { 
        database::RandomWalkMonteCarloRunPointer randomWalkMonteCarlorun 
            = getRandomWalkMonteCarloRun( );
    }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk Monte Carlo run failed.
    BOOST_CHECK( isError );
}



//! Test comparison of RandomWalkMonteCarloRun pointers using overloaded == operator.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunEqualComparison )
{
    // Create RandomWalkMonteCarloRun pointers.
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarloRun 
        = getRandomWalkMonteCarloRun( );
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarloRun2 
        = getRandomWalkMonteCarloRun( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkMonteCarloRun == *randomWalkMonteCarloRun2 );
}

//! Test comparison of RandomWalkMonteCarloRun pointers using overloaded != operator.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunNonEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkMonteCarloRun pointers.
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarloRun 
        = getRandomWalkMonteCarloRun( );

    monteCarloRun = 2;
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarloRun2 
        = getRandomWalkMonteCarloRun( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkMonteCarloRun != *randomWalkMonteCarloRun2 );
}

//! Test comparison of RandomWalkMonteCarloRun pointers using overloaded < operator.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunLessThanComparison )
{
    // Create RandomWalkMonteCarloRun pointers.
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarloRun 
        = getRandomWalkMonteCarloRun( );

    monteCarloRun = 2;
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarloRun2 
        = getRandomWalkMonteCarloRun( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkMonteCarloRun < *randomWalkMonteCarloRun2 );
}

//! Test comparison of RandomWalkMonteCarloRun pointers using overloaded > operator.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunGreaterThanComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkMonteCarloRun pointers.
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarloRun 
        = getRandomWalkMonteCarloRun( );

    monteCarloRun = 2;
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarloRun2 
        = getRandomWalkMonteCarloRun( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkMonteCarloRun2 > *randomWalkMonteCarloRun );
}

//! Test comparison of RandomWalkMonteCarloRun pointers using overloaded <= operator.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunLessThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkMonteCarloRun pointers.
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarloRun 
        = getRandomWalkMonteCarloRun( );

    monteCarloRun = 2;
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarloRun2 
        = getRandomWalkMonteCarloRun( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkMonteCarloRun <= *randomWalkMonteCarloRun2 );
}

//! Test comparison of RandomWalkMonteCarloRun pointers using overloaded >= operator.
BOOST_AUTO_TEST_CASE( testRandomWalkMonteCarloRunGreaterThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkMonteCarloRun pointers.
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarloRun 
        = getRandomWalkMonteCarloRun( );

    monteCarloRun = 2;
    database::RandomWalkMonteCarloRunPointer randomWalkMonteCarloRun2 
        = getRandomWalkMonteCarloRun( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkMonteCarloRun2 >= *randomWalkMonteCarloRun );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stochastic_migration      