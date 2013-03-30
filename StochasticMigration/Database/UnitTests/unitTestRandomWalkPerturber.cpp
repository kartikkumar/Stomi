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
 *      130329    K. Kumar          Added test to check that test particle input table is sorted
 *                                  correctly; added typedef check.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/make_shared.hpp>
#include <boost/ptr_container/ptr_set.hpp>
#include <boost/test/unit_test.hpp>

#include <Assist/Basics/operatorOverloadFunctions.h>

#include "StochasticMigration/Database/randomWalkPerturber.h" 

namespace stochastic_migration
{
namespace unit_tests
{

//! Test fixture used to test the RandomWalkPerturber struct.
struct RandomWalkPerturberFixture
{
public:

    //! Constructor initializing valid parameters.
    RandomWalkPerturberFixture( )
        : monteCarloRun( 1 ),
          testParticleSimulationNumber( 1 ),
          massRatio( 0.2 )
    { }

    //! Declare parameters of random walk perturber.
    //! Monte Carlo run.
    int monteCarloRun;

    //! Test particle simulation number.
    int testParticleSimulationNumber;

    //! Mass ratio between body receiving kick and body causing kick [-].
    double massRatio;

    //! Get random walk perturber created from specified parameters.
    database::RandomWalkPerturberPointer getRandomWalkPerturber( )
    {
        return boost::make_shared< database::RandomWalkPerturber >( 
            monteCarloRun, testParticleSimulationNumber, massRatio );
    }

protected:

private:
};    

BOOST_FIXTURE_TEST_SUITE( test_random_walk_perturber, RandomWalkPerturberFixture )    

//! Test definition of typedef for random walk perturber database struct.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberTypedef )
{
    BOOST_CHECK( typeid( database::RandomWalkPerturberTable )
                 == typeid( boost::ptr_set< database::RandomWalkPerturber > ) );
}

//! Test correct construction of random walk perturber.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberStructContruction )
{
    // Create random walk perturber.
    database::RandomWalkPerturberPointer randomWalkPerturber = getRandomWalkPerturber( );

    // Check that the random walk perturber created contains all the data as required.
    BOOST_CHECK_EQUAL( randomWalkPerturber->monteCarloRun, monteCarloRun );
    BOOST_CHECK_EQUAL( randomWalkPerturber->testParticleSimulationNumber, 
                       testParticleSimulationNumber );
    BOOST_CHECK_EQUAL( randomWalkPerturber->massRatio, massRatio );
}

//! Test initialization of random walk perturber with non-positive Monte Carlo run.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberNonPositiveMonteCarloRunError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set Monte Carlo run to invalid (non-positive) number.
    monteCarloRun = -1;

    // Try to create random walk perturber.
    try { database::RandomWalkPerturberPointer testParticleCase = getRandomWalkPerturber( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk perturber failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk perturber with non-positive test particle simulation number.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberNonPositiveTestParticleSimulationNumberError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set test particle simulation number to invalid (non-positive) number.
    testParticleSimulationNumber = -1;

    // Try to create random walk perturber.
    try { database::RandomWalkPerturberPointer testParticleCase = getRandomWalkPerturber( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk perturber failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk perturber with negative mass factor.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberNegativeMassRatioError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set mass factor to invalid (negative) number.
    massRatio = -1.0;

    // Try to create random walk perturber.
    try { database::RandomWalkPerturberPointer randomWalkPerturber = getRandomWalkPerturber( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk perturber failed.
    BOOST_CHECK( isError );
}

//! Test initialization of random walk perturber with >1.0 mass factor.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberGreaterThanOneMassRatioError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set mass factor to invalid (greater than 1.0) number.
    massRatio = 1.1;

    // Try to create random walk perturber.
    try { database::RandomWalkPerturberPointer randomWalkPerturber = getRandomWalkPerturber( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of random walk perturber failed.
    BOOST_CHECK( isError );
}

//! Test comparison of RandomWalkPerturber pointers using overloaded == operator.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberEqualComparison )
{
    // Create RandomWalkPerturber pointers.
    database::RandomWalkPerturberPointer randomWalkPerturber = getRandomWalkPerturber( );
    database::RandomWalkPerturberPointer randomWalkPerturber2 = getRandomWalkPerturber( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkPerturber == *randomWalkPerturber2 );
}

//! Test comparison of RandomWalkPerturber pointers using overloaded != operator.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberNonEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkPerturber pointers.
    database::RandomWalkPerturberPointer randomWalkPerturber = getRandomWalkPerturber( );

    monteCarloRun = 2;
    database::RandomWalkPerturberPointer randomWalkPerturber2 = getRandomWalkPerturber( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkPerturber != *randomWalkPerturber2 );
}

//! Test comparison of RandomWalkPerturber pointers using overloaded < operator.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberLessThanComparison )
{
    // Create RandomWalkPerturber pointers.
    database::RandomWalkPerturberPointer randomWalkPerturber = getRandomWalkPerturber( );

    monteCarloRun = 2;
    database::RandomWalkPerturberPointer randomWalkPerturber2 = getRandomWalkPerturber( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkPerturber < *randomWalkPerturber2 );
}

//! Test comparison of RandomWalkPerturber pointers using overloaded > operator.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberGreaterThanComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkPerturber pointers.
    database::RandomWalkPerturberPointer randomWalkPerturber = getRandomWalkPerturber( );

    monteCarloRun = 2;
    database::RandomWalkPerturberPointer randomWalkPerturber2 = getRandomWalkPerturber( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkPerturber2 > *randomWalkPerturber );
}

//! Test comparison of RandomWalkPerturber pointers using overloaded <= operator.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberLessThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkPerturber pointers.
    database::RandomWalkPerturberPointer randomWalkPerturber = getRandomWalkPerturber( );

    monteCarloRun = 2;
    database::RandomWalkPerturberPointer randomWalkPerturber2 = getRandomWalkPerturber( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkPerturber <= *randomWalkPerturber2 );
}

//! Test comparison of RandomWalkPerturber pointers using overloaded >= operator.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberGreaterThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create RandomWalkPerturber pointers.
    database::RandomWalkPerturberPointer randomWalkPerturber = getRandomWalkPerturber( );

    monteCarloRun = 2;
    database::RandomWalkPerturberPointer randomWalkPerturber2 = getRandomWalkPerturber( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *randomWalkPerturber2 >= *randomWalkPerturber );
}

//! Test sorting of RandomWalkPerturberTable.
BOOST_AUTO_TEST_CASE( testRandomWalkPerturberTableSorting )
{
    // Add RandomWalkPerturber objects to table.
    database::RandomWalkPerturberTable perturberTable;

    // Insert random walk perturber struct.
    perturberTable.insert( new database::RandomWalkPerturber( *getRandomWalkPerturber( ) ) );

    // Insert random walk perturber struct with another Monte Carlo run.
    monteCarloRun = 3;
    perturberTable.insert( new database::RandomWalkPerturber( *getRandomWalkPerturber( ) ) );

    // Insert random walk perturber struct with another Monte Carlo run.
    monteCarloRun = 2;
    perturberTable.insert( new database::RandomWalkPerturber( *getRandomWalkPerturber( ) ) );

    // Check that the table is sorted according to Monte Carlo run.
    database::RandomWalkPerturberTable::iterator iteratorPerturberTable = perturberTable.begin( );
    BOOST_CHECK_EQUAL( iteratorPerturberTable->monteCarloRun, 1 );

    iteratorPerturberTable++;
    BOOST_CHECK_EQUAL( iteratorPerturberTable->monteCarloRun, 2 );

    iteratorPerturberTable++;
    BOOST_CHECK_EQUAL( iteratorPerturberTable->monteCarloRun, 3 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stochastic_migration  
