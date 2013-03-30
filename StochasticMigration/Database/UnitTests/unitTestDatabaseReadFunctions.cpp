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
 *      130214    K. Kumar          Added unit test to check that run-time error is thrown if
 *                                  cases table contains multiple rows; added unit test for
 *                                  input_data table get-functions().
 *      130217    K. Kumar          Optimized and completed suite of unit tests for all database
 *                                  read functions; updated "mab simulations" references to
 *                                  "stochastic migration".
 *      130218    K. Kumar          Updated "encounter" to "conjunction".
 *      130329    K. Kumar          Updated unit tests to use pointers to data structs.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <stdexcept>
#include <string>
#include <vector>

#include <boost/assign/list_of.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Basics/testMacros.h>
#include <TudatCore/InputOutput/matrixTextFileReader.h>

#include "StochasticMigration/Basics/basics.h"

#include "StochasticMigration/Database/databaseReadFunctions.h"
// #include "StochasticMigration/Database/randomWalkMonteCarloRun.h"
#include "StochasticMigration/Database/testParticleCase.h"
// #include "StochasticMigration/Database/testParticleInput.h"
// #include "StochasticMigration/Database/testParticleKick.h"


namespace stochastic_migration
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_database_read_functions )

//! Test implementation of function to get test particle case data from SQLite3 database.
BOOST_AUTO_TEST_CASE( testGetTestParticleCaseFunction )
{
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace basics;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase = getStochasticMigrationRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleCase.db";

    // Retrieve test particle case data.
    const TestParticleCasePointer testParticleCase 
        = getTestParticleCase( absolutePathToTestDatabase );

    // Check that the values read from the database are correct.
    BOOST_CHECK_EQUAL( testParticleCase->caseNumber, 1 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyGravitationalParameter, 724466.0 );
    BOOST_CHECK_EQUAL( testParticleCase->randomWalkSimulationDuration, 1577880000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->synodicPeriodLimit, 1577880000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->outputInterval, 14400.0 );
    BOOST_CHECK_EQUAL( testParticleCase->startUpIntegrationDuration, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->centralBodyGravitationalParameter, 5.79397e15 );
    BOOST_CHECK_EQUAL( testParticleCase->centralBodyEquatorialRadius, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->centralBodyJ2GravityCoefficient, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->semiMajorAxisLimit, 1000000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->meanEccentricity, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->fullWidthHalfMaxmimumEccentricityDistribution, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->meanInclination, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->fullWidthHalfMaxmimumInclinationDistribution, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                           semiMajorAxisIndex ), 97736000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                           eccentricityIndex ), 0.00254 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                           inclinationIndex ), 0.00244346 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                           argumentOfPeriapsisIndex ), 0.330904 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                           longitudeOfAscendingNodeIndex ), 4.39704 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                           trueAnomalyIndex ), 6.18747 );
    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorType, "DOPRI853" );
    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorRelativeTolerance, 1.0e-12 );
    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorAbsoluteTolerance, 1.0e-15 );
    BOOST_CHECK_EQUAL( testParticleCase->initialStepSize, 60.0 ); 
}

//! Test run-time error in case of multiple test particle cases defined in SQLite3 database.
BOOST_AUTO_TEST_CASE( testGetTestParticleCaseFunctionExtraRow )
{
    using namespace basics;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStochasticMigrationRootPath( )
            + "/Database/UnitTests/testDatabaseMultipleTestParticleCases.db";

    // Try to retrieve case data.
    bool isExtraRowPresent = false;

    try
    {
        // Retrieve test particle case data.
        const TestParticleCasePointer testParticleCase
                = getTestParticleCase( absolutePathToTestDatabase );
    }

    // Catch expected run-time error.
    catch( std::runtime_error& )
    {
        isExtraRowPresent = true;
    }

    // Check that the expected run-time error was thrown.
    BOOST_CHECK( isExtraRowPresent );
}

//! Test implementation of function to get incomplete simulations from test particle input table in
//! SQLite3 database.
BOOST_AUTO_TEST_CASE( testGetTestParticleInputTableFunctionIncompleteSimulations )
{
    using tudat::input_output::readMatrixFromFile;

    using namespace basics;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStochasticMigrationRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleInputTable.db";

    // Retrieve table of input data for test particle simulations.
    const TestParticleInputTable testParticleInputTable
            = getTestParticleInputTable( absolutePathToTestDatabase );

    // Read in table of test particle input data from test data file.
    const Eigen::Matrix< double, 10, 8 > testDataTestParticleInputTable
            = readMatrixFromFile( getStochasticMigrationRootPath( )
                                  + "/Database/UnitTests/testDataTestParticleInputTable.csv" );

    // Check that the input data table retrieved matches the test data.
    unsigned int i = 0;

    for ( TestParticleInputTable::iterator iteratorInputTable = testParticleInputTable.begin( );
          iteratorInputTable != testParticleInputTable.end( ); iteratorInputTable++ )
    {
        BOOST_CHECK_EQUAL( iteratorInputTable->simulationNumber,
                           testDataTestParticleInputTable( i, 0 ) );
        BOOST_CHECK_EQUAL( iteratorInputTable->isCompleted, false );

        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        iteratorInputTable->initialStateInKeplerianElements,
                        testDataTestParticleInputTable.block( i, 2, 1, 6 ).transpose( ),
                        1.0e-14 );
        }

        i++;
    }
}

//! Test run-time error in case of no fetched test particle input data from SQLite3 database,
//! when requesting data for completed simulations.
BOOST_AUTO_TEST_CASE( testGetTestParticleInputTableFunctionNoRows )
{
    using namespace basics;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStochasticMigrationRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleInputTable.db";

    // Try to retrieve test particle input data.
    bool isNoRowPresent = false;

    try
    {
        // Retrieve all input data for test particle simulations that are complete.
        const TestParticleInputTable testParticleInputTable
                = getTestParticleInputTable( absolutePathToTestDatabase, true );
    }

    // Catch expected run-time error.
    catch( std::runtime_error& )
    {
        isNoRowPresent = true;
    }

    // Check that the expected run-time error was thrown.
    BOOST_CHECK( isNoRowPresent );
}

//! Test implementation of function to get selected test particle input data from SQLite3 database.
BOOST_AUTO_TEST_CASE( testGetTestParticleInputTableFunctionSpecificSimulations )
{
    using tudat::input_output::readMatrixFromFile;

    using namespace basics;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStochasticMigrationRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleInputTable.db";

    // Set string of selected test particle simulation numbers.
    const std::string testParticleSimulationNumbers = "1 3 5 9";

    // Set vector of selected test particle simulation numbers.
    const std::vector< int > testParticleSimulationNumbersVector
            = boost::assign::list_of( 1 )( 3 )( 5 )( 9 );

    // Retrieve table of input data for test particle simulations.
    const TestParticleInputTable testParticleInputTable = getTestParticleInputTable(
                absolutePathToTestDatabase, testParticleSimulationNumbers );

    // Read in table of test particle input data from test data file.
    const Eigen::Matrix< double, 4, 8 > testDataTestParticleInputTable = readMatrixFromFile(
                getStochasticMigrationRootPath( )
                + "/Database/UnitTests/testDataSelectedTestParticleInputTable.csv" ); 

    // Check that the input data table retrieved matches the test data.
    unsigned int i = 0;

    for ( TestParticleInputTable::iterator iteratorInputTable = testParticleInputTable.begin( );
          iteratorInputTable != testParticleInputTable.end( ); iteratorInputTable++ )
    {
        BOOST_CHECK_EQUAL( iteratorInputTable->simulationNumber,
                           testDataTestParticleInputTable( i, 0 ) );
        BOOST_CHECK_EQUAL( iteratorInputTable->isCompleted, false );

        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        iteratorInputTable->initialStateInKeplerianElements,
                        testDataTestParticleInputTable.block( i, 2, 1, 6 ).transpose( ),
                        1.0e-14 );
        }

        i++;
    }
}

//! Test run-time error in case of no fetched test particle input data from SQLite3 database,
//! when requesting non-existent test particle simulations.
BOOST_AUTO_TEST_CASE( testGetTestParticleInputTableFunctionNonExistentSimulation )
{
    using namespace basics;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStochasticMigrationRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleInputTable.db";

    // Set string of selected test particle simulation numbers.
    const std::string testParticleSimulationNumbers = "1 999999";

    // Try to retrieve table of test particle input data, and catch expected error.
    bool isRunTimeErrorThrown = false;

    try
    {
        // Retrieve all input data for test particle simulations.
        const TestParticleInputTable testParticleInputTable
                = getTestParticleInputTable( absolutePathToTestDatabase,
                                             testParticleSimulationNumbers );
    }

    catch( std::runtime_error& )
    {
        isRunTimeErrorThrown = true;
    }

    // Check that expected run-time error was thrown due to request for non-existent test particle
    // simulation number.
    BOOST_CHECK( isRunTimeErrorThrown );
}

//! Test implementation of function to get selected test particle kick data from SQLite3 database.
BOOST_AUTO_TEST_CASE( testGetTestParticleKickTableFunctionSpecificSimulations )
{
    using boost::assign::map_list_of;

    using tudat::input_output::readMatrixFromFile;

    using namespace basics;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase = getStochasticMigrationRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleKickTable.db";

    // Set random walk simulation duration [s].
    const double randomWalkSimulationDuration = 1577880000.0;

    // Set vector of selected test particle simulation numbers.
    const TestParticleSimulationNumbersAndMassRatios testParticleSimulationNumbersAndMassRatios
            = map_list_of( 1, 0.1 )( 5, 0.2 )( 9, 0.05 );

    // Retrieve table of test particle kick data.
    const TestParticleKickTable testParticleKickTable = getTestParticleKickTable(
                absolutePathToTestDatabase, randomWalkSimulationDuration,
                testParticleSimulationNumbersAndMassRatios );

    // Read in table of test particle kick data from test data file.
    const Eigen::MatrixXd testDataTestParticleKickTable
            = readMatrixFromFile( getStochasticMigrationRootPath( )
                                  + "/Database/UnitTests/testDataTestParticleKickTable.csv" );

    // Set iterator to map of test particle simulation numbers and mass ratios to beginning.
    TestParticleSimulationNumbersAndMassRatios::const_iterator
            iteratorTestParticleSimulationNumbersAndMassRatios;

    // Check that the test particle kick table retrieved matches the test data.
    unsigned int i = 0;

    for ( TestParticleKickTable::iterator iteratorKickTable = testParticleKickTable.begin( );
          iteratorKickTable != testParticleKickTable.end( ); iteratorKickTable++ )
    {
        BOOST_CHECK_EQUAL( iteratorKickTable->simulationNumber,
                           testDataTestParticleKickTable( i, 1 ) );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->conjunctionEpoch,
                                    testDataTestParticleKickTable( i, 2 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->conjunctionDistance,
                                    testDataTestParticleKickTable( i, 3 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->conjunctionDuration,
                                    testDataTestParticleKickTable( i, 4 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionEpoch,
                                    testDataTestParticleKickTable( i, 5 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionDistance,
                                    testDataTestParticleKickTable( i, 6 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionSemiMajorAxis,
                                    testDataTestParticleKickTable( i, 7 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionEccentricity,
                                    testDataTestParticleKickTable( i, 8 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionInclination,
                                    testDataTestParticleKickTable( i, 9 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionEpoch,
                                    testDataTestParticleKickTable( i, 10 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionDistance,
                                    testDataTestParticleKickTable( i, 11 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionSemiMajorAxis,
                                    testDataTestParticleKickTable( i, 12 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionEccentricity,
                                    testDataTestParticleKickTable( i, 13 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionInclination,
                                    testDataTestParticleKickTable( i, 14 ), 1.0e-14 );

        // Find the test particle simulation number corresponding to the current row in map.
        iteratorTestParticleSimulationNumbersAndMassRatios
                = testParticleSimulationNumbersAndMassRatios.find(
                    iteratorKickTable->simulationNumber );

        BOOST_CHECK_EQUAL( iteratorKickTable->massRatio,
                           iteratorTestParticleSimulationNumbersAndMassRatios->second );

        i++;
    }
}

// //! Test expected run-time error thrown by function to get selected test particle kick data from
// //! SQLite3 database.
// BOOST_AUTO_TEST_CASE( testGetTestParticleKickTableFunctionNonExistentSimulation )
// {
//     using boost::assign::map_list_of;

//     using namespace basics;
//     using namespace database;

//     // Set absolute path to test database.
//     const std::string absolutePathToTestDatabase = getStochasticMigrationRootPath( )
//             + "/Database/UnitTests/testDatabaseTestParticleKickTable.db";

//     // Set random walk simulation duration [s].
//     const double randomWalkSimulationDuration = 1577880000.0;

//     // Set vector of selected test particle simulation numbers.
//     const TestParticleSimulationNumbersAndMassRatios testParticleSimulationNumbersAndMassRatios
//             = map_list_of( 1, 0.5 )( 999999, 0.1 );

//     // Try to retrieve table of test particle kick data, and catch expected error.
//     bool isRunTimeErrorThrown = false;

//     try
//     {
//         // Retrieve kick data for test particle simulations.
//         const TestParticleKickTable testParticleKickTable = getTestParticleKickTable(
//                     absolutePathToTestDatabase, randomWalkSimulationDuration,
//                     testParticleSimulationNumbersAndMassRatios );
//     }

//     catch( std::runtime_error& )
//     {
//         isRunTimeErrorThrown = true;
//     }

//     // Check that expected run-time error was thrown due to request for non-existent test particle
//     // simulation number.
//     BOOST_CHECK( isRunTimeErrorThrown );
// }

// //! Test implementation of function to get selected random walk Monte Carlo run data from SQLite3
// //! database.
// BOOST_AUTO_TEST_CASE( testGetRandomWalkMonteCarloRunTableFunctionSpecificRuns )
// {
//     using boost::assign::list_of;

//     using namespace basics;
//     using namespace database;

//     // Set absolute path to test database.
//     const std::string absolutePathToTestDatabase = getStochasticMigrationRootPath( )
//             + "/Database/UnitTests/testDatabaseRandomWalkMonteCarloRunTable.db";

//     // Set vector of selected Monte Carlo runs.
//     const std::vector< unsigned int > monteCarloRuns = list_of( 1 )( 5 )( 1001 )( 1005 );

//     // Retrieve table of random walk Monte Carlo runs.
//     const RandomWalkMonteCarloRunTable randomWalkMonteCarloRunTable
//             = getRandomWalkMonteCarloRunsTable( absolutePathToTestDatabase, monteCarloRuns );

//     // Set expected random walk Monte Carlo run table.
//     const RandomWalkMonteCarloRunTable expectedRandomWalkMonteCarloRunTable
//             = list_of(
//                 RandomWalkMonteCarloRun( 1, 100, "EQUAL", 0.001, 0.0, 94672800.0, 7776000.0, 4 ) )(
//                 RandomWalkMonteCarloRun( 5, 100, "EQUAL", 0.001, 0.0, 94672800.0, 7776000.0, 4 ) )(
//                 RandomWalkMonteCarloRun( 1001, 1000, "EQUAL", 0.0001,
//                                          0.0, 94672800.0, 7776000.0, 4 ) )(
//                 RandomWalkMonteCarloRun( 1005, 1000, "EQUAL", 0.0001,
//                                          0.0, 94672800.0, 7776000.0, 4 ) );

//     // Check that the random walk Monte Carlo run table retrieved matches the test data.
//     for ( unsigned int i = 0; i < expectedRandomWalkMonteCarloRunTable.size( ); i++ )
//     {
//         BOOST_CHECK_EQUAL( randomWalkMonteCarloRunTable.at( i ).monteCarloRun,
//                            expectedRandomWalkMonteCarloRunTable.at( i ).monteCarloRun );
//         BOOST_CHECK_EQUAL( randomWalkMonteCarloRunTable.at( i ).perturberPopulation,
//                            expectedRandomWalkMonteCarloRunTable.at( i ).perturberPopulation );
//         BOOST_CHECK_EQUAL( randomWalkMonteCarloRunTable.at( i ).massDistributionType,
//                            expectedRandomWalkMonteCarloRunTable.at( i ).massDistributionType );
//         BOOST_CHECK_EQUAL( randomWalkMonteCarloRunTable.at( i ).massDistributionParameter1,
//                            expectedRandomWalkMonteCarloRunTable.at(
//                                i ).massDistributionParameter1 );
//         BOOST_CHECK_EQUAL( randomWalkMonteCarloRunTable.at( i ).massDistributionParameter2,
//                            expectedRandomWalkMonteCarloRunTable.at(
//                                i ).massDistributionParameter2 );
//         BOOST_CHECK_EQUAL( randomWalkMonteCarloRunTable.at( i ).observationPeriod,
//                            expectedRandomWalkMonteCarloRunTable.at( i ).observationPeriod );
//         BOOST_CHECK_EQUAL( randomWalkMonteCarloRunTable.at( i ).epochWindowSize,
//                            expectedRandomWalkMonteCarloRunTable.at( i ).epochWindowSize );
//         BOOST_CHECK_EQUAL( randomWalkMonteCarloRunTable.at( i ).numberOfEpochWindows,
//                            expectedRandomWalkMonteCarloRunTable.at( i ).numberOfEpochWindows );
//     }
// }

// //! Test expected run-time error thrown by function to get selected random walk Monte Carlo data
// //! from SQLite3 database.
// BOOST_AUTO_TEST_CASE( testGetRandomWalkMonteCarloRunTableFunctionNonExistentRun )
// {
//     using boost::assign::list_of;

//     using namespace basics;
//     using namespace database;

//     // Set absolute path to test database.
//     const std::string absolutePathToTestDatabase = getStochasticMigrationRootPath( )
//             + "/Database/UnitTests/testDatabaseRandomWalkMonteCarloRunTable.db";

//     // Set vector of selected Monte Carlo runs.
//     const std::vector< unsigned int > monteCarloRuns = list_of( 1 )( 99999 );

//     // Try to retrieve table of random walk Monte Carlo run data, and catch expected error.
//     bool isRunTimeErrorThrown = false;

//     try
//     {
//         // Retrieve table of random walk Monte Carlo runs.
//         const RandomWalkMonteCarloRunTable randomWalkMonteCarloRunTable
//                 = getRandomWalkMonteCarloRunsTable( absolutePathToTestDatabase, monteCarloRuns );
//     }

//     catch( std::runtime_error& )
//     {
//         isRunTimeErrorThrown = true;
//     }

//     // Check that expected run-time error was thrown due to request for non-existent Monte Carlo
//     // run.
//     BOOST_CHECK( isRunTimeErrorThrown );
// }

// //! Test implementation of function to get selected random walk perturber data from SQLite3
// //! database.
// BOOST_AUTO_TEST_CASE( testGetRandomWalkPerturberTableFunctionSpecificRun )
// {
//     using tudat::input_output::readMatrixFromFile;

//     using namespace basics;
//     using namespace database;

//     // Set absolute path to test database.
//     const std::string absolutePathToTestDatabase = getStochasticMigrationRootPath( )
//             + "/Database/UnitTests/testDatabaseRandomWalkPerturberTable.db";

//     // Set Monte Carlo run to retrieve data for.
//     const unsigned int monteCarloRun = 1;

//     // Retrieve table of perturbers for random walk.
//     const RandomWalkPerturberTable randomWalkPerturberTable = getRandomWalkPerturberTable(
//                 absolutePathToTestDatabase, monteCarloRun );

//     // Read in table of random walk perturber data from test data file.
//     const Eigen::MatrixXd testRandomWalkPerturberData = readMatrixFromFile(
//                 getStochasticMigrationRootPath( )
//                 + "/Database/UnitTests/testDataRandomWalkPerturberTable.csv" );

//     // Check that the random walk perturber table retrieved matches the test data.
//     for ( int i = 0; i < testRandomWalkPerturberData.rows( ); i++ )
//     {
//         BOOST_CHECK_EQUAL( randomWalkPerturberTable.at( i ).monteCarloRun,
//                            testRandomWalkPerturberData( i, 1 ) );
//         BOOST_CHECK_EQUAL( randomWalkPerturberTable.at( i ).testParticleSimulationNumber,
//                            testRandomWalkPerturberData( i, 2 ) );
//         BOOST_CHECK_EQUAL( randomWalkPerturberTable.at( i ).massFactor,
//                            testRandomWalkPerturberData( i, 3 ) );
//     }
// }

// //! Test expected run-time error thrown by function to get random walk perturber data from SQLite3
// //! database.
// BOOST_AUTO_TEST_CASE( testGetRandomWalkPerturberTableFunctionNonExistentMonteCarloRun )
// {
//     using boost::assign::list_of;

//     using namespace basics;
//     using namespace database;

//     // Set absolute path to test database.
//     const std::string absolutePathToTestDatabase = getStochasticMigrationRootPath( )
//             + "/Database/UnitTests/testDatabaseRandomWalkPerturberTable.db";

//     // Set selected Monte Carlo run.
//     const unsigned int monteCarloRun = 999999;

//     // Try to retrieve table of random walk perturber data, and catch expected error.
//     bool isRunTimeErrorThrown = false;

//     try
//     {
//         // Retrieve table of perturbers for random walk.
//         const RandomWalkPerturberTable randomWalkPerturberTable = getRandomWalkPerturberTable(
//                     absolutePathToTestDatabase, monteCarloRun );
//     }

//     catch( std::runtime_error& )
//     {
//         isRunTimeErrorThrown = true;
//     }

//     // Check that expected run-time error was thrown due to request for non-existent Monte Carlo
//     // run.
//     BOOST_CHECK( isRunTimeErrorThrown );
// }

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stochastic_migration
