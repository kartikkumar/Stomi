/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
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
 *      130704    K. Kumar          Updated tests based on revised case and input table schemas.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

// #include <stdexcept>
#include <string>
#include <vector>

#include <boost/assign/list_of.hpp>
// #include <boost/assign/ptr_list_inserter.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// #include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Basics/testMacros.h>
#include <TudatCore/InputOutput/matrixTextFileReader.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "StochasticMigration/Basics/basics.h"
#include "StochasticMigration/Database/databaseReadFunctions.h"
// #include "StochasticMigration/Database/randomWalkMonteCarloRun.h"
#include "StochasticMigration/Database/testParticleCase.h"
#include "StochasticMigration/Database/testParticleInput.h"
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
    using namespace tudat::basic_mathematics::mathematical_constants;
    using namespace basics;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase = getStochasticMigrationRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleCase.sqlite";

    // Retrieve test particle case data.
    const TestParticleCasePointer testParticleCase 
        = getTestParticleCase( absolutePathToTestDatabase, "test_case", "test_particle_case" );

    // Check that the values read from the database are correct.
    BOOST_CHECK_EQUAL( testParticleCase->caseId, 1 );
    BOOST_CHECK_EQUAL( testParticleCase->caseName, "test_case" );
    BOOST_CHECK_EQUAL( testParticleCase->randomWalkSimulationPeriod, 1577880000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->centralBodyGravitationalParameter, 5.793966e15 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyRadius, 12000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyBulkDensity, 1500.0 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                           semiMajorAxisIndex ), 97736000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                           eccentricityIndex ), 0.00254 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                           inclinationIndex ), 0.00244346095279206 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                           argumentOfPeriapsisIndex ), 0.330903954202613 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                           longitudeOfAscendingNodeIndex ), 4.39704289113435 );
    BOOST_CHECK_EQUAL( testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                           trueAnomalyIndex ), 6.18747145100022 );
    BOOST_CHECK_EQUAL( testParticleCase->semiMajorAxisDistributionLimit, 1016598.4241361 );
    BOOST_CHECK_EQUAL( testParticleCase->synodicPeriodMaximum, 1577880000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->startUpIntegrationPeriod, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->centralBodyJ2GravityCoefficient, 0.0 );    
    BOOST_CHECK_EQUAL( testParticleCase->centralBodyEquatorialRadius, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->conjunctionEventDetectionDistance, 48868000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->oppositionEventDetectionDistance, 146604000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->eccentricityDistributionMean, 0.0 );
    BOOST_CHECK_SMALL( testParticleCase->eccentricityDistributionAngle - ( PI / 4.0 ), 1.0e-15 );    
    BOOST_CHECK_EQUAL( testParticleCase->eccentricityDistributionFullWidthHalfMaximum, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->inclinationDistributionMean, 0.0 );
    BOOST_CHECK_SMALL( testParticleCase->inclinationDistributionAngle - ( PI / 4.0 ), 1.0e-15 );    
    BOOST_CHECK_EQUAL( testParticleCase->inclinationDistributionFullWidthHalfMaximum, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorType, DOPRI853 );
    BOOST_CHECK_EQUAL( testParticleCase->initialStepSize, 60.0 ); 
    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorRelativeTolerance, 1.0e-12 );
    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorAbsoluteTolerance, 1.0e-15 );
}

//! Test run-time error in case of multiple identical test particle cases defined in SQLite3 
//! database.
BOOST_AUTO_TEST_CASE( testGetTestParticleCaseFunctionExtraRow )
{
    using namespace basics;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStochasticMigrationRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleCaseMultipleError.sqlite";

    // Try to retrieve case data.
    bool isExtraRowPresent = false;

    try
    {
        // Retrieve test particle case data.
        const TestParticleCasePointer testParticleCase
                = getTestParticleCase( absolutePathToTestDatabase,
                                       "test_case", "test_particle_case" );
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
//! SQLite3 database for a given case ID.
BOOST_AUTO_TEST_CASE( testGetTestParticleInputTableFunctionIncompleteSimulations )
{
    using tudat::input_output::readMatrixFromFile;
    using namespace basics;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStochasticMigrationRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleInputTable.sqlite";

    // Set requested case ID.
    const int caseId = 1;

    // Retrieve table of input data for test particle simulations.
    const TestParticleInputTable testParticleInputTable
            = getCompleteTestParticleInputTable( 
                absolutePathToTestDatabase, caseId, "test_particle_input" );

    // Read in table of test particle input data from test data file.
    const Eigen::Matrix< double, 10, 9 > testDataTestParticleInputTable
            = readMatrixFromFile( getStochasticMigrationRootPath( )
                                  + "/Database/UnitTests/testDataTestParticleInputTable.csv" );

    // Check that the input data table retrieved matches the test data.
    unsigned int i = 0;

    for ( TestParticleInputTable::iterator iteratorInputTable = testParticleInputTable.begin( );
          iteratorInputTable != testParticleInputTable.end( ); iteratorInputTable++ )
    {
        BOOST_CHECK_EQUAL( iteratorInputTable->simulationId,
                           testDataTestParticleInputTable( i, 0 ) );
        BOOST_CHECK_EQUAL( iteratorInputTable->caseId, caseId );
        BOOST_CHECK_EQUAL( iteratorInputTable->isCompleted, false );

        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        iteratorInputTable->initialStateInKeplerianElements,
                        testDataTestParticleInputTable.block( i, 3, 1, 6 ).transpose( ),
                        1.0e-15 );
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
            + "/Database/UnitTests/testDatabaseTestParticleInputTable.sqlite";

    // Try to retrieve test particle input data.
    bool isNoRowPresent = false;

    try
    {
        // Retrieve all input data for test particle simulations that are complete.
        const TestParticleInputTable testParticleInputTable
                = getCompleteTestParticleInputTable( 
                    absolutePathToTestDatabase, 1, "test_particle_input", true );
    }

    // Catch expected run-time error.
    catch( std::runtime_error& )
    {
        isNoRowPresent = true;
    }

    // Check that the expected run-time error was thrown.
    BOOST_CHECK( isNoRowPresent );
}

//! Test implementation of function to get selected test particle input data from SQLite3 database
//! for given case ID.
BOOST_AUTO_TEST_CASE( testGetTestParticleInputTableFunctionSpecificSimulations )
{
    using tudat::input_output::readMatrixFromFile;
    using namespace basics;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStochasticMigrationRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleInputTable.sqlite";

    // Set string of selected test particle simulation numbers.
    const std::string testParticlesimulationIds = "1 3 5 9";

    // Set vector of selected test particle simulation numbers.
    const std::vector< int > testParticlesimulationIdsVector
            = boost::assign::list_of( 1 )( 3 )( 5 )( 9 );

    // Set requested case ID.
    const int caseId = 1;

    // Retrieve table of input data for selected test particle simulations.
    const TestParticleInputTable testParticleInputTable = getSelectedTestParticleInputTable(
                absolutePathToTestDatabase, caseId, 
                testParticlesimulationIds, "test_particle_input" );

    // Read in table of test particle input data from test data file.
    const Eigen::Matrix< double, 4, 9 > testDataTestParticleInputTable = readMatrixFromFile(
                getStochasticMigrationRootPath( )
                + "/Database/UnitTests/testDataSelectedTestParticleInputTable.csv" ); 

    // Check that the input data table retrieved matches the test data.
    unsigned int i = 0;

    for ( TestParticleInputTable::iterator iteratorInputTable = testParticleInputTable.begin( );
          iteratorInputTable != testParticleInputTable.end( ); iteratorInputTable++ )
    {
        BOOST_CHECK_EQUAL( iteratorInputTable->simulationId,
                           testDataTestParticleInputTable( i, 0 ) );
        BOOST_CHECK_EQUAL( iteratorInputTable->caseId, caseId );
        BOOST_CHECK_EQUAL( iteratorInputTable->isCompleted, false );

        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        iteratorInputTable->initialStateInKeplerianElements,
                        testDataTestParticleInputTable.block( i, 3, 1, 6 ).transpose( ),
                        1.0e-15 );
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
            + "/Database/UnitTests/testDatabaseTestParticleInputTable.sqlite";

    // Set string of selected test particle simulation numbers.
    const std::string testParticlesimulationIds = "1 999999";

    // Set requested case ID.
    const int caseId = 1;

    // Try to retrieve table of test particle input data, and catch expected error.
    bool isRunTimeErrorThrown = false;

    try
    {
        // Retrieve all input data for selected test particle simulations.
        const TestParticleInputTable testParticleInputTable
                = getSelectedTestParticleInputTable( 
                    absolutePathToTestDatabase, caseId,
                    testParticlesimulationIds, "test_particle_input" );
    }

    catch( std::runtime_error& )
    {
        isRunTimeErrorThrown = true;
    }

    // Check that expected run-time error was thrown due to request for non-existent test particle
    // simulation number.
    BOOST_CHECK( isRunTimeErrorThrown );
}

// //! Test implementation of function to get selected test particle kick data from SQLite3 database.
// BOOST_AUTO_TEST_CASE( testGetTestParticleKickTableFunctionSpecificSimulations )
// {
//     using boost::assign::map_list_of;

//     using tudat::input_output::readMatrixFromFile;

//     using namespace basics;
//     using namespace database;

//     // Set absolute path to test database.
//     const std::string absolutePathToTestDatabase = getStochasticMigrationRootPath( )
//             + "/Database/UnitTests/testDatabaseTestParticleKickTable.sqlite";

//     // Set random walk simulation period [s].
//     const double randomWalkSimulationPeriod = 1577880000.0;

//     // Set vector of selected test particle simulation numbers.
//     const TestParticlesimulationIdsAndMassRatios testParticlesimulationIdsAndMassRatios
//             = map_list_of( 1, 0.1 )( 5, 0.2 )( 9, 0.05 );

//     // Retrieve table of test particle kick data.
//     const TestParticleKickTable testParticleKickTable = getTestParticleKickTable(
//                 absolutePathToTestDatabase, randomWalkSimulationPeriod,
//                 testParticlesimulationIdsAndMassRatios );

//     // Read in table of test particle kick data from test data file.
//     const Eigen::MatrixXd testDataTestParticleKickTable
//             = readMatrixFromFile( getStochasticMigrationRootPath( )
//                                   + "/Database/UnitTests/testDataTestParticleKickTable.csv" );

//     // Set iterator to map of test particle simulation numbers and mass ratios to beginning.
//     TestParticlesimulationIdsAndMassRatios::const_iterator
//             iteratorTestParticlesimulationIdsAndMassRatios;

//     // Check that the test particle kick table retrieved matches the test data.
//     unsigned int i = 0;

//     for ( TestParticleKickTable::iterator iteratorKickTable = testParticleKickTable.begin( );
//           iteratorKickTable != testParticleKickTable.end( ); iteratorKickTable++ )
//     {
//         BOOST_CHECK_EQUAL( iteratorKickTable->simulationId,
//                            testDataTestParticleKickTable( i, 1 ) );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->conjunctionEpoch,
//                                     testDataTestParticleKickTable( i, 2 ), 1.0e-14 );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->conjunctionDistance,
//                                     testDataTestParticleKickTable( i, 3 ), 1.0e-14 );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->conjunctionPeriod,
//                                     testDataTestParticleKickTable( i, 4 ), 1.0e-14 );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionEpoch,
//                                     testDataTestParticleKickTable( i, 5 ), 1.0e-14 );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionDistance,
//                                     testDataTestParticleKickTable( i, 6 ), 1.0e-14 );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionSemiMajorAxis,
//                                     testDataTestParticleKickTable( i, 7 ), 1.0e-14 );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionEccentricity,
//                                     testDataTestParticleKickTable( i, 8 ), 1.0e-14 );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionInclination,
//                                     testDataTestParticleKickTable( i, 9 ), 1.0e-14 );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionEpoch,
//                                     testDataTestParticleKickTable( i, 10 ), 1.0e-14 );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionDistance,
//                                     testDataTestParticleKickTable( i, 11 ), 1.0e-14 );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionSemiMajorAxis,
//                                     testDataTestParticleKickTable( i, 12 ), 1.0e-14 );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionEccentricity,
//                                     testDataTestParticleKickTable( i, 13 ), 1.0e-14 );
//         BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionInclination,
//                                     testDataTestParticleKickTable( i, 14 ), 1.0e-14 );

//         // Find the test particle simulation number corresponding to the current row in map.
//         iteratorTestParticlesimulationIdsAndMassRatios
//                 = testParticlesimulationIdsAndMassRatios.find(
//                     iteratorKickTable->simulationId );

//         BOOST_CHECK_EQUAL( iteratorKickTable->massRatio,
//                            iteratorTestParticlesimulationIdsAndMassRatios->second );

//         i++;
//     }
// }

// //! Test expected run-time error thrown by function to get selected test particle kick data from
// //! SQLite3 database.
// BOOST_AUTO_TEST_CASE( testGetTestParticleKickTableFunctionNonExistentSimulation )
// {
//     using boost::assign::map_list_of;

//     using namespace basics;
//     using namespace database;

//     // Set absolute path to test database.
//     const std::string absolutePathToTestDatabase = getStochasticMigrationRootPath( )
//             + "/Database/UnitTests/testDatabaseTestParticleKickTable.sqlite";

//     // Set random walk simulation period [s].
//     const double randomWalkSimulationPeriod = 1577880000.0;

//     // Set vector of selected test particle simulation numbers.
//     const TestParticlesimulationIdsAndMassRatios testParticlesimulationIdsAndMassRatios
//             = map_list_of( 1, 0.5 )( 999999, 0.1 );

//     // Try to retrieve table of test particle kick data, and catch expected error.
//     bool isRunTimeErrorThrown = false;

//     try
//     {
//         // Retrieve kick data for test particle simulations.
//         const TestParticleKickTable testParticleKickTable = getTestParticleKickTable(
//                     absolutePathToTestDatabase, randomWalkSimulationPeriod,
//                     testParticlesimulationIdsAndMassRatios );
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
//     using boost::assign::ptr_insert;

//     using namespace basics;
//     using namespace database;

//     // Set absolute path to test database.
//     const std::string absolutePathToTestDatabase = getStochasticMigrationRootPath( )
//             + "/Database/UnitTests/testDatabaseRandomWalkMonteCarloRunTable.sqlite";

//     // Set vector of selected Monte Carlo runs.
//     const std::vector< unsigned int > monteCarloRuns = list_of( 1 )( 5 )( 1001 )( 1005 );

//     // Retrieve table of random walk Monte Carlo runs.
//     const RandomWalkMonteCarloRunTable randomWalkMonteCarloRunTable
//             = getRandomWalkMonteCarloRunsTable( absolutePathToTestDatabase, monteCarloRuns );

//     // Set expected random walk Monte Carlo run table.
//     RandomWalkMonteCarloRunTable expectedRandomWalkMonteCarloRunTable;
//     expectedRandomWalkMonteCarloRunTable.insert( 
//         new RandomWalkMonteCarloRun( 1, 100, "EQUAL", 0.001, 0.0, 94672800.0, 7776000.0, 4 ) );
//     expectedRandomWalkMonteCarloRunTable.insert( 
//         new RandomWalkMonteCarloRun( 5, 100, "EQUAL", 0.001, 0.0, 94672800.0, 7776000.0, 4 ) );
//     expectedRandomWalkMonteCarloRunTable.insert( 
//         new RandomWalkMonteCarloRun( 1001, 1000, "EQUAL", 0.0001, 0.0, 
//                                      94672800.0, 7776000.0, 4 ) );
//     expectedRandomWalkMonteCarloRunTable.insert( 
//         new RandomWalkMonteCarloRun( 1005, 1000, "EQUAL", 0.0001, 0.0, 
//                                      94672800.0, 7776000.0, 4 ) );

//     // Check that the random walk Monte Carlo run table retrieved matches the test data.
//     RandomWalkMonteCarloRunTable::iterator iteratorExpectedMonteCarloRunTable 
//         = expectedRandomWalkMonteCarloRunTable.begin( );

//     for ( RandomWalkMonteCarloRunTable::iterator iteratorMonteCarloRunTable 
//           = randomWalkMonteCarloRunTable.begin( );
//           iteratorMonteCarloRunTable != randomWalkMonteCarloRunTable.end( );
//           iteratorMonteCarloRunTable++ )
//     {
//         BOOST_CHECK_EQUAL( iteratorMonteCarloRunTable->monteCarloRun,
//                            iteratorExpectedMonteCarloRunTable->monteCarloRun );
//         BOOST_CHECK_EQUAL( iteratorMonteCarloRunTable->perturberPopulation,
//                            iteratorExpectedMonteCarloRunTable->perturberPopulation );
//         BOOST_CHECK_EQUAL( iteratorMonteCarloRunTable->massDistributionType,
//                            iteratorExpectedMonteCarloRunTable->massDistributionType );
//         BOOST_CHECK_EQUAL( iteratorMonteCarloRunTable->massDistributionParameter1,
//                            iteratorExpectedMonteCarloRunTable->massDistributionParameter1 );
//         BOOST_CHECK_EQUAL( iteratorMonteCarloRunTable->massDistributionParameter2,
//                            iteratorExpectedMonteCarloRunTable->massDistributionParameter2 );
//         BOOST_CHECK_EQUAL( iteratorMonteCarloRunTable->observationPeriod,
//                            iteratorExpectedMonteCarloRunTable->observationPeriod );
//         BOOST_CHECK_EQUAL( iteratorMonteCarloRunTable->epochWindowSize,
//                            iteratorExpectedMonteCarloRunTable->epochWindowSize );
//         BOOST_CHECK_EQUAL( iteratorMonteCarloRunTable->numberOfEpochWindows,
//                            iteratorExpectedMonteCarloRunTable->numberOfEpochWindows );

//         iteratorExpectedMonteCarloRunTable++;
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
//             + "/Database/UnitTests/testDatabaseRandomWalkMonteCarloRunTable.sqlite";

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
//             + "/Database/UnitTests/testDatabaseRandomWalkPerturberTable.sqlite";

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
//     unsigned int i = 0;

//     for ( RandomWalkPerturberTable::iterator iteratorPerturberTable 
//           = randomWalkPerturberTable.begin( );
//           iteratorPerturberTable != randomWalkPerturberTable.end( );
//           iteratorPerturberTable++ )
//     {
//         BOOST_CHECK_EQUAL( iteratorPerturberTable->monteCarloRun,
//                            testRandomWalkPerturberData( i, 1 ) );
//         BOOST_CHECK_EQUAL( iteratorPerturberTable->testParticlesimulationId,
//                            testRandomWalkPerturberData( i, 2 ) );
//         BOOST_CHECK_EQUAL( iteratorPerturberTable->massRatio,
//                            testRandomWalkPerturberData( i, 3 ) );

//         i++;
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
//             + "/Database/UnitTests/testDatabaseRandomWalkPerturberTable.sqlite";

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
