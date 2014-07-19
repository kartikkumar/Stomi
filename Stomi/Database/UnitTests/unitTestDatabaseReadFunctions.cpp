/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

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
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Stomi/Database/databaseReadFunctions.h"
#include "Stomi/Database/randomWalkCase.h"
#include "Stomi/Database/testParticleCase.h"
#include "Stomi/Database/testParticleInput.h"
#include "Stomi/Database/testParticleKick.h"
#include "Stomi/InputOutput/rootPath.h"

namespace stomi
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_database_read_functions )

//! Test implementation of function to get test particle case data from SQLite3 database.
BOOST_AUTO_TEST_CASE( testGetTestParticleCaseFunction )
{
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_mathematics::mathematical_constants;
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase = getStomiRootPath( )
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
    BOOST_CHECK_EQUAL( testParticleCase->semiMajorAxisDistributionLimit, 508299.21206805098 );
    BOOST_CHECK_EQUAL( testParticleCase->synodicPeriodMaximum, 1577880000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->startUpIntegrationPeriod, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->centralBodyJ2GravityCoefficient, 0.0 );    
    BOOST_CHECK_EQUAL( testParticleCase->centralBodyEquatorialRadius, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->conjunctionEventDetectionDistance, 48868000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->oppositionEventDetectionDistance, 146604000.0 );
    BOOST_CHECK_EQUAL( testParticleCase->eccentricityDistributionMean, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->eccentricityDistributionFullWidthHalfMaximum, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->inclinationDistributionMean, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->inclinationDistributionFullWidthHalfMaximum, 0.0 );
    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorType, DOPRI853 );
    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorInitialStepSize, 60.0 ); 
    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorRelativeTolerance, 1.0e-12 );
    BOOST_CHECK_EQUAL( testParticleCase->numericalIntegratorAbsoluteTolerance, 1.0e-7 );
}

//! Test run-time error in case of multiple identical test particle cases defined in SQLite3 
//! database.
BOOST_AUTO_TEST_CASE( testGetTestParticleCaseFunctionExtraRow )
{
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStomiRootPath( )
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
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleInputTable.sqlite";

    // Set requested test particle case ID.
    const int testParticleCaseId = 1;

    // Retrieve table of input data for test particle simulations.
    const TestParticleInputTable testParticleInputTable
            = getCompleteTestParticleInputTable( 
                absolutePathToTestDatabase, testParticleCaseId, "test_particle_input" );

    // Read in table of test particle input data from test data file.
    const Eigen::Matrix< double, 10, 9 > testDataTestParticleInputTable
            = readMatrixFromFile( getStomiRootPath( )
                                  + "/Database/UnitTests/testDataTestParticleInputTable.csv" );

    // Check that the input data table retrieved matches the test data.
    unsigned int i = 0;

    for ( TestParticleInputTable::iterator iteratorInputTable = testParticleInputTable.begin( );
          iteratorInputTable != testParticleInputTable.end( ); iteratorInputTable++ )
    {
        BOOST_CHECK_EQUAL( iteratorInputTable->simulationId,
                           testDataTestParticleInputTable( i, 0 ) );
        BOOST_CHECK_EQUAL( iteratorInputTable->testParticleCaseId, testParticleCaseId );
        BOOST_CHECK_EQUAL( iteratorInputTable->isCompleted, false );

        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        iteratorInputTable->initialStateInKeplerianElements,
                        testDataTestParticleInputTable.block( i, 3, 1, 6 ).transpose( ),
                        1.0e-14 );
        }

        i++;
    }
}

//! Test run-time error in case of no fetched test particle input data from SQLite3 database,
//! when requesting data for completed simulations.
BOOST_AUTO_TEST_CASE( testGetTestParticleInputTableFunctionNoRows )
{
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStomiRootPath( )
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
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleInputTable.sqlite";

    // Set string of selected test particle simulation IDs.
    const std::string testParticleSimulationIds = "1 3 5 9";

    // Set requested test particle case ID.
    const int testParticleCaseId = 1;

    // Retrieve table of input data for selected test particle simulations.
    const TestParticleInputTable testParticleInputTable = getSelectedTestParticleInputTable(
                absolutePathToTestDatabase, testParticleCaseId, 
                testParticleSimulationIds, "test_particle_input" );

    // Read in table of test particle input data from test data file.
    const Eigen::Matrix< double, 4, 9 > testDataTestParticleInputTable = readMatrixFromFile(
                getStomiRootPath( )
                + "/Database/UnitTests/testDataSelectedTestParticleInputTable.csv" ); 

    // Check that the input data table retrieved matches the test data.
    unsigned int i = 0;

    for ( TestParticleInputTable::iterator iteratorInputTable = testParticleInputTable.begin( );
          iteratorInputTable != testParticleInputTable.end( ); iteratorInputTable++ )
    {
        BOOST_CHECK_EQUAL( iteratorInputTable->simulationId,
                           testDataTestParticleInputTable( i, 0 ) );
        BOOST_CHECK_EQUAL( iteratorInputTable->testParticleCaseId, testParticleCaseId );
        BOOST_CHECK_EQUAL( iteratorInputTable->isCompleted, false );

        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        iteratorInputTable->initialStateInKeplerianElements,
                        testDataTestParticleInputTable.block( i, 3, 1, 6 ).transpose( ),
                        1.0e-14 );
        }

        i++;
    }
}

//! Test run-time error in case of no fetched test particle input data from SQLite3 database,
//! when requesting non-existent test particle simulations.
BOOST_AUTO_TEST_CASE( testGetTestParticleInputTableFunctionNonExistentSimulation )
{
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleInputTable.sqlite";

    // Set string of selected test particle simulation numbers.
    const std::string testParticleSimulationIds = "1 999999";

    // Set requested test particle case ID.
    const int testParticleCaseId = 1;

    // Try to retrieve table of test particle input data, and catch expected error.
    bool isRunTimeErrorThrown = false;

    try
    {
        // Retrieve all input data for selected test particle simulations.
        const TestParticleInputTable testParticleInputTable
                = getSelectedTestParticleInputTable( 
                    absolutePathToTestDatabase, testParticleCaseId,
                    testParticleSimulationIds, "test_particle_input" );
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
    using boost::assign::list_of;

    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using tudat::input_output::readMatrixFromFile;

    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleKickTable.sqlite";

    // Set random walk simulation period [s].
    const double randomWalkSimulationPeriod = 1577880000.0;

    // Set vector of selected test particle simulation numbers.
    const std::vector< int > selectedSimulationIds = list_of( 1 )( 5 )( 9 );

    // Retrieve table of test particle kick data.
    const TestParticleKickTable testParticleKickTable = getTestParticleKickTable(
                absolutePathToTestDatabase, randomWalkSimulationPeriod,
                selectedSimulationIds, "test_particle_kicks" );

    // Read in table of test particle kick data from test data file.
    const Eigen::MatrixXd testDataTestParticleKickTable
            = readMatrixFromFile( 
                getStomiRootPath( )
                + "/Database/UnitTests/testDataSelectedTestParticleKickTable.csv" );

    // Check that the kick data table retrieved matches the test data.
    unsigned int i = 0; 

    for ( TestParticleKickTable::iterator iteratorKickTable = testParticleKickTable.begin( );
          iteratorKickTable != testParticleKickTable.end( ); iteratorKickTable++ )
    {
        BOOST_CHECK_EQUAL( iteratorKickTable->testParticleSimulationId,
                           testDataTestParticleKickTable( i, 1 ) );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->conjunctionEpoch,
                                    testDataTestParticleKickTable( i, 2 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->conjunctionDistance,
                                    testDataTestParticleKickTable( i, 3 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionEpoch,
                                    testDataTestParticleKickTable( i, 4 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionDistance,
                                    testDataTestParticleKickTable( i, 5 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionStateInKeplerianElements( 
                                        semiMajorAxisIndex ),
                                    testDataTestParticleKickTable( i, 6 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionStateInKeplerianElements( 
                                        eccentricityIndex ),
                                    testDataTestParticleKickTable( i, 7 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionStateInKeplerianElements( 
                                        inclinationIndex ),
                                    testDataTestParticleKickTable( i, 8 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionStateInKeplerianElements( 
                                        argumentOfPeriapsisIndex ),
                                    testDataTestParticleKickTable( i, 9 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionStateInKeplerianElements( 
                                        longitudeOfAscendingNodeIndex ),
                                    testDataTestParticleKickTable( i, 10 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->preConjunctionStateInKeplerianElements( 
                                        trueAnomalyIndex ),
                                    testDataTestParticleKickTable( i, 11 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionEpoch,
                                    testDataTestParticleKickTable( i, 12 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionDistance,
                                    testDataTestParticleKickTable( i, 13 ), 1.0e-14 );        

        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionStateInKeplerianElements( 
                                        semiMajorAxisIndex ),
                                    testDataTestParticleKickTable( i, 14 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionStateInKeplerianElements( 
                                        eccentricityIndex ),
                                    testDataTestParticleKickTable( i, 15 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionStateInKeplerianElements( 
                                        inclinationIndex ),
                                    testDataTestParticleKickTable( i, 16 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionStateInKeplerianElements( 
                                        argumentOfPeriapsisIndex ),
                                    testDataTestParticleKickTable( i, 17 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionStateInKeplerianElements( 
                                        longitudeOfAscendingNodeIndex ),
                                    testDataTestParticleKickTable( i, 18 ), 1.0e-14 );
        BOOST_CHECK_CLOSE_FRACTION( iteratorKickTable->postConjunctionStateInKeplerianElements( 
                                        trueAnomalyIndex ),
                                    testDataTestParticleKickTable( i, 19 ), 1.0e-14 );

        i++;
    }
}

//! Test expected run-time error thrown by function to get selected test particle kick data from
//! SQLite3 database.
BOOST_AUTO_TEST_CASE( testGetTestParticleKickTableFunctionNonExistentSimulation )
{
    using boost::assign::list_of;

    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseTestParticleKickTable.sqlite";

    // Set random walk simulation period [s].
    const double randomWalkSimulationPeriod = 1577880000.0;

    // Set vector of selected test particle simulation numbers.
    const std::vector< int > selectedSimulationIds = list_of( 1 )( 9999 );

    // Try to retrieve table of test particle kick data, and catch expected error.
    bool isRunTimeErrorThrown = false;

    try
    {
        // Retrieve table of test particle kick data.
        const TestParticleKickTable testParticleKickTable = getTestParticleKickTable(
                    absolutePathToTestDatabase, randomWalkSimulationPeriod,
                    selectedSimulationIds, "test_particle_kicks" );
    }

    catch( std::runtime_error& )
    {
        isRunTimeErrorThrown = true;
    }

    // Check that expected run-time error was thrown due to request for non-existent test particle
    // simulation number.
    BOOST_CHECK( isRunTimeErrorThrown );
}

//! Test implementation of function to get random walk case data from SQLite3 database.
BOOST_AUTO_TEST_CASE( testGetRandomWalkCaseFunction )
{
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseRandomWalkCase.sqlite";

    // Retrieve random walk case data.
    const RandomWalkCasePointer randomWalkCase 
        = getRandomWalkCase( absolutePathToTestDatabase, 
                             "circular_equatorial_nominal", 
                             "random_walk_case" );

    // Check that the values read from the database are correct.
    BOOST_CHECK_EQUAL( randomWalkCase->caseId, 1 );
    BOOST_CHECK_EQUAL( randomWalkCase->caseName, "circular_equatorial_nominal" );
    BOOST_CHECK_EQUAL( randomWalkCase->testParticleCaseId, 1 );
    BOOST_CHECK_EQUAL( randomWalkCase->perturberDensity, 1.0 );
    BOOST_CHECK_EQUAL( randomWalkCase->perturberRingMass, 1.0 );  
    BOOST_CHECK_EQUAL( randomWalkCase->observationPeriod, 94672800.0 );    
    BOOST_CHECK_EQUAL( randomWalkCase->numberOfEpochWindows, 4 );    
    BOOST_CHECK_EQUAL( randomWalkCase->epochWindowSize, 7776000.0 );    
}

//! Test run-time error in case of multiple identical random walk cases defined in SQLite3 
//! database.
BOOST_AUTO_TEST_CASE( testGetRandomWalkCaseFunctionExtraRow )
{
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseRandomWalkCaseMultipleError.sqlite";

    // Try to retrieve case data.
    bool isExtraRowPresent = false;

    try
    {
        // Retrieve test particle case data.
        const RandomWalkCasePointer randomWalkCase
                = getRandomWalkCase( absolutePathToTestDatabase, 
                                     "circular_equatorial_nominal", 
                                     "random_walk_case" );
    }

    // Catch expected run-time error.
    catch( std::runtime_error& )
    {
        isExtraRowPresent = true;
    }

    // Check that the expected run-time error was thrown.
    BOOST_CHECK( isExtraRowPresent );
}

//! Test implementation of function to get incomplete simulations from random walk input table in
//! SQLite3 database for a given case ID.
BOOST_AUTO_TEST_CASE( testGetRandomWalkInputTableFunctionIncompleteSimulations )
{
    using tudat::input_output::readMatrixFromFile;
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseRandomWalkInputTable.sqlite";

    // Set requested random walk case ID.
    const int randomWalkCaseId = 1;

    // // Retrieve table of input data for test particle simulations.
    const RandomWalkInputTable randomWalkInputTable 
        = getCompleteRandomWalkInputTable(
            absolutePathToTestDatabase, randomWalkCaseId,
            "random_walk_input", "random_walk_perturbers" );

    // Read in table of random walk input data from test data file.
    const Eigen::Matrix< double, 10, 4 > testDataRandomWalkInputTable
            = readMatrixFromFile( getStomiRootPath( )
                                  + "/Database/UnitTests/testDataRandomWalkInputTable.csv" );

    // Read in table of random walk perturber data from test data file.
    const Eigen::Matrix< double, 1000, 3 > testDataRandomWalkPerturberTable
            = readMatrixFromFile( getStomiRootPath( )
                                  + "/Database/UnitTests/testDataRandomWalkPerturberTable.csv" ); 

    // Check that the input data table retrieved matches the test data.
    unsigned int i = 0;
    unsigned int j = 0;    

    for ( RandomWalkInputTable::iterator iteratorInputTable = randomWalkInputTable.begin( );
          iteratorInputTable != randomWalkInputTable.end( ); iteratorInputTable++ )
    {
        BOOST_CHECK_EQUAL( iteratorInputTable->monteCarloRunId,
                           testDataRandomWalkInputTable( i, 0 ) );
        BOOST_CHECK_EQUAL( iteratorInputTable->randomWalkCaseId, 
                           testDataRandomWalkInputTable( i, 1 ) );
        BOOST_CHECK_EQUAL( iteratorInputTable->isCompleted, false );
        BOOST_CHECK_CLOSE_FRACTION( 
                           iteratorInputTable->observationPeriodStartEpoch,
                           testDataRandomWalkInputTable( i, 3 ),
                           1.0e-14 );

        for ( std::vector< int >::const_iterator iteratorTestParticleSimulationIds
                = iteratorInputTable->testParticleSimulationIds.begin( );
              iteratorTestParticleSimulationIds 
                != iteratorInputTable->testParticleSimulationIds.end( );
              iteratorTestParticleSimulationIds++  )
        {           
            BOOST_CHECK_EQUAL( *iteratorTestParticleSimulationIds,
                               testDataRandomWalkPerturberTable( j, 2 ) );  

            j++;                                      
        }

        i++;
    }
}

//! Test run-time error in case of no fetched random wlak input data from SQLite3 database,
//! when requesting data for completed simulations.
BOOST_AUTO_TEST_CASE( testGetRandomWalkInputTableFunctionNoRows )
{
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseRandomWalkInputTable.sqlite";


    // Try to retrieve random walk input data.
    bool isNoRowPresent = false;

    try
    {
        // Retrieve all input data for random walk simulations that are complete.
        const RandomWalkInputTable randomWalkInputTable
                = getCompleteRandomWalkInputTable( 
                    absolutePathToTestDatabase, 1, 
                    "random_walk_input", "random_walk_perturbers", true );
    }

    // Catch expected run-time error.
    catch( std::runtime_error& )
    {
        isNoRowPresent = true;
    }

    // Check that the expected run-time error was thrown.
    BOOST_CHECK( isNoRowPresent );
}

//! Test implementation of function to get selected random walk Monte Carlo run data from SQLite3
//! database.
BOOST_AUTO_TEST_CASE( testGetRandomWalkInputTableFunctionSpecificRuns )
{
    using tudat::input_output::readMatrixFromFile;
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseRandomWalkInputTable.sqlite";

    // Set vector of selected Monte Carlo run IDs.
    const std::string monteCarloRunIds = "1 3 5 9";

    // Set requested random walk case ID.
    const int randomWalkCaseId = 1;

    // Retrieve table of random walk Monte Carlo runs.
    const RandomWalkInputTable randomWalkInputTable
            = getSelectedRandomWalkInputTable( 
                absolutePathToTestDatabase, randomWalkCaseId, monteCarloRunIds,
                "random_walk_input", "random_walk_perturbers" );         

    // Read in table of random walk input data from test data file.
    const Eigen::Matrix< double, 4, 4 > testDataRandomWalkInputTable = readMatrixFromFile(
                getStomiRootPath( )
                + "/Database/UnitTests/testDataSelectedRandomWalkInputTable.csv" ); 

    // Read in table of random walk perturber data from test data file.
    const Eigen::Matrix< double, 400, 3 > testDataRandomWalkPerturberTable
            = readMatrixFromFile( 
                getStomiRootPath( )
                + "/Database/UnitTests/testDataSelectedRandomWalkPerturberTable.csv" ); 

    // Check that the input data table retrieved matches the test data.
    unsigned int i = 0;
    unsigned int j = 0;    

    for ( RandomWalkInputTable::iterator iteratorInputTable = randomWalkInputTable.begin( );
          iteratorInputTable != randomWalkInputTable.end( ); iteratorInputTable++ )
    {
        BOOST_CHECK_EQUAL( iteratorInputTable->monteCarloRunId,
                           testDataRandomWalkInputTable( i, 0 ) );
        BOOST_CHECK_EQUAL( iteratorInputTable->randomWalkCaseId, 
                           testDataRandomWalkInputTable( i, 1 ) );
        BOOST_CHECK_EQUAL( iteratorInputTable->isCompleted, false );
        BOOST_CHECK_CLOSE_FRACTION( 
                           iteratorInputTable->observationPeriodStartEpoch,
                           testDataRandomWalkInputTable( i, 3 ),
                           1.0e-14 );

        for ( std::vector< int >::const_iterator iteratorTestParticleSimulationIds
                = iteratorInputTable->testParticleSimulationIds.begin( );
              iteratorTestParticleSimulationIds 
                != iteratorInputTable->testParticleSimulationIds.end( );
              iteratorTestParticleSimulationIds++  )
        {           
            BOOST_CHECK_EQUAL( *iteratorTestParticleSimulationIds,
                               testDataRandomWalkPerturberTable( j, 2 ) );  

            j++;                                      
        }

        i++;
    }    
}

//! Test run-time error in case of no fetched random walk input data from SQLite3 database,
//! when requesting non-existent random walk Monte Carlo runs.
BOOST_AUTO_TEST_CASE( testGetRandomWalkInputTableFunctionSpecificNonExistentRun )
{
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseRandomWalkInputTable.sqlite";

    // Set vector of selected Monte Carlo run IDs.
    const std::string monteCarloRunIds = "1 99999";

    // Set requested random walk case ID.
    const int randomWalkCaseId = 1;

    // Try to retrieve table of random walk Monte Carlo run data, and catch expected error.
    bool isRunTimeErrorThrown = false;

    try
    {
        const RandomWalkInputTable randomWalkInputTable
                = getSelectedRandomWalkInputTable( 
                    absolutePathToTestDatabase, randomWalkCaseId, monteCarloRunIds,
                    "random_walk_input", "random_walk_perturbers" );    
    }

    catch( std::runtime_error& )
    {
        isRunTimeErrorThrown = true;
    }

    // Check that expected run-time error was thrown due to request for non-existent Monte Carlo
    // run.
    BOOST_CHECK( isRunTimeErrorThrown );
}

//! Test implementation of function to get selected random walk perturbers from SQLite3
//! database.
BOOST_AUTO_TEST_CASE( testGetRandomWalkPerturbersFunctionSpecificRun )
{
    using tudat::input_output::readMatrixFromFile;
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseRandomWalkInputTable.sqlite";

    // Set Monte Carlo run to retrieve data for.
    const unsigned int monteCarloRun = 1;

    // Retrieve list of perturbers, expressed as vector of test particle simulation IDs for random 
    // walk.
    const std::vector< int > randomWalkPerturberList 
        = getRandomWalkPerturberList(
                absolutePathToTestDatabase, monteCarloRun, "random_walk_perturbers" );

    // Read in list of random walk perturbers from test data file.
    const Eigen::MatrixXd testRandomWalkPerturberData
        = readMatrixFromFile(
                getStomiRootPath( )
                + "/Database/UnitTests/testDataRandomWalkSinglePerturberList.csv" );

    // Check that the random walk perturber list retrieved matches the test data.
    for ( unsigned int i = 0; i < randomWalkPerturberList.size( ); i++ )
    {
        BOOST_CHECK_EQUAL( randomWalkPerturberList.at( i ), testRandomWalkPerturberData( i, 2 ) );
    }
}

//! Test expected run-time error thrown by function to get random walk perturber data from SQLite3
//! database.
BOOST_AUTO_TEST_CASE( testGetRandomWalkPerturbersFunctionNonExistentRun )
{
    using namespace input_output;
    using namespace database;

    // Set absolute path to test database.
    const std::string absolutePathToTestDatabase
            = getStomiRootPath( )
            + "/Database/UnitTests/testDatabaseRandomWalkInputTable.sqlite";

    // Set selected Monte Carlo run ID.
    const unsigned int monteCarloRun = 999999;

    // Try to retrieve list of random walk perturbers, and catch expected error.
    bool isRunTimeErrorThrown = false;

    try
    {
        // Retrieve list of perturbers for Monte Carlo run.
        const std::vector< int > randomWalkPerturberList 
            = getRandomWalkPerturberList(
                    absolutePathToTestDatabase, monteCarloRun, "random_walk_perturbers" );        
    }

    catch( std::runtime_error& )
    {
        isRunTimeErrorThrown = true;
    }

    // Check that expected run-time error was thrown due to request for non-existent Monte Carlo
    // run.
    BOOST_CHECK( isRunTimeErrorThrown );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stomi
