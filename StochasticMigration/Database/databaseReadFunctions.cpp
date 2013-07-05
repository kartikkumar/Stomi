/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120402    K. Kumar          File created from old mabSimulationDatabaseFunctions.cpp.
 *      130212    K. Kumar          Added planetary_rings namespace; rewrote getCaseData() function
 *                                  using SQLite3DatabaseConnector.
 *      130214    K. Kumar          Split file to contain only read-functions (write-functions
 *                                  ported to new file); Rewrote getInputDataTable() function using
 *                                  SQLite3DatabaseConnector.
 *      130217    K. Kumar          Optimized and standardized implementation of database read
 *                                  functions; added implementation of auxilliary functions;
 *                                  updated "mab simulations" references to "stochastic migration".
 *      130218    K. Kumar          Updated "encounter" to "conjunction".
 *
 *    References
 *
 *    Notes
 *
 */

// #include <algorithm>
#include <sstream>
// #include <vector>

#include <boost/make_shared.hpp>
// #include <boost/algorithm/string.hpp>
// #include <boost/algorithm/string/case_conv.hpp>
// #include <boost/lexical_cast.hpp>

#include <Eigen/Core>

#include <SQLiteC++.h> 

#include "StochasticMigration/Database/databaseReadFunctions.h"

namespace stochastic_migration
{
namespace database
{

// using namespace assist::database;

//! Get test particle case data.
TestParticleCasePointer getTestParticleCase( const std::string& databaseAbsolutePath,
                                             const std::string& testParticleCaseTableName )
{
    // Set stream with database query.
    std::ostringstream testParticleCaseQuery;
    testParticleCaseQuery << "SELECT * FROM " <<  testParticleCaseTableName << ";";

    // Oopen database in read-only mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READONLY );

    // Set up database query.
    SQLite::Statement query( database, testParticleCaseQuery.str( ).c_str( ));

    // Get row of case data.
    query.executeStep( );

    // Store data in test particle case object.
    const TestParticleCasePointer testParticleCase = boost::make_shared< TestParticleCase >(
        TestParticleCase( query.getColumn( 0 ), query.getColumn( 1 ), query.getColumn( 2 ),  
                          query.getColumn( 3 ), query.getColumn( 4 ), query.getColumn( 5 ),
                          query.getColumn( 6 ), query.getColumn( 7 ), query.getColumn( 8 ), 
                          query.getColumn( 9 ), query.getColumn( 10 ), query.getColumn( 11 ),
                          query.getColumn( 12 ), query.getColumn( 13 ), query.getColumn( 14 ), 
                          query.getColumn( 15 ), query.getColumn( 16 ), query.getColumn( 17 ),
                          query.getColumn( 18 ),
                          ( Eigen::VectorXd( 6 ) << 
                                query.getColumn( 19 ), query.getColumn( 20 ),
                                query.getColumn( 21 ), query.getColumn( 22 ), 
                                query.getColumn( 23 ), query.getColumn( 24 ) ).finished( ),
                          query.getColumn( 25 ), query.getColumn( 26 ), query.getColumn( 27 ), 
                          query.getColumn( 28 ) ) );

    // Throw an error if there are multiple rows present in the table.
    if ( query.executeStep( ) )
    {
        throw std::runtime_error( "Multiple cases in table!" );
    }

    return testParticleCase;
}

// //! Get test particle input table.
// TestParticleInputTable getTestParticleInputTable( const std::string& databaseAbsolutePath,
//                                                   bool isCompleted,
//                                                   const std::string& testParticleInputTableName )
// {
//     // Set stream with query.
//     std::ostringstream testParticleInputQuery;
//     testParticleInputQuery << "SELECT * FROM " << testParticleInputTableName
//                            << " WHERE \"completed\" = " << isCompleted << ";";

//     // Open connection to database.
//     SQLite::Database database( databaseAbsolutePath.c_str( ) );

//     // Compile a SQL query.
//     SQLite::Statement query( database, testParticleInputQuery.str( ).c_str( ) );

//     // Declare test particle input table.
//     TestParticleInputTable testParticleInputTable;

//     // Loop through the table retrieved from the database, step-by-step.
//     while ( query.executeStep( ) )
//     {
//         // Store fetched row in test particle input struct.
//         testParticleInputTable.insert(
//             new TestParticleInput(
//                             query.getColumn( 0 ),
//                             boost::lexical_cast< bool >( query.getColumn( 1 ) ),
//                             ( Eigen::VectorXd( 6 ) << query.getColumn( 2 ),
//                                 query.getColumn( 3 ),
//                                 query.getColumn( 4 ),
//                                 query.getColumn( 5 ),
//                                 query.getColumn( 6 ),
//                                 query.getColumn( 7 ) ).finished( ) ) );
//     }

//     // Check if input table is empty.
//     if ( testParticleInputTable.size( ) == 0 )
//     {
//         // Throw run-time error.
//         throw std::runtime_error( "Test particle input table is empty!" );
//     }

//     // Return test particle input table.
//     return testParticleInputTable;
// }

// //! Get test particle input table.
// TestParticleInputTable getTestParticleInputTable(
//         const std::string& databaseAbsolutePath, const std::string& testParticleSimulationNumbers,
//         const std::string& testParticleInputTableName )
// {
//     // Initiate database connector.
//     Sqlite3DatabaseConnectorPointer databaseConnector
//             = initiateDatabaseConnector( databaseAbsolutePath );

//     // Cast simulation numbers to vector of string tokens.
//     std::vector< std::string > testParticleSimulationNumberTokens;
//     boost::split( testParticleSimulationNumberTokens, testParticleSimulationNumbers,
//                   boost::is_any_of( " " ), boost::token_compress_on );

//     // Set up query statement.
//     std::ostringstream testParticleInputQuery;
//     testParticleInputQuery << "SELECT * FROM " << testParticleInputTableName
//                            << " WHERE \"testParticleSimulation\" IN ("
//                            << testParticleSimulationNumberTokens.at( 0 );

//     // Loop over test particle simulation numbers and construct SQLite query.
//     for ( unsigned int i = 1; i < testParticleSimulationNumberTokens.size( ); i++ )
//     {
//         testParticleInputQuery << ", " << testParticleSimulationNumberTokens.at( i );
//     }

//     testParticleInputQuery << ");";

//     // Populate vector of test particle simulation numbers
//     std::vector< unsigned int > testParticleSimulationNumbersVector;

//     for ( unsigned int i = 0; i < testParticleSimulationNumberTokens.size( ); i++ )
//     {
//         testParticleSimulationNumbersVector.push_back(
//                     boost::lexical_cast< unsigned int >(
//                         testParticleSimulationNumberTokens.at( i ) ) );
//     }

//     // Prepare database query.
//     databaseConnector->prepare_v2( testParticleInputQuery.str( ) );

//     // Declare test particle input table.
//     TestParticleInputTable testParticleInputTable;

//     // Declare database handler status.
//     unsigned int databaseStatus = 0;

//     // Loop through the table retrieved from the database, step-by-step.
//     while ( ( databaseStatus = databaseConnector->step( ) ) == SQLITE_ROW )
//     {
//         // Store fetched row in test particle input struct.
//         testParticleInputTable.insert(
//             new TestParticleInput(
//                         databaseConnector->fetchInteger( 0 ),
//                         boost::lexical_cast< bool >( databaseConnector->fetchString( 1 ) ),
//                         ( Eigen::VectorXd( 6 ) << databaseConnector->fetchDouble( 2 ),
//                           databaseConnector->fetchDouble( 3 ),
//                           databaseConnector->fetchDouble( 4 ),
//                           databaseConnector->fetchDouble( 5 ),
//                           databaseConnector->fetchDouble( 6 ),
//                           databaseConnector->fetchDouble( 7 ) ).finished( ) ) );

//         // Delete the test particle simulation number if found in the STL vector.
//         std::vector< unsigned int >::iterator iteratorSimulationNumber
//                 = std::find( testParticleSimulationNumbersVector.begin( ),
//                              testParticleSimulationNumbersVector.end( ),
//                              testParticleInputTable.rbegin( )->simulationNumber );

//         if ( iteratorSimulationNumber != testParticleSimulationNumbersVector.end( ) )
//         {
//             testParticleSimulationNumbersVector.erase( iteratorSimulationNumber );
//         }
//     }

//     // Check if the end of the table has been reached, and whether all simulations have been found.
//     if ( databaseStatus != SQLITE_DONE
//          || ( databaseStatus == SQLITE_DONE && testParticleSimulationNumbersVector.size( ) > 0 ) )
//     {
//         // Throw run-time error.
//         throwDatabaseError( databaseConnector, databaseStatus );
//     }

//     // Terminate database connector cleanly.
//     terminateDatabaseConnector( databaseConnector );

//     // Return case simulation table.
//     return testParticleInputTable;
// }

// //! Get test particle kick table.
// TestParticleKickTable getTestParticleKickTable(
//         const std::string& databaseAbsolutePath, const double randomWalkDuration,
//         const TestParticleSimulationNumbersAndMassRatios&
//         testParticleSimulationNumbersAndMassRatios,
//         const std::string& testParticleKickTableName )
// {
//     // Initiate database connector.
//     Sqlite3DatabaseConnectorPointer databaseConnector
//             = initiateDatabaseConnector( databaseAbsolutePath );

//     // Set up query statement.
//     std::ostringstream testParticleKickQuery;
//     testParticleKickQuery << "SELECT * FROM " << testParticleKickTableName
//                           << " WHERE \"simulation\" IN (";

//     // Loop over test particle simulation numbers and construct SQLite query.
//     TestParticleSimulationNumbersAndMassRatios::const_iterator
//             testParticleSimulationsOneButLastIterator
//             = testParticleSimulationNumbersAndMassRatios.end( );
//     testParticleSimulationsOneButLastIterator--;

//     for ( TestParticleSimulationNumbersAndMassRatios::const_iterator
//           testParticleSimulationsIterator = testParticleSimulationNumbersAndMassRatios.begin( );
//           testParticleSimulationsIterator != testParticleSimulationsOneButLastIterator;
//           testParticleSimulationsIterator++ )
//     {
//         testParticleKickQuery << testParticleSimulationsIterator->first << ",";
//     }

//     testParticleKickQuery << testParticleSimulationsOneButLastIterator->first
//                          << ") AND \"conjunctionEpoch\" > 0.0 AND \"conjunctionEpoch\" <= "
//                          << randomWalkDuration << ";";

//     // Populate vector of test particle simulation numbers.
//     std::vector< unsigned int > testParticleSimulationNumbers;

//     for ( TestParticleSimulationNumbersAndMassRatios::const_iterator
//           testParticleSimulationsIterator = testParticleSimulationNumbersAndMassRatios.begin( );
//           testParticleSimulationsIterator != testParticleSimulationNumbersAndMassRatios.end( );
//           testParticleSimulationsIterator++ )
//     {
//         testParticleSimulationNumbers.push_back( testParticleSimulationsIterator->first );
//     }

//     // Prepare database query.
//     databaseConnector->prepare_v2( testParticleKickQuery.str( ) );

//     // Declare test particle kick table.
//     TestParticleKickTable testParticleKickTable;

//     // Declare database handler status.
//     unsigned int databaseStatus = 0;

//     // Loop through the table retrieved from the database, step-by-step.
//     while ( ( databaseStatus = databaseConnector->step( ) ) == SQLITE_ROW )
//     {
//         // Store simulation number.
//         int simulationNumber = databaseConnector->fetchInteger( 1 );

//         // Store fetched row in test particle kick struct.
//         testParticleKickTable.insert(
//             new TestParticleKick(
//                         simulationNumber, databaseConnector->fetchDouble( 2 ),
//                         databaseConnector->fetchDouble( 3 ), databaseConnector->fetchDouble( 4 ),
//                         databaseConnector->fetchDouble( 5 ), databaseConnector->fetchDouble( 6 ),
//                         databaseConnector->fetchDouble( 7 ), databaseConnector->fetchDouble( 8 ),
//                         databaseConnector->fetchDouble( 9 ), databaseConnector->fetchDouble( 10 ),
//                         databaseConnector->fetchDouble( 11 ), databaseConnector->fetchDouble( 12 ),
//                         databaseConnector->fetchDouble( 13 ), databaseConnector->fetchDouble( 14 ),
//                         testParticleSimulationNumbersAndMassRatios.find(
//                             databaseConnector->fetchInteger( 1 ) )->second ) );

//         // Delete the test particle simulation number if found in the STL vector.
//         std::vector< unsigned int >::iterator iteratorSimulationNumber
//                 = std::find( testParticleSimulationNumbers.begin( ),
//                              testParticleSimulationNumbers.end( ),
//                              simulationNumber );

//         if ( iteratorSimulationNumber != testParticleSimulationNumbers.end( ) )
//         {
//             testParticleSimulationNumbers.erase( iteratorSimulationNumber );
//         }
//     }

//     // Check if the end of the table has been reached, and whether all simulations have been found.
//     if ( databaseStatus != SQLITE_DONE
//          || ( databaseStatus == SQLITE_DONE && testParticleSimulationNumbers.size( ) > 0 ) )
//     {
//         // Throw run-time error.
//         throwDatabaseError( databaseConnector, databaseStatus );
//     }

//     // Terminate database connector cleanly.
//     terminateDatabaseConnector( databaseConnector );

//     // Return aggregated test particle kick table.
//     return testParticleKickTable;
// }

// //! Get table of random walk Monte Carlo runs.
// RandomWalkMonteCarloRunTable getRandomWalkMonteCarloRunsTable(
//         const std::string& databaseAbsolutePath, const std::vector< unsigned int >& monteCarloRuns,
//         const std::string& randomWalkMonteCarloRunTableName )
// {
//     // Initiate database connector.
//     Sqlite3DatabaseConnectorPointer databaseConnector
//             = initiateDatabaseConnector( databaseAbsolutePath );

//     // Set stream with query.
//     std::ostringstream randomWalkMonteCarloRunQuery;
//     randomWalkMonteCarloRunQuery << "SELECT * FROM " << randomWalkMonteCarloRunTableName
//                                  << " WHERE \"run\" in (" << monteCarloRuns.at( 0 );

//     for ( unsigned int i = 1; i < monteCarloRuns.size( ); i++ )
//     {
//         randomWalkMonteCarloRunQuery << ", " << monteCarloRuns.at( i );
//     }

//     randomWalkMonteCarloRunQuery << ");";

//     // Copy vector of Monte Carlo runs.
//     std::vector< unsigned int > monteCarloRunsCopy = monteCarloRuns;

//     // Prepare database query.
//     databaseConnector->prepare_v2( randomWalkMonteCarloRunQuery.str( ) );

//     // Declare random walk Monte Carlo run table.
//     RandomWalkMonteCarloRunTable randomWalkMonteCarloRunTable;

//     // Declare database handler status.
//     unsigned int databaseStatus = 0;

//     // Loop through the table retrieved from the database, step-by-step.
//     while ( ( databaseStatus = databaseConnector->step( ) ) == SQLITE_ROW )
//     {
//         // Store Monte Carlo run.
//         int monteCarloRun = databaseConnector->fetchInteger( 0 );

//         // Store fetched row in random walk Monte Carlo run struct.
//         randomWalkMonteCarloRunTable.insert(
//             new RandomWalkMonteCarloRun(
//                 monteCarloRun, databaseConnector->fetchInteger( 1 ),
//                 databaseConnector->fetchString( 2 ),  databaseConnector->fetchDouble( 3 ),
//                 databaseConnector->fetchDouble( 4 ), databaseConnector->fetchDouble( 5 ),
//                 databaseConnector->fetchDouble( 6 ),
//                 databaseConnector->fetchInteger( 7 ) ) );

//         // Delete the Monte Carlo run number if found in the STL vector.
//         std::vector< unsigned int >::iterator iteratorMonteCarloRunNumber
//                 = std::find( monteCarloRunsCopy.begin( ), monteCarloRunsCopy.end( ),
//                              monteCarloRun );

//         if ( iteratorMonteCarloRunNumber != monteCarloRunsCopy.end( ) )
//         {
//             monteCarloRunsCopy.erase( iteratorMonteCarloRunNumber );
//         }
//     }

//     // Check if the end of the table has been reached, and whether all Monte Carlo runs have been
//     // found.
//     if ( databaseStatus != SQLITE_DONE
//          || ( databaseStatus == SQLITE_DONE && monteCarloRunsCopy.size( ) > 0 ) )
//     {
//         // Throw run-time error.
//         throwDatabaseError( databaseConnector, databaseStatus );
//     }

//     // Terminate database connector cleanly.
//     terminateDatabaseConnector( databaseConnector );

//     // Return random walk Monte Carlo run table.
//     return randomWalkMonteCarloRunTable;
// }

// //! Get table of selected perturbers for random walk Monte Carlo run.
// RandomWalkPerturberTable getRandomWalkPerturberTable(
//         const std::string& databaseAbsolutePath, const unsigned int monteCarloRun,
//         const std::string& randomWalkPerturberTableName )
// {
//     // Initiate database connector.
//     Sqlite3DatabaseConnectorPointer databaseConnector
//             = initiateDatabaseConnector( databaseAbsolutePath );

//     // Set stream with query.
//     std::ostringstream randomWalkPerturberTableQuery;
//     randomWalkPerturberTableQuery << "SELECT * FROM " << randomWalkPerturberTableName
//                                   << " WHERE \"run\" = " << monteCarloRun << ";";

//     // Prepare database query.
//     databaseConnector->prepare_v2( randomWalkPerturberTableQuery.str( ) );

//     // Declare random walk perturber table.
//     RandomWalkPerturberTable randomWalkPerturberTable;

//     // Declare database handler status.
//     unsigned int databaseStatus = 0;

//     // Loop through the table retrieved from the database, step-by-step.
//     while ( ( databaseStatus = databaseConnector->step( ) ) == SQLITE_ROW )
//     {
//         // Store fetched row in random walk perturber struct.
//         randomWalkPerturberTable.insert(
//             new RandomWalkPerturber( databaseConnector->fetchInteger( 1 ),
//                                      databaseConnector->fetchInteger( 2 ),
//                                      databaseConnector->fetchDouble( 3 ) ) );
//     }

//     // Check if the end of the table has been reached, and whether any perturbers have been found.
//     if ( databaseStatus != SQLITE_DONE
//          || ( databaseStatus == SQLITE_DONE && randomWalkPerturberTable.size( ) == 0 ) )
//     {
//         // Throw run-time error.
//         throwDatabaseError( databaseConnector, databaseStatus );
//     }

//     // Terminate database connector cleanly.
//     terminateDatabaseConnector( databaseConnector );

//     // Return random walk perturber table for specified Monte Carlo run.
//     return randomWalkPerturberTable;
// }

} // namespace database
} // namespace stochastic_migration
