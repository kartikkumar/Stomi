/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <stdexcept>
#include <sstream>
#include <vector>

#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/lexical_cast.hpp>

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
                                             const std::string& caseName,
                                             const std::string& testParticleCaseTableName )
{
    // Set stream with database query.
    std::ostringstream testParticleCaseQuery;
    testParticleCaseQuery << "SELECT * FROM " <<  testParticleCaseTableName 
                          << " WHERE \"caseName\" = \"" << caseName << "\";";

    // Open database in read-only mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READONLY );

    // Set up database query.
    SQLite::Statement query( database, testParticleCaseQuery.str( ).c_str( ));

    // Get row of case data.
    query.executeStep( );

    // Store data in test particle case object.
    const TestParticleCasePointer testParticleCase = boost::make_shared< TestParticleCase >(
        TestParticleCase( query.getColumn( 0 ), query.getColumn( 1 ), query.getColumn( 2 ),  
                          query.getColumn( 3 ), query.getColumn( 4 ), query.getColumn( 5 ),
                          ( Eigen::VectorXd( 6 ) << 
                                query.getColumn( 6 ), query.getColumn( 7 ),
                                query.getColumn( 8 ), query.getColumn( 9 ), 
                                query.getColumn( 10 ), query.getColumn( 11 ) ).finished( ),                          
                          query.getColumn( 12 ), query.getColumn( 13 ), query.getColumn( 14 ), 
                          query.getColumn( 15 ), query.getColumn( 16 ), query.getColumn( 17 ),
                          query.getColumn( 18 ), query.getColumn( 19 ), query.getColumn( 20 ), 
                          query.getColumn( 21 ), query.getColumn( 22 ), query.getColumn( 23 ),
                          query.getColumn( 24 ), query.getColumn( 25 ), query.getColumn( 26 ) ) );

    // Throw an error if there are multiple rows present in the table.
    if ( query.executeStep( ) )
    {
        throw std::runtime_error( "Multiple rows for case in table!" );
    }

    return testParticleCase;
}

// //! Get test particle input table.
TestParticleInputTable getCompleteTestParticleInputTable(
        const std::string& databaseAbsolutePath, const int caseId,
        const std::string& testParticleInputTableName, bool isCompleted )
{
    // Set stream with query.
    std::ostringstream testParticleInputQuery;
    testParticleInputQuery << "SELECT * FROM " << testParticleInputTableName
                           << " WHERE \"testParticleCaseId\" = " << caseId 
                           << " AND \"completed\" = " << isCompleted << ";";

    // Open database in read-only mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READONLY );

    // Set up database query.
    SQLite::Statement query( database, testParticleInputQuery.str( ).c_str( ));

    // Declare test particle input table.
    TestParticleInputTable testParticleInputTable;

    // Loop through the table retrieved from the database, step-by-step.
    while ( query.executeStep( ) )
    {
        // Store fetched row in test particle input struct.
        testParticleInputTable.insert(
            new TestParticleInput(
                            query.getColumn( 0 ), query.getColumn( 1 ),
                            boost::lexical_cast< bool >( query.getColumn( 2 ) ),
                            ( Eigen::VectorXd( 6 ) << query.getColumn( 3 ),
                                query.getColumn( 4 ),
                                query.getColumn( 5 ),
                                query.getColumn( 6 ),
                                query.getColumn( 7 ),
                                query.getColumn( 8 ) ).finished( ) ) );
    }

    // Check if input table is empty.
    if ( testParticleInputTable.size( ) == 0 )
    {
        // Throw run-time error.
        throw std::runtime_error( "Test particle input table is empty!" );
    }

    // Return test particle input table.
    return testParticleInputTable;
}

//! Get test particle input table.
TestParticleInputTable getSelectedTestParticleInputTable(
        const std::string& databaseAbsolutePath, const int caseId,
        const std::string& testParticleSimulationIds,
        const std::string& testParticleInputTableName )
{ 
    // Cast simulation numbers to vector of string tokens.
    std::vector< std::string > testParticleSimulationIdTokens;
    boost::split( testParticleSimulationIdTokens, testParticleSimulationIds,
                  boost::is_any_of( " " ), boost::token_compress_on );

    // Populate vector of test particle simulation numbers
    std::vector< int > testParticularSimulationIdsVector;

    for ( unsigned int i = 0; i < testParticleSimulationIdTokens.size( ); i++ )
    {
        testParticularSimulationIdsVector.push_back(
                    boost::lexical_cast< int >( testParticleSimulationIdTokens.at( i ) ) );
    }

    // Open database in read-only mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READONLY );

    // Open transaction to query database.
    SQLite::Transaction transaction( database );

    // Set up query statement.
    std::ostringstream testParticleInputQuery;
    testParticleInputQuery << "SELECT * FROM " << testParticleInputTableName
                           << " WHERE \"testParticleCaseId\" == " << caseId
                           << " AND \"simulationId\" == :simulationId;";

    // Compile a SQL query.
    SQLite::Statement query( database, testParticleInputQuery.str( ).c_str( ) );

    // Declare test particle input table.
    TestParticleInputTable testParticleInputTable;

    // Loop through the table retrieved from the database, step-by-step.
    for ( unsigned int i = 0; i < testParticularSimulationIdsVector.size( ); i++ )
    {
        // Bind simulation number to query.
        query.bind( ":simulationId", testParticularSimulationIdsVector.at( i ) );

        // Execute select query.
        // A run-time error will be thrown if the requested simulation ID can't be found.
        query.executeStep( );

        // Store fetched row in test particle input struct.
        testParticleInputTable.insert(
            new TestParticleInput(
                        query.getColumn( 0 ), query.getColumn( 1 ),
                        boost::lexical_cast< bool >( query.getColumn( 2 ) ),
                        ( Eigen::VectorXd( 6 ) << query.getColumn( 3 ), query.getColumn( 4 ), 
                          query.getColumn( 5 ), query.getColumn( 6 ), query.getColumn( 7 ),
                          query.getColumn( 8 ) ).finished( ) ) );

        // Reset query.
        query.reset( );
    }

    // Commit database transaction.
    transaction.commit( );

    // Return case simulation table.
    return testParticleInputTable;
}

//! Get test particle kick table.
TestParticleKickTable getTestParticleKickTable(
        const std::string& databaseAbsolutePath, const double randomWalkSimulationPeriod, 
        const std::vector< int >& selectedSimulationIds, 
        const std::string& testParticleKickTableName )
{
    // Open database in read-only mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READONLY );

    // Open transaction to query database.
    SQLite::Transaction transaction( database );

        // Set up query statement.
    std::ostringstream testParticleKickQuery;
    testParticleKickQuery << "SELECT * FROM " << testParticleKickTableName
                          << " WHERE \"testParticleSimulationId\" == :simulationId"
                          << " AND \"conjunctionEpoch\" > 0.0"
                          << " AND \"conjunctionEpoch\" < " << randomWalkSimulationPeriod << ";";

    // Compile a SQL query.
    SQLite::Statement query( database, testParticleKickQuery.str( ).c_str( ) );

    // Declare test particle kick table.
    TestParticleKickTable testParticleKickTable;

    // Loop through the table retrieved from the database, step-by-step.
    for ( unsigned int i = 0; i < selectedSimulationIds.size( ); i++ )
    {
        // Bind simulation number to query.
        query.bind( ":simulationId", selectedSimulationIds.at( i ) );

        // Set flag to indicate if simulation ID was found.
        bool isSimulationFound = false;

        // Execute select query.
        // A run-time error will be thrown if the requested simulation ID can't be found.
        while( query.executeStep( ) == true )
        {
            // Set flag to indicate that simulation ID was found to true.
            isSimulationFound = true;

            // Store fetched row in test particle input struct.
            testParticleKickTable.insert(
                new TestParticleKick(
                        query.getColumn( 0 ), query.getColumn( 1 ), query.getColumn( 2 ),
                        query.getColumn( 3 ), query.getColumn( 4 ), query.getColumn( 5 ),
                        ( Eigen::VectorXd( 6 ) << query.getColumn( 6 ), query.getColumn( 7 ), 
                          query.getColumn( 8 ), query.getColumn( 9 ), query.getColumn( 10 ),
                          query.getColumn( 11 ) ).finished( ),
                        query.getColumn( 12 ), query.getColumn( 13 ),
                        ( Eigen::VectorXd( 6 ) << query.getColumn( 14 ), query.getColumn( 15 ), 
                          query.getColumn( 16 ), query.getColumn( 17 ), query.getColumn( 18 ),
                          query.getColumn( 19 ) ).finished( ) ) );   
        }

        // Check if simulation ID was not found (by checking flag).
        if ( !isSimulationFound )
        {
            std::ostringstream errorMessage;
            errorMessage << "ERROR Simulation ID " << selectedSimulationIds.at( i ) << " "
                         << "was not found in the test particle kicks table!";
            throw std::runtime_error( errorMessage.str( ).c_str( ) );
        }

        // Reset query.
        query.reset( );
    }

    // Commit database transaction.
    transaction.commit( );    

    // Return aggregated test particle kick table.
    return testParticleKickTable;
}

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
//                 databaseConnector->fetchString( 2 ),  query.getColumn( 3 ),
//                 query.getColumn( 4 ), query.getColumn( 5 ),
//                 query.getColumn( 6 ),
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
//                                      query.getColumn( 3 ) ) );
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
