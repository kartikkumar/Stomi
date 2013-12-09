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

#include "StoMi/Database/databaseReadFunctions.h"

namespace stomi
{
namespace database
{


//! Get test particle case data.
TestParticleCasePointer getTestParticleCase( const std::string& databaseAbsolutePath, 
                                             const int caseId,
                                             const std::string& testParticleCaseTableName )
{
    // Set stream with database query.
    std::ostringstream testParticleCaseQuery;
    testParticleCaseQuery << "SELECT * FROM " <<  testParticleCaseTableName 
                          << " WHERE \"caseId\" = \"" << caseId << "\";";

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
        throw std::runtime_error( "ERROR: Multiple rows for case in table!" );
    }

    return testParticleCase;
}

//! Get test particle case data.
TestParticleCasePointer getTestParticleCase( const std::string& databaseAbsolutePath, 
                                             const std::string& caseName,
                                             const std::string& testParticleCaseTableName )
{
    // Set stream with database query.
    std::ostringstream testParticleCaseIdQuery;
    testParticleCaseIdQuery << "SELECT caseId FROM " << testParticleCaseTableName 
                            << " WHERE \"caseName\" = \"" << caseName << "\";";

    // Open database in read-only mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READONLY ); 

    // Set up database query.
    SQLite::Statement query( database, testParticleCaseIdQuery.str( ).c_str( ));

    // Get row of case data.
    query.executeStep( );

    const int caseId = query.getColumn( 0 );

    // Throw an error if there are multiple rows present in the table.
    if ( query.executeStep( ) )
    {
        throw std::runtime_error( "ERROR: Multiple rows for case in table!" );
    }
    
    // Get test particle case from database and return data object.
    return getTestParticleCase( databaseAbsolutePath, caseId, testParticleCaseTableName );   
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
        throw std::runtime_error( "ERROR: Test particle input table is empty!" );
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
    // Cast simulation IDs to vector of string tokens.
    std::vector< std::string > testParticleSimulationIdTokens;
    boost::split( testParticleSimulationIdTokens, testParticleSimulationIds,
                  boost::is_any_of( " " ), boost::token_compress_on );

    // Populate vector of test particle simulation IDs.
    std::vector< int > testParticleSimulationIdsVector;

    for ( unsigned int i = 0; i < testParticleSimulationIdTokens.size( ); i++ )
    {
        testParticleSimulationIdsVector.push_back(
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
    for ( unsigned int i = 0; i < testParticleSimulationIdsVector.size( ); i++ )
    {
        // Bind simulation ID to query.
        query.bind( ":simulationId", testParticleSimulationIdsVector.at( i ) );

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

    // Return input table.
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
                          << " WHERE \"testParticleSimulationId\" IN (";

    // Loop through the list of test particle simulation IDs.
    for ( unsigned int i = 0; i < selectedSimulationIds.size( ) - 1; i++ )
    {
        testParticleKickQuery << selectedSimulationIds.at( i ) << ",";
    }

    testParticleKickQuery << selectedSimulationIds.back( ) << ")"
                          << " AND \"conjunctionEpoch\" > 0.0"
                          << " AND \"conjunctionEpoch\" < " << randomWalkSimulationPeriod << ";";
    // Compile a SQL query.
    SQLite::Statement query( database, testParticleKickQuery.str( ).c_str( ) );

    // Declare test particle kick table.
    TestParticleKickTable testParticleKickTable;

    // // Set flag to indicate if simulation ID was found.
    // bool isSimulationFound = false;

    // Execute select query.
    // A run-time error will be thrown if the requested simulation ID can't be found.
    while( query.executeStep( ) )
    {
        // // Set flag to indicate that simulation ID was found to true.
        // isSimulationFound = true;

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

    // // Check if simulation ID was not found (by checking flag).
    // if ( !isSimulationFound )
    // {
    //     std::ostringstream errorMessage;
    //     errorMessage << "ERROR: Simulation ID " << selectedSimulationIds.at( i ) << " "
    //                  << "was not found in the test particle kicks table!";
    //     throw std::runtime_error( errorMessage.str( ).c_str( ) );
    // }

    // Commit database transaction.
    transaction.commit( );    

    // Return aggregated test particle kick table.
    return testParticleKickTable;
}

//! Get random walk case data.
RandomWalkCasePointer getRandomWalkCase( const std::string& databaseAbsolutePath, 
                                         const std::string& caseName,
                                         const std::string& randomWalkCaseTableName )
{
    // Set stream with database query.
    std::ostringstream randomWalkCaseQuery;
    randomWalkCaseQuery << "SELECT * FROM " <<  randomWalkCaseTableName 
                        << " WHERE \"caseName\" = \"" << caseName << "\";";

    // Open database in read-only mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READONLY );

    // Set up database query.
    SQLite::Statement query( database, randomWalkCaseQuery.str( ).c_str( ));

    // Get row of case data.
    query.executeStep( );

    // Store data in random walk case object.
    const RandomWalkCasePointer randomWalkCase = boost::make_shared< RandomWalkCase >(
        query.getColumn( 0 ), query.getColumn( 1 ), query.getColumn( 2 ),  
        query.getColumn( 3 ), query.getColumn( 4 ), query.getColumn( 5 ),
        query.getColumn( 6 ), query.getColumn( 7 ) );

    // Throw an error if there are multiple rows present in the table.
    if ( query.executeStep( ) )
    {
        throw std::runtime_error( "ERROR: Multiple rows for case in table!" );
    }

    return randomWalkCase;    
}

//! Get complete random walk input table.
RandomWalkInputTable getCompleteRandomWalkInputTable(
        const std::string& databaseAbsolutePath, const int caseId,
        const std::string& randomWalkInputTableName, 
        const std::string& randomWalkPerturberTableName, bool isCompleted )
{
    // Set stream with input table query.
    std::ostringstream randomWalkInputQuery;
    randomWalkInputQuery << "SELECT * FROM " << randomWalkInputTableName
                         << " WHERE \"randomWalkCaseId\" = " << caseId 
                         << " AND \"completed\" = " << isCompleted << ";";

    // Open database in read-only mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READONLY );

    // Open transaction to query database.
    SQLite::Transaction transaction( database );

    // Set up database query.
    SQLite::Statement query( database, randomWalkInputQuery.str( ).c_str( ) );

    // Declare random walk input table.
    RandomWalkInputTable randomWalkInputTable;

    // Loop through the table retrieved from the database, step-by-step.
    while ( query.executeStep( ) )
    {
        // Store Monte Carlo run ID.
        const int monteCarloRunId = query.getColumn( 0 );

        // Fetch list of perturbers.
        const std::vector< int > perturbers = getRandomWalkPerturberList( 
            databaseAbsolutePath, monteCarloRunId, randomWalkPerturberTableName );

        // Store fetched row in random walk input struct.
        randomWalkInputTable.insert(
            new RandomWalkInput( query.getColumn( 0 ), query.getColumn( 1 ),
                                 boost::lexical_cast< bool >( query.getColumn( 2 ) ), 
                                 query.getColumn( 3 ), perturbers ) );
    }

    // Commit database transaction.
    transaction.commit( );    

    // Check if input table is empty.
    if ( randomWalkInputTable.size( ) == 0 )
    {
        // Throw run-time error.
        throw std::runtime_error( "ERROR: Random walk input table is empty!" );
    }

    // Return random walk input table.
    return randomWalkInputTable;
}

//! Get selected random walk input table.
RandomWalkInputTable getSelectedRandomWalkInputTable(
        const std::string& databaseAbsolutePath, const int caseId,
        const std::string& monteCarloRunIds,
        const std::string& randomWalkInputTableName, 
        const std::string& randomWalkPerturberTableName )
{
    // Cast Monte Carlo run IDs to vector of string tokens.
    std::vector< std::string > monteCarloRunIdTokens;
    boost::split( monteCarloRunIdTokens, monteCarloRunIds,
                  boost::is_any_of( " " ), boost::token_compress_on );

    // Populate vector of random walk monte Carlo Run IDs.
    std::vector< int > monteCarloRunIdsVector;

    for ( unsigned int i = 0; i < monteCarloRunIdTokens.size( ); i++ )
    {
        monteCarloRunIdsVector.push_back(
                    boost::lexical_cast< int >( monteCarloRunIdTokens.at( i ) ) );
    }

    // Open database in read-only mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READONLY );

    // Open transaction to query database.
    SQLite::Transaction transaction( database );

    // Set up query statement.
    std::ostringstream randomWalkInputQuery;
    randomWalkInputQuery << "SELECT * FROM " << randomWalkInputTableName
                         << " WHERE \"randomWalkCaseId\" = " << caseId
                         << " AND \"monteCarloRunId\" = :monteCarloRunId;";         

    // Compile a SQL query.
    SQLite::Statement query( database, randomWalkInputQuery.str( ).c_str( ) );

    // Declare random walk input table.
    RandomWalkInputTable randomWalkInputTable;

    // Loop through the table retrieved from the database, step-by-step.
    for ( unsigned int i = 0; i < monteCarloRunIdsVector.size( ); i++ )
    {
        // Bind Monte Carlo run ID to query.
        query.bind( ":monteCarloRunId", monteCarloRunIdsVector.at( i ) );

        // Execute select query.
        // A run-time error will be thrown if the requested simulation ID can't be found.
        query.executeStep( );

        // Fetch list of perturbers.
        const std::vector< int > perturbers = getRandomWalkPerturberList( 
            databaseAbsolutePath, monteCarloRunIdsVector.at( i ), 
            randomWalkPerturberTableName );

        // Store fetched row in random walk input struct.
        randomWalkInputTable.insert(
            new RandomWalkInput( query.getColumn( 0 ), query.getColumn( 1 ),
                                 boost::lexical_cast< bool >( query.getColumn( 2 ) ), 
                                 query.getColumn( 3 ), perturbers ) );

        // Reset query.
        query.reset( );
    }

    // Commit database transaction.
    transaction.commit( );

    // Return input table.
    return randomWalkInputTable;
}

//! Get list of selected test particle simulation IDs for random walk Monte Carlo run.
std::vector< int > getRandomWalkPerturberList(
        const std::string& databaseAbsolutePath, const unsigned int monteCarloRunId,
        const std::string& randomWalkPerturberTableName )
{
    // Set stream with perturber table query.
    std::ostringstream randomWalkPerturberQuery;
    randomWalkPerturberQuery << "SELECT * FROM " << randomWalkPerturberTableName
                             << " WHERE \"monteCarloRunId\" = " << monteCarloRunId << ";";

    // Open database in read-only mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READONLY );

    // Open transaction to query database.
    SQLite::Transaction transaction( database );

    // Set up database query.
    SQLite::Statement query( database, randomWalkPerturberQuery.str( ).c_str( ));          
    
    // Declare list of perturbers.
    std::vector< int > perturbers;

    // Loop through the table retrieved from the database, step-by-step.
    while ( query.executeStep( ) )
    {
        perturbers.push_back( query.getColumn( 2 ) );
    }       

    // Commit database transaction.
    transaction.commit( );   

    // Check if list of perturbers is empty.
    if ( perturbers.size( ) == 0 )
    {
        // Throw run-time error.
        throw std::runtime_error( "ERROR: List of perturbers is empty!" );
    }                    

    // Return perturber list.
    return perturbers;    
}

} // namespace database
} // namespace stomi
