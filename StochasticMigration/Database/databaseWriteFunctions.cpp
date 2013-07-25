/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120402    K. Kumar          File created from old databaseFunctions.cpp.
 *      130217    K. Kumar          updated "mab simulations" references to "stochastic migration".
 *      130218    K. Kumar          Updated "encounter" to "conjunction".
 *      130330    K. Kumar          Updated function to populate kick table.
 *      130725    K. Kumar          Updated function to populate kick table using SQLiteCpp 
 *                                  library. Updated function to reflect database schema.
 *
 *    References
 *
 *    Notes
 *
 */

#include <iostream>
#include <sstream>
 
#include <SQLiteC++.h> 

#include "StochasticMigration/Database/databaseWriteFunctions.h"

namespace stochastic_migration
{
namespace database
{

// using namespace assist::database;

//! Populate test particle kick table.
void populateTestParticleKickTable( const std::string& databaseAbsolutePath, 
                                    const int simulationNumber,
                                    const TestParticleKickTable& kickTable,
                                    const std::string& testParticleKickTableName, 
                                    const std::string& testParticleInputTableName )
{  
    // Set completed status variable to 1 (completed), unless table is empty (then set it to -1).
    const int isCompleted = ( kickTable.size( ) > 0 ) ? 1 : -1;

    // Set up update statement for test particle input table.
    std::ostringstream inputTableUpdate;
    inputTableUpdate << "UPDATE " << testParticleInputTableName << " SET \"completed\" = " 
                     << isCompleted << " WHERE \"simulationId\" = " << simulationNumber << ";" 
                     << std::endl;

    // Open database in write mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READWRITE );    

    // Execute command to update input table.
    database.exec( inputTableUpdate.str( ).c_str( ) );

    // If completed status is -1, return the function handler so the rest is not executed, since
    // the kick table is empty. Emit warning message.
    if ( isCompleted == -1 ) 
    {   
        std::cout << "WARNING: Kick table is empty, so updating input table and ";
        std::cout << "skipping population of kick table in database ..." << std::endl;
        return;
    }

    // Else, populate kick table in the database.
    else if ( isCompleted == 1 )
    {
        // Set up database transaction.
        SQLite::Transaction kicktTableTransaction( database );

        // Set up insert statement for test particle kick table.
        std::ostringstream kickTableInsert;
        kickTableInsert << "INSERT INTO " << testParticleKickTableName << " "
                        << "VALUES (NULL, :simulationId, :conjunctionEpoch, :conjunctionDistance, "
                        << ":preConjunctionEpoch, :preConjunctionDistance, "
                        << ":preConjunctionSemiMajorAxis, :preConjunctionEccentricity, "
                        << ":preConjunctionInclination, :postConjunctionEpoch, "
                        << ":postConjunctionDistance, :postConjunctionSemiMajorAxis, "
                        << ":postConjunctionEccentricity, :postConjunctionInclination);";

        // Compile a SQL query.
        SQLite::Statement kickTableInsertQuery( database, kickTableInsert.str( ).c_str( ) );

        // Loop over kick table and add data rows to database table.
        for ( TestParticleKickTable::iterator iteratorKickTable = kickTable.begin( );
              iteratorKickTable != kickTable.end( ); 
              iteratorKickTable++ )
        {
            // Bind values to prepared SQLite statement.
            kickTableInsertQuery.bind( ":simulationId", iteratorKickTable->simulationNumber );
            kickTableInsertQuery.bind( ":conjunctionEpoch", iteratorKickTable->conjunctionEpoch );
            kickTableInsertQuery.bind( ":conjunctionDistance", 
                iteratorKickTable->conjunctionDistance );
            kickTableInsertQuery.bind( ":preConjunctionEpoch", 
                iteratorKickTable->preConjunctionEpoch );
            kickTableInsertQuery.bind( ":preConjunctionDistance", 
                iteratorKickTable->preConjunctionDistance );
            kickTableInsertQuery.bind( ":preConjunctionSemiMajorAxis", 
                iteratorKickTable->preConjunctionSemiMajorAxis );
            kickTableInsertQuery.bind( ":preConjunctionEccentricity", 
                iteratorKickTable->preConjunctionEccentricity );
            kickTableInsertQuery.bind( ":preConjunctionInclination", 
                iteratorKickTable->preConjunctionInclination );
            kickTableInsertQuery.bind( ":postConjunctionEpoch", 
                iteratorKickTable->postConjunctionEpoch );
            kickTableInsertQuery.bind( ":postConjunctionDistance", 
                iteratorKickTable->postConjunctionDistance );
            kickTableInsertQuery.bind( ":postConjunctionSemiMajorAxis", 
                iteratorKickTable->postConjunctionSemiMajorAxis );
            kickTableInsertQuery.bind( ":postConjunctionEccentricity", 
                iteratorKickTable->postConjunctionEccentricity );
            kickTableInsertQuery.bind( ":postConjunctionInclination", 
                iteratorKickTable->postConjunctionInclination );

            // Execute insert query.
            kickTableInsertQuery.exec( );

            // Reset query.
            kickTableInsertQuery.reset( );
        }

        // Commit transaction.
        kicktTableTransaction.commit( );
    }
}

// //! Populate random walk Monte Carlo run and output tables.
// void populateRandomWalkRunAndOutputTables(
//         const std::string& databaseAbsolutePath,
//         const TestParticleSimulationNumbersAndMassFactors&
//         testParticleSimulationNumbersAndMassFactors,
//         const std::string& massDistributionType,
//         const std::vector< double > massDistributionParameters, const double observationPeriod,
//         const double epochWindowSize, const int numberOfEpochWindows,
//         const double maximumEccentricityChange, const double maximumLongitudeResidualChange,
//         const double maximumInclinationChange )
// {
//     // Initiate database connector.
//     Sqlite3DatabaseConnectorPointer databaseConnector
//             = initiateDatabaseConnector( databaseAbsolutePath );

//     // Query number of random walk runs already present in database.
//     // Set up query statement.
//     std::ostringstream randomWalkRunsQuery;
//     randomWalkRunsQuery << "SELECT COUNT(*) FROM random_walk_runs;";

// //    // Declare sqlite3 statement structure.
// //    sqlite3_stmt* fetchedRandomWalkRunsRowCount;

// //    // Declare pointer to unused part of SQL statement.
// //    const char** unusedPartOfStatement_ = 0;

// //    // Prepare database query.
// //    if ( sqlite3_prepare_v2( database, randomWalkRunsQuery.str( ).c_str( ),
// //                             randomWalkRunsQuery.str( ).size( ),
// //                             &fetchedRandomWalkRunsRowCount, unusedPartOfStatement_ ) != SQLITE_OK )
// //    {
// //        sqlite3_close( database );
// //        boost::throw_exception(
// //                    boost::enable_error_info(
// //                        std::runtime_error( "No data fetched!" ) ) );
// //    }

// //    // Loop through fetched results.
// //    int rowCountRandomWalkRuns = 0;
// //    while ( sqlite3_step( fetchedRandomWalkRunsRowCount ) == SQLITE_ROW )
// //    {
// //        rowCountRandomWalkRuns = boost::lexical_cast< int >(
// //                    sqlite3_column_text( fetchedRandomWalkRunsRowCount, 0 ) );
// //    }
// //    rowCountRandomWalkRuns++;

// //    // Insert random walk run in table.
// //    std::stringstream randomWalkRunInsert;
// //    randomWalkRunInsert << "INSERT INTO random_walk_runs VALUES (\""
// //                        << rowCountRandomWalkRuns << "\",\""
// //                        << selectedSimulationNumbers.size( )
// //                        << "\",\"" << massDistributionType << "\",\""
// //                        << massDistributionParameters.at( 0 );

// //    if ( massDistributionParameters.size( ) == 2 )
// //    {
// //        randomWalkRunInsert << "\",\"" << massDistributionParameters.at( 1 );
// //    }

// //    else
// //    {
// //        randomWalkRunInsert << "\",\"0.0\"";
// //    }

// //    randomWalkRunInsert << ",\"" << observationPeriod << "\","
// //                        << "\"" << epochWindowSize << "\","
// //                        << "\"" << numberOfEpochWindows << "\"";

// //    randomWalkRunInsert << ");" << std::endl;

// //    // Execute insert command.
// //    sqlite3_exec( database, randomWalkRunInsert.str( ).c_str( ), NULL, 0, &databaseErrorMessage );

// //    // Check if sql statement exited with an error message and throw runtime error if necessary.
// //    if ( databaseErrorMessage )
// //    {
// //        sqlite3_close( database );
// //        std::cout << randomWalkRunInsert.str( ) << std::endl;
// //        boost::throw_exception(
// //                    boost::enable_error_info(
// //                        std::runtime_error( "Updating random_walk_runs failed!" ) ) );
// //    }

// //    // Set up insert statement.
// //    std::ostringstream randomWalkSelectionInsert;
// //    randomWalkSelectionInsert << "INSERT INTO random_walk_selection VALUES (NULL, ?1, ?2, ?3);";

// //    // Declare sqlite3 statement structure.
// //    sqlite3_stmt* randomWalkSelectionInsertStatement;

// //    // Prepare sql statement.
// //    sqlite3_prepare_v2( database, randomWalkSelectionInsert.str( ).c_str( ),
// //                        randomWalkSelectionInsert.str( ).size( ),
// //                        &randomWalkSelectionInsertStatement, unusedPartOfStatement_ );

// //    // Check if any part of the sql statement is unused and throw runtime error if necessary.
// //    if ( unusedPartOfStatement_ )
// //    {
// //        sqlite3_close( database );
// //        boost::throw_exception(
// //                    boost::enable_error_info(
// //                        std::runtime_error( "Part of sql statement is unused!" ) ) );
// //    }

// //    // Insert simulation numbers and mass factors into database.
// //    // Loop over selected simulation numbers and mass factors.
// //    for ( common_typedefs::MassFactors::const_iterator iteratorMassFactors = massFactors.begin( );
// //          iteratorMassFactors != massFactors.end( ); iteratorMassFactors++ )
// //    {
// //        sqlite3_bind_int( randomWalkSelectionInsertStatement, 1, rowCountRandomWalkRuns );
// //        sqlite3_bind_int( randomWalkSelectionInsertStatement, 2, iteratorMassFactors->first );
// //        sqlite3_bind_double( randomWalkSelectionInsertStatement, 3, iteratorMassFactors->second );

// //        // Insert simulation data row.
// //        if ( sqlite3_step( randomWalkSelectionInsertStatement ) != SQLITE_DONE )
// //        {
// //            sqlite3_close( database );
// //            boost::throw_exception(
// //                        boost::enable_error_info(
// //                            std::runtime_error( "Random walk selection could not be inserted!" ) ) );
// //        }

// //        // Reset values in insert statement.
// //        sqlite3_clear_bindings( randomWalkSelectionInsertStatement );

// //        // Reset insert statement.
// //        sqlite3_reset( randomWalkSelectionInsertStatement );
// //    }

// //    // Finalize sql statement to prevent memory leaks.
// //    sqlite3_finalize( randomWalkSelectionInsertStatement );

// //    // Insert row data into database.
// //    // Set up insert statement.
// //    std::ostringstream randomWalkRowInsert;
// //    randomWalkRowInsert << "INSERT INTO random_walk VALUES (NULL, "
// //                        << rowCountRandomWalkRuns << ", "
// //                        << maximumEccentricityChange << ", "
// //                        << maximumLongitudeResidualChange << ", "
// //                        << maximumInclinationChange << ");";

// //    // Declare sqlite3 statement structure.
// //    sqlite3_stmt* randomWalkRowInsertStatement;

// //    // Prepare sql statement.
// //    sqlite3_prepare_v2( database, randomWalkRowInsert.str( ).c_str( ),
// //                        randomWalkRowInsert.str( ).size( ),
// //                        &randomWalkRowInsertStatement, unusedPartOfStatement_ );

// //    // Check if any part of the sql statement is unused and throw runtime error if necessary.
// //    if ( unusedPartOfStatement_ )
// //    {
// //        sqlite3_close( database );
// //        boost::throw_exception(
// //                    boost::enable_error_info(
// //                        std::runtime_error( "Part of sql statement is unused!" ) ) );
// //    }

// //    // Insert simulation data row.
// //    if ( sqlite3_step( randomWalkRowInsertStatement ) != SQLITE_DONE )
// //    {
// //        sqlite3_close( database );
// //        boost::throw_exception(
// //                    boost::enable_error_info(
// //                        std::runtime_error( "Random walk row could not be inserted!" ) ) );
// //    }

// //    // Finalize sql statement to prevent memory leaks.
// //    sqlite3_finalize( randomWalkRowInsertStatement );

// //    // End database transaction.
// //    sqlite3_exec( database, "END;", 0, 0, &databaseErrorMessage );

// //    // Disconnect from SQLite database.
// //    sqlite3_close( database );
// }

} // namespace database
} // namespace stochastic_migration
