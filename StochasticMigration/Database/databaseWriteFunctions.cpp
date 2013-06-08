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
 *      120402    K. Kumar          File created from old databaseFunctions.cpp.
 *      130217    K. Kumar          updated "mab simulations" references to "stochastic migration".
 *      130218    K. Kumar          Updated "encounter" to "conjunction".
 *      130330    K. Kumar          Updated function to populate kick table.
 *
 *    References
 *
 *    Notes
 *
 */

#include <sstream>

#include <sqlite3.h>

#include <boost/make_shared.hpp>

#include <Assist/Database/sqlite3DatabaseConnector.h>

#include "StochasticMigration/Database/databaseHelpFunctions.h"
#include "StochasticMigration/Database/databaseWriteFunctions.h"

namespace stochastic_migration
{
namespace database
{

using namespace assist::database;

//! Populate test particle kick table.
void populateTestParticleKickTable(
        const std::string& databaseAbsolutePath, const int testParticleSimulationNumber,
        const TestParticleKickTable& testParticleKickTable, const double perturbedBodyEnergyError,
        const double perturbedBodyAngularMomentumError,
        const std::string& testParticleKickTableName,
        const std::string& testParticleInputTableName )
{  
    // // Set completed status variable to 1 (completed), unless table is empty (then set it to -1).
    // int isCompleted = ( testParticleKickTable.size( ) > 0 ) ? 1 : -1;

    // // Set up update statement for test particle input table.
    // std::ostringstream testParticleInputTableUpdate;
    // testParticleInputTableUpdate
    //         << "UPDATE " << testParticleInputTableName << " SET \"completed\" = " << isCompleted
    //         << ", \"perturbedBodyEnergyError\" = "  << perturbedBodyEnergyError
    //         << ", \"perturbedBodyAngularMomentumError\" = " << perturbedBodyAngularMomentumError
    //         << " WHERE \"simulation\" = " << testParticleSimulationNumber << ";" << std::endl;

    // // Initiate database connector.
    // Sqlite3DatabaseConnectorPointer databaseConnector
    //         = initiateDatabaseConnector( databaseAbsolutePath );        

    // Execute completed command.
    // unsigned int databaseStatus = 0; 

    //     databaseConnector->prepare_v2( 
    //     "UPDATE \"test_particle_input\" SET \"completed\" = 1 WHERE \"simulation\" == 9" );
    //     if ( ( databaseStatus = databaseConnector->step( ) ) != SQLITE_DONE )
    //     {
    //         // Throw run-time error.
    //         throwDatabaseError( databaseConnector, databaseStatus );
    //     }
    //         // Terminate database connector cleanly.
    // terminateDatabaseConnector( databaseConnector );

    // Sqlite3DatabaseConnectorPointer databaseConnector  
    //         = boost::make_shared< Sqlite3DatabaseConnector >( databaseAbsolutePath );
    // databaseConnector->execute( "UPDATE test_particle_input SET completed = 1 WHERE simulation == 9;" );    

    // // Terminate database connector cleanly.
    // databaseConnector->closeDatabase( );   

    // sqlite3* database;
    // sqlite3_open( databaseAbsolutePath.c_str( ), &database );

    // // Declare database error message.
    // char* databaseErrorMessage;

    // sqlite3_exec( database, "UPDATE test_particle_input SET completed = 1 WHERE simulation == 9;", NULL, 0, &databaseErrorMessage );

    // // Disconnect from SQLite database.
    // sqlite3_close( database );

    exit( 0 );   





    // // If completed status is -1, return the function handler so the rest is not executed.
    // if ( isCompleted == -1 ) { return; }

    // // Initiate database connector.
    // Sqlite3DatabaseConnectorPointer databaseConnector2
    //         = initiateDatabaseConnector( databaseAbsolutePath );

    // // Set up insert statement for test particle kick table.
    // std::ostringstream testParticleKickInsert;
    // testParticleKickInsert << "INSERT INTO " << testParticleKickTableName
    //                        << " VALUES (NULL, ?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8,"
    //                        << "?9, ?10, ?11, ?12, ?13, ?14);";

    // // Prepare database query.
    // databaseConnector->prepare_v2( testParticleKickInsert.str( ) );

    // // Declare database handler status.
    // unsigned int databaseStatus = 0;

    // // Loop over kick table and add data rows to database table.
    // for ( TestParticleKickTable::iterator iteratorKickTable = testParticleKickTable.begin( );
    //       iteratorKickTable != testParticleKickTable.end( ); iteratorKickTable++ )
    // {
    //     // Bind values to prepared SQLite statement.
    //     databaseConnector->bind( testParticleSimulationNumber, 1 );
    //     databaseConnector->bind( iteratorKickTable->conjunctionEpoch, 2 );
    //     databaseConnector->bind( iteratorKickTable->conjunctionDistance, 3 );
    //     databaseConnector->bind( iteratorKickTable->conjunctionDuration, 4 );
    //     databaseConnector->bind( iteratorKickTable->preConjunctionEpoch, 5 );
    //     databaseConnector->bind( iteratorKickTable->preConjunctionDistance, 6 );
    //     databaseConnector->bind( iteratorKickTable->preConjunctionSemiMajorAxis, 7 );
    //     databaseConnector->bind( iteratorKickTable->preConjunctionEccentricity, 8 );
    //     databaseConnector->bind( iteratorKickTable->preConjunctionInclination, 9 );
    //     databaseConnector->bind( iteratorKickTable->postConjunctionEpoch, 10 );
    //     databaseConnector->bind( iteratorKickTable->postConjunctionDistance, 11 );
    //     databaseConnector->bind(
    //                 iteratorKickTable->postConjunctionSemiMajorAxis, 12 );
    //     databaseConnector->bind( iteratorKickTable->postConjunctionEccentricity, 13 );
    //     databaseConnector->bind( iteratorKickTable->postConjunctionInclination, 14 );

    //     // Insert simulation data row.
    //     if ( ( databaseStatus = databaseConnector->step( ) ) != SQLITE_OK )
    //     {
    //         // Throw run-time error.
    //         throwDatabaseError( databaseConnector, databaseStatus );
    //     }

    //     // Reset values in select statement.
    //     databaseConnector->clearBindings( );

    //     // Reset select statement.
    //     databaseConnector->resetStatement( );
    // }

    // // Terminate database connector cleanly.
    // terminateDatabaseConnector( databaseConnector );
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
