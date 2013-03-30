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
#include <boost/lexical_cast.hpp>

#include <Eigen/Core>

#include <sqlite3.h>

// #include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>

#include <Assist/Database/sqlite3DatabaseConnector.h>

#include "StochasticMigration/Database/databaseHelpFunctions.h"
#include "StochasticMigration/Database/databaseReadFunctions.h"

namespace stochastic_migration
{
namespace database
{

// //! Using statements.
// using namespace tudat::basic_astrodynamics::orbital_element_conversions;
using namespace assist::database;

// //! Get test particle case data.
// TestParticleCasePointer getTestParticleCase( const std::string& databaseAbsolutePath,
//                                              const std::string& testParticleCaseTableName )
// {
//     // Initiate database connector.
//     Sqlite3DatabaseConnectorPointer databaseConnector
//             = initiateDatabaseConnector( databaseAbsolutePath );

//     // Set stream with query.
//     std::ostringstream testParticleCaseQuery;
//     testParticleCaseQuery << "SELECT * FROM " <<  testParticleCaseTableName << ";";

//     // Prepare database query.
//     databaseConnector->prepare_v2( testParticleCaseQuery.str( ) );

//     // Declare database handler status.
//     unsigned int databaseStatus = 0;

//     // Check if status indicates that there is a row of data to follow.
//     if ( ( databaseStatus = databaseConnector->step( ) ) != SQLITE_ROW )
//     {
//         // Throw run-time error.
//         throwDatabaseError( databaseConnector, databaseStatus );
//     } 

//     // Store fetched row in test particle case struct.
//     const TestParticleCasePointer testParticleCase =
//         boost::make_shared< TestParticleCase >
//             ( TestParticleCase(
//                 databaseConnector->fetchInteger( 0 ), databaseConnector->fetchDouble( 1 ),
//                 databaseConnector->fetchDouble( 2 ),  databaseConnector->fetchDouble( 3 ),
//                 databaseConnector->fetchDouble( 4 ), databaseConnector->fetchDouble( 5 ),
//                 databaseConnector->fetchDouble( 6 ), databaseConnector->fetchDouble( 7 ),
//                 databaseConnector->fetchDouble( 8 ), databaseConnector->fetchDouble( 9 ),
//                 databaseConnector->fetchDouble( 10 ), databaseConnector->fetchDouble( 11 ),
//                 databaseConnector->fetchDouble( 12 ), databaseConnector->fetchDouble( 13 ),
//                 databaseConnector->fetchDouble( 14 ),  databaseConnector->fetchDouble( 15 ),
//                 ( Eigen::VectorXd( 6 ) << databaseConnector->fetchDouble( 16 ),
//                   databaseConnector->fetchDouble( 17 ), databaseConnector->fetchDouble( 18 ),
//                   databaseConnector->fetchDouble( 19 ), databaseConnector->fetchDouble( 20 ),
//                   databaseConnector->fetchDouble( 21 ) ).finished( ),
//                 databaseConnector->fetchString( 22 ), databaseConnector->fetchDouble( 23 ),
//                 databaseConnector->fetchDouble( 24 ), databaseConnector->fetchDouble( 25 ) ) );

//     // Attempt to step through to next row.
//     databaseStatus = databaseConnector->step( );

//     // Throw an error if the end of the query table is not reached.
//     if ( databaseStatus != SQLITE_DONE )
//     {
//         // Throw run-time error.
//         throwDatabaseError( databaseConnector, databaseStatus );
//     }

//     // Terminate database connector cleanly.
//     terminateDatabaseConnector( databaseConnector );

//     // Return test particle case.
//     return testParticleCase;
// }

// //! Get test particle input table.
// TestParticleInputTable getTestParticleInputTable( const std::string& databaseAbsolutePath,
//                                                   bool isCompleted,
//                                                   const std::string& testParticleInputTableName )
// {
//      // Initiate database connector.
//      Sqlite3DatabaseConnectorPointer databaseConnector
//              = initiateDatabaseConnector( databaseAbsolutePath );

//      // Set stream with query.
//      std::ostringstream testParticleInputQuery;
//      testParticleInputQuery << "SELECT * FROM " << testParticleInputTableName
//                             << " WHERE \"completed\" = " << isCompleted << ";";

//     // Prepare database query.
//     databaseConnector->prepare_v2( testParticleInputQuery.str( ) );

//     // Declare database handler status.
//     unsigned int databaseStatus = 0;

//     // Declare test particle input table.
//     TestParticleInputTable testParticleInputTable;

//     // Loop through the table retrieved from the database, step-by-step.
//     while ( ( databaseStatus = databaseConnector->step( ) ) == SQLITE_ROW )
//     {
//         // Store fetched row in test particle input struct.
//         testParticleInputTable.insert(
//                     TestParticleInputPointer(
//                         boost::make_shared< TestParticleInput >(
//                             databaseConnector->fetchInteger( 0 ),
//                             boost::lexical_cast< bool >( databaseConnector->fetchString( 1 ) ),
//                             ( Eigen::VectorXd( 6 ) << databaseConnector->fetchDouble( 2 ),
//                             databaseConnector->fetchDouble( 3 ),
//                             databaseConnector->fetchDouble( 4 ),
//                             databaseConnector->fetchDouble( 5 ),
//                             databaseConnector->fetchDouble( 6 ),
//                             databaseConnector->fetchDouble( 7 ) ).finished( ) ) ) );
//     }

//     // Check if the end of the query table has been reached, and whether all simulations have been
//     // found.
//     if ( databaseStatus != SQLITE_DONE
//          || ( databaseStatus == SQLITE_DONE && testParticleInputTable.size( ) == 0 ) )
//     {
//         // Throw run-time error.
//         throwDatabaseError( databaseConnector, databaseStatus );
//     }

//     // Terminate database connector cleanly.
//     terminateDatabaseConnector( databaseConnector );

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
//         testParticleInputTable.push_back(
//                     TestParticleInput(
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
//                              testParticleInputTable.back( ).simulationNumber );

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
//         const TestParticleSimulationNumbersAndMassFactors&
//         testParticleSimulationNumbersAndMassFactors, const std::string& testParticleKickTableName )
// {
//     // Initiate database connector.
//     Sqlite3DatabaseConnectorPointer databaseConnector
//             = initiateDatabaseConnector( databaseAbsolutePath );

//     // Set up query statement.
//     std::ostringstream testParticleKickQuery;
//     testParticleKickQuery << "SELECT * FROM " << testParticleKickTableName
//                           << " WHERE \"simulation\" IN (";

//     // Loop over test particle simulation numbers and construct SQLite query.
//     TestParticleSimulationNumbersAndMassFactors::const_iterator
//             testParticleSimulationsOneButLastIterator
//             = testParticleSimulationNumbersAndMassFactors.end( );
//     testParticleSimulationsOneButLastIterator--;

//     for ( TestParticleSimulationNumbersAndMassFactors::const_iterator
//           testParticleSimulationsIterator = testParticleSimulationNumbersAndMassFactors.begin( );
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

//     for ( TestParticleSimulationNumbersAndMassFactors::const_iterator
//           testParticleSimulationsIterator = testParticleSimulationNumbersAndMassFactors.begin( );
//           testParticleSimulationsIterator != testParticleSimulationNumbersAndMassFactors.end( );
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
//         // Store fetched row in test particle kick struct.
//         testParticleKickTable.push_back(
//                     TestParticleKick(
//                         databaseConnector->fetchInteger( 1 ), databaseConnector->fetchDouble( 2 ),
//                         databaseConnector->fetchDouble( 3 ), databaseConnector->fetchDouble( 4 ),
//                         databaseConnector->fetchDouble( 5 ), databaseConnector->fetchDouble( 6 ),
//                         databaseConnector->fetchDouble( 7 ), databaseConnector->fetchDouble( 8 ),
//                         databaseConnector->fetchDouble( 9 ), databaseConnector->fetchDouble( 10 ),
//                         databaseConnector->fetchDouble( 11 ), databaseConnector->fetchDouble( 12 ),
//                         databaseConnector->fetchDouble( 13 ), databaseConnector->fetchDouble( 14 ),
//                         testParticleSimulationNumbersAndMassFactors.find(
//                             databaseConnector->fetchInteger( 1 ) )->second ) );

//         // Delete the test particle simulation number if found in the STL vector.
//         std::vector< unsigned int >::iterator iteratorSimulationNumber
//                 = std::find( testParticleSimulationNumbers.begin( ),
//                              testParticleSimulationNumbers.end( ),
//                              testParticleKickTable.back( ).simulationNumber );

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
//         // Store fetched row in random walk Monte Carlo run struct.
//         randomWalkMonteCarloRunTable.push_back(
//                     RandomWalkMonteCarloRun(
//                         databaseConnector->fetchInteger( 0 ), databaseConnector->fetchInteger( 1 ),
//                         databaseConnector->fetchString( 2 ),  databaseConnector->fetchDouble( 3 ),
//                         databaseConnector->fetchDouble( 4 ), databaseConnector->fetchDouble( 5 ),
//                         databaseConnector->fetchDouble( 6 ),
//                         databaseConnector->fetchInteger( 7 ) ) );

//         // Delete the Monte Carlo run number if found in the STL vector.
//         std::vector< unsigned int >::iterator iteratorMonteCarloRunNumber
//                 = std::find( monteCarloRunsCopy.begin( ), monteCarloRunsCopy.end( ),
//                              randomWalkMonteCarloRunTable.back( ).monteCarloRun );

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
//         randomWalkPerturberTable.push_back(
//                     RandomWalkPerturber( databaseConnector->fetchInteger( 1 ),
//                                          databaseConnector->fetchInteger( 2 ),
//                                          databaseConnector->fetchDouble( 3 ) ) );
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
