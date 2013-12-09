/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <iostream>
#include <sstream>

#include <boost/math/special_functions/fpclassify.hpp>

#include <SQLiteC++.h> 

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>

#include "StochasticMigration/Database/databaseWriteFunctions.h"

namespace stochastic_migration
{
namespace database
{

using namespace tudat::basic_astrodynamics::orbital_element_conversions;

//! Populate test particle kick table.
void populateTestParticleKickTable( const std::string& databaseAbsolutePath, 
                                    const int simulationId,
                                    const TestParticleKickTable& kickTable,
                                    const std::string& testParticleKickTableName, 
                                    const std::string& testParticleInputTableName )
{  
    // Set completed status variable to 1 (completed), unless table is empty (then set it to -1).
    const int isCompleted = ( kickTable.size( ) > 0 ) ? 1 : -1;

    // Set up update statement for test particle input table.
    std::ostringstream inputTableUpdate;
    inputTableUpdate << "UPDATE " << testParticleInputTableName << " SET \"completed\" = " 
                     << isCompleted << " WHERE \"simulationId\" = " << simulationId << ";" 
                     << std::endl;

    // Open database in write mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READWRITE );    

    // Execute command to update input table.
    database.exec( inputTableUpdate.str( ).c_str( ) );

    // If completed status is -1, return the function handler so the rest is not executed, since
    // the kick table is empty. Emit warning message.
    if ( isCompleted == -1 ) 
    {   
        std::cerr << "WARNING: Kick table is empty, so updating input table and ";
        std::cerr << "skipping population of kick table in database ..." << std::endl;
        return;
    }

    // Else, populate kick table in the database.
    else if ( isCompleted == 1 )
    {
        // Set up database transaction.
        SQLite::Transaction kickTableTransaction( database );

        // Set up insert statement for test particle kick table.
        std::ostringstream kickTableInsert;
        kickTableInsert << "INSERT INTO " << testParticleKickTableName << " "
                        << "VALUES (NULL, :simulationId, :conjunctionEpoch, :conjunctionDistance, "
                        << ":preConjunctionEpoch, :preConjunctionDistance, "
                        << ":preConjunctionSemiMajorAxis, :preConjunctionEccentricity, "
                        << ":preConjunctionInclination, :preConjunctionArgumentOfPeriapsis, "
                        << ":preConjunctionLongitudeOfAscendingNode, :preConjunctionTrueAnomaly, "
                        << ":postConjunctionEpoch, :postConjunctionDistance, "
                        << ":postConjunctionSemiMajorAxis, :postConjunctionEccentricity, "
                        << ":postConjunctionInclination, :postConjunctionArgumentOfPeriapsis, "
                        << ":postConjunctionLongitudeOfAscendingNode, :postConjunctionTrueAnomaly);";

        // Compile a SQL query.
        SQLite::Statement kickTableInsertQuery( database, kickTableInsert.str( ).c_str( ) );

        // Loop over kick table and add data rows to database table.
        for ( TestParticleKickTable::iterator iteratorKickTable = kickTable.begin( );
              iteratorKickTable != kickTable.end( ); 
              iteratorKickTable++ )
        {
            // Bind values to prepared SQLite statement.
            kickTableInsertQuery.bind( ":simulationId", 
                                       iteratorKickTable->testParticleSimulationId );
            kickTableInsertQuery.bind( ":conjunctionEpoch", iteratorKickTable->conjunctionEpoch );
            kickTableInsertQuery.bind( ":conjunctionDistance", 
                iteratorKickTable->conjunctionDistance );
            kickTableInsertQuery.bind( ":preConjunctionEpoch", 
                iteratorKickTable->preConjunctionEpoch );
            kickTableInsertQuery.bind( ":preConjunctionDistance", 
                iteratorKickTable->preConjunctionDistance );
            kickTableInsertQuery.bind( ":preConjunctionSemiMajorAxis", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( semiMajorAxisIndex ) );
            kickTableInsertQuery.bind( ":preConjunctionEccentricity", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( eccentricityIndex ) );
            kickTableInsertQuery.bind( ":preConjunctionInclination", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( inclinationIndex ) );
            kickTableInsertQuery.bind( ":preConjunctionArgumentOfPeriapsis", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( 
                    argumentOfPeriapsisIndex ) );
            kickTableInsertQuery.bind( ":preConjunctionLongitudeOfAscendingNode", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( 
                    longitudeOfAscendingNodeIndex ) );
            kickTableInsertQuery.bind( ":preConjunctionTrueAnomaly", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( trueAnomalyIndex ) );                                            
            kickTableInsertQuery.bind( ":postConjunctionEpoch", 
                iteratorKickTable->postConjunctionEpoch );
            kickTableInsertQuery.bind( ":postConjunctionDistance", 
                iteratorKickTable->postConjunctionDistance );
            kickTableInsertQuery.bind( ":postConjunctionSemiMajorAxis", 
                iteratorKickTable->postConjunctionStateInKeplerianElements( semiMajorAxisIndex ) );
            kickTableInsertQuery.bind( ":postConjunctionEccentricity", 
                iteratorKickTable->postConjunctionStateInKeplerianElements( eccentricityIndex ) );
            kickTableInsertQuery.bind( ":postConjunctionInclination", 
                iteratorKickTable->postConjunctionStateInKeplerianElements( inclinationIndex ) );
            kickTableInsertQuery.bind( ":postConjunctionArgumentOfPeriapsis", 
                iteratorKickTable->postConjunctionStateInKeplerianElements( 
                    argumentOfPeriapsisIndex ) );
            kickTableInsertQuery.bind( ":postConjunctionLongitudeOfAscendingNode", 
                iteratorKickTable->postConjunctionStateInKeplerianElements( 
                    longitudeOfAscendingNodeIndex ) );
            kickTableInsertQuery.bind( ":postConjunctionTrueAnomaly", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( trueAnomalyIndex ) );     

            // Execute insert query.
            kickTableInsertQuery.exec( );

            // Reset query.
            kickTableInsertQuery.reset( );
        }

        // Commit transaction.
        kickTableTransaction.commit( );
    }
}

//! Populate random walk output tables.
void populateRandomWalkOutputTable(
        const std::string& databaseAbsolutePath, const int monteCarloRunId,
        const double averageLongitudeResidual, const double maximumLongitudeResidualChange,
        const double averageEccentricity, const double maximumEccentricityChange,
        const double averageInclination, const double maximumInclinationChange,
        const std::string& randomWalkOutputTableName, const std::string& randomWalkInputTableName )
{
    using boost::math::isnan;

    // Set completed status variable to 1 (completed), unless any values are NaN.
    const int isCompleted = ( isnan( averageLongitudeResidual ) 
                              || isnan( maximumLongitudeResidualChange )
                              || isnan( averageEccentricity ) 
                              || isnan( maximumEccentricityChange ) 
                              || isnan( averageInclination ) 
                              || isnan( maximumInclinationChange ) ) ? -1 : 1;

    // Set up update statement for random walk input table.
    std::ostringstream inputTableUpdate;
    inputTableUpdate << "UPDATE " << randomWalkInputTableName << " SET \"completed\" = " 
                     << isCompleted << " WHERE \"monteCarloRunId\" = " << monteCarloRunId << ";" 
                     << std::endl;

    // Open database in write mode.          
    SQLite::Database database( databaseAbsolutePath.c_str( ), SQLITE_OPEN_READWRITE );     

   // Set up database transaction.
    SQLite::Transaction outputTableTransaction( database );

    // Execute command to update input table.
    database.exec( inputTableUpdate.str( ).c_str( ) );

    // If completed status is -1, return the function handler so the rest is not executed, since
    // the output is incomplete. Emit warning message.
    if ( isCompleted == -1 ) 
    {   
        std::cerr << "WARNING: Output for Monte Carlo run is incomplete so updating input table ";
        std::cerr << "and skipping population of output table in database ..." << std::endl;
        return;
    }

    // Else, populate output table in the database.
    else if ( isCompleted == 1 )
    {
        // Set up insert statement for test particle kick table.
        std::ostringstream outputTableInsert;
        outputTableInsert << "INSERT INTO " << randomWalkOutputTableName << " "
                          << "VALUES (NULL, " << monteCarloRunId << ", "
                          << averageLongitudeResidual << ", " 
                          << maximumLongitudeResidualChange << ", "
                          << averageEccentricity << ", "
                          << maximumEccentricityChange << ", "
                          << averageInclination << ", "
                          << maximumInclinationChange << ");";

        // Execute insert query.
        database.exec( outputTableInsert.str( ).c_str( ) );
    }

    // Commit transaction.
    outputTableTransaction.commit( );    
}

} // namespace database
} // namespace stochastic_migration
