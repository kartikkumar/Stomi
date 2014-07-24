/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <iostream>
#include <sstream>

#include <boost/math/special_functions/fpclassify.hpp>

#include <SQLiteCpp/SQLiteCpp.h> 

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>

#include "Stomi/Database/databaseWriteFunctions.h"

namespace stomi
{
namespace database
{

using namespace tudat::basic_astrodynamics::orbital_element_conversions;

//! Populate test particle kick table.
void populateTestParticleKickTable( const std::string& databaseAbsolutePath,
                                    const int testParticleSimulationId,
                                    const TestParticleKickTable& testParticleKickTable,
                                    const std::string& testParticleKickTableName, 
                                    const std::string& testParticleInputTableName )
{  
    // Set completed status variable to 1 (completed), unless table is empty (then set it to -1).
    const int isCompleted = ( testParticleKickTable.size( ) > 0 ) ? 1 : -1;

    // Set up update statement for test particle input table.
    std::ostringstream inputTableUpdate;
    inputTableUpdate << "UPDATE " << testParticleInputTableName << " SET \"completed\" = " 
                     << isCompleted << " WHERE \"testParticleSimulationId\" = " 
                     << testParticleSimulationId << ";" << std::endl;

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
        SQLite::Transaction testParticleKickTableTransaction( database );

        // Set up insert statement for test particle kick table.
        std::ostringstream testParticleKickTableInsert;
        testParticleKickTableInsert << "INSERT INTO " << testParticleKickTableName << " "
                        << "VALUES (NULL, :testParticleSimulationId, :conjunctionEpoch, "
                        << ":conjunctionDistance, "
                        << ":preConjunctionEpoch, :preConjunctionDistance, "
                        << ":preConjunctionSemiMajorAxis, :preConjunctionEccentricity, "
                        << ":preConjunctionInclination, :preConjunctionArgumentOfPeriapsis, "
                        << ":preConjunctionLongitudeOfAscendingNode, :preConjunctionTrueAnomaly, "
                        << ":postConjunctionEpoch, :postConjunctionDistance, "
                        << ":postConjunctionSemiMajorAxis, :postConjunctionEccentricity, "
                        << ":postConjunctionInclination, :postConjunctionArgumentOfPeriapsis, "
                        << ":postConjunctionLongitudeOfAscendingNode, "
                        << ":postConjunctionTrueAnomaly);";

        // Compile a SQL query.
        SQLite::Statement testParticleKickTableInsertQuery( 
            database, testParticleKickTableInsert.str( ).c_str( ) );

        // Loop over kick table and add data rows to database table.
        for ( TestParticleKickTable::iterator iteratorKickTable = testParticleKickTable.begin( );
              iteratorKickTable != testParticleKickTable.end( ); 
              iteratorKickTable++ )
        {
            // Bind values to prepared SQLite statement.
            testParticleKickTableInsertQuery.bind( ":testParticleSimulationId", 
                                       iteratorKickTable->testParticleSimulationId );
            testParticleKickTableInsertQuery.bind( ":conjunctionEpoch", 
                iteratorKickTable->conjunctionEpoch );
            testParticleKickTableInsertQuery.bind( ":conjunctionDistance", 
                iteratorKickTable->conjunctionDistance );
            testParticleKickTableInsertQuery.bind( ":preConjunctionEpoch", 
                iteratorKickTable->preConjunctionEpoch );
            testParticleKickTableInsertQuery.bind( ":preConjunctionDistance", 
                iteratorKickTable->preConjunctionDistance );
            testParticleKickTableInsertQuery.bind( ":preConjunctionSemiMajorAxis", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( semiMajorAxisIndex ) );
            testParticleKickTableInsertQuery.bind( ":preConjunctionEccentricity", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( eccentricityIndex ) );
            testParticleKickTableInsertQuery.bind( ":preConjunctionInclination", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( inclinationIndex ) );
            testParticleKickTableInsertQuery.bind( ":preConjunctionArgumentOfPeriapsis", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( 
                    argumentOfPeriapsisIndex ) );
            testParticleKickTableInsertQuery.bind( ":preConjunctionLongitudeOfAscendingNode", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( 
                    longitudeOfAscendingNodeIndex ) );
            testParticleKickTableInsertQuery.bind( ":preConjunctionTrueAnomaly", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( 
                    trueAnomalyIndex ) );                                            
            testParticleKickTableInsertQuery.bind( ":postConjunctionEpoch", 
                iteratorKickTable->postConjunctionEpoch );
            testParticleKickTableInsertQuery.bind( ":postConjunctionDistance", 
                iteratorKickTable->postConjunctionDistance );
            testParticleKickTableInsertQuery.bind( ":postConjunctionSemiMajorAxis", 
                iteratorKickTable->postConjunctionStateInKeplerianElements( semiMajorAxisIndex ) );
            testParticleKickTableInsertQuery.bind( ":postConjunctionEccentricity", 
                iteratorKickTable->postConjunctionStateInKeplerianElements( eccentricityIndex ) );
            testParticleKickTableInsertQuery.bind( ":postConjunctionInclination", 
                iteratorKickTable->postConjunctionStateInKeplerianElements( inclinationIndex ) );
            testParticleKickTableInsertQuery.bind( ":postConjunctionArgumentOfPeriapsis", 
                iteratorKickTable->postConjunctionStateInKeplerianElements( 
                    argumentOfPeriapsisIndex ) );
            testParticleKickTableInsertQuery.bind( ":postConjunctionLongitudeOfAscendingNode", 
                iteratorKickTable->postConjunctionStateInKeplerianElements( 
                    longitudeOfAscendingNodeIndex ) );
            testParticleKickTableInsertQuery.bind( ":postConjunctionTrueAnomaly", 
                iteratorKickTable->preConjunctionStateInKeplerianElements( 
                    trueAnomalyIndex ) );     

            // Execute insert query.
            testParticleKickTableInsertQuery.exec( );

            // Reset query.
            testParticleKickTableInsertQuery.reset( );
        }

        // Commit transaction.
        testParticleKickTableTransaction.commit( );
    }
}

//! Populate random walk output tables.
void populateRandomWalkOutputTable(
        const std::string& databaseAbsolutePath, const int randomWalkSimulationId,
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
                     << isCompleted << " WHERE \"randomWalkSimulationId\" = " 
                     << randomWalkSimulationId << ";" << std::endl;

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
                          << "VALUES (NULL, " << randomWalkSimulationId << ", "
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
} // namespace stomi
