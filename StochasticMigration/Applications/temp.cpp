#include <SQLiteC++.h>
#include <vector>
#include <iostream>
#include <sstream>

int main( )
{
    // Create/open database.          
    SQLite::Database database( "/Users/kartikkumar/Desktop/stochasticMigrationResults.sqlite", 
        SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE );

    std::vector< double > case1SemiMajorAxes;

    // Set up database query.
    SQLite::Statement query( database, "SELECT \"semiMajorAxis\" FROM test_particle_input WHERE \"caseId\" == 1;");

    while ( query.executeStep( ) )
    {
        case1SemiMajorAxes.push_back( query.getColumn( 0 ) );
    }

    // Update case 2 semi-major axes.
    for ( unsigned int simulationId = 10001; simulationId < 20001; simulationId++ )
    {
        std::ostringstream update;
        update << "UPDATE test_particle_input SET \"semiMajorAxis\" == " 
               << case1SemiMajorAxes.at( simulationId - 10001 ) 
               << " WHERE \"simulationId\" == " 
               << simulationId 
               << ";";
        database.exec( update.str( ).c_str( ) );
    }

    // Update case 3 semi-major axes.
    for ( unsigned int simulationId = 20001; simulationId < 30001; simulationId++ )
    {
        std::ostringstream update;
        update << "UPDATE test_particle_input SET \"semiMajorAxis\" == " 
               << case1SemiMajorAxes.at( simulationId - 20001 ) 
               << " WHERE \"simulationId\" == " 
               << simulationId 
               << ";";
        database.exec( update.str( ).c_str( ) );
    }

    // Update case 4 semi-major axes.
    for ( unsigned int simulationId = 30001; simulationId < 40001; simulationId++ )
    {
        std::ostringstream update;
        update << "UPDATE test_particle_input SET \"semiMajorAxis\" == " 
               << case1SemiMajorAxes.at( simulationId - 30001 ) 
               << " WHERE \"simulationId\" == " 
               << simulationId 
               << ";";
        database.exec( update.str( ).c_str( ) );
    }    

    return 0;
}