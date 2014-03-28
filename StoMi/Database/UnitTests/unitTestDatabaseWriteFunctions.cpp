/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#define BOOST_TEST_MAIN

// #include <iostream>

#include <string>

#include <boost/filesystem.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>  

#include <Eigen/Core>

#include <SQLiteCpp/SQLiteCpp.h>  

#include <TudatCore/InputOutput/matrixTextFileReader.h>

#include "StoMi/Basics/basics.h"

// #include "StoMi/Database/databaseWriteFunctions.h"
// #include "StoMi/Database/testParticleKick.h"

namespace stomi
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_database_read_functions )

//! Test implementation of function to write test particle kick data to SQLite3 database.
BOOST_AUTO_TEST_CASE( testWriteTestParticleKickFunction )
{
    using tudat::input_output::readMatrixFromFile; 

    using namespace basics;
    // using namespace database; 

    // Set absolute path to empty test database.
    const std::string absolutePathToEmptyTestDatabase = getStoMiRootPath( )
            + "/Database/UnitTests/testDatabaseEmptyTestParticleKickTable.db";

    // Copy empty test database to temporary file.
    const std::string absolutePathToTestDatabase = getStoMiRootPath( )
        + "/Database/UnitTests/testDatabaseWritableTestParticleKickTable.db";

    boost::filesystem::copy_file( absolutePathToEmptyTestDatabase, absolutePathToTestDatabase,
                                  boost::filesystem::copy_option::overwrite_if_exists );

    SQLite::Database database( absolutePathToTestDatabase.c_str( ) );

    // Compile a SQL query.
    SQLite::Statement selectInputTableQuery( database, "SELECT * FROM test_particle_input" );

    // Loop to execute the query step by step.
    while ( selectInputTableQuery.executeStep( ) )
    {
        Eigen::VectorXd test = Eigen::VectorXd( 6 );

    //     // Demonstrate how to get some typed column value
    //     int         id          = query.getColumn(0);
    //     int         completed   = query.getColumn(1);
    //     test( 0 )               = query.getColumn(2);
    //     test( 1 )               = query.getColumn(3);
    //     test( 2 )               = query.getColumn(4);
    //     test( 3 )               = query.getColumn(5);
    //     test( 4 )               = query.getColumn(6);
    //     test( 5 )               = query.getColumn(7);
    }

    // // Read in table of test particle kick data from test data file.
    // const Eigen::MatrixXd testDataTestParticleKickTable
    //         = readMatrixFromFile( getStoMiRootPath( )
    //                               + "/Database/UnitTests/testDataWriteKickTableToDatabase.csv" );

    // // Create table of test particle kicks from test data matrix.
    // // Mass ratio used is arbitrary.
    // TestParticleKickTable kickTable;
    // const double massRatio = 0.1;

    // for ( int i = 0; i < testDataTestParticleKickTable.rows( ); i++ ) 
    // {
    //     kickTable.insert(
    //         new TestParticleKick( testDataTestParticleKickTable( i, 1 ),
    //                               testDataTestParticleKickTable( i, 2 ),
    //                               testDataTestParticleKickTable( i, 3 ),
    //                               testDataTestParticleKickTable( i, 4 ),
    //                               testDataTestParticleKickTable( i, 5 ),
    //                               testDataTestParticleKickTable( i, 6 ),
    //                               testDataTestParticleKickTable( i, 7 ),
    //                               testDataTestParticleKickTable( i, 8 ),
    //                               testDataTestParticleKickTable( i, 9 ),
    //                               testDataTestParticleKickTable( i, 10 ),
    //                               testDataTestParticleKickTable( i, 11 ),
    //                               testDataTestParticleKickTable( i, 12 ),
    //                               testDataTestParticleKickTable( i, 13 ),
    //                               testDataTestParticleKickTable( i, 14 ),
    //                               massRatio ) );
    // }        
 
    // // Write kick table to test database.
    // populateTestParticleKickTable( absolutePathToTestDatabase,
    //                                testDataTestParticleKickTable( 0, 1 ), kickTable ); 

    // // Remove test database used for this unit test.
    // boost::filesystem::remove( absolutePathToTestDatabase ); 
 


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stomi
