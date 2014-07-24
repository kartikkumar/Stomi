/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <string>

#include <boost/filesystem.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>  

#include <Eigen/Core>

#include <SQLiteCpp/SQLiteCpp.h>  

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/InputOutput/matrixTextFileReader.h>

#include "Stomi/InputOutput/rootPath.h"

#include "Stomi/Database/databaseWriteFunctions.h"
#include "Stomi/Database/testParticleKick.h"

namespace stomi
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_database_read_functions )

//! Test implementation of function to write test particle kick data to SQLite3 database.
BOOST_AUTO_TEST_CASE( testWriteTestParticleKickFunction )
{
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;   
    using tudat::input_output::readMatrixFromFile; 
    using namespace database; 
    using namespace input_output;

    // Set absolute path to empty test database.
    const std::string absolutePathToEmptyTestDatabase = getStomiRootPath( )
        + "/Database/UnitTests/testDatabaseEmptyTestParticleKickTable.sqlite";

    // Copy empty test database to temporary file.
    const std::string absolutePathToTestDatabase = getStomiRootPath( )
        + "/Database/UnitTests/testDatabaseWritableTestParticleKickTable.sqlite";

    boost::filesystem::copy_file( absolutePathToEmptyTestDatabase, absolutePathToTestDatabase,
                                  boost::filesystem::copy_option::overwrite_if_exists );

    SQLite::Database database( absolutePathToTestDatabase.c_str( ) );

    // Read in table of test particle kick data from test data file.
    const Eigen::Matrix< double, 266, 20 > testDataTestParticleKickTable
            = readMatrixFromFile( getStomiRootPath( )
                                  + "/Database/UnitTests/testDataWriteTestParticleKickTable.csv" );

    // Set test particle simulation ID.
    const unsigned int testParticleSimulationId = 1;

    // Create table of test particle kicks from test data matrix.
    TestParticleKickTable testParticleKickTable;

    for ( int i = 0; i < testDataTestParticleKickTable.rows( ); i++ ) 
    {
        testParticleKickTable.insert(
            new TestParticleKick( 
                testDataTestParticleKickTable( i, 0 ),
                testDataTestParticleKickTable( i, 1 ),
                testDataTestParticleKickTable( i, 2 ),
                testDataTestParticleKickTable( i, 3 ),
                testDataTestParticleKickTable( i, 4 ),
                testDataTestParticleKickTable( i, 5 ),
                ( Eigen::VectorXd( 6 ) 
                    << testDataTestParticleKickTable( i, 6 ),
                       testDataTestParticleKickTable( i, 7 ),
                       testDataTestParticleKickTable( i, 8 ),
                       testDataTestParticleKickTable( i, 9 ),
                       testDataTestParticleKickTable( i, 10 ),
                       testDataTestParticleKickTable( i, 11 ) ).finished( ),
                testDataTestParticleKickTable( i, 12 ),
                testDataTestParticleKickTable( i, 13 ),
                ( Eigen::VectorXd( 6 ) 
                    << testDataTestParticleKickTable( i, 14 ),
                       testDataTestParticleKickTable( i, 15 ),
                       testDataTestParticleKickTable( i, 16 ),
                       testDataTestParticleKickTable( i, 17 ),
                       testDataTestParticleKickTable( i, 18 ),
                       testDataTestParticleKickTable( i, 19 ) ).finished( ) ) );
    }        
 
    // Write kick table to test database.
    populateTestParticleKickTable( absolutePathToTestDatabase,
                                   testParticleSimulationId,
                                   testParticleKickTable,
                                   "test_particle_kicks",
                                   "test_particle_input" );

    // Read data written to kick table.

    // Compile a SQL query.
    SQLite::Statement query( database, "SELECT * FROM \"test_particle_kicks\"" );    

    // Loop through the kick table and compare to the data retrieved from the database.
    for ( TestParticleKickTable::iterator iteratorKickTable = testParticleKickTable.begin( );
          iteratorKickTable != testParticleKickTable.end( ); 
          iteratorKickTable++ )
    {
        // Step through data retrieved.
        query.executeStep( );

        // Check that data retrieved matches kick table.
        BOOST_CHECK_EQUAL( iteratorKickTable->testParticleKickId, 
                           static_cast< int >( query.getColumn( 0 ) ) );
        BOOST_CHECK_EQUAL( iteratorKickTable->testParticleSimulationId, 
                           static_cast< int >( query.getColumn( 1 ) ) );
        BOOST_CHECK_EQUAL( iteratorKickTable->conjunctionEpoch, 
                           static_cast< double >( query.getColumn( 2 ) ) );
        BOOST_CHECK_EQUAL( iteratorKickTable->conjunctionDistance, 
                           static_cast< double >( query.getColumn( 3 ) ) );
        BOOST_CHECK_EQUAL( iteratorKickTable->preConjunctionEpoch, 
                           static_cast< double >( query.getColumn( 4 ) ) );
        BOOST_CHECK_EQUAL( iteratorKickTable->preConjunctionDistance, 
                           static_cast< double >( query.getColumn( 5 ) ) );
        BOOST_CHECK_EQUAL( 
          iteratorKickTable->preConjunctionStateInKeplerianElements( semiMajorAxisIndex ), 
          static_cast< double >( query.getColumn( 6 ) ) );
        BOOST_CHECK_EQUAL( 
          iteratorKickTable->preConjunctionStateInKeplerianElements( eccentricityIndex ), 
          static_cast< double >( query.getColumn( 7 ) ) );
        BOOST_CHECK_EQUAL( 
          iteratorKickTable->preConjunctionStateInKeplerianElements( inclinationIndex ), 
          static_cast< double >( query.getColumn( 8 ) ) );
        BOOST_CHECK_EQUAL( 
          iteratorKickTable->preConjunctionStateInKeplerianElements( argumentOfPeriapsisIndex ), 
          static_cast< double >( query.getColumn( 9 ) ) );
        BOOST_CHECK_EQUAL( 
          iteratorKickTable->preConjunctionStateInKeplerianElements( 
            longitudeOfAscendingNodeIndex ), 
          static_cast< double >( query.getColumn( 10 ) ) );
        BOOST_CHECK_EQUAL( 
          iteratorKickTable->preConjunctionStateInKeplerianElements( trueAnomalyIndex ), 
          static_cast< double >( query.getColumn( 11 ) ) );                    
        BOOST_CHECK_EQUAL( iteratorKickTable->postConjunctionEpoch, 
                           static_cast< double >( query.getColumn( 12 ) ) );
        BOOST_CHECK_EQUAL( iteratorKickTable->postConjunctionDistance, 
                           static_cast< double >( query.getColumn( 13 ) ) );
        BOOST_CHECK_EQUAL( 
          iteratorKickTable->postConjunctionStateInKeplerianElements( semiMajorAxisIndex ), 
          static_cast< double >( query.getColumn( 14 ) ) );
        BOOST_CHECK_EQUAL( 
          iteratorKickTable->postConjunctionStateInKeplerianElements( eccentricityIndex ), 
          static_cast< double >( query.getColumn( 15 ) ) );
        BOOST_CHECK_EQUAL( 
          iteratorKickTable->postConjunctionStateInKeplerianElements( inclinationIndex ), 
          static_cast< double >( query.getColumn( 16 ) ) );
        BOOST_CHECK_EQUAL( 
          iteratorKickTable->postConjunctionStateInKeplerianElements( argumentOfPeriapsisIndex ), 
          static_cast< double >( query.getColumn( 17 ) ) );
        BOOST_CHECK_EQUAL( 
          iteratorKickTable->postConjunctionStateInKeplerianElements( 
            longitudeOfAscendingNodeIndex ), 
          static_cast< double >( query.getColumn( 18 ) ) );
        BOOST_CHECK_EQUAL( 
          iteratorKickTable->postConjunctionStateInKeplerianElements( trueAnomalyIndex ), 
          static_cast< double >( query.getColumn( 19 ) ) );         
    }

    // Remove test database used for this unit test.
    boost::filesystem::remove( absolutePathToTestDatabase ); 
}

//! Test implementation of function to write random walk output data to SQLite3 database.
BOOST_AUTO_TEST_CASE( testWriteRandomWalkOutputFunction )
{
    using tudat::input_output::readMatrixFromFile;   
    using namespace database; 
    using namespace input_output;

    // Set absolute path to empty test database.
    const std::string absolutePathToEmptyTestDatabase = getStomiRootPath( )
        + "/Database/UnitTests/testDatabaseEmptyRandomWalkOutputTable.sqlite";

    // Copy empty test database to temporary file.
    const std::string absolutePathToTestDatabase = getStomiRootPath( )
        + "/Database/UnitTests/testDatabaseWritableRandomWalkOutputTable.sqlite";

    boost::filesystem::copy_file( absolutePathToEmptyTestDatabase, absolutePathToTestDatabase,
                                  boost::filesystem::copy_option::overwrite_if_exists );

    SQLite::Database database( absolutePathToTestDatabase.c_str( ) );

    // Read in table of random walk output data from test data file.
    const Eigen::Matrix< double, 1, 8 > testDataRandomWalkOutputTable
            = readMatrixFromFile( getStomiRootPath( )
                                  + "/Database/UnitTests/testDataWriteRandomWalkOutputTable.csv" );

    // Set random walk simulation ID.
    const unsigned int randomWalkSimulationId = 1;

    // Write output table to test database.
    populateRandomWalkOutputTable( absolutePathToTestDatabase, 
                                   randomWalkSimulationId,
                                   testDataRandomWalkOutputTable( 0, 2 ),
                                   testDataRandomWalkOutputTable( 0, 3 ),
                                   testDataRandomWalkOutputTable( 0, 4 ),
                                   testDataRandomWalkOutputTable( 0, 5 ),
                                   testDataRandomWalkOutputTable( 0, 6 ),
                                   testDataRandomWalkOutputTable( 0, 7 ),
                                   "random_walk_output",
                                   "random_walk_input" );

    // Read data written to output table.

    // Compile a SQL query.
    SQLite::Statement query( database, "SELECT * FROM \"random_walk_output\"" );

    // Step through data retrieved.
    query.executeStep( );   
    
    // Check that data retrieved matches kick table.    
    BOOST_CHECK_EQUAL( testDataRandomWalkOutputTable( 0, 0 ), 
                       static_cast< int >( query.getColumn( 0 ) ) );   
    BOOST_CHECK_EQUAL( randomWalkSimulationId, static_cast< int >( query.getColumn( 1 ) ) ); 
    BOOST_CHECK_EQUAL( testDataRandomWalkOutputTable( 0, 2 ),
                       static_cast< double >( query.getColumn( 2 ) ) );   
    BOOST_CHECK_EQUAL( testDataRandomWalkOutputTable( 0, 3 ),
                       static_cast< double >( query.getColumn( 3 ) ) );  
    BOOST_CHECK_EQUAL( testDataRandomWalkOutputTable( 0, 4 ),
                       static_cast< double >( query.getColumn( 4 ) ) );  
    BOOST_CHECK_EQUAL( testDataRandomWalkOutputTable( 0, 5 ),
                       static_cast< double >( query.getColumn( 5 ) ) );  
    BOOST_CHECK_EQUAL( testDataRandomWalkOutputTable( 0, 6 ),
                       static_cast< double >( query.getColumn( 6 ) ) );  
    BOOST_CHECK_EQUAL( testDataRandomWalkOutputTable( 0, 7 ),
                       static_cast< double >( query.getColumn( 7 ) ) );      

    // Remove test database used for this unit test.
    boost::filesystem::remove( absolutePathToTestDatabase );    
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stomi
