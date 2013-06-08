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
 *      130217    K. Kumar          File created.
 *
 *    References
 *
 *    Notes 
 *
 */

#define BOOST_TEST_MAIN

#include <iostream>

#include <string>

#include <boost/filesystem.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>  

#include <Eigen/Core>

#include <SQLiteC++.h>  

#include <TudatCore/InputOutput/matrixTextFileReader.h>

#include "StochasticMigration/Basics/basics.h"

#include "StochasticMigration/Database/databaseWriteFunctions.h"
#include "StochasticMigration/Database/testParticleKick.h"

namespace stochastic_migration
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_database_read_functions )

//! Test implementation of function to write test particle kick data to SQLite3 database.
BOOST_AUTO_TEST_CASE( testWriteTestParticleKickFunction )
{
    using tudat::input_output::readMatrixFromFile; 

    using namespace basics;
    using namespace database; 

    // Set absolute path to empty test database.
    const std::string absolutePathToEmptyTestDatabase = getStochasticMigrationRootPath( )
            + "/Database/UnitTests/testDatabaseEmptyTestParticleKickTable.db";

    // Copy empty test database to temporary file.
    const std::string absolutePathToTestDatabase = getStochasticMigrationRootPath( )
        + "/Database/UnitTests/testDatabaseWritableTestParticleKickTable.db";

    boost::filesystem::copy_file( absolutePathToEmptyTestDatabase, absolutePathToTestDatabase,
                                  boost::filesystem::copy_option::overwrite_if_exists );

    SQLite::Database    db( absolutePathToTestDatabase.c_str( ) );

        // Compile a SQL query, containing one parameter (index 1)
    SQLite::Statement   query(db, "SELECT * FROM test_particle_input");

        // Loop to execute the query step by step, to get rows of result
    while (query.executeStep())
    {
        Eigen::VectorXd test = Eigen::VectorXd( 6 );

        // Demonstrate how to get some typed column value
        int         id          = query.getColumn(0);
        int         completed   = query.getColumn(1);
        test( 0 )               = query.getColumn(2);
        test( 1 )               = query.getColumn(3);
        test( 2 )               = query.getColumn(4);
        test( 3 )               = query.getColumn(5);
        test( 4 )               = query.getColumn(6);
        test( 5 )               = query.getColumn(7);
        
        std::cout << id << ", " << completed << ", "
                  << test( 0 ) << "," << test( 1 ) << "," << test( 2 ) << ","
                  << test( 3 ) << "," << test( 4 ) << "," << test( 5 ) << "," << std::endl;
    }

    // // Read in table of test particle kick data from test data file.
    // const Eigen::MatrixXd testDataTestParticleKickTable
    //         = readMatrixFromFile( getStochasticMigrationRootPath( )
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
    // // Energy and angular momenta errors are chosen arbitrarily.
    // const double perturbedBodyEnergyError = 1.0e-5;
    // const double perturbedBodyAngularMomentumError = 1.0e-5;

    // populateTestParticleKickTable( absolutePathToTestDatabase,
    //                                testDataTestParticleKickTable( 0, 1 ), kickTable,  
    //                                perturbedBodyEnergyError, perturbedBodyAngularMomentumError ); 

    // // Remove test database used for this unit test.
    // boost::filesystem::remove( absolutePathToTestDatabase ); 
 


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stochastic_migration
