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
 *      130212    K. Kumar          File created.
 *      130217    K. Kumar          Updated "mab simulations" references to "stochastic migration".
 *      130218    K. Kumar          Updated "encounter" to "conjunction".
 *      130329    K. Kumar          Updated unit test to new kick table definition; updated test 
 *                                  data to generate valid test particle kicks.
 *
 *    References
 *
 *    Notes
 *      Mass factors are not written to the output file, so they are not tested for in the unit
 *      test.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/InputOutput/matrixTextFileReader.h>

#include "StochasticMigration/Basics/basics.h"

#include "StochasticMigration/Database/testParticleKick.h"

#include "StochasticMigration/InputOutput/basicInputOutput.h"

namespace stochastic_migration
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_basic_input_output )

//! Test implementation of write-to-file function for test particle kick table.
BOOST_AUTO_TEST_CASE( testWritingTestParticleKickPointerTableToFile )
{
    using namespace database;

    // Set test particle kick table with test data.
    TestParticleKickPointerTable testParticleKickTable;
    testParticleKickTable.insert(
            TestParticleKickPointer( 
                boost::make_shared< TestParticleKick > (
                    TestParticleKick(
                        1, 1.23456789012345678, 2.34567890123456789, 2.12345678901234567,
                        1.23456789012345678, 2.34567890123456789, 3.45678901234567890,
                        0.12345678901234567, 2.34567890123456789, 3.45678901234567890,
                        3.45678901234567890, 2.34567890123456789, 0.01234567890123456,
                        1.23456789012345678, 0.23456789012345678 ) ) ) );
    testParticleKickTable.insert(
            TestParticleKickPointer( 
                boost::make_shared< TestParticleKick > (
                    TestParticleKick(
                        2, 1.23456789012345678, 2.34567890123456789, 2.12345678901234567,
                        1.23456789012345678, 2.34567890123456789, 3.45678901234567890,
                        0.12345678901234567, 2.34567890123456789, 3.45678901234567890,
                        3.45678901234567890, 2.34567890123456789, 0.01234567890123456,
                        1.23456789012345678, 0.23456789012345678 ) ) ) );

    // Set output file name.
    const std::string outputFilename = "testParticleKickTableTestOutputFile.txt";

    // Set output directory.
    const std::string outputDirectory
           = basics::getStochasticMigrationRootPath( ) + "/InputOutput/UnitTests/";

    // Write test data to output file.
    input_output::writeTestParticleKickPointerTableToFile(
               testParticleKickTable, outputFilename, outputDirectory );

    // Read in test particle kick table output file and store in matrix.
    Eigen::Matrix< double, 2, 13 > testParticleKickTableFromFile
           = tudat::input_output::readMatrixFromFile( outputDirectory + outputFilename );

    // Remove output file.
    boost::filesystem::remove( outputDirectory + outputFilename );

    // Check that the file read in from the output file corresponds with the initial input data.
    unsigned int i = 0;

    for ( TestParticleKickPointerTable::iterator iteratorKickTable 
          = testParticleKickTable.begin( );
          iteratorKickTable != testParticleKickTable.end( ); iteratorKickTable++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->conjunctionEpoch,
                                    testParticleKickTableFromFile( i, 0 ),  1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->conjunctionDistance,
                                    testParticleKickTableFromFile( i, 1 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->conjunctionDuration,
                                    testParticleKickTableFromFile( i, 2 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->preConjunctionEpoch,
                                    testParticleKickTableFromFile( i, 3 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->preConjunctionDistance,
                                    testParticleKickTableFromFile( i, 4 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->preConjunctionSemiMajorAxis,
                                    testParticleKickTableFromFile( i, 5 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->preConjunctionEccentricity,
                                    testParticleKickTableFromFile( i, 6 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->preConjunctionInclination,
                                    testParticleKickTableFromFile( i, 7 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->postConjunctionEpoch,
                                    testParticleKickTableFromFile( i, 8 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->postConjunctionDistance,
                                    testParticleKickTableFromFile( i, 9 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->postConjunctionSemiMajorAxis,
                                    testParticleKickTableFromFile( i, 10 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->postConjunctionEccentricity,
                                    testParticleKickTableFromFile( i, 11 ), 1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( ( *iteratorKickTable )->postConjunctionInclination,
                                    testParticleKickTableFromFile( i, 12 ), 1.0e-14 );

        i++;
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stochastic_migration
