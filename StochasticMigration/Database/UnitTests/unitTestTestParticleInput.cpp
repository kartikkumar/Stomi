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
 *      130310    K. Kumar          File created.
 *      130328    K. Kumar          Updated unit tests to use boiler plate operator overload 
 *                                  functions in Assist.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <set>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Basics/testMacros.h>

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include <Assist/Basics/operatorOverloadFunctions.h>

#include "StochasticMigration/Database/testParticleInput.h"

namespace stochastic_migration
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics::orbital_element_conversions;

//! Test fixture used to test the TestParticleInput struct.
struct TestParticleInputFixture
{
public:

    // Set Eigen macro to correctly align class with fixed-size vectorizable types.
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor initializing valid parameters.
    TestParticleInputFixture( )
        : simulationNumber( 1 ),
          isCompleted( true ),
          initialStateInKeplerianElements(
              ( Eigen::VectorXd( 6 ) << 1234.5, 0.55, 0.975, 1.234, 4.567, 3.456 ).finished( ) )
    { }

    //! Declare parameters of test particle input.
    int simulationNumber;
    bool isCompleted;
    tudat::basic_mathematics::Vector6d initialStateInKeplerianElements;

    //! Get test particle input created from specified parameters.
    database::TestParticleInputPointer getTestParticleInput( )
    {
        return boost::make_shared< database::TestParticleInput >(
                    simulationNumber, isCompleted, initialStateInKeplerianElements );
    }

protected:

private:
};

BOOST_FIXTURE_TEST_SUITE( test_test_particle_input, TestParticleInputFixture )

//! Test definition of typedef for test particle input database struct.
BOOST_AUTO_TEST_CASE( testTestParticleInputPointerTypedef )
{
    BOOST_CHECK( typeid( database::TestParticleInputTable )
                 == typeid( std::set< database::TestParticleInputPointer > ) );
}

//! Test correct construction of test particle input.
BOOST_AUTO_TEST_CASE( testTestParticleInputStructContruction )
{
    // Create test particle input.
    database::TestParticleInputPointer testParticleInput = getTestParticleInput( );

    // Check that the test particle input created contains all the data as required.
    BOOST_CHECK_EQUAL( testParticleInput->simulationNumber, simulationNumber );
    BOOST_CHECK_EQUAL( testParticleInput->isCompleted, isCompleted );

    {
        TUDAT_CHECK_MATRIX_BASE( testParticleInput->initialStateInKeplerianElements,
                                 initialStateInKeplerianElements )
                BOOST_CHECK_EQUAL(
                    testParticleInput->initialStateInKeplerianElements.coeff( row, col ),
                    initialStateInKeplerianElements.coeff( row, col ) );
    }
}

//! Test initialization of test particle input with non-positive simulation number.
BOOST_AUTO_TEST_CASE( testTestParticleInputNonPositiveSimulationNumberError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set simulation number to invalid (non-positive) number.
    simulationNumber = -1;

    // Try to create test particle input.
    try { database::TestParticleInputPointer testParticleInput = getTestParticleInput( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle input failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle input with negative initial semi-major axis.
BOOST_AUTO_TEST_CASE( testTestParticleInputNegativeInitialSemiMajorAxisError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set initial semi-major axis to invalid (negative) number.
    initialStateInKeplerianElements( semiMajorAxisIndex ) = -1.0;

    // Try to create test particle input.
    try { database::TestParticleInputPointer testParticleInput = getTestParticleInput( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle input failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle input with negative initial eccentricity.
BOOST_AUTO_TEST_CASE( testTestParticleInputNegativeInitialEccentricityError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set initial eccentricity to invalid (negative) number.
    initialStateInKeplerianElements( eccentricityIndex ) = -1.0;

    // Try to create test particle input.
    try { database::TestParticleInputPointer testParticleInput = getTestParticleInput( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle input failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle input with >1.0 initial eccentricity.
BOOST_AUTO_TEST_CASE( testTestParticleInputGreaterThanOneInitialEccentricityError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set initial eccentricity to invalid (greater than 1.0) number.
    initialStateInKeplerianElements( eccentricityIndex ) = 1.1;

    // Try to create test particle input.
    try { database::TestParticleInputPointer testParticleInput = getTestParticleInput( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle input failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative initial inclination.
BOOST_AUTO_TEST_CASE( testTestParticleInputNegativeInitialInclinationError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set initial eccentricity to invalid (negative) number.
    initialStateInKeplerianElements( inclinationIndex ) = -1.0;

    // Try to create test particle input.
    try { database::TestParticleInputPointer testParticleInput = getTestParticleInput( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle input failed.
    BOOST_CHECK( isError );
}

//! Test comparison of TestParticleInput pointers using overloaded == operator.
BOOST_AUTO_TEST_CASE( testTestParticleInputEqualComparison )
{
    // Create TestParticleInput pointers.
    database::TestParticleInputPointer testParticleInput = getTestParticleInput( );
    database::TestParticleInputPointer testParticleInput2 = getTestParticleInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleInput == *testParticleInput2 );
}

//! Test comparison of TestParticleInput pointers using overloaded != operator.
BOOST_AUTO_TEST_CASE( testTestParticleInputNonEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleInput pointers.
    database::TestParticleInputPointer testParticleInput = getTestParticleInput( );

    simulationNumber = 2;
    database::TestParticleInputPointer testParticleInput2 = getTestParticleInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleInput != *testParticleInput2 );
}

//! Test comparison of TestParticleInput pointers using overloaded < operator.
BOOST_AUTO_TEST_CASE( testTestParticleInputLessThanComparison )
{
    // Create TestParticleInput pointers.
    database::TestParticleInputPointer testParticleInput = getTestParticleInput( );

    simulationNumber = 2;
    database::TestParticleInputPointer testParticleInput2 = getTestParticleInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleInput < *testParticleInput2 );
}

//! Test comparison of TestParticleInput pointers using overloaded > operator.
BOOST_AUTO_TEST_CASE( testTestParticleInputGreaterThanComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleInput pointers.
    database::TestParticleInputPointer testParticleInput = getTestParticleInput( );

    simulationNumber = 2;
    database::TestParticleInputPointer testParticleInput2 = getTestParticleInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleInput2 > *testParticleInput );
}

//! Test comparison of TestParticleInput pointers using overloaded <= operator.
BOOST_AUTO_TEST_CASE( testTestParticleInputLessThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleInput pointers.
    database::TestParticleInputPointer testParticleInput = getTestParticleInput( );

    simulationNumber = 2;
    database::TestParticleInputPointer testParticleInput2 = getTestParticleInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleInput <= *testParticleInput2 );
}

//! Test comparison of TestParticleInput pointers using overloaded >= operator.
BOOST_AUTO_TEST_CASE( testTestParticleInputGreaterThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleInput pointers.
    database::TestParticleInputPointer testParticleInput = getTestParticleInput( );

    simulationNumber = 2;
    database::TestParticleInputPointer testParticleInput2 = getTestParticleInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleInput2 >= *testParticleInput );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stochastic_migration
