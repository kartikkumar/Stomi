/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130310    K. Kumar          File created.
 *      130328    K. Kumar          Updated unit tests to use boiler plate operator overload 
 *                                  functions in Assist.
 *      130329    K. Kumar          Added test to check that test particle input table is sorted
 *                                  correctly.
 *      130706    K. Kumar          Added test for caseId variable.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <typeinfo>

#include <boost/make_shared.hpp>
#include <boost/ptr_container/ptr_set.hpp>
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
        : simulationId( 1 ),
          caseId( 1 ),
          isCompleted( true ),
          initialStateInKeplerianElements(
              ( Eigen::VectorXd( 6 ) << 1234.5, 0.55, 0.975, 1.234, 4.567, 3.456 ).finished( ) )
    { }

    //! Simulation number.
    int simulationId;

    //! Case ID.
    int caseId;

    //! Flag indicating if simulation has been completed/executed already.
    bool isCompleted;

    //! Initial state of test particle in Keplerian elements.
    tudat::basic_mathematics::Vector6d initialStateInKeplerianElements;

    //! Get test particle input created from specified parameters.
    database::TestParticleInputPointer getTestParticleInput( )
    {
        return boost::make_shared< database::TestParticleInput >(
                    simulationId, caseId, isCompleted, initialStateInKeplerianElements );
    }

protected:

private:
};

BOOST_FIXTURE_TEST_SUITE( test_test_particle_input, TestParticleInputFixture )

//! Test definition of typedef for test particle input database struct.
BOOST_AUTO_TEST_CASE( testTestParticleInputPointerTypedef )
{
    BOOST_CHECK( typeid( database::TestParticleInputTable )
                 == typeid( boost::ptr_set< database::TestParticleInput > ) );
}

//! Test correct construction of test particle input.
BOOST_AUTO_TEST_CASE( testTestParticleInputStructContruction )
{
    // Create test particle input.
    database::TestParticleInputPointer testParticleInput = getTestParticleInput( );

    // Check that the test particle input created contains all the data as required.
    BOOST_CHECK_EQUAL( testParticleInput->simulationId, simulationId );
    BOOST_CHECK_EQUAL( testParticleInput->caseId, caseId );
    BOOST_CHECK_EQUAL( testParticleInput->isCompleted, isCompleted );

    {
        TUDAT_CHECK_MATRIX_BASE( testParticleInput->initialStateInKeplerianElements,
                                 initialStateInKeplerianElements )
                BOOST_CHECK_EQUAL(
                    testParticleInput->initialStateInKeplerianElements.coeff( row, col ),
                    initialStateInKeplerianElements.coeff( row, col ) );
    }
}

//! Test initialization of test particle input with non-positive simulation ID.
BOOST_AUTO_TEST_CASE( testTestParticleInputNonPositiveSimulationIdError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set simulation ID to invalid (non-positive) number.
    simulationId = -1;

    // Try to create test particle input.
    try { database::TestParticleInputPointer testParticleInput = getTestParticleInput( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle input failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle input with non-positive case ID.
BOOST_AUTO_TEST_CASE( testTestParticleInputNonPositiveCaseIdError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set case ID to invalid (non-positive) number.
    caseId = -1;

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

    simulationId = 2;
    database::TestParticleInputPointer testParticleInput2 = getTestParticleInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleInput != *testParticleInput2 );
}

//! Test comparison of TestParticleInput pointers using overloaded < operator.
BOOST_AUTO_TEST_CASE( testTestParticleInputLessThanComparison )
{
    // Create TestParticleInput pointers.
    database::TestParticleInputPointer testParticleInput = getTestParticleInput( );

    simulationId = 2;
    database::TestParticleInputPointer testParticleInput2 = getTestParticleInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleInput < *testParticleInput2 );
}

//! Test comparison of TestParticleInput pointers using overloaded < operator (pointers).
BOOST_AUTO_TEST_CASE( testTestParticleInputLessThanComparisonPointers )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleInput pointers.
    const database::TestParticleInputPointer testParticleInput = getTestParticleInput( );

    simulationId = 2;
    const database::TestParticleInputPointer testParticleInput2 = getTestParticleInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( testParticleInput < testParticleInput2 );
}

//! Test comparison of TestParticleInput pointers using overloaded > operator.
BOOST_AUTO_TEST_CASE( testTestParticleInputGreaterThanComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleInput pointers.
    database::TestParticleInputPointer testParticleInput = getTestParticleInput( );

    simulationId = 2;
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

    simulationId = 2;
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

    simulationId = 2;
    database::TestParticleInputPointer testParticleInput2 = getTestParticleInput( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleInput2 >= *testParticleInput );
}

//! Test sorting of TestParticleInputTable.
BOOST_AUTO_TEST_CASE( testTestParticleInputTableSorting )
{
    // Add TestParticleInput objects to table.
    database::TestParticleInputTable inputTable;

    // Insert input struct.
    inputTable.insert( new database::TestParticleInput( *getTestParticleInput( ) ) );

    // Insert input struct with another simulation number.
    simulationId = 3;
    inputTable.insert( new database::TestParticleInput( *getTestParticleInput( ) ) );

    // Insert input struct with another simulation number.
    simulationId = 2;
    inputTable.insert( new database::TestParticleInput( *getTestParticleInput( ) ) );

    // Check that the table is sorted according to simulation number.
    database::TestParticleInputTable::iterator iteratorInputTable = inputTable.begin( );
    BOOST_CHECK_EQUAL( iteratorInputTable->simulationId, 1 );

    iteratorInputTable++;
    BOOST_CHECK_EQUAL( iteratorInputTable->simulationId, 2 );

    iteratorInputTable++;
    BOOST_CHECK_EQUAL( iteratorInputTable->simulationId, 3 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stochastic_migration
