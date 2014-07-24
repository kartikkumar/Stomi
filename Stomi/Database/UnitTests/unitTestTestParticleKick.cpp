/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <stdexcept>
#include <typeinfo>

#include <boost/make_shared.hpp>
#include <boost/ptr_container/ptr_set.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Basics/testMacros.h>

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include <Assist/Basics/operatorOverloadFunctions.h>

#include "Stomi/Database/testParticleKick.h"

namespace stomi
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics::orbital_element_conversions;
using tudat::basic_mathematics::Vector6d;

//! Test fixture used to test the TestParticleKick struct.
struct TestParticleKickFixture
{
public:

    //! Constructor initializing valid parameters.
    TestParticleKickFixture( )
        : testParticleKickId( 1 ),
          testParticleSimulationId( 1 ),
          conjunctionEpoch( 1.23e4 ),
          conjunctionDistance( 2.0e7 ),
          preConjunctionEpoch( 5.67e2 ),
          preConjunctionDistance( 789.54e8 ),
          preConjunctionStateInKeplerianElements( 
           ( Eigen::VectorXd( 6 ) << 112233.44, 0.456, 1.234, 2.345, 3.456, 4.567 ).finished( ) ),
          postConjunctionEpoch( 9.87e5 ),
          postConjunctionDistance( 987.65e8 ),
          postConjunctionStateInKeplerianElements( 
           ( Eigen::VectorXd( 6 ) << 112244.55, 0.567, 1.345, 5.432, 6.543, 1.654 ).finished( ) )
    { }

    //! Declare parameters of test particle kick.

    //! Unique ID for test particle kick in database.
    int testParticleKickId;

    //! Test particle simulation ID.
    int testParticleSimulationId;

    //! Conjunction epoch [s].
    double conjunctionEpoch;

    //! Conjunction distance [m].
    double conjunctionDistance;

    //! Pre-Conjunction (opposition) epoch [s].
    double preConjunctionEpoch;

    //! Pre-Conjunction (opposition) distance [m].
    double preConjunctionDistance;

    //! Pre-Conjunction state of test particle in Keplerian elements.
    Vector6d preConjunctionStateInKeplerianElements;

    //! Post-Conjunction (opposition) epoch [s].
    double postConjunctionEpoch;

    //! Post-Conjunction (opposition) distance [m].
    double postConjunctionDistance;

    //! Post-Conjunction state of test particle in Keplerian elements.
    Vector6d postConjunctionStateInKeplerianElements;

    //! Get test particle kick created from specified parameters.
    database::TestParticleKickPointer getTestParticleKick( )
    {
        return boost::make_shared< database::TestParticleKick >(
                    database::TestParticleKick(
                        testParticleKickId, testParticleSimulationId, 
                        conjunctionEpoch, conjunctionDistance,
                        preConjunctionEpoch, preConjunctionDistance, 
                        preConjunctionStateInKeplerianElements,
                        postConjunctionEpoch, postConjunctionDistance,
                        postConjunctionStateInKeplerianElements ) );
    }

protected:

private:
};

BOOST_FIXTURE_TEST_SUITE( test_test_particle_kick, TestParticleKickFixture )

//! Test definition of typedef for test particle kick database struct.
BOOST_AUTO_TEST_CASE( testTestParticleKickTypedef )
{
    BOOST_CHECK( typeid( database::TestParticleKickTable )
                 == typeid( boost::ptr_set< database::TestParticleKick > ) );
}

//! Test correct construction of test particle kick.
BOOST_AUTO_TEST_CASE( testTestParticleKickStructContruction )
{
    // Create test particle kick.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );

    // Check that the test particle kick created contains all the data as required.
    BOOST_CHECK_EQUAL( testParticleKick->testParticleKickId, testParticleKickId );    
    BOOST_CHECK_EQUAL( testParticleKick->testParticleSimulationId, testParticleSimulationId );
    BOOST_CHECK_EQUAL( testParticleKick->conjunctionEpoch, conjunctionEpoch );
    BOOST_CHECK_EQUAL( testParticleKick->conjunctionDistance, conjunctionDistance );
    BOOST_CHECK_EQUAL( testParticleKick->preConjunctionEpoch, preConjunctionEpoch );
    BOOST_CHECK_EQUAL( testParticleKick->preConjunctionDistance, preConjunctionDistance );
    {
        TUDAT_CHECK_MATRIX_BASE( testParticleKick->preConjunctionStateInKeplerianElements,
                                 preConjunctionStateInKeplerianElements )
                BOOST_CHECK_EQUAL(
                    testParticleKick->preConjunctionStateInKeplerianElements.coeff( row, col ),
                    preConjunctionStateInKeplerianElements.coeff( row, col ) );
    }
    BOOST_CHECK_EQUAL( testParticleKick->postConjunctionEpoch, postConjunctionEpoch );
    BOOST_CHECK_EQUAL( testParticleKick->postConjunctionDistance, postConjunctionDistance );
    {
        TUDAT_CHECK_MATRIX_BASE( testParticleKick->postConjunctionStateInKeplerianElements,
                                 postConjunctionStateInKeplerianElements )
                BOOST_CHECK_EQUAL(
                    testParticleKick->postConjunctionStateInKeplerianElements.coeff( row, col ),
                    postConjunctionStateInKeplerianElements.coeff( row, col ) );
    }                       
}

//! Test initialization of test particle kick with non-positive test particle kick ID.
BOOST_AUTO_TEST_CASE( testTestParticleKickNonPositiveTestParticleKickIdError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set test particle kick ID to invalid (non-positive) number.
    testParticleKickId = -1;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with non-positive simulation ID.
BOOST_AUTO_TEST_CASE( testTestParticleKickNonPositiveTestParticleSimulationIdError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set simulation ID to invalid (non-positive) number.
    testParticleSimulationId = -1;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative conjunction epoch.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativeConjunctionEpochError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set conjunction epoch to invalid (negative) number.
    conjunctionEpoch = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative conjunction distance.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativeConjunctionDistanceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set conjunction distance to invalid (negative) number.
    conjunctionDistance = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with pre-conjunction epoch after conjunction epoch.
BOOST_AUTO_TEST_CASE( testTestParticleKickPreConjunctionAfterConjunctionEpochError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction epoch to value after conjunction epoch.
    preConjunctionEpoch = 2.34e5;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative pre-conjunction distance.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePreConjunctionDistanceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction distance to invalid (negative) number.
    preConjunctionDistance = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with pre-conjunction distance less than conjunction
//! distance.
BOOST_AUTO_TEST_CASE( testTestParticleKickPreConjunctionLessThanConjunctionDistanceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction distance to value less than conjunction distance.
    preConjunctionDistance = 1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative pre-conjunction semi-major axis.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePreConjunctionSemiMajorAxisError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction semi-major axis to invalid (negative) number.
    preConjunctionStateInKeplerianElements( semiMajorAxisIndex ) = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative pre-conjunction eccentricity.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePreConjunctionEccentricityError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction eccentricity to invalid (negative) number.
    preConjunctionStateInKeplerianElements( eccentricityIndex ) = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with >1.0 pre-conjunction eccentricity.
BOOST_AUTO_TEST_CASE( testTestParticleKickGreaterThanOnePreConjunctionEccentricityError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction eccentricity to invalid (greater than 1.0) number.
    preConjunctionStateInKeplerianElements( eccentricityIndex ) = 1.1;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative pre-conjunction inclination.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePreConjunctionInclinationError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set pre-conjunction eccentricity to invalid (negative) number.
    preConjunctionStateInKeplerianElements( inclinationIndex ) = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with conjunction epoch after post-conjunction epoch.
BOOST_AUTO_TEST_CASE( testTestParticleKickConjunctionAfterPostConjunctionEpochError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction epoch to value before conjunction epoch.
    postConjunctionEpoch = 0.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with post-conjunction distance less than conjunction
//! distance.
BOOST_AUTO_TEST_CASE( testTestParticleKickPostConjunctionLessThanConjunctionDistanceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction distance to value less than conjunction distance.
    postConjunctionDistance = 1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative post-conjunction distance.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePostConjunctionDistanceError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction distance to invalid (negative) number.
    postConjunctionDistance = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative post-conjunction semi-major axis.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePostConjunctionSemiMajorAxisError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction semi-major axis to invalid (negative) number.
    postConjunctionStateInKeplerianElements( semiMajorAxisIndex ) = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative post-conjunction eccentricity.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePostConjunctionEccentricityError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction eccentricity to invalid (negative) number.
    postConjunctionStateInKeplerianElements( eccentricityIndex ) = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with >1.0 post-conjunction eccentricity.
BOOST_AUTO_TEST_CASE( testTestParticleKickGreaterThanOnePostConjunctionEccentricityError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction eccentricity to invalid (greater than 1.0) number.
    postConjunctionStateInKeplerianElements( eccentricityIndex ) = 1.1;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test initialization of test particle kick with negative post-conjunction inclination.
BOOST_AUTO_TEST_CASE( testTestParticleKickNegativePostConjunctionInclinationError )
{
    // Set flag to indicate if error is thrown to false.
    bool isError = false;

    // Set post-conjunction eccentricity to invalid (negative) number.
    postConjunctionStateInKeplerianElements( inclinationIndex ) = -1.0;

    // Try to create test particle kick.
    try { database::TestParticleKickPointer testParticleKick = getTestParticleKick( ); }

    // Catch expected run-time error.
    catch ( std::runtime_error& error ) { isError = true; }

    // Check that construction of test particle kick failed.
    BOOST_CHECK( isError );
}

//! Test comparison of TestParticleKick pointers using overloaded == operator.
BOOST_AUTO_TEST_CASE( testTestParticleKickEqualComparison )
{
    // Create TestParticleKick pointers.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );
    database::TestParticleKickPointer testParticleKick2 = getTestParticleKick( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleKick == *testParticleKick2 );
}

//! Test comparison of TestParticleKick pointers using overloaded != operator.
BOOST_AUTO_TEST_CASE( testTestParticleKickNonEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleKick pointers.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );

    testParticleSimulationId = 2;
    database::TestParticleKickPointer testParticleKick2 = getTestParticleKick( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleKick != *testParticleKick2 );
}

//! Test comparison of TestParticleKick pointers using overloaded < operator.
BOOST_AUTO_TEST_CASE( testTestParticleKickLessThanComparison )
{
    // Create TestParticleKick pointers.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );

    conjunctionEpoch = 1.24e4;
    database::TestParticleKickPointer testParticleKick2 = getTestParticleKick( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleKick < *testParticleKick2 );
}

//! Test comparison of TestParticleKick pointers using overloaded > operator.
BOOST_AUTO_TEST_CASE( testTestParticleKickGreaterThanComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleKick pointers.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );

    conjunctionEpoch = 1.24e4;
    database::TestParticleKickPointer testParticleKick2 = getTestParticleKick( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleKick2 > *testParticleKick );
}

//! Test comparison of TestParticleKick pointers using overloaded <= operator.
BOOST_AUTO_TEST_CASE( testTestParticleKickLessThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleKick pointers.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );

    conjunctionEpoch = 1.24e4;
    database::TestParticleKickPointer testParticleKick2 = getTestParticleKick( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleKick <= *testParticleKick2 );
}

//! Test comparison of TestParticleKick pointers using overloaded >= operator.
BOOST_AUTO_TEST_CASE( testTestParticleKickGreaterThanOrEqualComparison )
{
    using namespace assist::basics::operator_overload_functions;

    // Create TestParticleKick pointers.
    database::TestParticleKickPointer testParticleKick = getTestParticleKick( );

    conjunctionEpoch = 1.24e4;
    database::TestParticleKickPointer testParticleKick2 = getTestParticleKick( );

    // Check that operator is overloaded correctly.
    BOOST_CHECK( *testParticleKick2 >= *testParticleKick );
}

//! Test sorting of TestParticleKickTable.
BOOST_AUTO_TEST_CASE( testTestParticleKickTableSorting )
{
    // Add TestParticleKick objects to table.
    database::TestParticleKickTable kickTable;

    // Insert kick with conjunction epoch = 1.23e4. Set test particle kick ID.
    testParticleKickId = 2;
    kickTable.insert( new database::TestParticleKick( *getTestParticleKick( ) ) );

    // Insert kick with later conjunction epoch. Set test particle kick ID.
    testParticleKickId = 3;
    conjunctionEpoch = 1.24e4;
    kickTable.insert( new database::TestParticleKick( *getTestParticleKick( ) ) );

    // Insert kick with earlier conjunction epoch. Set test particle kick ID.
    testParticleKickId = 1;
    conjunctionEpoch = 1.2e3;
    kickTable.insert( new database::TestParticleKick( *getTestParticleKick( ) ) );

    // Check that the table is sorted according to conjunction epoch.
    database::TestParticleKickTable::iterator iteratorKickTable = kickTable.begin( );
    BOOST_CHECK_EQUAL( iteratorKickTable->testParticleKickId, 1 );
    BOOST_CHECK_EQUAL( iteratorKickTable->conjunctionEpoch, 1.2e3 );

    iteratorKickTable++;
    BOOST_CHECK_EQUAL( iteratorKickTable->testParticleKickId, 2 );
    BOOST_CHECK_EQUAL( iteratorKickTable->conjunctionEpoch, 1.23e4 );

    iteratorKickTable++;
    BOOST_CHECK_EQUAL( iteratorKickTable->testParticleKickId, 3 );    
    BOOST_CHECK_EQUAL( iteratorKickTable->conjunctionEpoch, 1.24e4 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stomi
