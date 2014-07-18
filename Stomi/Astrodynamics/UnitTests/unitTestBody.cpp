/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <string>

#include <boost/make_shared.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <TudatCore/Basics/testMacros.h>

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include "Stomi/Astrodynamics/body.h"

namespace stomi
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_body )

//! Test implementation of  body containing 3D state (position, velocity) and time.
BOOST_AUTO_TEST_CASE( testBodyCreationAndUpdatingWithTime )
{
    // Set name of body.
    std::string bodyName = "testBody";

    // Set initial state of body.
    tudat::basic_mathematics::Vector6d initialState
            = ( Eigen::VectorXd( 6 ) << 1.1, 2.2, 3.3, -4.4, -5.5, -6.6 ).finished( );

    // Set updated state of body.
    tudat::basic_mathematics::Vector6d updatedState
            = ( Eigen::VectorXd( 6 ) << 0.2, -0.4, 0.6, 6.5, 4.3, -12.1 ).finished( );

    // Set initial time.
    const double initialTime = 2.45;

    // Set updated time.
    const double updatedTime = 4.6;

    // Create Body object.
    const astrodynamics::BodyPointer testBody
            = boost::make_shared< astrodynamics::Body >( bodyName, initialState, initialTime );

    // Check that Body object is initialized as expected.
    BOOST_CHECK_EQUAL( testBody->getBodyName( ), bodyName );

    {
        TUDAT_CHECK_MATRIX_BASE( testBody->getCurrentState( ),  initialState )
                BOOST_CHECK_EQUAL( testBody->getCurrentState( ).coeff( row, col ),
                                   initialState.coeff( row, col ) );
    }

    {
        TUDAT_CHECK_MATRIX_BASE( testBody->getCurrentPosition( ), initialState.segment( 0, 3 ) )
                BOOST_CHECK_EQUAL( testBody->getCurrentPosition( ).coeff( row, col ),
                                   initialState.segment( 0, 3 ).coeff( row, col ) );
    }

    {
        TUDAT_CHECK_MATRIX_BASE( testBody->getCurrentVelocity( ), initialState.segment( 3, 3 ) )
                BOOST_CHECK_EQUAL( testBody->getCurrentVelocity( ).coeff( row, col ),
                                   initialState.segment( 3, 3 ).coeff( row, col ) );
    }

    BOOST_CHECK_EQUAL( testBody->getCurrentTime( ), initialTime );

    // Update state and time stored in Body object.
    testBody->setCurrentTimeAndState( updatedState, updatedTime );

    // Check that Body object is updated as expected.
    {
        TUDAT_CHECK_MATRIX_BASE( testBody->getCurrentState( ), updatedState )
                BOOST_CHECK_EQUAL( testBody->getCurrentState( ).coeff( row, col ),
                                   updatedState.coeff( row, col ) );
    }

    {
        TUDAT_CHECK_MATRIX_BASE( testBody->getCurrentPosition( ), updatedState.segment( 0, 3 ) )
                BOOST_CHECK_EQUAL( testBody->getCurrentPosition( ).coeff( row, col ),
                                   updatedState.segment( 0, 3 ).coeff( row, col ) );
    }

    {
        TUDAT_CHECK_MATRIX_BASE( testBody->getCurrentVelocity( ), updatedState.segment( 3, 3 ) )
                BOOST_CHECK_EQUAL( testBody->getCurrentVelocity( ).coeff( row, col ),
                                   updatedState.segment( 3, 3 ).coeff( row, col ) );
    }

    BOOST_CHECK_EQUAL( testBody->getCurrentTime( ), updatedTime );
}

//! Test implementation of  body containing 3D state (position, velocity) and no time.
BOOST_AUTO_TEST_CASE( testBodyCreationAndUpdatingWithoutTime )
{
    // Set name of body.
    std::string bodyName = "testBody";

    // Set initial state of body.
    tudat::basic_mathematics::Vector6d initialState
            = ( Eigen::VectorXd( 6 ) << 1.1, 2.2, 3.3, -4.4, -5.5, -6.6 ).finished( );

    // Set updated state of body.
    tudat::basic_mathematics::Vector6d updatedState
            = ( Eigen::VectorXd( 6 ) << 0.2, -0.4, 0.6, 6.5, 4.3, -12.1 ).finished( );

    // Create Body object.
    const astrodynamics::BodyPointer testBody
            = boost::make_shared< astrodynamics::Body >( bodyName, initialState );

    // Check that Body object is initialized as expected.
    BOOST_CHECK_EQUAL( testBody->getBodyName( ), bodyName );

    {
        TUDAT_CHECK_MATRIX_BASE( testBody->getCurrentState( ),  initialState )
                BOOST_CHECK_EQUAL( testBody->getCurrentState( ).coeff( row, col ),
                                   initialState.coeff( row, col ) );
    }

    {
        TUDAT_CHECK_MATRIX_BASE( testBody->getCurrentPosition( ), initialState.segment( 0, 3 ) )
                BOOST_CHECK_EQUAL( testBody->getCurrentPosition( ).coeff( row, col ),
                                   initialState.segment( 0, 3 ).coeff( row, col ) );
    }

    {
        TUDAT_CHECK_MATRIX_BASE( testBody->getCurrentVelocity( ), initialState.segment( 3, 3 ) )
                BOOST_CHECK_EQUAL( testBody->getCurrentVelocity( ).coeff( row, col ),
                                   initialState.segment( 3, 3 ).coeff( row, col ) );
    }

    BOOST_CHECK_EQUAL( testBody->getCurrentTime( ), 0.0 );

    // Update state and time stored in Body object (time should remain as 0.0).
    testBody->setCurrentTimeAndState( updatedState );

    // Check that Body object is updated as expected.
    {
        TUDAT_CHECK_MATRIX_BASE( testBody->getCurrentState( ), updatedState )
                BOOST_CHECK_EQUAL( testBody->getCurrentState( ).coeff( row, col ),
                                   updatedState.coeff( row, col ) );
    }

    {
        TUDAT_CHECK_MATRIX_BASE( testBody->getCurrentPosition( ), updatedState.segment( 0, 3 ) )
                BOOST_CHECK_EQUAL( testBody->getCurrentPosition( ).coeff( row, col ),
                                   updatedState.segment( 0, 3 ).coeff( row, col ) );
    }

    {
        TUDAT_CHECK_MATRIX_BASE( testBody->getCurrentVelocity( ), updatedState.segment( 3, 3 ) )
                BOOST_CHECK_EQUAL( testBody->getCurrentVelocity( ).coeff( row, col ),
                                   updatedState.segment( 3, 3 ).coeff( row, col ) );
    }

    BOOST_CHECK_EQUAL( testBody->getCurrentTime( ), 0.0 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stomi
