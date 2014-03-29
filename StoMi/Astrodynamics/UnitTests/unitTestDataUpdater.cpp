/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <boost/make_shared.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>

#include "StoMi/Astrodynamics/body.h"
#include "StoMi/Astrodynamics/dataUpdater.h"

namespace stomi
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_data_updater )

//! Test implementation of data updater for two bodies.
BOOST_AUTO_TEST_CASE( testDataUpdater )
{
    // Declare bodies.
    astrodynamics::BodyPointer body1 = boost::make_shared< astrodynamics::Body >(
                "body1",
                ( Eigen::VectorXd( 6 ) << 1.2, 3.4, -9.7, 4.6, -5.1, -0.1 ).finished( ),
                0.0 );

    astrodynamics::BodyPointer body2  = boost::make_shared< astrodynamics::Body >(
                "body2",
                ( Eigen::VectorXd( 6 ) << -4.5, 4.4, 7.8, -0.2, 6.7, -6.9 ).finished( ),
                0.0 );

    // Declare data updater.
    astrodynamics::DataUpdaterPointer dataUpdater
            = boost::make_shared< astrodynamics::DataUpdater >( body1, body2 );

    // Set updated composite state (body1 state first, followed by body2).
    const Eigen::Matrix< double, 12, 1 > updatedCompositeState
            = ( Eigen::VectorXd( 12 ) << 6.0, -1.3, 5.8, 4.3, 2.0, -0.9,
                1.1, -2.5, -9.8, 1.2, -0.4, 1.3 ).finished( );

    // Set updated time.
    const double updatedTime = 2.4;

    // Update states and time.
    dataUpdater->updateTimeAndCompositeState( updatedTime, updatedCompositeState );

    // Check that states and time were updated correctly.
    {
        TUDAT_CHECK_MATRIX_BASE( body1->getCurrentState( ), updatedCompositeState.segment( 0, 6 ) )
                BOOST_CHECK_EQUAL( body1->getCurrentState( ).coeff( row, col ),
                                   updatedCompositeState.segment( 0, 6 ).coeff( row, col ) );
    }

    {
        TUDAT_CHECK_MATRIX_BASE( body2->getCurrentState( ), updatedCompositeState.segment( 6, 6 ) )
                BOOST_CHECK_EQUAL( body2->getCurrentState( ).coeff( row, col ),
                                   updatedCompositeState.segment( 6, 6 ).coeff( row, col ) );
    }

    BOOST_CHECK_EQUAL( body1->getCurrentTime( ), updatedTime );
    BOOST_CHECK_EQUAL( body2->getCurrentTime( ), updatedTime );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stomi

/*
 *  The unit tests in this file are not strictly "unit tests", but rather system tests, since
 *  they depend on the code implemented in StoMi/Astrodynamics/body.h.
 */
