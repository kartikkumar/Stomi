/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130920   K. Kumar          File created.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>

#include <StochasticMigration/Astrodynamics/randomWalkFunctions.h> 
#include <StochasticMigration/Database/testParticleKick.h> 

namespace stochastic_migration
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( testRandomWalkFunctions )

//! Test execution of zero-kick on perturbed body.
BOOST_AUTO_TEST_CASE( testExecuteZeroKick )
{
    using namespace stochastic_migration::astrodynamics;
    using namespace stochastic_migration::database;

    // Set up kick table with single entry. 
    // Pre- and post-conjunction Keplerian elements for perturber are equal.
    TestParticleKickTable kickTable;
    kickTable.insert( new TestParticleKick( 
                      1, 1, 1.0, 0.1, 0.5, 1.9,
                      ( Eigen::VectorXd( 6 ) << 1.0, 0.1, 1.3, 2.3, 3.4, 4.5 ).finished( ),
                      1.5, 1.8,
                      ( Eigen::VectorXd( 6 ) << 1.0, 0.1, 1.3, 2.3, 3.4, 4.5 ).finished( ) ) );

    // Set perturbed body state before kick.
    const Eigen::Vector3d perturbedBodyStateBeforeKick( 1.0, 0.5, 1.2 );

    // Execute zero-kick on perturbed body and check that 
    const Eigen::Vector3d perturbedBodyStateAfterKick 
        = executeKick( perturbedBodyStateBeforeKick, kickTable.begin( ), 1.0 );

    // Check that perturbed body state before and after kick are equal.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( perturbedBodyStateBeforeKick, 
                                       perturbedBodyStateAfterKick, 
                                       1.0e-15 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stochastic_migration    
