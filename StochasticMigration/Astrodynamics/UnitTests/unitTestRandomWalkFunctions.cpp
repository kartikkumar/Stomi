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

#include <iterator>

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

    // Execute zero-kick on perturbed body.
    const Eigen::Vector3d perturbedBodyStateAfterKick 
        = executeKick( perturbedBodyStateBeforeKick, kickTable.begin( ), 1.0 );

    // Check that perturbed body state before and after kick are equal.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( perturbedBodyStateBeforeKick, 
                                       perturbedBodyStateAfterKick, 
                                       1.0e-15 );
}

//! Test execution of two kicks that cancel each other out on perturbed body.
BOOST_AUTO_TEST_CASE( testExecuteCancellingKicks )
{
    using namespace stochastic_migration::astrodynamics;
    using namespace stochastic_migration::database;

    // Set mass ratio between perturber and perturbed body.
    const double massRatio = 0.01;

    // Set pre-conjunction state of perturber.
    Eigen::VectorXd preConjunctionState 
      = ( Eigen::VectorXd( 6 ) << 1.0, 0.1, 1.3, 2.3, 3.4, 4.5 ).finished( );

    // Set post-conjunction state of perturber.
    Eigen::VectorXd postConjunctionState 
      = ( Eigen::VectorXd( 6 ) << 2.0, 0.2, 1.5, 1.2, 3.6, 2.9 ).finished( );

    // Set up kick table with two equal and oppositive kicks. 
    TestParticleKickTable kickTable;
    kickTable.insert( new TestParticleKick( 1, 1, 1.0, 0.1, 0.5, 1.9, preConjunctionState,
                                            1.5, 1.8, postConjunctionState ) );
    kickTable.insert( new TestParticleKick( 2, 1, 1.3, 0.1, 0.8, 1.9, postConjunctionState,
                                            1.8, 1.8, preConjunctionState ) );
    // Set perturbed body state before kick.
    const Eigen::Vector3d perturbedBodyInitialState( 1.0, 0.5, 1.2 );

    // Execute kicks on perturbed body.
    TestParticleKickTable::iterator iteratorKickTable = kickTable.begin( );
    const Eigen::Vector3d perturbedBodyStateAfterFirstKick 
        = executeKick( perturbedBodyInitialState, iteratorKickTable, massRatio );

    std::advance( iteratorKickTable, 1 );    
    const Eigen::Vector3d perturbedBodyStateAfterSecondKick 
        = executeKick( perturbedBodyStateAfterFirstKick, iteratorKickTable, massRatio );

    // Check that perturbed body state before and after kick are equal.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( perturbedBodyInitialState, 
                                       perturbedBodyStateAfterSecondKick, 
                                       1.0e-15 );
}

//! Test execution of arbitrary kick on perturbed body.
BOOST_AUTO_TEST_CASE( testExecuteGeneralKick )
{
    using namespace stochastic_migration::astrodynamics;
    using namespace stochastic_migration::database;

    // Set mass ratio between perturber and perturbed body.
    const double massRatio = 0.01;

    // Set pre-conjunction state of perturber.
    Eigen::VectorXd preConjunctionState 
      = ( Eigen::VectorXd( 6 ) << 9.877700000000000e+07, 
                                  2.000000000000000e-04, 
                                  1.919862177193763e-03, 
                                  0.0, 0.0, 0.0 ).finished( );

    // Set post-conjunction state of perturber.
    Eigen::VectorXd postConjunctionState 
      = ( Eigen::VectorXd( 6 ) << 9.912300000000000e+07, 8.800000000000000e-04, 1.570796326794896e-03, 
                                  1.2, 3.6, 2.9 ).finished( );

    // Set up kick table with two equal and oppositive kicks. 
    TestParticleKickTable kickTable;
    kickTable.insert( new TestParticleKick( 1, 1, 1.0, 0.1, 0.5, 1.9, preConjunctionState,
                                            1.5, 1.8, postConjunctionState ) );

    // Set perturbed body state before kick.
    const Eigen::Vector3d perturbedBodyInitialState( 
      9.773600000000000e+07, 2.540000000000000e-03, 2.443460952792062e-03 );

    // Set expected perturbed body state after kick.
    const Eigen::Vector3d expectedPerturbedBodyFinalState( 
      9.773262448570548e+07, 2.662823514764514e-03, 2.445978887439132e-03 );

    // Execute kicks on perturbed body.
    const Eigen::Vector3d perturbedBodyStateAfterKick 
        = executeKick( perturbedBodyInitialState, kickTable.begin( ), massRatio );

    // Check that perturbed body state before and after kick are equal.
    // This test fails for a tolerance lower than 1.0e-10 for the eccentricity and inclination.
    // This is due to the precision of the Python data and the relative error computation.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( perturbedBodyStateAfterKick, 
                                       expectedPerturbedBodyFinalState, 
                                       1.0e-10 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace stochastic_migration    
