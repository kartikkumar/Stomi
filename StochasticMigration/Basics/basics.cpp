/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130705    K. Kumar          File created.
 *      130717    K. Kumar          Updated getRungeKuttaCoefficients() to use enum for coefficient 
 *                                  set.
 *
 *    References
 *
 *    Notes
 *      A unit test is needed for the getRungeKuttaCoefficients() function.
 *
 */

#include <iostream>
#include <stdexcept>
 
#include "StochasticMigration/Basics/basics.h"

namespace stochastic_migration
{
namespace basics
{

using namespace tudat::numerical_integrators;
using namespace database;

//! Get Runge-Kutta integrator coefficient set.
RungeKuttaCoefficients getRungeKuttaCoefficients( const NumericalIntegratorType coefficientSet )
{
    // Declare Runge-Kutta coefficients.
    RungeKuttaCoefficients rungeKuttaCoefficients;

    // Set output message.
    std::cout << "Runge-Kutta coefficient set                               ";

    if ( coefficientSet == DOPRI853 )
    {
        rungeKuttaCoefficients = RungeKuttaCoefficients::get(
                    RungeKuttaCoefficients::rungeKutta87DormandPrince );
        std::cout << "Dormand Prince 8(7)" << std::endl;
    }

    else if ( coefficientSet == RKF78 )
    {
        rungeKuttaCoefficients = RungeKuttaCoefficients::get(
                    RungeKuttaCoefficients::rungeKuttaFehlberg78 );
        std::cout << "Runge-Kutta-Fehlberg 7(8)" << std::endl;
    }

    else
    {
        throw std::runtime_error( "Invalid coefficient set specified for numerical integrator." );
    }

    return rungeKuttaCoefficients;
}

} // namespace basics
} // namespace stochastic_migration

