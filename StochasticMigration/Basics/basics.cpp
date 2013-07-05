/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130705    K. Kumar          File created.
 *
 *    References
 *
 *    Notes
 *      A unit test is needed for the getRungeKuttaCoefficients() function.
 *
 */

#include <iostream>
#include <stdexcept>
 
#include <boost/algorithm/string/predicate.hpp>

#include "StochasticMigration/Basics/basics.h"

namespace stochastic_migration
{
namespace basics
{

using namespace tudat::numerical_integrators;

//! Get Runge-Kutta integrator coefficient set.
RungeKuttaCoefficients getRungeKuttaCoefficients( const std::string& coefficientSet )
{
    // Declare Runge-Kutta coefficients.
    RungeKuttaCoefficients rungeKuttaCoefficients;

    // Set output message.
    std::cout << "Runge-Kutta coefficient set                               ";

    if ( boost::iequals( coefficientSet, "DOPRI853" ) )
    {
        rungeKuttaCoefficients = RungeKuttaCoefficients::get(
                    RungeKuttaCoefficients::rungeKutta87DormandPrince );
        std::cout << "Dormand Prince 8(7)" << std::endl;
    }

    else if ( boost::iequals( coefficientSet, "RKF78" ) )
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

