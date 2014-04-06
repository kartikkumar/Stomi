/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <iostream>
#include <stdexcept>
 
#include "StoMi/Mathematics/basicMathematics.h"

namespace stomi
{
namespace mathematics
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
        throw std::runtime_error( 
            "ERROR: Invalid coefficient set specified for numerical integrator!" );
    }

    return rungeKuttaCoefficients;
}

} // namespace mathematics
} // namespace stomi
