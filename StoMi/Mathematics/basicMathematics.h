/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOMI_BASICS_H
#define STOMI_BASICS_H

#include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h>

#include "StoMi/Database/testParticleCase.h"

namespace stomi
{
namespace mathematics
{

//! Get Runge-Kutta integrator coefficient set.
/*!
 * Returns Runge-Kutta integrator coefficient set based on string-name provided as input. If the 
 * string-name provided cannot be located, a run-time error is thrown.
 * \param coefficientSet Runge-Kutta integrator coefficient set. Currently, the options available 
 *          are: DOPRI853, RKF78. The enum for this variable is defined in testParticleCase.h
 * \return RungeKuttaCoefficients object with coefficient set loaded.
 */
tudat::numerical_integrators::RungeKuttaCoefficients getRungeKuttaCoefficients( 
    const database::NumericalIntegratorType coefficientSet );

} // namespace mathematics
} // namespace stomi

#endif // STOMI_BASICS_H

/*
 *    Unit tests are needed for the getRungeKuttaCoefficients() function.
 */
