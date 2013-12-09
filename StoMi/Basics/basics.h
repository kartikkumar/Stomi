/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOCHASTIC_MIGRATION_BASICS_H
#define STOCHASTIC_MIGRATION_BASICS_H

#include <string>

#include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h>

#include "StoMi/Database/testParticleCase.h"

namespace stochastic_migration
{
namespace basics
{

//! Get root-path for StochasticMigration directory.
/*!
 * Returns root-path corresponding with root-directory of StochasticMigration as a string with
 * trailing slash included.
 * \return StochasticMigration root-path.
 */
static inline std::string getStochasticMigrationRootPath( )
{
#ifdef STOCHASTIC_MIGRATION_CUSTOM_ROOT_PATH
    return std::string( STOCHASTIC_MIGRATION_CUSTOM_ROOT_PATH );
#else
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path in the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return filePath_.substr( 0, filePath_.length( ) -
                                std::string( "Basics/basics.h" ).length( ) );
#endif
}

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

} // namespace basics
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_BASICS_H

/*
 *    Unit tests are needed for the getStochasticMigrationRootPath() and 
 *    getRungeKuttaCoefficients() functions.
 */
