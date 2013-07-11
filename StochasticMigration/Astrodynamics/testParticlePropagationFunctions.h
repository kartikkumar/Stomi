/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120817    K. Kumar          File created.
 *      130218    K. Kumar          Updated "perturbedBody simulations" references to "stochastic migration";
 *                                  renamed file.
 *
 *    References
 *
 *    Notes
 *      The functions implemented in this file need to be unit tested.
 *
 */

#ifndef STOCHASTIC_MIGRATION_TEST_PARTICLE_PROPAGATION_FUNCTIONS_H
#define STOCHASTIC_MIGRATION_TEST_PARTICLE_PROPAGATION_FUNCTIONS_H

#include <map>
#include <utility>

#include <Assist/Astrodynamics/body.h>

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h>

#include "StochasticMigration/Database/testParticleCase.h"
#include "StochasticMigration/Database/testParticleKick.h"

namespace stochastic_migration
{
namespace astrodynamics
{

inline bool checkIfMaximum( const double currentValue, const double previousValue, const double nextValue )
{
    if ( currentValue > previousValue && currentValue > nextValue )        
    {
        return true;
    }

    return false;
}

inline bool checkIfMinimum( const double currentValue, const double previousValue, const double nextValue )
{
    if ( currentValue < previousValue && currentValue < nextValue )        
    {
        return true;
    }

    return false;
}

//! Propagate test-particle-perturbed-body system and generate test particle kick table.
/*!
* Propagates system of test particle, perturbed body, orbiting the specified central body forward
* in time using the numerical integrator provided.
* \param perturbedBody Shared-pointer to perturbed body.
* \param testParticle Shared-pointer to perturbed test particle.
* \param testParticleCase TestParticleCase object containing all of the simulation case data.
* \param numericalIntegrator Shared-pointer to numerical integrator used to propagate system.
* \return Table of test particle kicks, stored as a set of pointers.
*/
database::TestParticleKickTable propagateSystemAndGenerateKickTable(
        const assist::astrodynamics::BodyPointer perturbedBody,
        const assist::astrodynamics::BodyPointer testParticle,
        const database::TestParticleCasePointer testParticleCase,
        const tudat::numerical_integrators::
                RungeKuttaVariableStepSizeIntegratorXdPointer numericalIntegrator );

} // namespace astrodynamics
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_TEST_PARTICLE_PROPAGATION_FUNCTIONS_H
