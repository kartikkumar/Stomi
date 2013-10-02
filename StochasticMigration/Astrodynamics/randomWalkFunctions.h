/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */


#ifndef STOCHASTIC_MIGRATION_RANDOM_WALK_FUNCTIONS_H
#define STOCHASTIC_MIGRATION_RANDOM_WALK_FUNCTIONS_H

// #include <string>
// #include <vector>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h> 

#include <Assist/Basics/commonTypedefs.h> 
 
#include "StochasticMigration/Database/testParticleKick.h"

namespace stochastic_migration
{
namespace astrodynamics
{

//! Execute kick.
/*!
 * Executes kick on perturbed body, specified by row data from an aggregate kick table and returns
 * perturbed body's state after the kick. Only the change in the action variables is computed 
 * (semi-major axis, eccentricity, inclination).
 * \param stateInKeplerianElementsBeforeKick Perturbed body's state in Keplerian elements before 
 *          kick (only action variables).
 * \param kick Row data from aggregate kick table.
 * \param perturberMassRatio Mass ratio between perturber and perturbed body.
 * \return Perturbed body's state in Keplerian elements after kick (only action variables).
 */
Eigen::Vector3d executeKick( const Eigen::Vector3d& stateInKeplerianElementsBeforeKick,
                             const database::TestParticleKickTable::iterator kick, 
                             const double perturberMassRatio );

// //! Compute longitude history.
// DoubleKeyDoubleValueMap computeLongitudeHistory(
//         const ActionPropagationHistory& keplerianActionElementsHistory,
//         const double uranusGravitationalParameter );

// //! Reduce longitude history to observation period epoch windows.
// DoubleKeyDoubleValueMap reduceLongitudeHistory( const DoubleKeyDoubleValueMap& longitudeHistory,
//                                                 const double observationPeriodEpoch,
//                                                 const double epochWindowSpacing,
//                                                 const double epochWindowSize,
//                                                 const unsigned int numberOfEpochWindows );

} // namespace astrodynamics
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_RANDOM_WALK_FUNCTIONS_H
