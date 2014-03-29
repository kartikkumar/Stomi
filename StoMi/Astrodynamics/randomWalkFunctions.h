/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOMI_RANDOM_WALK_FUNCTIONS_H
#define STOMI_RANDOM_WALK_FUNCTIONS_H

// #include <string>
// #include <vector>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h> 

#include <Assist/Basics/commonTypedefs.h> 
 
#include "StoMi/Database/testParticleKick.h"

namespace stomi
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

//! Compute longitude history.
/*!
 * Computes longitude history derived from time-history of semi-major axis of perturbed body. The
 * semi-major axis of the perturbed body is used to compute its mean motion, which is integrated to
 * yield the longitude. 
 * \param semiMajorAxisHistory Map containing the semi-major axis propagation history for the 
 *          perturbed body.
 * \param centralBodyGravitationalParameter Gravitational parameter of central body that perturbed 
 *          body is in orbit about [kg m^3 s^-2].
 * \return Map containing time-history of longitude of perturbed body.
 */
assist::basics::DoubleKeyDoubleValueMap computeLongitudeHistory(
        const assist::basics::DoubleKeyDoubleValueMap& semiMajorAxisHistory,
        const double centralBodyGravitationalParameter );

//! Reduce longitude history to observation period epoch windows.
/*!
 * Reduces longitude history data of perturbed body to observation period epoch windows only. This 
 * data is used for further analysis, to derive longitude residuals.
 * \param longitudeHistory Time-history of longitude of perturbed body.
 * \param observationPeriodStartEpoch Epoch at start of observation period [s].
 * \param epochWindowSpacing Spacing between windows centered on epochs within observation 
 *          period [s].
 * \param epochWindowSize Size of epoch windows [s].
 * \param numberOfEpochWindows Number of epoch windows in observation period.
 * \return Reduced (to epoch windows) time-history of longitude of perturbed body.
 */
assist::basics::DoubleKeyDoubleValueMap reduceLongitudeHistory( 
        const assist::basics::DoubleKeyDoubleValueMap& longitudeHistory,
        const double observationPeriodStartEpoch,
        const double epochWindowSpacing,
        const double epochWindowSize,
        const unsigned int numberOfEpochWindows );

} // namespace astrodynamics
} // namespace stomi

#endif // STOMI_RANDOM_WALK_FUNCTIONS_H
