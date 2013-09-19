/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120402    K. Kumar          File created from old randomWalkFunctions.h.
 *      120522    K. Kumar          Added general Keplerian element averaging function; added
 *                                  wrapper for inclination.
 *    References
 *
 *    Notes
 *
 */

#ifndef STOCHASTIC_MIGRATION_RANDOM_WALK_FUNCTIONS_H
#define STOCHASTIC_MIGRATION_RANDOM_WALK_FUNCTIONS_H

// #include <string>
// #include <vector>

// #include <Eigen/Core>

// #include "MabSimulations/AuxilliaryFiles/commonTypedefs.h"
// #include "MabSimulations/AuxilliaryFiles/kickTableRow.h"

namespace stochastic_migration
{
namespace astrodynamics
{

// typedef common_typedefs::ActionPropagationHistory ActionPropagationHistory;
// typedef common_typedefs::DoubleKeyDoubleValueMap DoubleKeyDoubleValueMap;
// typedef tudat::orbital_element_conversions::KeplerianElementVectorIndices
// KeplerianElementVectorIndices;

// //! Execute kick.
// /*!
//  * Executes kick on Mab, specified by row data from an aggregate kick table and returns Mab's state
//  * after the kick. Only the change in the action variables is computed (semi-major axis,
//  * eccentricity, inclination.
//  * \param stateInKeplerianElementsBeforeKick Mab's state in Keplerian elements before kick (only
//  *          action variables).
//  * \param aggregateKickTableRowData Row data from aggregate kick table.
//  * \return Mab's state in Keplerian elements after kick (only action variables).
//  */
// Eigen::Vector3d executeKick( const Eigen::Vector3d& stateInKeplerianElementsBeforeKick,
//                              const kick_table::KickTableRow& aggregateKickTableRowData );

// //! Compare elements from a DoubleKeyDoubleValue map.
// /*!
//  * Compares elements stored in a DoubleKeyDoubleValue map and returns true if the first element
//  * value is smaller than the second.
//  * \param element1Value Value of first element.
//  * \param element2Value Value of second element.
//  */
// bool compareDoubleKeyDoubleValueElements(
//         const DoubleKeyDoubleValueMap::value_type& element1Value,
//         const DoubleKeyDoubleValueMap::value_type& element2Value );

// //! Compute average longitude residual for epoch window.
// /*!
//  * Computes average longitude residual for a given epoch window.
//  * \param propagationHistoryInKeplerianElements Mab's propagation history.
//  * \param epochWindowStart Epoch at start of window.
//  * \param epochWindowEnd Epoch at end of window.
//  */
// double computeAverageLongitudeResidualForEpochWindow(
//         const DoubleKeyDoubleValueMap& longitudeResidualsHistory,
//         const double epochWindowStart, const double epochWindowEnd );

// //! Compute average Keplerian element for epoch window.
// !
//  * Computes average Keplerian element for a given epoch window.
//  * \param keplerianActionElementsHistory Mab's Keplerian action elements history.
//  * \param epochWindowStart Epoch at start of window.
//  * \param epochWindowEnd Epoch at end of window.
//  * \param keplerianElementIndex Index in vector of desired Keplerian element.
 
// double computeAverageKeplerianElementForEpochWindow(
//         const ActionPropagationHistory& keplerianActionElementsHistory,
//         const double epochWindowStart, const double epochWindowEnd,
//         KeplerianElementVectorIndices keplerianElementIndex );

// //! Compute average eccentricity for epoch window.
// /*!
//  * Computes average eccentricity for a given epoch window. This calls the general function with
//  * eccentricity as argument.
//  * \param keplerianActionElementsHistory Mab's Keplerian action elements history.
//  * \param epochWindowStart Epoch at start of window.
//  * \param epochWindowEnd Epoch at end of window.
//  * \sa computeAverageKeplerianElementForEpochWindow().
//  */
// double computeAverageEccentricityForEpochWindow(
//         const ActionPropagationHistory& keplerianActionElementsHistory,
//         const double epochWindowStart, const double epochWindowEnd );

// //! Compute average inclination for epoch window.
// /*!
//  * Computes average inclination for a given epoch window. This calls the general function with
//  * eccentricity as argument.
//  * \param keplerianActionElementsHistory Mab's Keplerian action elements history.
//  * \param epochWindowStart Epoch at start of window.
//  * \param epochWindowEnd Epoch at end of window.
//  * \sa computeAverageKeplerianElementForEpochWindow().
//  */
// double computeAverageInclinationForEpochWindow(
//         const ActionPropagationHistory& keplerianActionElementsHistory,
//         const double epochWindowStart, const double epochWindowEnd );

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
