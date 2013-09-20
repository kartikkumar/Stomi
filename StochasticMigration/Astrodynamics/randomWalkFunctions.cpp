/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120402    K. Kumar          File created from old randomWalkFunctions.cpp.
 *      120522    K. Kumar          Added inclination kick computation. Added general Keplerian
 *                                  element averaging function; added wrapper for inclination.
 *      130206    K. Kumar          Fixed error in inclination kick computation.
 *      130919    K. Kumar          Uncommented and updated executeKick() function.
 *
 *    References
 *
 *    Notes
 *
 */

#include <cmath>
#include <stdexcept>
// #include <iomanip>
// #include <iostream>
// #include <iterator>
// #include <limits>
// #include <numeric>
// #include <sstream>

#include <boost/math/special_functions/fpclassify.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
// #include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "StochasticMigration/Astrodynamics/randomWalkFunctions.h"

namespace stochastic_migration
{
namespace astrodynamics
{

using namespace tudat::astrodynamics;
using namespace tudat::orbital_element_conversions;

//! Execute kick.
Eigen::Vector3d executeKick( const Eigen::Vector3d& stateInKeplerianElementsBeforeKick,
                             const database::TestParticleKickTable::iterator kick, 
                             const double perturberMassRatio )
{
    // Declare state in Keplerian elements after kick.
    Eigen::Vector3d stateInKeplerianElementsAfterKick = Eigen::Vector3d::Zero( );

    // Declare and set pre- and post-conjunction semi-major axes of perturber.
    const double preConjunctionSemiMajorAxis 
        = kick->preConjunctionStateInKeplerianElements( semiMajorAxisIndex );
    const double postConjunctionSemiMajorAxis 
        = kick->postConjunctionStateInKeplerianElements( semiMajorAxisIndex );

    // Declare and set pre- and post-conjunction eccentricities of perturber.
    const double preConjunctionEccentricity 
        = kick->preConjunctionStateInKeplerianElements( eccentricityIndex );
    const double postConjunctionEccentricity 
        = kick->postConjunctionStateInKeplerianElements( eccentricityIndex );

    // Declare and set pre- and post-conjunction inclinations of perturber.
    const double preConjunctionInclination 
        = kick->preConjunctionStateInKeplerianElements( inclinationIndex );
    const double postConjunctionInclination 
        = kick->postConjunctionStateInKeplerianElements( inclinationIndex );

    // Compute semi-major axis after kick by using conservation of energy.
    stateInKeplerianElementsAfterKick( semiMajorAxisIndex )
            = 1.0 / ( 1.0 / stateInKeplerianElementsBeforeKick( semiMajorAxisIndex )
                      + perturberMassRatio
                      * ( 1.0 / preConjunctionSemiMajorAxis 
                          - 1.0 / postConjunctionSemiMajorAxis ) );

    // Check that semi-major axis after kick is valid (has to be positive for bounded orbit).
    if ( stateInKeplerianElementsAfterKick( semiMajorAxisIndex ) < 0.0 )
    {
        throw std::runtime_error( "Semi-major axis of perturbed body after kick is erroneous!" );
    }

    // Compute angular momentum of perturbed body before conjunction.
    const double angularMomentumOfPerturbedBodyBeforeConjunction = computeKeplerAngularMomentum(
                stateInKeplerianElementsBeforeKick( semiMajorAxisIndex ),
                stateInKeplerianElementsBeforeKick( eccentricityIndex ), 1.0, 1.0 );

    // Compute angular momentum of perturber before conjunction.
    const double angularMomentumOfPerturberBeforeConjunction
            = computeKeplerAngularMomentum( 
                preConjunctionSemiMajorAxis, preConjunctionEccentricity, 1.0, perturberMassRatio );

    // Compute angular momentum of perturber after conjunction.
    const double angularMomentumOfPerturberAfterConjunction
            = computeKeplerAngularMomentum(
                postConjunctionSemiMajorAxis, postConjunctionEccentricity, 
                1.0, perturberMassRatio );

    // Compute eccentricity after kick by using conservation of angular momentum.
    stateInKeplerianElementsAfterKick( eccentricityIndex )
            = std::sqrt( 1.0 - 1.0 / stateInKeplerianElementsAfterKick( semiMajorAxisIndex )
                         * ( angularMomentumOfPerturbedBodyBeforeConjunction
                             + angularMomentumOfPerturberBeforeConjunction
                             - angularMomentumOfPerturberAfterConjunction )
                         * ( angularMomentumOfPerturbedBodyBeforeConjunction
                             + angularMomentumOfPerturberBeforeConjunction
                             - angularMomentumOfPerturberAfterConjunction ) );

    // Check that eccentricity after kick is valid (has to be 0.0 < e_new < 1.0).
    if ( stateInKeplerianElementsAfterKick( eccentricityIndex ) < 0.0
         || stateInKeplerianElementsAfterKick( eccentricityIndex ) > 1.0
         || boost::math::isnan( stateInKeplerianElementsAfterKick( eccentricityIndex ) ) )
    {
        throw std::runtime_error( "Eccentricity of perturbed body after kick is erroneous!" );
    }            

    // Compute angular momentum of perturbed body after kick.
    const double angularMomentumOfPerturbedBodyAfterConjunction = computeKeplerAngularMomentum(
                stateInKeplerianElementsAfterKick( semiMajorAxisIndex ),
                stateInKeplerianElementsAfterKick( eccentricityIndex ), 1.0, 1.0 );

    // Compute inclination after kick by using conservation of z-component of angular momentum.
    stateInKeplerianElementsAfterKick( inclinationIndex )
            = std::acos( 1.0 / angularMomentumOfPerturbedBodyAfterConjunction
                         * ( angularMomentumOfPerturbedBodyBeforeConjunction
                           * cos( stateInKeplerianElementsBeforeKick( inclinationIndex ) )
                           + angularMomentumOfPerturberBeforeConjunction
                           * cos( preConjunctionInclination )
                           - angularMomentumOfPerturberAfterConjunction
                           * cos( postConjunctionInclination ) ) );

    // Check that inclination after kick is valid (must not be NaN).
    if ( boost::math::isnan( stateInKeplerianElementsAfterKick( inclinationIndex ) ) )
    {
        throw std::runtime_error( "Inclination of perturbed body after kick is erroneous!" );
    }

    // Return state in Keplerian elements after kick.
    return stateInKeplerianElementsAfterKick;
}

// //! Compare elements from a DoubleKeyDoubleValue map.
// bool compareDoubleKeyDoubleValueElements(
//         const common_typedefs::DoubleKeyDoubleValueMap::value_type& element1Value,
//         const common_typedefs::DoubleKeyDoubleValueMap::value_type& element2Value )
// {
//     return element1Value.second < element2Value.second;
// }

// //! Compute average longitude residual for epoch window.
// double computeAverageLongitudeResidualForEpochWindow(
//         const DoubleKeyDoubleValueMap& longitudeResidualsHistory,
//         const double epochWindowStart, const double epochWindowEnd )
// {
//     // Set iterator to start of epoch window.
//     DoubleKeyDoubleValueMap::const_iterator iteratorStartOfEpochWindow
//             = longitudeResidualsHistory.lower_bound( epochWindowStart );

//     // Set iterator to end of epoch window.
//     DoubleKeyDoubleValueMap::const_iterator iteratorEndOfEpochWindow
//             = longitudeResidualsHistory.lower_bound( epochWindowEnd );
//     std::advance( iteratorEndOfEpochWindow, -1 );

//     // Declare value of first moment.
//     double firstMomentOfLongitudeResidualHistory = 0.0;

//     // Set iterator to one element before start of epoch window.
//     DoubleKeyDoubleValueMap::const_iterator iteratorOneBeforeFirstElementEpochWindow
//             = iteratorStartOfEpochWindow;
//     if ( iteratorStartOfEpochWindow != longitudeResidualsHistory.begin( ) )
//     {
//         std::advance( iteratorOneBeforeFirstElementEpochWindow, -1 );
//     }

//     // Check if start and end of epoch window are reversed in order. This means that the epoch
//     // window occurs within one step in the propagation history.
//     if ( iteratorStartOfEpochWindow->first > iteratorEndOfEpochWindow->first )
//     {
//         // Add contribution to first moment.
//         firstMomentOfLongitudeResidualHistory = ( epochWindowEnd - epochWindowStart )
//                 * iteratorOneBeforeFirstElementEpochWindow->second;
//     }

//     else
//     {
//         // Check that start of the epoch window is not the first element in the propagation
//         // history.
//         if ( iteratorStartOfEpochWindow != longitudeResidualsHistory.begin( ) )
//         {
//             firstMomentOfLongitudeResidualHistory
//                     += ( iteratorStartOfEpochWindow->first - epochWindowStart )
//                     * iteratorOneBeforeFirstElementEpochWindow->second;
//         }

//         // Declare iterator to next element in epoch window.
//         DoubleKeyDoubleValueMap::const_iterator iteratorNextElementEpochWindow;

//         // Loop through the epoch window to compute the first moment.
//         for ( DoubleKeyDoubleValueMap::const_iterator iteratorEpochWindow
//               = iteratorStartOfEpochWindow;
//               iteratorEpochWindow != iteratorEndOfEpochWindow; iteratorEpochWindow++ )
//         {
//             iteratorNextElementEpochWindow = iteratorEpochWindow;
//             std::advance( iteratorNextElementEpochWindow, 1 );

//             firstMomentOfLongitudeResidualHistory
//                     += ( iteratorNextElementEpochWindow->first - iteratorEpochWindow->first )
//                     * iteratorEpochWindow->second;
//         }

//         // Compute contribution to the first moment for the end of the epoch window.
//         firstMomentOfLongitudeResidualHistory
//                 += ( epochWindowEnd - iteratorEndOfEpochWindow->first )
//                 * iteratorEndOfEpochWindow->second;
//     }

//     // Return the average longitude residual over the epoch window.
//     return firstMomentOfLongitudeResidualHistory / ( epochWindowEnd - epochWindowStart );
// }

// //! Compute average Keplerian element for epoch window.
// double computeAverageKeplerianElementForEpochWindow(
//         const ActionPropagationHistory& keplerianActionElementsHistory,
//         const double epochWindowStart, const double epochWindowEnd,
//         KeplerianElementVectorIndices keplerianElementIndex )
// {
//     // Set iterator to start of epoch window.
//     ActionPropagationHistory::const_iterator iteratorStartOfEpochWindow
//             = keplerianActionElementsHistory.lower_bound( epochWindowStart );

//     // Set iterator to end of epoch window.
//     ActionPropagationHistory::const_iterator iteratorEndOfEpochWindow
//             = keplerianActionElementsHistory.lower_bound( epochWindowEnd );
//     std::advance( iteratorEndOfEpochWindow, -1 );

//     // Declare first moment.
//     double firstMoment = 0.0;

//     // Set iterator to one element before start of epoch window.
//     ActionPropagationHistory::const_iterator iteratorOneBeforeFirstElementEpochWindow
//             = iteratorStartOfEpochWindow;
//     if ( iteratorStartOfEpochWindow != keplerianActionElementsHistory.begin( ) )
//     {
//         std::advance( iteratorOneBeforeFirstElementEpochWindow, -1 );
//     }

//     // Check if start and end of epoch window are reversed in order. This means that the epoch
//     // window occurs within one step in the propagation history.
//     if ( iteratorStartOfEpochWindow->first > iteratorEndOfEpochWindow->first )
//     {
//         // Add contribution to first moment.
//         firstMoment = ( epochWindowEnd - epochWindowStart )
//                 * iteratorOneBeforeFirstElementEpochWindow->second( keplerianElementIndex );
//     }

//     else
//     {
//         // Check that start of the epoch window is not the first element in the propagation
//         // history.
//         if ( iteratorStartOfEpochWindow != keplerianActionElementsHistory.begin( ) )
//         {
//             // Add contribution to first moment.
//             firstMoment += ( iteratorStartOfEpochWindow->first - epochWindowStart )
//                     * iteratorOneBeforeFirstElementEpochWindow->second( keplerianElementIndex );
//         }

//         // Declare iterator to next element in epoch window.
//         ActionPropagationHistory::const_iterator iteratorNextElementEpochWindow;

//         // Loop through the epoch window to compute the first moment.
//         for ( ActionPropagationHistory::const_iterator iteratorEpochWindow
//               = iteratorStartOfEpochWindow;
//               iteratorEpochWindow != iteratorEndOfEpochWindow; iteratorEpochWindow++ )
//         {
//             // Set iterator to next element in epoch window.
//             iteratorNextElementEpochWindow = iteratorEpochWindow;
//             std::advance( iteratorNextElementEpochWindow, 1 );

//             // Add contribution to first moment.
//             firstMoment += ( iteratorNextElementEpochWindow->first - iteratorEpochWindow->first )
//                     * iteratorEpochWindow->second( keplerianElementIndex );
//         }

//         // Compute contribution to the first moment for the end of the epoch window.
//         firstMoment += ( epochWindowEnd - iteratorEndOfEpochWindow->first )
//                 * iteratorEndOfEpochWindow->second( keplerianElementIndex );
//     }

//     // Return the average Keplerian element over the epoch window.
//     return firstMoment / ( epochWindowEnd - epochWindowStart );
// }

// //! Compute average eccentricity for epoch window.
// double computeAverageEccentricityForEpochWindow(
//         const ActionPropagationHistory& keplerianActionElementsHistory,
//         const double epochWindowStart, const double epochWindowEnd )
// {
//     // Call general function.
//     return computeAverageKeplerianElementForEpochWindow(
//                 keplerianActionElementsHistory,
//                 epochWindowStart, epochWindowEnd, eccentricityIndex );
// }

// //! Compute average inclination for epoch window.
// double computeAverageInclinationForEpochWindow(
//         const ActionPropagationHistory& keplerianActionElementsHistory,
//         const double epochWindowStart, const double epochWindowEnd )
// {
//     // Call general function.
//     return computeAverageKeplerianElementForEpochWindow(
//                 keplerianActionElementsHistory,
//                 epochWindowStart, epochWindowEnd, inclinationIndex );
// }

// //! Compute longitude history.
// DoubleKeyDoubleValueMap computeLongitudeHistory(
//         const ActionPropagationHistory& keplerianActionElementsHistory,
//         const double uranusGravitationalParameter )
// {
//     // Declare longitude history.
//     DoubleKeyDoubleValueMap longitudeHistory;

//     // Set initial longitude to zero.
//     longitudeHistory[ keplerianActionElementsHistory.begin( )->first ] = 0.0;

//     // Declare and set action elements history iterator to one element after start.
//     ActionPropagationHistory::const_iterator iteratorActionElementsHistory
//             = keplerianActionElementsHistory.begin( );
//     std::advance( iteratorActionElementsHistory, 1 );

//     // Declare iterator to previous action elements history element.
//     ActionPropagationHistory::const_iterator iteratorPreviousActionElementsHistory;
//     iteratorPreviousActionElementsHistory = keplerianActionElementsHistory.begin( );

//     // Declare iterator to previous longitude history element.
//     DoubleKeyDoubleValueMap::iterator iteratorPreviousLongitudeHistory;
//     iteratorPreviousLongitudeHistory = longitudeHistory.begin( );

//     // Loop over propagation history.
//     for ( ; iteratorActionElementsHistory != keplerianActionElementsHistory.end( );
//           iteratorActionElementsHistory++ )
//     {
//         longitudeHistory[ iteratorActionElementsHistory->first ]
//                 = iteratorPreviousLongitudeHistory->second
//                 + std::sqrt( uranusGravitationalParameter
//                              / ( iteratorPreviousActionElementsHistory->second(
//                                      semiMajorAxisIndex )
//                                  * iteratorPreviousActionElementsHistory->second(
//                                      semiMajorAxisIndex )
//                                  * iteratorPreviousActionElementsHistory->second(
//                                      semiMajorAxisIndex ) ) )
//                 * ( iteratorActionElementsHistory->first
//                     - iteratorPreviousActionElementsHistory->first );

//         iteratorPreviousActionElementsHistory++;
//         iteratorPreviousLongitudeHistory++;
//     }

//     // Return longitude history.
//     return longitudeHistory;
// }

// //! Reduce longitude history to observation period epoch windows.
// DoubleKeyDoubleValueMap reduceLongitudeHistory( const DoubleKeyDoubleValueMap& longitudeHistory,
//                                                 const double observationPeriodEpoch,
//                                                 const double epochWindowSpacing,
//                                                 const double epochWindowSize,
//                                                 const unsigned int numberOfEpochWindows )
// {
//     // Declare reduced (non-shifted) longitude history.
//     DoubleKeyDoubleValueMap reducedNonShiftedLongitudeHistory;

//     // Declare longitude history iterators.
//     DoubleKeyDoubleValueMap::const_iterator iteratorEpochWindowStart;
//     DoubleKeyDoubleValueMap::const_iterator iteratorEpochWindowEnd;

//     // Loop over observation period and extract longitude epoch windows
//     for ( unsigned int windowNumber = 0; windowNumber < numberOfEpochWindows; windowNumber++ )
//     {
//         const double epochWindowCenter = observationPeriodEpoch
//                 + windowNumber * epochWindowSpacing;

//         iteratorEpochWindowStart = longitudeHistory.lower_bound(
//                     epochWindowCenter - epochWindowSize / 2.0 );

//         if ( iteratorEpochWindowStart != longitudeHistory.begin( ) )
//         {
//             iteratorEpochWindowStart--;
//         }

//         iteratorEpochWindowEnd = longitudeHistory.lower_bound(
//                     epochWindowCenter + epochWindowSize / 2.0 );

//         for ( DoubleKeyDoubleValueMap::const_iterator iteratorEpochWindow = iteratorEpochWindowStart;
//               iteratorEpochWindow != iteratorEpochWindowEnd;
//               iteratorEpochWindow++ )
//         {
//             reducedNonShiftedLongitudeHistory[ iteratorEpochWindow->first ]
//                     = iteratorEpochWindow->second;
//         }
//     }

//     // Declare reduced (shifted) longitude history.
//     DoubleKeyDoubleValueMap reducedLongitudeHistory;

//     // Shift reduced longitude history to start at origin.
//     for ( DoubleKeyDoubleValueMap::iterator iteratorReducedLongitudeHistory
//           = reducedNonShiftedLongitudeHistory.begin( );
//           iteratorReducedLongitudeHistory != reducedNonShiftedLongitudeHistory.end( );
//           iteratorReducedLongitudeHistory++ )
//     {
//         reducedLongitudeHistory[ iteratorReducedLongitudeHistory->first
//                 - reducedNonShiftedLongitudeHistory.begin( )->first ]
//                 = iteratorReducedLongitudeHistory->second
//                 - reducedNonShiftedLongitudeHistory.begin( )->second;
//     }

//     // Return reduced longitude history.
//     return reducedLongitudeHistory;
// }

} // namespace astrodynamics
} // namespace stochastic_migration
