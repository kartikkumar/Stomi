/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */


#include <cmath>
#include <stdexcept>
// #include <iomanip>
// #include <iostream>
#include <iterator>
// #include <limits>
// #include <numeric>
// #include <sstream>

#include <boost/math/special_functions/fpclassify.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
// #include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "StoMi/Astrodynamics/randomWalkFunctions.h"

namespace stomi
{
namespace astrodynamics
{

using namespace assist::basics;

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

//! Compute longitude history.
DoubleKeyDoubleValueMap computeLongitudeHistory(
        const DoubleKeyDoubleValueMap& semiMajorAxisHistory,
        const double centralBodyGravitationalParameter )
{
    // Declare longitude history.
    DoubleKeyDoubleValueMap longitudeHistory;

    // Set initial longitude to zero.
    longitudeHistory[ semiMajorAxisHistory.begin( )->first ] = 0.0;

    // Declare and set semi-major axis history iterator to one element after start.
    DoubleKeyDoubleValueMap::const_iterator iteratorSemiMajorAxisHistory
            = semiMajorAxisHistory.begin( );
    std::advance( iteratorSemiMajorAxisHistory, 1 );

    // Declare iterator to previous semi-major axis history element.
    DoubleKeyDoubleValueMap::const_iterator iteratorPreviousSemiMajorAxisHistory;
    iteratorPreviousSemiMajorAxisHistory = semiMajorAxisHistory.begin( );

    // Declare iterator to previous longitude history element.
    DoubleKeyDoubleValueMap::iterator iteratorPreviousLongitudeHistory;
    iteratorPreviousLongitudeHistory = longitudeHistory.begin( );

    // Loop over propagation history.
    for ( ; iteratorSemiMajorAxisHistory != semiMajorAxisHistory.end( );
          iteratorSemiMajorAxisHistory++ )
    {
        longitudeHistory[ iteratorSemiMajorAxisHistory->first ]
                = iteratorPreviousLongitudeHistory->second
                + computeKeplerMeanMotion( iteratorPreviousSemiMajorAxisHistory->second,
                                           centralBodyGravitationalParameter )
                * ( iteratorSemiMajorAxisHistory->first
                    - iteratorPreviousSemiMajorAxisHistory->first );

        iteratorPreviousSemiMajorAxisHistory++;
        iteratorPreviousLongitudeHistory++;
    }

    // Return longitude history.
    return longitudeHistory;
}

//! Reduce longitude history to observation period epoch windows.
DoubleKeyDoubleValueMap reduceLongitudeHistory( 
        const DoubleKeyDoubleValueMap& longitudeHistory,
        const double observationPeriodStartEpoch,
        const double epochWindowSpacing,
        const double epochWindowSize,
        const unsigned int numberOfEpochWindows )
{
    // Declare reduced longitude history.
    DoubleKeyDoubleValueMap reducedLongitudeHistory;

    // Declare longitude history iterators.
    DoubleKeyDoubleValueMap::const_iterator iteratorEpochWindowStart;
    DoubleKeyDoubleValueMap::const_iterator iteratorEpochWindowEnd;

    // Loop over observation period and extract longitude epoch windows.
    for ( unsigned int windowNumber = 0; windowNumber < numberOfEpochWindows; windowNumber++ )
    {
        // Compute the epoch at the center of the current window.
        const double epochWindowCenter = observationPeriodStartEpoch
                + windowNumber * epochWindowSpacing;

        // Set start iterator to element before start of the epoch window.
        iteratorEpochWindowStart = longitudeHistory.lower_bound(
            epochWindowCenter - 0.5 * epochWindowSize );

        if ( iteratorEpochWindowStart != longitudeHistory.begin( ) )
        {
            iteratorEpochWindowStart--;
        }

        // Set end iterator to element after start of the epoch window.
        iteratorEpochWindowEnd = longitudeHistory.lower_bound(
                    epochWindowCenter + 0.5 * epochWindowSize );

        // Loop through the epoch window and save the data.
        for ( DoubleKeyDoubleValueMap::const_iterator iteratorEpochWindow 
                = iteratorEpochWindowStart;
              iteratorEpochWindow != iteratorEpochWindowEnd;
              iteratorEpochWindow++ )
        {
            reducedLongitudeHistory[ iteratorEpochWindow->first ] = iteratorEpochWindow->second;
        }
    }

    // Return reduced longitude history.
    return reducedLongitudeHistory;
}

} // namespace astrodynamics
} // namespace stomi
