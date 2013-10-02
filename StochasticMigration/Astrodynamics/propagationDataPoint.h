/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOCHASTIC_MIGRATION_PROPAGATION_DATA_POINT_H
#define STOCHASTIC_MIGRATION_PROPAGATION_DATA_POINT_H

#include <iostream>
#include <set> 

#include <boost/shared_ptr.hpp> 

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

namespace stochastic_migration
{
namespace astrodynamics
{

struct PropagationDataPoint
{
    //! Overload default constructor, initializing member variables.
    PropagationDataPoint( 
        const double anEpoch,
        const double aMutualDistance,
        const tudat::basic_mathematics::Vector6d& aTestParticleStateInKeplerianElements,
        const tudat::basic_mathematics::Vector6d& aPerturbedBodyStateInKeplerianElements )
            : epoch( anEpoch ),
              mutualDistance( aMutualDistance ),
              testParticleStateInKeplerianElements( aTestParticleStateInKeplerianElements ),
              perturbedBodyStateInKeplerianElements( aPerturbedBodyStateInKeplerianElements ) 
    { }

    //! Overload copy constructor.
    PropagationDataPoint( const PropagationDataPoint& sourceDataPoint )
            : epoch( sourceDataPoint.epoch ),
              mutualDistance( sourceDataPoint.mutualDistance ),
              testParticleStateInKeplerianElements( 
                sourceDataPoint.testParticleStateInKeplerianElements ),
              perturbedBodyStateInKeplerianElements( 
                sourceDataPoint.perturbedBodyStateInKeplerianElements ) 
    { }

    //! Overload assignment-operator.
    PropagationDataPoint& operator=( const PropagationDataPoint& sourceDataPoint );

    //! Wpoch.
    double epoch;

    //! Mutual distance.
    double mutualDistance;

    //! Test particle state in Keplerian elements.
    tudat::basic_mathematics::Vector6d testParticleStateInKeplerianElements;

    //! Perturbed body state in Keplerian elements.
    tudat::basic_mathematics::Vector6d perturbedBodyStateInKeplerianElements;
};

//! Typedef for shared-pointer to PropagationDataPoint object.
typedef boost::shared_ptr< PropagationDataPoint > PropagationDataPointPointer;

//! Typedef for table of propgation data points test particle simulations.
typedef std::set< PropagationDataPoint > PropagationDataPointTable;

// // Define all of the operator overloads as non-member functions (sbi, 2010).

//! Overload == operator.
bool operator==( const PropagationDataPoint& propagationDataPoint1,
                 const PropagationDataPoint& propagationDataPoint2 );

//! Overload < operator.
bool operator<( const PropagationDataPoint& PropagationDataPoint1,
                const PropagationDataPoint& propagationDataPoint2 );

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, 
                          const PropagationDataPoint& testParticleEvent );

inline bool compareMutualDistances( const PropagationDataPoint& propagationDataPoint1,
                                    const PropagationDataPoint& propagationDataPoint2 )
{
  return propagationDataPoint1.mutualDistance < propagationDataPoint2.mutualDistance;
}

} // namespace astrodynamics
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_PROPAGATION_DATA_POINT_H

/*    
 * The PropagationDataPoint data struct needs to be unit tested.
 */