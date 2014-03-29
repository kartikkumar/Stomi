/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOMI_BODY_H
#define STOMI_BODY_H

#include <string>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

namespace stomi
{
namespace astrodynamics
{

//! Body class.
/*!
 * This class functions as a container for body properties. At the moment, the state of a body can
 * be stored in this class, and accompanying set- and get-functions are provided to use in
 * conjunction with acceleration models and the Cartesian state derivative model in Tudat.
 * \sa AccelerationModel, CartesianStateDerivativeModel.
 */
class Body
{
private:

    //! Typedef for a double-returning function.
    typedef boost::function< double( ) > DoubleReturningFunction;

public:

    // Set Eigen macro to correctly align class with fixed-size vectorizable types.
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor taking the body's name, current state vector, and associated current time.
    /*!
     * Constructor taking the body's name, current state vector, and associated current time
     * (which has a default of 0.0).
     * \param aBodyName Body name.
     * \param aCurrentState Current state vector of body.
     * \param aCurrentTime Current time associated with current state vector of body.
     */
    Body( const std::string& aBodyName,
          const tudat::basic_mathematics::Vector6d& aCurrentState,
          const double aCurrentTime = 0.0 )
        : bodyName( aBodyName ),
          currentPosition( aCurrentState.segment( 0, 3 ) ),
          currentState( aCurrentState ),
          currentVelocity( aCurrentState.segment( 3, 3 ) ),
          currentTime( aCurrentTime )
    { }

    //! Get body name.
    std::string getBodyName( ) { return bodyName; }

    //! Set current state and time.
    void setCurrentTimeAndState( const tudat::basic_mathematics::Vector6d& aCurrentState,
                                 const double aCurrentTime = 0.0 )
    {
        currentTime = aCurrentTime;
        currentState = aCurrentState;
        currentPosition = aCurrentState.segment( 0, 3 );
        currentVelocity = aCurrentState.segment( 3, 3 );
    }

    //! Get current time.
    double getCurrentTime( ) { return currentTime; }

    //! Get current position.
    Eigen::Vector3d getCurrentPosition( ) { return currentPosition; }

    //! Get current position.
    Eigen::Vector3d getCurrentVelocity( ) { return currentVelocity; }

    //! Get current state.
    tudat::basic_mathematics::Vector6d getCurrentState( ) { return currentState; }

protected:

private:

    //! Body name.
    const std::string bodyName;

    //! Current position vector of body.
    Eigen::Vector3d currentPosition;

    //! Current state vector of body.
    tudat::basic_mathematics::Vector6d currentState;

    //! Current velocity vector of body.
    Eigen::Vector3d currentVelocity;

    //! Current time associated with current position vector of body.
    double currentTime;
};

//! Typedef for a shared-pointer to a Body object.
typedef boost::shared_ptr< Body > BodyPointer;

} // namespace astrodynamics
} // namespace stomi

#endif // STOMI_BODY_H
