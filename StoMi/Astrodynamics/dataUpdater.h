/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120817    K. Kumar          File created; migrated code from
 *                                  MabSimulations/Astrodynamics/body.h
 *    References
 *
 *    Notes
 *
 */

#ifndef STOMI_DATA_UPDATER_H
#define STOMI_DATA_UPDATER_H

#include <boost/shared_ptr.hpp>

#include "Assist/Basics/commonTypedefs.h"

#include "StoMi/Astrodynamics/body.h"

namespace stomi
{
namespace astrodynamics
{

//! Data updater class.
/*!
 * Essentially, this class serves to provide a wrapper for an update()-function that is called by
 * the composite state derivative model in Tudat. It provides an interface to update a composition
 * of two states, carried by Body objects. The current time associated with each state is also
 * updated.
 * \sa CompositeStateDerivativeModel, Body.
 */
class DataUpdater
{
public:

    //! Constructor taking shared-pointers to bodies.
    /*!
     * Constructor taking shared-pointers to bodies, which are used internally to update the times
     * and states stored in these bodies. This is used in conjunction with the composite state
     * derivative model in Tudat. The order in which the bodies are passed must be respected!
     * \param aBody1 Shared-pointer to body1.
     * \param aBody2 Shared-pointer to body2.
     */
    DataUpdater( const BodyPointer aBody1, const BodyPointer aBody2 )
        : body1( aBody1 ),
          body2( aBody2 )
    { }

    //! Update time and composite state.
    /*!
     * Updates the stored time and composite state information in body1 and body1 and body2.
     * Essentially, this function serves as a means to disassemble the composite state
     * and accordingly update the current state of the two bodies. For this to work, the
     * order of the states in the composite state MUST be respected where ever the composite
     * state is used.
     * \param time Current time.
     * \param compositeState Current composite state (represented as vector of 12 elements).
     */
    void updateTimeAndCompositeState( const double time, 
                                      const assist::basics::Vector12d& compositeState )
    {
        body1->setCurrentTimeAndState( compositeState.segment( 0, 6 ), time  );
        body2->setCurrentTimeAndState( compositeState.segment( 6, 6 ), time  );
    }

protected:
private:

    //! Shared-pointer to body1.
    const BodyPointer body1;

    //! Shared-pointer to body2.
    const BodyPointer body2;
};

//! Typedef for a shared-pointer to a DataUpdater object.
typedef boost::shared_ptr< DataUpdater > DataUpdaterPointer;

} // namespace astrodynamics
} // namespace stomi

#endif // STOMI_DATA_UPDATER_H
