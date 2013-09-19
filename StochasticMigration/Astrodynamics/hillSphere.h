/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130714    K. Kumar          File created.
 *      130715    K. Kumar          Added missing Doxygen documentation, changed naming from
 *                                  kilometers to meters.
 *      130919    K. Kumar          Added missing include guard.
 *
 *    References
 *
 *    Notes
 *      The ConvertHillRadiiToMeters needs to be unit tested.
 *      Should be migrated to the Assist project.
 *
 */

#ifndef STOCHASTIC_MIGRATION_HILL_SPHERE_H
#define STOCHASTIC_MIGRATION_HILL_SPHERE_H

#include <iostream>
#include <cmath>

namespace stochastic_migration
{
namespace astrodynamics
{

//! Functor to convert Hill radii to meters.
class ConvertHillRadiiToMeters
{
public:

    //! Declare default constructor taking parameters to convert Hill radii to meters.
    ConvertHillRadiiToMeters( const double aCentralBodyGravitationalParameter,
                              const double anOrbitingBodyGravitationalParameter,
                              const double aSemiMajorAxis )
        : centralBodyGravitationalParameter( aCentralBodyGravitationalParameter ),
          orbitingBodyGravitationalParameter( anOrbitingBodyGravitationalParameter ),
          semiMajorAxis( aSemiMajorAxis )
    { }

    //! Overload ()-operator to convert Hill radii to meters.
    /*!
     * Overloads ()-operator so that when this function is called with the input
     * given in terms of number of Hill radii, the output returned is the value in
     * meters. This function internally uses the parameters set through the constructor.
     * \param numberOfHillRadii Number of Hill radii to be converted to meters.
     * \return Number of Hill radii expressed in meters.
     */
    double operator( )( const double numberOfHillRadii )
    {
        return numberOfHillRadii * semiMajorAxis 
            * std::pow( orbitingBodyGravitationalParameter 
                / ( 3.0 * centralBodyGravitationalParameter ), 1.0 / 3.0 );
    }

protected:

private:

    //! Gravitational parameter of central body [m^3 s^-2].
    const double centralBodyGravitationalParameter;

    //! Gravitational parameter of orbiting body [m^3 s^-2].
    const double orbitingBodyGravitationalParameter;

    //! Semi-major axis of orbiting body [m].
    const double semiMajorAxis;
};

} // namespace astrodynamics
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_HILL_SPHERE_H
