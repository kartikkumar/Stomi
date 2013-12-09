/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */


#ifndef STOMI_HILL_SPHERE_H
#define STOMI_HILL_SPHERE_H

#include <cmath>

namespace stomi
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
} // namespace stomi

#endif // STOMI_HILL_SPHERE_H

/*
 *    The ConvertHillRadiiToMeters needs to be unit tested. Should be migrated to Assist.
 */
