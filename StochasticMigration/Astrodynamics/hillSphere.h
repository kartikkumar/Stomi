/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130714    K. Kumar          File created.
 *
 *    References
 *
 *    Notes
 *
 */

#include <iostream>
#include <cmath>

namespace stochastic_migration
{
namespace astrodynamics
{

//! Functor to convert Hill radii to kilometers.
class ConvertHillRadiiToKilometers
{
public:

    //! Declare default constructor taking parameters to compute Hill sphere radius.
    ConvertHillRadiiToKilometers( const double aCentralBodyGravitationalParameter,
                                  const double anOrbitingBodyGravitationalParameter,
                                  const double aSemiMajorAxis )
        : centralBodyGravitationalParameter( aCentralBodyGravitationalParameter ),
          orbitingBodyGravitationalParameter( anOrbitingBodyGravitationalParameter ),
          semiMajorAxis( aSemiMajorAxis )
    { }

    //! Overload ()-operator to convert Hill radii to kilometers.
    double operator( )( const double numberOfHillRadii )
    {
        std::cout << "Yaay!" << std::endl;
        return numberOfHillRadii * semiMajorAxis * std::pow( orbitingBodyGravitationalParameter 
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