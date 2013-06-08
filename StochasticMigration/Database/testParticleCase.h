/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120402    K. Kumar          File created from old caseDataRow.h.
 *      130212    K. Kumar          Added Doxygen comments and a note. Added planetary_rings
 *                                  namespace; renamed file and struct.
 *      130217    K. Kumar          Updated "mab simulations" references to "stochastic migration";
 *                                  renamed file.
 *
 *    References
 *      sbi. C++ Operator Overloading, Stack Overflow,
 *          http://stackoverflow.com/questions/4421706/operator-overloading, 2010, last accessed:
 *          9th March, 2013.
 *
 *    Notes
 *
 */

#ifndef STOCHASTIC_MIGRATION_TEST_PARTICLE_CASE_H
#define STOCHASTIC_MIGRATION_TEST_PARTICLE_CASE_H

#include <iostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <SQLiteC++.h> 

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

namespace stochastic_migration
{
namespace database
{

//! Data struct that contains all of the case information for a set of test particle simulations.
/*!
 * This data struct contains all of the case information for a set of test particle simulations,
 * stored in an SQLite3 database. The data stored is, in essence, metadata for test particle
 * simulations.
 */
struct TestParticleCase
{
public:

    // Set Eigen macro to correctly align class with fixed-size vectorizable types.
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor taking all case data as input.
    TestParticleCase(
            const int aCaseNumber,
            const double aRandomWalkSimulationDuration,
            const double aSynodicPeriodLimit,
            const double anOutputInterval,
            const double aStartUpIntegrationDuration,
            const double aConjunctionEventDetectionDistance,
            const double anOppositionEventDetectionDistance,
            const double aCentralBodyGravitationalParameter,
            const double aCentralBodyJ2GravityCoefficient,
            const double aCentralBodyEquatorialRadius,
            const double aSemiMajorAxisLimit,
            const double aMeanEccentricity,
            const double aFullWidthHalfMaxmimumEccentricityDistribution,
            const double aMeanInclination,
            const double aFullWidthHalfMaxmimumInclinationDistribution,
            const double aPerturbedBodyGravitationalParameter,
            const tudat::basic_mathematics::Vector6d& aPerturbedBodyStateInKeplerianElementsAtT0,
            const std::string& aNumericalIntegratorType,
            const double aNumericalIntegratorRelativeTolerance,
            const double aNumericalIntegratorAbsoluteTolerance,
            const double anInitialStepSize );

    //! Case number.
    const int caseNumber;

    //! Random walk simulation duration [s].
    const  double randomWalkSimulationDuration;

    //! Synodic period limit [s].
    const double synodicPeriodLimit;

    //! Output interval, determining output frequency for data files [s].
    const double outputInterval;

    //! Startup integration duration [s].
    const double startUpIntegrationDuration;

    //! Mutual distance used to detect start and end of conjunction events [m].
    const double conjunctionEventDetectionDistance;

    //! Distance used to detect start and end of opposition events [m].
    const double oppositionEventDetectionDistance;

    //! Central body gravitational parameter [m^3 s^-2].
    const double centralBodyGravitationalParameter;

    //! Central body J2 gravity field coefficient.
    const double centralBodyJ2GravityCoefficient;

    //! Central body equatorial radius [m].
    const double centralBodyEquatorialRadius;

    //! Limits on maximum semi-major axis values wrt perturbed body [m].
    const double semiMajorAxisLimit;

    //! Mean eccentricity value for distribution.
    const double meanEccentricity;

    //! FWHM eccentricity value for distribution.
    const double fullWidthHalfMaxmimumEccentricityDistribution;

    //! Mean inclination value for distribution.
    const double meanInclination;

    //! FWHM inclination value for distribution.
    const double fullWidthHalfMaxmimumInclinationDistribution;

    //! Perturbed gravitational parameter [m^3 s^-2].
    const double perturbedBodyGravitationalParameter;

    //! Perturbed body state in Keplerian elements at T0.
    const tudat::basic_mathematics::Vector6d perturbedBodyStateInKeplerianElementsAtT0;

    //! Numerical integrator type.
    const std::string numericalIntegratorType;

    //! Relative tolerance for numerical integrator.
    const double numericalIntegratorRelativeTolerance;

    //! Absolute tolerance for numerical integrator.
    const double numericalIntegratorAbsoluteTolerance;

    //! Initial step size for numerical integrator.
    const double initialStepSize;

protected:
private:
};

//! Typedef for shared-pointer to TestParticleCase object.
typedef boost::shared_ptr< TestParticleCase > TestParticleCasePointer;

// Define all of the operator overloads as non-member functions (sbi, 2010).

//! Overload == operator.
bool operator==( const TestParticleCase& testParticleCase1,
                 const TestParticleCase& testParticleCase2 );

//! Overload < operator.
bool operator<( const TestParticleCase& testParticleCase1,
                const TestParticleCase& testParticleCase2 );

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, const TestParticleCase& testParticleCase );

} // namespace database
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_TEST_PARTICLE_CASE_H
