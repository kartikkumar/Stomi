/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOMI_TEST_PARTICLE_CASE_H
#define STOMI_TEST_PARTICLE_CASE_H

#include <iostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <SQLiteC++.h> 

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

namespace stomi
{
namespace database
{

 enum NumericalIntegratorType { DOPRI853, RKF78 };   

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
            // Required parameters.
            const int aCaseId,
            const std::string& aCaseName,
            const double aRandomWalkSimulationPeriod,
            const double aCentralBodyGravitationalParameter,
            const double aPerturbedBodyRadius,
            const double aPerturbedBodyBulkDensity,
            const tudat::basic_mathematics::Vector6d& aPerturbedBodyStateInKeplerianElementsAtT0,
            const double aSemiMajorAxisDistributionLimit,
            // Optional parameters.
            const double aSynodicPeriodMaximum,
            const double aStartUpIntegrationPeriod,
            const double aCentralBodyJ2GravityCoefficient,
            const double aCentralBodyEquatorialRadius,            
            const double aConjunctionEventDetectionDistance,
            const double anOppositionEventDetectionDistance,
            const double anEccentricityDistributionMean,
            const double anEccentricityDistributionFullWidthHalfMaximum,
            const double anInclinationDistributionMean,
            const double anInclinationDistributionFullWidthHalfMaximum,
            const std::string& aNumericalIntegratorType,
            const double aNumericalIntegratorInitialStepSize,
            const double aNumericalIntegratorRelativeTolerance,
            const double aNumericalIntegratorAbsoluteTolerance );

    // Required parameters.

    //! Case ID.
    const int caseId;

    //! Case name.
    const std::string caseName;

    //! Random walk simulation period [s].
    const double randomWalkSimulationPeriod;

    //! Central body gravitational parameter [m^3 s^-2].
    const double centralBodyGravitationalParameter;

    //! Perturbed body radius [m].
    const double perturbedBodyRadius;

    //! Perturbed body bulk density [kg m^-3].
    const double perturbedBodyBulkDensity;

    //! Perturbed body state in Keplerian elements at T0.
    const tudat::basic_mathematics::Vector6d perturbedBodyStateInKeplerianElementsAtT0;

    //! Limits on maximum semi-major axis values wrt perturbed body [m].
    const double semiMajorAxisDistributionLimit;

    // Optional parameters.

    //! Maximum synodic period permitted [s].
    const double synodicPeriodMaximum;

    //! Startup integration period [s].
    const double startUpIntegrationPeriod;

    //! Central body J2 gravity field coefficient.
    const double centralBodyJ2GravityCoefficient;

    //! Central body equatorial radius [m].
    const double centralBodyEquatorialRadius;

    //! Mutual distance used to detect start and end of conjunction events [m].
    const double conjunctionEventDetectionDistance;

    //! Distance used to detect start and end of opposition events [m].
    const double oppositionEventDetectionDistance;

    //! Mean eccentricity value for distribution.
    const double eccentricityDistributionMean;

    //! FWHM eccentricity value for distribution.
    const double eccentricityDistributionFullWidthHalfMaximum;

    //! Mean inclination value for distribution [rad].
    const double inclinationDistributionMean;

    //! FWHM inclination value for distribution [rad].
    const double inclinationDistributionFullWidthHalfMaximum;

    //! Numerical integrator type.
    NumericalIntegratorType numericalIntegratorType;

    //! Initial step size for numerical integrator.
    const double numericalIntegratorInitialStepSize;

    //! Relative tolerance for numerical integrator.
    const double numericalIntegratorRelativeTolerance;

    //! Absolute tolerance for numerical integrator.
    const double numericalIntegratorAbsoluteTolerance;

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
} // namespace stomi

#endif // STOMI_TEST_PARTICLE_CASE_H

/*
 *    References
 *      sbi. C++ Operator Overloading, Stack Overflow,
 *          http://stackoverflow.com/questions/4421706/operator-overloading, 2010, last accessed:
 *          9th March, 2013.
 */
