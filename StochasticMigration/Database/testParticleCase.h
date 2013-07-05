/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120402    K. Kumar          File created from old caseDataRow.h.
 *      130212    K. Kumar          Added Doxygen comments and a note. Added planetary_rings
 *                                  namespace; renamed file and struct.
 *      130217    K. Kumar          Updated "mab simulations" references to "stochastic migration";
 *                                  renamed file.
 *      130704    K. Kumar          Updated class contents based on revised table schema.
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
            const int aCaseId,
            const std::string aCaseName,
            const double aRandomWalkSimulationDuration,
            const double aSynodicPeriodLimit,
            const double anOutputInterval,
            const double aStartUpIntegrationDuration,
            const double aConjunctionEventDetectionDistance,
            const double anOppositionEventDetectionDistance,
            const double aCentralBodyGravitationalParameter,
            const double aCentralBodyJ2GravityCoefficient,
            const double aCentralBodyEquatorialRadius,
            const double aSemiMajorAxisDistributionLimit,
            const double anEccentricityDistributionMean,
            const double anEccentricityDistributionAngle,
            const double anEccentricityDistributionFullWidthHalfMaximum,
            const double anInclinationDistributionMean,
            const double anInclinationDistributionAngle,
            const double anInclinationDistributionFullWidthHalfMaximum,
            const double aPerturbedBodyRadius,
            const double aPerturbedBodyBulkDensity,
            const tudat::basic_mathematics::Vector6d& aPerturbedBodyStateInKeplerianElementsAtT0,
            const std::string& aNumericalIntegratorType,
            const double anInitialStepSize,
            const double aNumericalIntegratorRelativeTolerance,
            const double aNumericalIntegratorAbsoluteTolerance );

    //! Case ID.
    const int caseId;

    //! Case name.
    const std::string caseName;

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
    const double semiMajorAxisDistributionLimit;

    //! Mean eccentricity value for distribution.
    const double eccentricityDistributionMean;

    //! Angle between vector components of eccentricity distribution [rad].
    const double eccentricityDistributionAngle;

    //! FWHM eccentricity value for distribution.
    const double eccentricityDistributionFullWidthHalfMaximum;

    //! Mean inclination value for distribution [rad].
    const double inclinationDistributionMean;

    //! Angle between vector components of inclination distribution [rad].
    const double inclinationDistributionAngle;

    //! FWHM inclination value for distribution [rad].
    const double inclinationDistributionFullWidthHalfMaximum;

    //! Perturbed body radius [m].
    const double perturbedBodyRadius;

    //! Perturbed body bulk density [kg m^-3].
    const double perturbedBodyBulkDensity;

    //! Perturbed body state in Keplerian elements at T0.
    const tudat::basic_mathematics::Vector6d perturbedBodyStateInKeplerianElementsAtT0;

    //! Numerical integrator type.
    const std::string numericalIntegratorType;

    //! Initial step size for numerical integrator.
    const double initialStepSize;

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
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_TEST_PARTICLE_CASE_H
