/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOCHASTIC_MIGRATION_RANDOM_WALK_CASE_H
#define STOCHASTIC_MIGRATION_RANDOM_WALK_CASE_H

#include <iostream>
#include <string>

#include <boost/shared_ptr.hpp>

namespace stochastic_migration
{
namespace database
{

//! Data struct that contains all of the case information for a set of random walk simulations.
/*!
 * This data struct contains all of the case information for a set of random walk simulations,
 * stored in an SQLite3 database. The data stored is, in essence, metadata for random walk
 * simulations.
 */
struct RandomWalkCase
{
public:

    //! Constructor taking all case data as input.
    RandomWalkCase(
            // Required parameters.
            const int aCaseId,
            const std::string& aCaseName, 
            const int aTestParticleCaseId,
            const double aPerturberDensity,
            const double aPerturberRingMass,
            const double anObservationPeriod,
            const int aNumberOfEpochWindows,
            const double anEpochWindowSize );

    // Required parameters.

    //! Case ID.
    const int caseId;

    //! Case name.
    const std::string caseName;

    //! Test particle case ID.
    const int testParticleCaseId;

    //! Perturber density [Number of perturber per Hill radius of perturbed body].
    const double perturberDensity;

    //! Perturber ring mass [M_perturbedBody].
    const double perturberRingMass;

    //! Observation period [s].
    const double observationPeriod;

    //! Number of epoch windows within observation period.
    const int numberOfEpochWindows;

    //! Epoch window size [s].
    const double epochWindowSize;

protected:
private:
};

//! Typedef for shared-pointer to RandomWalkCase object.
typedef boost::shared_ptr< RandomWalkCase > RandomWalkCasePointer;

// Define all of the operator overloads as non-member functions (sbi, 2010).

//! Overload == operator.
bool operator==( const RandomWalkCase& randomWalkCase1,
                 const RandomWalkCase& randomWalkCase2 );

//! Overload < operator.
bool operator<( const RandomWalkCase& randomWalkCase1,
                const RandomWalkCase& randomWalkCase2 );

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, const RandomWalkCase& randomWalkCase );

} // namespace database
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_RANDOM_WALK_CASE_H

/*
 *    References
 *      sbi. C++ Operator Overloading, Stack Overflow,
 *          http://stackoverflow.com/questions/4421706/operator-overloading, 2010, last accessed:
 *          9th March, 2013.
 */
