/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOMI_RANDOM_WALK_RUN_H
#define STOMI_RANDOM_WALK_RUN_H

#include <iostream>
#include <string>

#include <boost/shared_ptr.hpp>

namespace stomi
{
namespace database
{

//! Struct that contains all of the information for a set of random walk runs.
/*!
 * This struct contains all of the information for a set of random walk runs, stored in an 
 * SQLite3 database. The data stored is, in essence, meta-data for random walk simulations.
 */
struct RandomWalkRun
{
public:

    //! Constructor taking all random walk run data as input.
    RandomWalkRun(
            // Required parameters.
            const int aRandomWalkRunId,
            const std::string& aRandomWalkRunName, 
            const int aTestParticleCaseId,
            const double aPerturberRingNumberDensity,
            const double aPerturberRingMass,
            const double anObservationPeriod,
            const int aNumberOfEpochWindows,
            const double anEpochWindowSize );

    // Required parameters.

    //! Random walk run ID.
    const int randomWalkRunId;

    //! Random walk run name.
    const std::string randomWalkRunName;

    //! Test particle case ID.
    const int testParticleCaseId;

    //! Perturber ring number density [Number of perturber per Hill radius of perturbed body].
    const double perturberRingNumberDensity;

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

//! Typedef for shared-pointer to RandomWalkRun object.
typedef boost::shared_ptr< RandomWalkRun > RandomWalkRunPointer;

// Define all of the operator overloads as non-member functions (sbi, 2010).

//! Overload == operator.
bool operator==( const RandomWalkRun& RandomWalkRun1,
                 const RandomWalkRun& RandomWalkRun2 );

//! Overload < operator.
bool operator<( const RandomWalkRun& RandomWalkRun1,
                const RandomWalkRun& RandomWalkRun2 );

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, const RandomWalkRun& RandomWalkRun );

} // namespace database
} // namespace stomi

#endif // STOMI_RANDOM_WALK_RUN_H

/*
 * References
 *   sbi. C++ Operator Overloading, Stack Overflow,
 *      http://stackoverflow.com/questions/4421706/operator-overloading, 2010, last accessed:
 *      9th March, 2013.
 */
