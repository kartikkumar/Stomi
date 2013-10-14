/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOCHSTIC_MIGRATION_RANDOM_WALK_INPUT_H
#define STOCHSTIC_MIGRATION_RANDOM_WALK_INPUT_H

#include <iostream>
#include <vector>

#include <boost/ptr_container/ptr_set.hpp>
#include <boost/shared_ptr.hpp>

#include "StochasticMigration/Database/testParticleKick.h" 

namespace stochastic_migration
{
namespace database
{

//! Data struct that contains data for a single random walk Monte Carlo run.
/*!
 * This data struct contains all of the data for a single random walk Monte Carlo run, stored in an
 * SQLite3 database. The data stored is metadata for the random walk simulations conducted.
 */
struct RandomWalkInput
{
public:

    //! Default constructor, initializing class members with specified values.
    RandomWalkInput( const int aMonteCarloRunId,
                     const int aRandomWalkCaseId,
                     const bool aCompletedFlag,
                     const double anObservationPeriodStartEpoch,
                     const std::vector< int >& someTestParticleSimulationIds );

    //! Monte Carlo run ID.
    const int monteCarloRunId;

    //! Random walk case ID.
    const int randomWalkCaseId;

    //! Flag indicating if simulation has been completed/executed and stored in database.
    const bool isCompleted;    

    //! Epoch at start of observation period [s].
    const double observationPeriodStartEpoch;

    //! List of test particle simulation IDs, used to generate perturbers.
    const std::vector< int > testParticleSimulationIds;

protected:
private:
};

//! Typedef for shared-pointer to RandomWalkInput object.
typedef boost::shared_ptr< RandomWalkInput > RandomWalkInputPointer;

//! Typedef for table of random walk Monte Carlo runs (pointers).
typedef boost::ptr_set< RandomWalkInput > RandomWalkInputTable;

// Define all of the operator overloads as non-member functions (sbi, 2010).

//! Overload == operator.
bool operator==( const RandomWalkInput& randomWalkInput1,
                 const RandomWalkInput& randomWalkInput2 );

//! Overload < operator.
bool operator<( const RandomWalkInput& randomWalkInput1,
                const RandomWalkInput& randomWalkInput2 );

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, 
                          const RandomWalkInput& randomWalkInput );

} // namespace database
} // namespace stochastic_migration

#endif // STOCHSTIC_MIGRATION_RANDOM_WALK_INPUT_H

/*
 *    References
 *      sbi. C++ Operator Overloading, Stack Overflow,
 *          http://stackoverflow.com/questions/4421706/operator-overloading, 2010, last accessed:
 *          9th March, 2013.
 */
