/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOMI_RANDOM_WALK_INPUT_H
#define STOMI_RANDOM_WALK_INPUT_H

#include <iostream>
#include <vector>

#include <boost/ptr_container/ptr_set.hpp>
#include <boost/shared_ptr.hpp>

#include "Stomi/Database/testParticleKick.h" 

namespace stomi
{
namespace database
{

//! Struct that contains data for a single random walk simulation.
/*!
 * This struct contains all of the data for a single random walk simulation, stored in an 
 * SQLite3 database. The data stored is meta-data for the random walk simulations conducted.
 */
struct RandomWalkInput
{
public:

    //! Default constructor, initializing class members with specified values.
    RandomWalkInput( const int aRandomWalkSimulationId,
                     const int aRandomWalkRunId,
                     const bool aCompletedFlag,
                     const double anObservationPeriodStartEpoch,
                     const std::vector< int >& someTestParticleSimulationIds );

    //! Random walk simulation ID.
    const int randomWalkSimulationId;

    //! Random walk run ID.
    const int randomWalkRunId;

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

//! Typedef for input table (pointers) for random walk simulations.
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
} // namespace stomi

#endif // STOMI_RANDOM_WALK_INPUT_H

/*
 * References
 *   sbi. C++ Operator Overloading, Stack Overflow,
 *      http://stackoverflow.com/questions/4421706/operator-overloading, 2010, last accessed:
 *      9th March, 2013.
 */
