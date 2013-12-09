/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOCHASTIC_MIGRATION_TEST_PARTICLE_INPUT_H
#define STOCHASTIC_MIGRATION_TEST_PARTICLE_INPUT_H

#include <iostream>
#include <limits>
#include <set>

#include <boost/ptr_container/ptr_set.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

namespace stochastic_migration
{
namespace database
{

//! Data struct that contains input data for one test particle simulation.
/*!
 * This data struct contains input data, stored in an SQLite3 database. The data stored is used to
 * set up a test particle simulation.
 */
struct TestParticleInput
{
public:

    // Set Eigen macro to correctly align class with fixed-size vectorizable types.
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Default constructor, initializing class members.
    TestParticleInput(
            const int aSimulationId,
            const int aTestParticleCaseId,
            const bool aCompletedFlag,
            const tudat::basic_mathematics::Vector6d& anInitialStateInKeplerianElements );

    //! Simulation ID.
    const int simulationId;

    //! Test particle case ID.
    const int testParticleCaseId;

    //! Flag indicating if simulation has been completed/executed and stored in database.
    const bool isCompleted;

    //! Initial state of test particle in Keplerian elements.
    const tudat::basic_mathematics::Vector6d initialStateInKeplerianElements;

protected:
private:
};

//! Typedef for shared-pointer to TestParticleInput object.
typedef boost::shared_ptr< TestParticleInput > TestParticleInputPointer;

//! Typedef for input table (pointers) for test particle simulations.
typedef boost::ptr_set< TestParticleInput > TestParticleInputTable;

// Define all of the operator overloads as non-member functions (sbi, 2010).

//! Overload == operator.
bool operator==( const TestParticleInput& testParticleInput1,
                 const TestParticleInput& testParticleInput2 );

//! Overload < operator.
bool operator<( const TestParticleInput& testParticleInput1,
                const TestParticleInput& testParticleInput2 );

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, const TestParticleInput& testParticleInput );

} // namespace database
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_TEST_PARTICLE_INPUT_H

/*
 *    References
 *      sbi. C++ Operator Overloading, Stack Overflow,
 *          http://stackoverflow.com/questions/4421706/operator-overloading, 2010, last accessed:
 *          9th March, 2013.
 */
 