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
 *      120523    K. Kumar          File created.
 *      130212    K. Kumar          Added Doxygen comments and a note. Added planetary_rings
 *                                  namespace.
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

#ifndef STOCHSTIC_MIGRATION_RANDOM_WALK_MONTE_CARLO_RUN_H
#define STOCHSTIC_MIGRATION_RANDOM_WALK_MONTE_CARLO_RUN_H

#include <iostream>
#include <set>
#include <string>

#include <boost/shared_ptr.hpp>

namespace stochastic_migration
{
namespace database
{

//! Data struct that contains data for a single random walk Monte Carlo run.
/*!
 * This data struct contains all of the data for a single random walk Monte Carlo run, stored in an
 * SQLite3 database. The data stored is metadata for the random walk simulations conducted.
 */
struct RandomWalkMonteCarloRun
{
public:

    //! Default constructor, initializing class members with specified values.
    RandomWalkMonteCarloRun( const int aMonteCarloRun,
                             const int aPerturberPopulation,
                             const std::string& aMassDistributionType,
                             const double aMassDistributionParameter1,
                             const double aMassDistributionParameter2,
                             const double anObservationPeriod,
                             const double anEpochWindowSize,
                             const int aNumberOfEpochWindows );

    //! Monte Carlo run.
    const int monteCarloRun;

    //! Perturber population.
    const int perturberPopulation;

    //! Mass distribution type used for random walk simulation (can be equal, uniform, linear, 
    //! or power-law).
    const std::string massDistributionType;

    //! Mass distribution parameter 1.
    const double massDistributionParameter1;

    //! Mass distribution parameter 2.
    const double massDistributionParameter2;

    //! Observation period [s].
    const double observationPeriod;

    //! Epoch window size used to average observation data [s].
    const double epochWindowSize;

    //! Number of epoch windows within prescribed observation period.
    const int numberOfEpochWindows;

protected:
private:
};

//! Typedef for shared-pointer to RandomWalkMonteCarloRun object.
typedef boost::shared_ptr< RandomWalkMonteCarloRun > RandomWalkMonteCarloRunPointer;

//! Typedef for table of random walk Monte Carlo runs (pointers).
typedef std::set< RandomWalkMonteCarloRunPointer > RandomWalkMonteCarloRunTable;

// Define all of the operator overloads as non-member functions (sbi, 2010).

//! Overload == operator.
bool operator==( const RandomWalkMonteCarloRun& randomWalkMonteCarloRun1,
                 const RandomWalkMonteCarloRun& randomWalkMonteCarloRun2 );

//! Overload < operator.
bool operator<( const RandomWalkMonteCarloRun& randomWalkMonteCarloRun1,
                const RandomWalkMonteCarloRun& randomWalkMonteCarloRun2 );

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, 
                          const RandomWalkMonteCarloRun& randomWalkMonteCarloRun );

} // namespace database
} // namespace stochastic_migration

#endif // STOCHSTIC_MIGRATION_RANDOM_WALK_MONTE_CARLO_RUN_H
