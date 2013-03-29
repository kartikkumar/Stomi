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
 *
 *    Notes
 *      It is unclear how a data struct as as this can be unit tested, other than by constructing
 *      a class object and checking that the typeid of the member variables is as required.
 *
 */

#ifndef STOCHASTIC_MIGRATION_RANDOM_WALK_PERTURBER_H
#define STOCHASTIC_MIGRATION_RANDOM_WALK_PERTURBER_H

#include <iostream>
#include <set>

#include <boost/shared_ptr.hpp>

namespace stochastic_migration
{
namespace database
{

//! Data struct that contains data for a perturber for a random walk Monte Carlo run.
/*!
 * This data struct contains all of the data for a perturber, used in a random walk Monte Carlo
 * run. The data stored is meta data for the random walk simulations conducted.
 */
struct RandomWalkPerturber
{
public:

    //! Default constructor, initializing class members with speficied values.
    RandomWalkPerturber( const int aMonteCarloRun,
                         const int aTestParticleSimulationNumber,
                         const double aMassRatio );

    //! Monte Carlo run.
    const int monteCarloRun;

    //! Test particle simulation number.
    const int testParticleSimulationNumber;

    //! Mass ratio between body receiving kick and body causing kick [-].
    const double massRatio;

protected:
private:
};

//! Typedef for shared-pointer to RandomWalkPerturber object.
typedef boost::shared_ptr< RandomWalkPerturber > RandomWalkPerturberPointer;

//! Typedef for table of perturbers (pointers) for a random walk Monte Carlo run.
typedef std::set< RandomWalkPerturberPointer > RandomWalkPerturberTable;

//! Overload == operator.
bool operator==( const RandomWalkPerturber& randomWalkPerturber1,
                 const RandomWalkPerturber& randomWalkPerturber2 );

//! Overload < operator.
bool operator<( const RandomWalkPerturber& randomWalkPerturber1,
                const RandomWalkPerturber& randomWalkPerturber2 );

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, 
                          const RandomWalkPerturber& randomWalkPerturber );


} // namespace database
} // namespace mab_simulations

#endif // STOCHASTIC_MIGRATION_RANDOM_WALK_PERTURBER_H
