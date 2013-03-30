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
 *      130218    K. Kumar          File created.
 *
 *    References
 *
 *    Notes
 *
 */

#include <Assist/Basics/comparisonFunctions.h>

#include "StochasticMigration/Database/randomWalkPerturber.h"

namespace stochastic_migration
{
namespace database
{

using namespace assist::basics;    

//! Default constructor, initializing class members with speficied values.
RandomWalkPerturber::RandomWalkPerturber( const int aMonteCarloRun,
                                          const int aTestParticleSimulationNumber,
                                          const double aMassRatio )
    : monteCarloRun( checkPositive( aMonteCarloRun, "Monte Carlo run" ) ),
      testParticleSimulationNumber( 
        checkPositive( aTestParticleSimulationNumber, "Test particle simulation number ") ),
      massRatio( checkInRange( aMassRatio, "Mass ratio [-]", 0.0, 1.0 ) )
{ }    

//! Overload == operator.
bool operator==( const RandomWalkPerturber& randomWalkPerturber1,
                 const RandomWalkPerturber& randomWalkPerturber2 )
{
    return ( randomWalkPerturber1.monteCarloRun == randomWalkPerturber2.monteCarloRun
             && randomWalkPerturber1.testParticleSimulationNumber
             == randomWalkPerturber2.testParticleSimulationNumber
             && randomWalkPerturber1.massRatio == randomWalkPerturber2.massRatio );
}

//! Overload < operator.
bool operator<( const RandomWalkPerturber& randomWalkPerturber1,
                const RandomWalkPerturber& randomWalkPerturber2 )
{
    if ( randomWalkPerturber1.monteCarloRun == randomWalkPerturber2.monteCarloRun )
    {
        return randomWalkPerturber1.testParticleSimulationNumber 
            < randomWalkPerturber2.testParticleSimulationNumber;
    }

    return randomWalkPerturber1.monteCarloRun < randomWalkPerturber2.monteCarloRun;
}

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream,
                          const RandomWalkPerturber& randomWalkPerturber )
{
    // Write contents of RandomWalkPerturber object to output stream.
    outputStream << "Monte Carlo run: " << randomWalkPerturber.monteCarloRun << std::endl;
    outputStream << "Test particle simulation: "
                 << randomWalkPerturber.testParticleSimulationNumber << std::endl;
    outputStream << "Mass factor [-]: "
                 << randomWalkPerturber.massRatio << std::endl;

    // Return output stream.
    return outputStream;
}

} // namespace database
} // namespace stochastic_migration

