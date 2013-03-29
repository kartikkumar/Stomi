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

#include <sstream>
#include <stdexcept>

#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/comparisonFunctions.h>

#include "StochasticMigration/Database/randomWalkMonteCarloRun.h"

namespace stochastic_migration
{
namespace database
{

using namespace assist::basics;

//! Default constructor, initializing class members with specified values.
RandomWalkMonteCarloRun::RandomWalkMonteCarloRun( const int aMonteCarloRun,
                                                  const int aPerturberPopulation,
                                                  const std::string& aMassDistributionType,
                                                  const double aMassDistributionParameter1,
                                                  const double aMassDistributionParameter2,
                                                  const double anObservationPeriod,
                                                  const double anEpochWindowSize,
                                                  const int aNumberOfEpochWindows )
: monteCarloRun( checkPositive( aMonteCarloRun, "Monte Carlo run" ) ),
  perturberPopulation( checkPositive( aPerturberPopulation, "Perturber population" ) ),
  massDistributionType( aMassDistributionType ),
  massDistributionParameter1( aMassDistributionParameter1 ),
  massDistributionParameter2( aMassDistributionParameter2 ),
  observationPeriod( checkPositive( anObservationPeriod, "Observation period [s]" ) ),
  epochWindowSize( checkPositive( anEpochWindowSize, "Epoch window size [s]" ) ),
  numberOfEpochWindows( checkPositive( aNumberOfEpochWindows, "Number of epoch windows" ) )
{
    if ( aMassDistributionType.empty( ) )
    {
        std::stringstream errorMessage;
        errorMessage << "Mass distribution type string is empty!";
        throw std::runtime_error( errorMessage.str( ) );
    }
}    

//! Overload == operator.
bool operator==( const RandomWalkMonteCarloRun& randomWalkMonteCarloRun1,
                 const RandomWalkMonteCarloRun& randomWalkMonteCarloRun2 )
{
    return ( randomWalkMonteCarloRun1.monteCarloRun == randomWalkMonteCarloRun2.monteCarloRun
         && randomWalkMonteCarloRun1.perturberPopulation
         == randomWalkMonteCarloRun2.perturberPopulation
         && !randomWalkMonteCarloRun1.massDistributionType.compare(
             randomWalkMonteCarloRun2.massDistributionType )
         && randomWalkMonteCarloRun1.massDistributionParameter1
         == randomWalkMonteCarloRun2.massDistributionParameter1
         && randomWalkMonteCarloRun1.massDistributionParameter2
         == randomWalkMonteCarloRun2.massDistributionParameter2
         && randomWalkMonteCarloRun1.observationPeriod
         == randomWalkMonteCarloRun2.observationPeriod
         && randomWalkMonteCarloRun1.epochWindowSize
         == randomWalkMonteCarloRun2.epochWindowSize
         && randomWalkMonteCarloRun1.numberOfEpochWindows
         == randomWalkMonteCarloRun2.numberOfEpochWindows );
}

//! Overload < operator.
bool operator<( const RandomWalkMonteCarloRun& randomWalkMonteCarloRun1,
                const RandomWalkMonteCarloRun& randomWalkMonteCarloRun2 )
{
    return randomWalkMonteCarloRun1.monteCarloRun < randomWalkMonteCarloRun2.monteCarloRun;    
}

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, 
                          const RandomWalkMonteCarloRun& randomWalkMonteCarloRun )
{
    using namespace assist::astrodynamics;

    // Write contents of RandomWalkMonteCarloRun object to output stream.
    outputStream << "Monte Carlo run: " << randomWalkMonteCarloRun.monteCarloRun << std::endl;
    outputStream << "Perturber population: " << randomWalkMonteCarloRun.perturberPopulation
                 << std::endl;
    outputStream << "Mass-distribution type: "
                 << randomWalkMonteCarloRun.massDistributionType << std::endl;
    outputStream << "Mass-distribution parameters: ("
                 << randomWalkMonteCarloRun.massDistributionParameter1
                 << randomWalkMonteCarloRun.massDistributionParameter2 << ")" << std::endl;
    outputStream << "Observation period [Jyr]: "
                 << convertSecondsToJulianYears( randomWalkMonteCarloRun.observationPeriod )
                 << std::endl;
    outputStream << "Epoch window size [Jdays]: "
                 << convertSecondsToJulianDays( randomWalkMonteCarloRun.epochWindowSize )
                 << std::endl;
    outputStream << "Number of epoch windows: "
                 << randomWalkMonteCarloRun.numberOfEpochWindows << std::endl;

    // Return output stream.
    return outputStream;
}

} // namespace database
} // namespace stochastic_migration
