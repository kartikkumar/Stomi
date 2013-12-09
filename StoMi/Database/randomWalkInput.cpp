/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <string>
 
#include <Assist/Basics/comparisonFunctions.h>

#include "StoMi/Database/randomWalkInput.h"

namespace stomi
{
namespace database
{

using namespace assist::basics;

//! Default constructor, initializing class members with specified values.
RandomWalkInput::RandomWalkInput( const int aMonteCarloRunId,
                                  const int aRandomWalkCaseId,
                                  const bool aCompletedFlag,
                                  const double anObservationPeriodStartEpoch,
                                  const std::vector< int >& someTestParticleSimulationIds )
  : monteCarloRunId( checkPositive( aMonteCarloRunId, "Monte Carlo run ID" ) ),
    randomWalkCaseId( checkPositive( aRandomWalkCaseId, "Random walk case ID" ) ),
    isCompleted( aCompletedFlag ),
    observationPeriodStartEpoch( checkPositive( anObservationPeriodStartEpoch, 
                                                "Observation period start epoch [s]" ) ),
    testParticleSimulationIds( someTestParticleSimulationIds )
{ }    

//! Overload == operator.
bool operator==( const RandomWalkInput& randomWalkInput1,
                 const RandomWalkInput& randomWalkInput2 )
{
    return ( randomWalkInput1.monteCarloRunId == randomWalkInput2.monteCarloRunId
             && randomWalkInput1.randomWalkCaseId == randomWalkInput2.randomWalkCaseId
             && randomWalkInput1.isCompleted == randomWalkInput2.isCompleted
             && randomWalkInput1.observationPeriodStartEpoch 
             == randomWalkInput2.observationPeriodStartEpoch 
             && randomWalkInput1.testParticleSimulationIds 
             == randomWalkInput2.testParticleSimulationIds );
}

//! Overload < operator.
bool operator<( const RandomWalkInput& randomWalkInput1,
                const RandomWalkInput& randomWalkInput2 )
{
    return randomWalkInput1.monteCarloRunId < randomWalkInput2.monteCarloRunId;    
}

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, 
                          const RandomWalkInput& randomWalkInput )
{ 
    using std::endl; 

    // Set string indicating whether random walk simulation has been completed or not.
    std::string completedStatus;

    if ( randomWalkInput.isCompleted == true )
    {
        completedStatus = "TRUE";
    }

    else
    {
        completedStatus = "FALSE";
    }    

    // Write contents of RandomWalkInput object to output stream.
    outputStream << "Monte Carlo run ID: " << randomWalkInput.monteCarloRunId << endl;
    outputStream << "Random walk case ID: " << randomWalkInput.randomWalkCaseId << endl;
    outputStream << "Random walk simulation completed?: " << completedStatus << endl;    
    outputStream << "Observation period start epoch [s]: " 
                 << randomWalkInput.observationPeriodStartEpoch << endl;

    outputStream << "Test particle simulation IDs: ";
    for ( unsigned int i = 0; i < randomWalkInput.testParticleSimulationIds.size( ) - 1; i++ )
    {
      outputStream << randomWalkInput.testParticleSimulationIds.at( i ) << ", ";
    }
    outputStream << randomWalkInput.testParticleSimulationIds.back( ) << endl;

    // Return output stream.
    return outputStream;
}

} // namespace database
} // namespace stomi
