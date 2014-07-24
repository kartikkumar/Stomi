/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <stdexcept>
#include <string>
 
#include <Assist/Basics/comparisonFunctions.h>

#include "Stomi/Database/randomWalkInput.h"

namespace stomi
{
namespace database
{

using namespace assist::basics;

//! Default constructor, initializing class members with specified values.
RandomWalkInput::RandomWalkInput( const int aRandomWalkSimulationId,
                                  const int aRandomWalkRunId,
                                  const bool aCompletedFlag,
                                  const double anObservationPeriodStartEpoch,
                                  const std::vector< int >& someTestParticleSimulationIds )
  : randomWalkSimulationId( 
      checkPositive( aRandomWalkSimulationId, "Rnadom walk simulation ID" ) ),
    randomWalkRunId( checkPositive( aRandomWalkRunId, "Random walk run ID" ) ),
    isCompleted( aCompletedFlag ),
    observationPeriodStartEpoch( checkPositive( anObservationPeriodStartEpoch, 
                                                "Observation period start epoch [s]" ) ),
    testParticleSimulationIds( someTestParticleSimulationIds )
{ 
  // Check that vector of test particule simulation IDs is not empty.
  if ( testParticleSimulationIds.size( ) == 0 )
  {
    throw std::runtime_error( 
      "ERROR: Vector of test particle simulation IDs for random walk perturbers is empty!" );
  }
}    

//! Overload == operator.
bool operator==( const RandomWalkInput& randomWalkInput1,
                 const RandomWalkInput& randomWalkInput2 )
{
    return ( randomWalkInput1.randomWalkSimulationId == randomWalkInput2.randomWalkSimulationId
             && randomWalkInput1.randomWalkRunId == randomWalkInput2.randomWalkRunId
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
    return randomWalkInput1.randomWalkSimulationId < randomWalkInput2.randomWalkSimulationId;    
}

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, const RandomWalkInput& randomWalkInput )
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
    outputStream << "Random walk simulation ID: " 
                 << randomWalkInput.randomWalkSimulationId << endl;
    outputStream << "Random walk run ID: " << randomWalkInput.randomWalkRunId << endl;
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
