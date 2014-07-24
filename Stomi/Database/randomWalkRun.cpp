/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <stdexcept>

#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/comparisonFunctions.h>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Stomi/Database/randomWalkRun.h"

namespace stomi
{
namespace database
{

using namespace assist::basics;

//! Constructor taking all random walk run data as input.
RandomWalkRun::RandomWalkRun(
    // Required parameters.
    const int aRandomWalkRunId,
    const std::string& aRandomWalkRunName, 
    const int aTestParticleCaseId,
    const double aPerturberRingNumberDensity,
    const double aPerturberRingMass,
    const double anObservationPeriod,
    const int aNumberOfEpochWindows,
    const double anEpochWindowSize )
        : randomWalkRunId( checkPositive( aRandomWalkRunId, "Random walk run ID" ) ),
          randomWalkRunName( aRandomWalkRunName ),
          testParticleCaseId( checkPositive( aTestParticleCaseId, "Test particle case ID" ) ),
          perturberRingNumberDensity( 
            checkPositive( aPerturberRingNumberDensity, "Perturber ring number density" ) ),
          perturberRingMass( checkPositive( aPerturberRingMass, "Perturber ring mass" ) ),
          observationPeriod( checkPositive( anObservationPeriod, "Observation period" ) ),
          numberOfEpochWindows( checkPositive( aNumberOfEpochWindows, 
                                               "Number of epoch windows" ) ),
          epochWindowSize( checkLessThan( 
                              checkPositive( anEpochWindowSize, "Epoch window size" ), 
                              "Epoch Window Size", anObservationPeriod ) )
{ 
    if ( randomWalkRunName.empty( ) )
    {
        throw std::runtime_error( "Random walk run name is empty!" );
    } 
    
    if ( numberOfEpochWindows * epochWindowSize > observationPeriod )
    {
        throw std::runtime_error( "ERROR: Observation period is less than sum of epoch windows!" );
    }  
}   

//! Overload == operator.
bool operator==( const RandomWalkRun& randomWalkRun1,
                 const RandomWalkRun& randomWalkRun2 )
{
    return ( randomWalkRun1.randomWalkRunId == randomWalkRun2.randomWalkRunId
             && randomWalkRun1.randomWalkRunName == randomWalkRun2.randomWalkRunName
             && randomWalkRun1.testParticleCaseId == randomWalkRun2.testParticleCaseId
             && randomWalkRun1.perturberRingNumberDensity 
             == randomWalkRun2.perturberRingNumberDensity
             && randomWalkRun1.perturberRingMass == randomWalkRun2.perturberRingMass
             && randomWalkRun1.observationPeriod == randomWalkRun2.observationPeriod
             && randomWalkRun1.numberOfEpochWindows == randomWalkRun2.numberOfEpochWindows
             && randomWalkRun1.epochWindowSize == randomWalkRun2.epochWindowSize );
}

//! Overload < operator.
bool operator<( const RandomWalkRun& randomWalkRun1,
                const RandomWalkRun& randomWalkRun2 )
{
    return randomWalkRun1.randomWalkRunId < randomWalkRun2.randomWalkRunId;
}

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, const RandomWalkRun& randomWalkRun )
{
    using std::endl;
    using namespace assist::astrodynamics;
    using namespace tudat::basic_astrodynamics::unit_conversions;
    
    // Write contents of RandomWalkRun object to output stream.

    // Write required parameters to stream.
    outputStream << endl;
    outputStream << "********************************************************************" << endl;
    outputStream << "Required parameters" << endl;
    outputStream << "********************************************************************" << endl;
    outputStream << endl;

    outputStream << "Random walk run ID: " << randomWalkRun.randomWalkRunId << endl;
    outputStream << "Random walk run: " << randomWalkRun.randomWalkRunName << endl; 
    outputStream << "Test particle case ID: " << randomWalkRun.testParticleCaseId << endl;
    outputStream << "Perturber ring number density [#/R_Hill]: " 
                 << randomWalkRun.perturberRingNumberDensity << endl;  
    outputStream << "Perturber ring mass [M_perturbed]: " 
                 << randomWalkRun.perturberRingMass << endl; 
    outputStream << "Observation period [Jyrs]: " 
                 << convertSecondsToJulianYears( randomWalkRun.observationPeriod ) << endl;    
    outputStream << "Number of epoch windows: " << randomWalkRun.numberOfEpochWindows << endl;
    outputStream << "Epoch window size [Jdays]: " 
                 << convertSecondsToJulianDays( randomWalkRun.epochWindowSize ) << endl;        

    // Return output stream.
    return outputStream;
}

} // namespace database
} // namespace stomi    
