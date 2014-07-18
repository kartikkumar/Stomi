/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <stdexcept>

#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/comparisonFunctions.h>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Stomi/Database/randomWalkCase.h"

namespace stomi
{
namespace database
{

using namespace assist::basics;

//! Constructor taking all case data as input.
RandomWalkCase::RandomWalkCase(
    // Required parameters.
    const int aCaseId,
    const std::string& aCaseName, 
    const int aTestParticleCaseId,
    const double aPerturberDensity,
    const double aPerturberRingMass,
    const double anObservationPeriod,
    const int aNumberOfEpochWindows,
    const double anEpochWindowSize )
        : caseId( checkPositive( aCaseId, "Case ID" ) ),
          caseName( aCaseName ),
          testParticleCaseId( checkPositive( aTestParticleCaseId, "Test particle case ID" ) ),
          perturberDensity( checkPositive( aPerturberDensity, "Perturber density" ) ),
          perturberRingMass( checkPositive( aPerturberRingMass, "Perturber ring mass" ) ),
          observationPeriod( checkPositive( anObservationPeriod, "Observation period" ) ),
          numberOfEpochWindows( checkPositive( aNumberOfEpochWindows, 
                                               "Number of epoch windows" ) ),
          epochWindowSize( checkLessThan( 
                              checkPositive( anEpochWindowSize, "Epoch window size" ), 
                              "Epoch Window Size", anObservationPeriod ) )
{ 
    if ( caseName.empty( ) )
    {
        throw std::runtime_error( "Case name is empty!" );
    } 
    
    if ( numberOfEpochWindows * epochWindowSize > observationPeriod )
    {
        throw std::runtime_error( "ERROR: Observation period is less than sum of epoch windows!" );
    }  
}   

//! Overload == operator.
bool operator==( const RandomWalkCase& randomWalkCase1,
                 const RandomWalkCase& randomWalkCase2 )
{
    return ( randomWalkCase1.caseId == randomWalkCase2.caseId
             && randomWalkCase1.caseName == randomWalkCase2.caseName
             && randomWalkCase1.testParticleCaseId == randomWalkCase2.testParticleCaseId
             && randomWalkCase1.perturberDensity == randomWalkCase2.perturberDensity
             && randomWalkCase1.perturberRingMass == randomWalkCase2.perturberRingMass
             && randomWalkCase1.observationPeriod == randomWalkCase2.observationPeriod
             && randomWalkCase1.numberOfEpochWindows == randomWalkCase2.numberOfEpochWindows
             && randomWalkCase1.epochWindowSize == randomWalkCase2.epochWindowSize );
}

//! Overload < operator.
bool operator<( const RandomWalkCase& randomWalkCase1,
                const RandomWalkCase& randomWalkCase2 )
{
    return randomWalkCase1.caseId < randomWalkCase2.caseId;
}

//! Overload << operator.
std::ostream& operator<<( std::ostream& outputStream, const RandomWalkCase& randomWalkCase )
{
    using std::endl;
    using namespace assist::astrodynamics;
    using namespace tudat::basic_astrodynamics::unit_conversions;
    
    // Write contents of RandomWalkCase object to output stream.

    // Write required parameters to stream.
    outputStream << endl;
    outputStream << "********************************************************************" << endl;
    outputStream << "Required parameters" << endl;
    outputStream << "********************************************************************" << endl;
    outputStream << endl;

    outputStream << "Case ID: " << randomWalkCase.caseId << endl;
    outputStream << "Case: " << randomWalkCase.caseName << endl; 
    outputStream << "Test particle case ID: " << randomWalkCase.testParticleCaseId << endl;
    outputStream << "Perturber density [#/R_Hill]: " << randomWalkCase.perturberDensity << endl;  
    outputStream << "Perturber ring mass [M_perturbed]: " 
                 << randomWalkCase.perturberRingMass << endl; 
    outputStream << "Observation period [Jyrs]: " 
                 << convertSecondsToJulianYears( randomWalkCase.observationPeriod ) << endl;    
    outputStream << "Number of epoch windows: " << randomWalkCase.numberOfEpochWindows << endl;
    outputStream << "Epoch window size [Jdays]: " 
                 << convertSecondsToJulianDays( randomWalkCase.epochWindowSize ) << endl;        

    // Return output stream.
    return outputStream;
}

} // namespace database
} // namespace stomi    