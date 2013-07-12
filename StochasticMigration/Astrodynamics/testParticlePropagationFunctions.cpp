/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120817    K. Kumar          File created.
 *      130218    K. Kumar          File renamed; namespaces updated; implementation of functions
 *                                  updated; updated "encounter" to "conjunction".
 *
 *    References
 *
 *    Notes
 *      The functions implemented in this file need to be unit tested.
 *
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits> 
#include <sstream> 


#include <algorithm>
//#include <cmath>
//#include <map>
//#include <sstream>
//#include <stdexcept>
#include <utility>

#include <Assist/Basics/commonTypedefs.h>
#include <Assist/Basics/comparisonFunctions.h>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include <Tudat/InputOutput/basicInputOutput.h>

#include "StochasticMigration/Astrodynamics/propagationDataPoint.h"
#include "StochasticMigration/Astrodynamics/testParticlePropagationFunctions.h"
//#include "StochasticMigration/Database/testParticleKick.h"

namespace stochastic_migration
{
namespace astrodynamics
{


//using std::runtime_error;

//using namespace tudat::numerical_integrators;

//using namespace database;

// database::TestParticleKickTable propagateSystemAndGenerateKickTable(
//         const assist::astrodynamics::BodyPointer perturbedBody,
//         const assist::astrodynamics::BodyPointer testParticle,
//         const database::TestParticleCasePointer testParticleCase,
//         const tudat::numerical_integrators::
//             RungeKuttaVariableStepSizeIntegratorXdPointer numericalIntegrator )
void propagateSystemAndGenerateKickTable(
        const assist::astrodynamics::BodyPointer perturbedBody,
        const assist::astrodynamics::BodyPointer testParticle,
        const database::TestParticleCasePointer testParticleCase,
        const tudat::numerical_integrators::
            RungeKuttaVariableStepSizeIntegratorXdPointer numericalIntegrator, 
        const double synodicPeriod, const double nextStepSize, const int i )
{
    // TEMP
    PropagationDataPointTable dataPointsAll;

    // Declare using-statements.
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;

    // Set time shift (this is the amount of time by which all the epochs have to be shifted
    // so that the kicks lie in the range [0.0, randomWalkSimulationDuration]).
    const double timeShift = -( numericalIntegrator->getCurrentIndependentVariable( ) 
                                + synodicPeriod );

    // Create a table of propagation data points.
    PropagationDataPointTable dataPoints;
    PropagationDataPointTable::iterator iteratorDataPoints = dataPoints.begin( );

    // Create a tables of opposition (global maxima) and conjunction events (global minima).
    PropagationDataPointTable oppositionEvents;
    PropagationDataPointTable conjunctionEvents;

    // Compute mutual distance between test particle and perturbed body [m].
    double mutualDistance = ( perturbedBody->getCurrentPosition( ) 
            - testParticle->getCurrentPosition( ) ).norm( );

    dataPointsAll.insert( PropagationDataPoint( 
        testParticle->getCurrentTime( ) + timeShift,
        mutualDistance,
        convertCartesianToKeplerianElements( 
            testParticle->getCurrentState( ),
            testParticleCase->centralBodyGravitationalParameter ),
        convertCartesianToKeplerianElements( 
            perturbedBody->getCurrentState( ),
            testParticleCase->centralBodyGravitationalParameter ) ) ); 
        
    // Integrate one step forwards.
    numericalIntegrator->performIntegrationStep( nextStepSize );

    // Integrate until the first opposition event is detected.
    while ( mutualDistance < testParticleCase->oppositionEventDetectionDistance )
    {
        dataPointsAll.insert( PropagationDataPoint( 
            testParticle->getCurrentTime( ) + timeShift,
            mutualDistance,
            convertCartesianToKeplerianElements( 
                testParticle->getCurrentState( ),
                testParticleCase->centralBodyGravitationalParameter ),
            convertCartesianToKeplerianElements( 
                perturbedBody->getCurrentState( ),
                testParticleCase->centralBodyGravitationalParameter ) ) ); 

        // Integrate one step forwards.
        numericalIntegrator->performIntegrationStep( numericalIntegrator->getNextStepSize( ) );

        // Compute mutual distance based on new test particle and perturbed body states.
        mutualDistance = ( perturbedBody->getCurrentPosition( ) 
            - testParticle->getCurrentPosition( ) ).norm( );                 
    }

    // Set flag indicating if an opposition event has been detected to true.
    bool isOppositionEventDetected = true;

    // Integrate until the end of the simulation period and detect local maxima and minima.
    while ( numericalIntegrator->getCurrentIndependentVariable( )
            < testParticleCase->startUpIntegrationDuration
            + testParticleCase->randomWalkSimulationDuration + 2.0 * synodicPeriod )
    {
        // Compute mutual distance between test particle and perturbed body.
        mutualDistance = ( perturbedBody->getCurrentPosition( ) 
            - testParticle->getCurrentPosition( ) ).norm( );

        // Check if the current event detected is an opposition event.
        if ( isOppositionEventDetected )
        {
            // Check if end of opposition event is detected (i.e., start of conjunction event).
            if ( mutualDistance < testParticleCase->conjunctionEventDetectionDistance )
            {
                // Set flag to false.
                isOppositionEventDetected = false;

                // Find data point at opposition event.
                iteratorDataPoints = std::max_element( 
                    dataPoints.begin( ), dataPoints.end( ), compareMutualDistances );

                // Add opposition event to table.
                oppositionEvents.insert( PropagationDataPoint( *iteratorDataPoints ) );

                // Clear table of data points.
                dataPoints.clear( );
            }
        }

        // Else, the current event must be a conjunction event.
        else
        {
            // Check if end of conjunction event is detected (i.e., start of opposition event).
            if ( mutualDistance > testParticleCase->conjunctionEventDetectionDistance )
            {
                // Set flag to true.
                isOppositionEventDetected = true;

                // Find data point at conjunctiion event.
                iteratorDataPoints = std::min_element( 
                    dataPoints.begin( ), dataPoints.end( ), compareMutualDistances );

                // Add conjunction event to table.
                conjunctionEvents.insert( PropagationDataPoint( *iteratorDataPoints ) );

                // Clear table of data points.
                dataPoints.clear( );
            }            
        }

        // Add data point to table.
        dataPoints.insert( PropagationDataPoint( 
            testParticle->getCurrentTime( ) + timeShift,
            mutualDistance,
            convertCartesianToKeplerianElements( 
                testParticle->getCurrentState( ),
                testParticleCase->centralBodyGravitationalParameter ),
            convertCartesianToKeplerianElements( 
                perturbedBody->getCurrentState( ),
                testParticleCase->centralBodyGravitationalParameter ) ) ); 

        dataPointsAll.insert( PropagationDataPoint( 
            testParticle->getCurrentTime( ) + timeShift,
            mutualDistance,
            convertCartesianToKeplerianElements( 
                testParticle->getCurrentState( ),
                testParticleCase->centralBodyGravitationalParameter ),
            convertCartesianToKeplerianElements( 
                perturbedBody->getCurrentState( ),
                testParticleCase->centralBodyGravitationalParameter ) ) );  

        // Integrate one step forward and store data in table.
        numericalIntegrator->performIntegrationStep( numericalIntegrator->getNextStepSize( ) );
    }

    if ( conjunctionEvents.size( ) == 0 || oppositionEvents.size( ) == 0 ) 
    {
         std::cout << "AHHHHH!" << std::endl; exit( 0 ); 
    }

    // Check if epoch of last conjunction event is greater than that of the last opposition event.
    // If so, remove from the table.
    if ( conjunctionEvents.rbegin( )->epoch > oppositionEvents.rbegin( )->epoch )
    {
        // Set iterator to last element in table of conjunction events.
        PropagationDataPointTable::iterator iteratorLastElement = conjunctionEvents.end( );
        iteratorLastElement--;

        // Erase last element.
        conjunctionEvents.erase( iteratorLastElement );
    }

    // Check if the epoch of the last conjunction event is beyond the random walk simulation 
    // window [0, randomWalkSimulationDuration].
    // If so, delete the last conjunction event and the last opposition event.
    // Repeat until the last epoch is within the window.
    while ( conjunctionEvents.rbegin( )->epoch > testParticleCase->randomWalkSimulationDuration )
    {
        // Set iterator to last element in table of conjunction events.
        PropagationDataPointTable::iterator iteratorLastElement = conjunctionEvents.end( );
        iteratorLastElement--;

        // Erase last element.
        conjunctionEvents.erase( iteratorLastElement );

        // Set iterator to last element in table of opposition events.
        iteratorLastElement = oppositionEvents.end( );
        iteratorLastElement--;

        // Erase last element.
        oppositionEvents.erase( iteratorLastElement );
    }

    // Check if the epoch of the first conjunction event is before the random walk simulation 
    // window [0.0, randomWalkSimulationDuration].
    // If so, delete the first conjunction event and the first opposition event.
    // Repeat until the first epoch is within the window.
    while ( conjunctionEvents.begin( )->epoch < 0.0 )
    {
        // Erase first element.
        conjunctionEvents.erase( conjunctionEvents.begin( ) );

        // Erase last element.
        oppositionEvents.erase( oppositionEvents.begin( ) );
    }    

    std::ostringstream fileName;
    fileName << "/Users/kartikkumar/Desktop/data" << i << ".csv";
    std::ofstream file( fileName.str( ).c_str( ) );

    file << "t,d" << std::endl;

    int counter = 0;

    for ( PropagationDataPointTable::iterator it = dataPointsAll.begin( );
          it != dataPointsAll.end( );
          it++ )
    {
        if ( it->epoch > testParticleCase->outputInterval * counter - synodicPeriod )
        {
            file << std::setprecision( std::numeric_limits< double >::digits10 )
                 << it->epoch << "," << it->mutualDistance << std::endl;
            counter++;
        }
    }

    file.close( );

    std::ostringstream fileName2;
    fileName2 << "/Users/kartikkumar/Desktop/opposition" << i << ".csv";
    std::ofstream file2( fileName2.str( ).c_str( ) );

    file2 << "t,d" << std::endl;

    for ( PropagationDataPointTable::iterator it = oppositionEvents.begin( );
          it != oppositionEvents.end( );
          it++ )
    {
        file2 << std::setprecision( std::numeric_limits< double >::digits10 )
              << it->epoch << "," << it->mutualDistance << std::endl;
    }

    file2.close( );

    std::ostringstream fileName3;
    fileName3 << "/Users/kartikkumar/Desktop/conjunction" << i << ".csv";
    std::ofstream file3( fileName3.str( ).c_str( ) );

    file3 << "t,d" << std::endl;

    for ( PropagationDataPointTable::iterator it = conjunctionEvents.begin( );
          it != conjunctionEvents.end( );
          it++ )
    {
        file3 << std::setprecision( std::numeric_limits< double >::digits10 )
              << it->epoch << "," << it->mutualDistance << std::endl;
    }

    file3.close( );

    std::cout << "Hi!" << std::endl;
}

} // namespace astrodynamics
} // namespace stochastic_migration
