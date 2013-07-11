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

#include <iomanip>
#include <iostream>
#include <limits> 


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

database::TestParticleKickTable propagateSystemAndGenerateKickTable(
        const assist::astrodynamics::BodyPointer perturbedBody,
        const assist::astrodynamics::BodyPointer testParticle,
        const database::TestParticleCasePointer testParticleCase,
        const tudat::numerical_integrators::
            RungeKuttaVariableStepSizeIntegratorXdPointer numericalIntegrator )
{
    // Declare using-statements.
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;

    // Set time shift (this is the amount of time by which all the epochs have to be shifted
    // so that the kicks lie in the range [0.0, randomWalkSimulationDuration]).
    const double timeShift = -( testParticleCase->startUpIntegrationDuration 
                                + testParticleCase->synodicPeriodLimit );

    // Create a table of propagation data points. This is used to store the previous, current and
    // next data points.
    PropagationDataPointTable dataPoints;

    // Create a tables of local maxima and minima.
    PropagationDataPointTable localMaxima;
    PropagationDataPointTable localMinima;

    // Create a tables of opposition (global maxima) and conjunction events (global minima).
    PropagationDataPointTable oppositionEvents;
    PropagationDataPointTable conjunctionEvents;

    // Extract local maxima and minima from propagation data. Numerical integration is done 
    // "on-the-fly", in the sense that after each integration step the check if performed for
    // local extrema.
    {
        // Add propagation data point to table based on initial data. The states of the test 
        // particle and perturbed body are converted to Keplerian elements.
        dataPoints.insert( new PropagationDataPoint( 
            testParticle->getCurrentTime( ) + timeShift,
            ( perturbedBody->getCurrentPosition( ) - testParticle->getCurrentPosition( ) ).norm( ),
            convertCartesianToKeplerianElements( 
                testParticle->getCurrentState( ),
                testParticleCase->centralBodyGravitationalParameter ),
            convertCartesianToKeplerianElements( 
                perturbedBody->getCurrentState( ),
                testParticleCase->centralBodyGravitationalParameter ) ) );        

        // Integrate one step forward and store data in table.
        numericalIntegrator->performIntegrationStep( numericalIntegrator->getNextStepSize( ) );

        dataPoints.insert( new PropagationDataPoint( 
            testParticle->getCurrentTime( ) + timeShift,
            ( perturbedBody->getCurrentPosition( ) - testParticle->getCurrentPosition( ) ).norm( ),
            convertCartesianToKeplerianElements( 
                testParticle->getCurrentState( ),
                testParticleCase->centralBodyGravitationalParameter ),
            convertCartesianToKeplerianElements( 
                perturbedBody->getCurrentState( ),
                testParticleCase->centralBodyGravitationalParameter ) ) );     

        // Integrate one step forward and store data in table.
        numericalIntegrator->performIntegrationStep( numericalIntegrator->getNextStepSize( ) );

        dataPoints.insert( new PropagationDataPoint( 
            testParticle->getCurrentTime( ) + timeShift,
            ( perturbedBody->getCurrentPosition( ) - testParticle->getCurrentPosition( ) ).norm( ),
            convertCartesianToKeplerianElements( 
                testParticle->getCurrentState( ),
                testParticleCase->centralBodyGravitationalParameter ),
            convertCartesianToKeplerianElements( 
                perturbedBody->getCurrentState( ),
                testParticleCase->centralBodyGravitationalParameter ) ) );     

        // Declare iterators to previous, current and next data points.
        PropagationDataPointTable::iterator iteratorDataPoint = dataPoints.begin( );
        const PropagationDataPointTable::iterator iteratorPreviousDataPoint = iteratorDataPoint;

        iteratorDataPoint++;
        const PropagationDataPointTable::iterator iteratorCurrentDataPoint = iteratorDataPoint;

        iteratorDataPoint++;
        const PropagationDataPointTable::iterator iteratorNextDataPoint = iteratorDataPoint;

        assist::basics::DoubleKeyDoubleValueMap distanceHistory;
        distanceHistory[ iteratorCurrentDataPoint->epoch ] = iteratorCurrentDataPoint->mutualDistance;

        int i = 0;

        // Integrate until the end of the simulation period and detect local maxima and minima.
        while ( numericalIntegrator->getCurrentIndependentVariable( )
                < testParticleCase->startUpIntegrationDuration
                + testParticleCase->randomWalkSimulationDuration 
                + 2.0 * testParticleCase->synodicPeriodLimit )
        {
            if ( i == 100000 ) break;

            // Check if the current data point is a local maximum.
            if ( checkIfMaximum( iteratorCurrentDataPoint->mutualDistance,
                                 iteratorPreviousDataPoint->mutualDistance,
                                 iteratorNextDataPoint->mutualDistance ) )
            {
                // Add the current data point to the table of local maxima.
                localMaxima.insert( new PropagationDataPoint( *iteratorCurrentDataPoint ) );
            }

            // Else, check if the current data point is a local minimum.
            else if ( checkIfMinimum( iteratorCurrentDataPoint->mutualDistance,
                                      iteratorPreviousDataPoint->mutualDistance,
                                      iteratorNextDataPoint->mutualDistance ) )
            {
                // Add the current data point to the table of local minima.
                localMinima.insert( new PropagationDataPoint( *iteratorCurrentDataPoint ) );
            }

            // Integrate one step forwards.
            numericalIntegrator->performIntegrationStep( numericalIntegrator->getNextStepSize( ) );

            // Advance previous and current data points. Store new data in next data point.
            iteratorPreviousDataPoint->epoch = iteratorCurrentDataPoint->epoch;
            iteratorPreviousDataPoint->mutualDistance = iteratorCurrentDataPoint->mutualDistance;
            iteratorPreviousDataPoint->testParticleStateInKeplerianElements 
                = iteratorCurrentDataPoint->testParticleStateInKeplerianElements;
            iteratorPreviousDataPoint->perturbedBodyStateInKeplerianElements
                = iteratorCurrentDataPoint->perturbedBodyStateInKeplerianElements;


            iteratorCurrentDataPoint->epoch = iteratorNextDataPoint->epoch;
            iteratorCurrentDataPoint->mutualDistance = iteratorNextDataPoint->mutualDistance;
            iteratorCurrentDataPoint->testParticleStateInKeplerianElements 
                = iteratorNextDataPoint->testParticleStateInKeplerianElements;
            iteratorCurrentDataPoint->perturbedBodyStateInKeplerianElements
                = iteratorNextDataPoint->perturbedBodyStateInKeplerianElements;

            iteratorNextDataPoint->epoch = testParticle->getCurrentTime( ) + timeShift;
            iteratorNextDataPoint->mutualDistance 
                = ( perturbedBody->getCurrentPosition( ) 
                    - testParticle->getCurrentPosition( ) ).norm( );
            iteratorNextDataPoint->testParticleStateInKeplerianElements
                = convertCartesianToKeplerianElements( 
                    testParticle->getCurrentState( ),
                    testParticleCase->centralBodyGravitationalParameter );
            iteratorNextDataPoint->perturbedBodyStateInKeplerianElements
                = convertCartesianToKeplerianElements( 
                    perturbedBody->getCurrentState( ),
                    testParticleCase->centralBodyGravitationalParameter );

            distanceHistory[ iteratorCurrentDataPoint->epoch ] 
                = iteratorCurrentDataPoint->mutualDistance;
            i++;
        }

    tudat::input_output::writeDataMapToTextFile( 
        distanceHistory.begin( ), distanceHistory.end( ),
        "distanceHistory.dat", "/Users/kartikkumar/Desktop", "Epoch,Distance\n",
        std::numeric_limits< double >::digits10, 
        std::numeric_limits< double >::digits10, "," );        
    }

    std::cout << std::endl;
    std::cout << localMaxima.size( ) << "detected!" << std::endl;
    std::cout << std::endl;

    // Loop through local maxima.
    for ( PropagationDataPointTable::iterator it = localMaxima.begin( );
          it != localMaxima.end( ); it++ )
    {
        std::cout << std::setprecision( std::numeric_limits< double >::digits10 )
                  << it->epoch << ", " << it->mutualDistance << std::endl;
    }

    std::cout << std::endl;

    std::cout << localMinima.size( ) << "detected!" << std::endl;

    for ( PropagationDataPointTable::iterator it = localMinima.begin( );
          it != localMinima.end( ); it++ )
    {
        std::cout << std::setprecision( std::numeric_limits< double >::digits10 )
                  << it->epoch << ", " << it->mutualDistance << std::endl;
    }




    // // Extract global maxima (opposition events) from table of local maxima.
    // {
    //     // Declare iterators to previous, current and next local maxima.
    //     PropagationDataPointTable::iterator iteratorLocalMaximum = localMaxima.begin( );
    //     PropagationDataPointTable::iterator iteratorPreviousLocalMaximum = iteratorLocalMaximum;

    //     iteratorLocalMaximum++;
    //     PropagationDataPointTable::iterator iteratorCurrentLocalMaximum = iteratorLocalMaximum;

    //     iteratorLocalMaximum++;
    //     PropagationDataPointTable::iterator iteratorNextLocalMaximum = iteratorLocalMaximum;

    //     // Loop through local maxima.
    //     for ( ; iteratorNextLocalMaximum != localMaxima.end( ); iteratorNextLocalMaximum++ )
    //     {
    //         // Check if the current local maximum is also an opposition event.
    //         if ( checkIfMaximum( iteratorCurrentLocalMaximum->mutualDistance,
    //                              iteratorPreviousLocalMaximum->mutualDistance,
    //                              iteratorNextLocalMaximum->mutualDistance ) )
    //         {
    //             // Add the current data point to the table of opposition events.
    //             oppositionEvents.insert( 
    //                 new PropagationDataPoint( *iteratorCurrentLocalMaximum ) );
    //         }

    //         // Advance iterators.
    //         iteratorPreviousLocalMaximum++;
    //         iteratorCurrentLocalMaximum++;
    //     }
    // }

    // // Extract global minima (conjunction events) from table of local minima.
    // {
    //     // Declare iterators to previous, current and next local minima.
    //     PropagationDataPointTable::iterator iteratorLocalMinimum = localMinima.begin( );
    //     PropagationDataPointTable::iterator iteratorPreviousLocalMinimum = iteratorLocalMinimum;

    //     iteratorLocalMinimum++;
    //     PropagationDataPointTable::iterator iteratorCurrentLocalMinimum = iteratorLocalMinimum;

    //     iteratorLocalMinimum++;
    //     PropagationDataPointTable::iterator iteratorNextLocalMinimum = iteratorLocalMinimum;

    //     // Loop through local minima.
    //     for ( ; iteratorNextLocalMinimum != localMinima.end( ); iteratorNextLocalMinimum++ )
    //     {
    //         // Check if the current local minimum is also a conjunction event.
    //         if ( checkIfMinimum( iteratorCurrentLocalMinimum->mutualDistance,
    //                              iteratorPreviousLocalMinimum->mutualDistance,
    //                              iteratorNextLocalMinimum->mutualDistance ) )
    //         {
    //             // Add the current data point to the table of conjunction events.
    //             conjunctionEvents.insert( 
    //                 new PropagationDataPoint( *iteratorCurrentLocalMinimum ) );
    //         }

    //         // Advance iterators.
    //         iteratorPreviousLocalMinimum++;
    //         iteratorCurrentLocalMinimum++;
    //     }
    // }

    // // Populate kick table.

    std::cout << "Hi! " << std::endl;

}

} // namespace astrodynamics
} // namespace stochastic_migration
