/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120402    K. Kumar          File created from old mabRandomWalk.cpp.
 *      120521    K. Kumar          Updated to use input file; updated to SQLite C API; fixed bug
 *                                  in sliding-window computations.
 *      120522    K. Kumar          Added inclination changes computation; changed filename to
 *                                  randomWalkSimulator.cpp.
 *      130729    K. Kumar          Ported over old cold to GitHub project and updated.
 *
 *    References
 *      Kumar, K., de Pater, I., Showalter, M.R. In prep, 2013.
 *
 *    Notes
 *
 */

// #include <algorithm>
// #include <ctime>
// #include <iomanip>
#include <iostream>
// #include <iterator>
// #include <limits>
// #include <stdexcept>
#include <string>
// #include <sstream>
// #include <vector>

// #include <omp.h>

// #include <boost/algorithm/string/classification.hpp>
// #include <boost/algorithm/string/split.hpp>
// #include <boost/exception/all.hpp>
// #include <boost/lexical_cast.hpp>
// #include <boost/random/mersenne_twister.hpp>
// #include <boost/random/exponential_distribution.hpp>
// #include <boost/random/uniform_real_distribution.hpp>
// #include <boost/random/uniform_int_distribution.hpp>
// #include <boost/random/variate_generator.hpp>

// #include <sqlite3.h>

#include <Assist/InputOutput/basicInputOutput.h>

// #include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/InputOutput/dictionaryTools.h>
#include <Tudat/InputOutput/fieldType.h>
#include <Tudat/InputOutput/separatedParser.h>
#include <Tudat/InputOutput/parsedDataVectorUtilities.h>
// #include <Tudat/Mathematics/Statistics/simpleLinearRegression.h>

// #include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>

// #include "AuxilliaryFiles/commonTypedefs.h"
// #include "AuxilliaryFiles/randomWalkFunctions.h"

// #include "Database/caseDataRow.h"
// #include "Database/databaseFunctions.h"

// #include "InputOutput/basicInputOutput.h"
#include "StochasticMigration/InputOutput/dictionaries.h"

//! Execute random walk simulations.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    using std::cout;
    using std::endl;
    using std::string;
//     using std::generate;
//     using std::runtime_error;

//     using namespace boost::random;
//     using boost::enable_error_info;
//     using boost::is_any_of;
//     using boost::lexical_cast;
//     using boost::mt19937;
//     using boost::split;
//     using boost::throw_exception;
//     using boost::variate_generator;

    using namespace assist::input_output;

//     using namespace tudat::basic_astrodynamics::physical_constants;
    using namespace tudat::input_output;
    using namespace tudat::input_output::dictionary;
    using namespace tudat::input_output::field_types::general;
    using namespace tudat::input_output::parsed_data_vector_utilities;
//     using namespace tudat::statistics;

//     using namespace stochastic_migration::common_typedefs;
//     using namespace stochastic_migration::database;
//     using namespace stochastic_migration::database_functions;
    using namespace stochastic_migration::input_output;
//     using namespace stochastic_migration::random_walk_functions;

    ///////////////////////////////////////////////////////////////////////////

    // Input deck.

    // Check input arguments.
    checkNumberOfInputArguments( numberOfInputs - 1 );

    // Get input parameter dictionary.
    const DictionaryPointer dictionary = getRandomWalkSimulatorDictionary( );

    // Read and filter input stream.
    string filteredInput = readAndFilterInputFile( inputArguments[ 1 ] );

    // Declare a separated parser.
    SeparatedParser parser( string( ": " ), 2, parameterName, parameterValue );

    // Parse filtered data.
    const ParsedDataVectorPtr parsedData = parser.parse( filteredInput );

    cout << endl;
    cout << "****************************************************************************" << endl;
    cout << "Input parameters" << endl;
    cout << "****************************************************************************" << endl;
    cout << endl;

    // Extract input parameters.

    // Extract required parameters. 
    const string caseName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "CASE" ) );
    cout << "Case                                                      " 
         << caseName << endl; 

    const string databasePath = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "DATABASE" ) );
    cout << "Database                                                  "
         << databasePath << endl;

    const unsigned int monteCarloPopulation = extractParameterValue< unsigned int >(
                parsedData->begin( ), parsedData->end( ), 
                findEntry( dictionary, "MONTECARLOPOPULATION" ) );
    cout << "Monte Carlo population                                    "
         << monteCarloPopulation << endl;     

    const double perturberDensity = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBERDENSITY" ) );    
    cout << "Perturber density                                         "
         << perturberDensity << endl;                         

    const double perturberRingMass = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBERRINGMASS" ) );    
    cout << "Perturber ring mass                                       "
         << perturberRingMass << endl;                       

    // Extract optional parameters. 
    const int numberOfThreads = extractParameterValue< int >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "NUMBEROFTHREADS" ), 1 );
    cout << "Number of threads                                         "
         << numberOfThreads << endl;  

    const string outputMode = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "OUTPUTMODE" ), "DATABASE" );
    cout << "Output mode                                               "
         << outputMode << endl;          

    const string fileOutputDirectory = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "FILEOUTPUTDIRECTORY" ), "" );
    cout << "File output directory                                     "
         << fileOutputDirectory << endl;

//     const double observationPeriod = extractParameterValue< double >(
//                 parsedData->begin( ), parsedData->end( ),
//                 findEntry( dictionary, "OBSERVATIONPERIOD" ),
//                 convertJulianYearsToSeconds( 3.0 ), &convertJulianYearsToSeconds );

//     const double epochWindowSize = extractParameterValue< double >(
//                 parsedData->begin( ), parsedData->end( ),
//                 findEntry( dictionary, "EPOCHWINDOWSIZE" ),
//                 convertJulianDaysToSeconds( 90.0 ), &convertJulianDaysToSeconds  );

//     const unsigned int numberOfEpochWindows = extractParameterValue< unsigned int >(
//                 parsedData->begin( ), parsedData->end( ),
//                 findEntry( dictionary, "NUMBEROFEPOCHWINDOWS" ), 4 );

//     // Compute epoch window spacing.
//     const double epochWindowSpacing = observationPeriod / ( numberOfEpochWindows - 1 );

// //    // DEBUG.
// //    std::cout << "Database: " << databasePath << std::endl;
// //    std::cout << "Threads: " << numberOfThreads << std::endl;
// //    std::cout << "Monte Carlo: "<< monteCarloPopulation << std::endl;
// //    std::cout << "Perturbers: " << perturberPopulation << std::endl;
// //    std::cout << "Mass distribution type: " << massDistributionType << std::endl;
// //    std::cout << "Mass distribution parameters: " << massDistributionParameters.at( 0 );
// //    if ( massDistributionParameters.size( ) == 2 )
// //        std::cout << ", " << massDistributionParameters.at( 1 ) << std::endl;
// //    else
// //        std::cout << std::endl;
// //    std::cout << "Observation period: " << observationPeriod / JULIAN_YEAR << std::endl;
// //    std::cout << "Epoch window size: " << epochWindowSize / JULIAN_DAY << std::endl;
// //    std::cout << "Number of epoch windows: " << numberOfEpochWindows << std::endl;
// //    std::cout << "Epoch window spacing: " << epochWindowSpacing / JULIAN_YEAR << std::endl;

//     // Retrieve and store case data.
//     const CaseDataRow caseData = getCaseData( databasePath );

// //    // DEBUG.
// //    std::cout << caseData << std::endl;

// //    // Set output data precision.
// //    const int outputDataPrecision = std::numeric_limits< double >::digits10;

    // Check that all required parameters have been set.
    checkRequiredParameters( dictionary );

    ///////////////////////////////////////////////////////////////////////////

//     ///////////////////////////////////////////////////////////////////////////

//     // Query completed test particle simulations in database.

//     // Get input table from database.
//     const InputTable inputTable = getInputTable( databasePath, true );

//     // Declare vector of completed simulation numbers.
//     SimulationNumbers completedSimulationNumbers;

//     // Store list of completed simulation numbers from input table.
//     for ( unsigned int i = 0; i < inputTable.size( ); i++ )
//     {
//         completedSimulationNumbers.push_back( inputTable.at( i ).simulationNumber );
//     }

// //    // DEBUG.
// //    for ( unsigned int i = 0; i < completedSimulationNumbers.size( ); i++ )
// //        std::cout << completedSimulationNumbers.at( i ) << std::endl;

//     ///////////////////////////////////////////////////////////////////////////

//     ///////////////////////////////////////////////////////////////////////////

//     // Create required random number generators.

//     // Define a random number generator and initialize it with a reproducible
//     // seed (current cpu time).
//     mt19937 randomNumbergenerator( std::time( 0 ) );

//     // Define an uniform random number distribution for simulation numbers.
//     uniform_int_distribution< > simulationNumberDistribution(
//                 1, completedSimulationNumbers.size( ) );

//     // Define variate generator for simulation numbers using the random number
//     // generator and uniform distribution of simulation numbers.
//     variate_generator< mt19937&, uniform_int_distribution< > > generatorSimulationNumber(
//                 randomNumbergenerator, simulationNumberDistribution );

//     // Define an exponential random number distribution for power-law scaling mass factors.
//     const double dohnanyiPowerLawIndex = 1.837;
//     exponential_distribution< > massFactorPowerLawDistribution( dohnanyiPowerLawIndex );

//     // Define variate generator for mass factors using the random number
//     // generator and uniform distribution of mass factors.
//     variate_generator< mt19937&, exponential_distribution< > > generatePowerLawMassFactors(
//                 randomNumbergenerator, massFactorPowerLawDistribution );

//     // Define an uniform random number distribution to select a random observation period during
//     // random walk simulation duration.
//     uniform_real_distribution< > observationPeriodDistribution(
//                 epochWindowSize / 2.0, caseData.randomWalkSimulationDuration - observationPeriod
//                 - epochWindowSize / 2.0 );

//     // Define variate generator for observation period start epoch.
//     variate_generator< mt19937&, uniform_real_distribution< > > generatorObservationPeriodEpoch(
//                 randomNumbergenerator, observationPeriodDistribution );

//     ///////////////////////////////////////////////////////////////////////////

//     ///////////////////////////////////////////////////////////////////////////

//     // Execute Monte Carlo simulation.

//     // Loop over Monte Carlo population size.
// #pragma omp parallel for num_threads( numberOfThreads ) schedule( static, 1 )
//     for ( unsigned int monteCarloIndividual = 0; monteCarloIndividual < monteCarloPopulation;
//           monteCarloIndividual++ )
//     {

// #pragma omp critical( outputToConsole )
//         {
//             std::cout << "Simulation " << monteCarloIndividual + 1 << " / " << monteCarloPopulation
//                       << " on thread " << omp_get_thread_num( ) + 1
//                       << " / " <<  omp_get_num_threads( ) << std::endl;
//         }

//         ///////////////////////////////////////////////////////////////////////////

//         // Declare selected simulation numbers.
//         SimulationNumbers selectedSimulationNumbers;

//         // If desired perturber population is greater than number of completed simulations, throw
//         // a runtime error.
//         if ( perturberPopulation > static_cast< int >( completedSimulationNumbers.size( ) ) )
//         {
//             throw_exception(
//                         enable_error_info(
//                             runtime_error(
//                                 "Perturber population > # completed simulations" ) ) );
//         }

//         // Else, if the desired perturber population is set to -1, select all completed simulations.
//         else if ( perturberPopulation == -1 )
//         {
//             selectedSimulationNumbers = completedSimulationNumbers;
//         }

//         // Else, populate the vector of selected simulation numbers, ensuring that all values are
//         // unique.
//         else
//         {
//             // Resize the vector of selected simulation numbers.
//             selectedSimulationNumbers.resize( perturberPopulation );

//             // Generate vector of randomly selected simulation numbers.
//             generate( selectedSimulationNumbers.begin( ), selectedSimulationNumbers.end( ),
//                       generatorSimulationNumber );

//             // Check if the simulation numbers are unique, and if not, generate new random number.
//             for ( unsigned int i = 0; i < selectedSimulationNumbers.size( ); i++ )
//             {
//                 for ( unsigned int j = 0; j < selectedSimulationNumbers.size( ); j++ )
//                 {
//                     if ( i == j )
//                     {
//                         continue;
//                     }

//                     else if ( selectedSimulationNumbers.at( j )
//                               == selectedSimulationNumbers.at( i ) )
//                     {
//                         selectedSimulationNumbers.at( j ) = generatorSimulationNumber( );
//                         i = 0;
//                         break;
//                     }
//                 }
//             }
//         }

//         // Declare map of mass factors.
//         MassFactors massFactors;

//         // Loop through if-statements to select mass factor type and populate map.
//         // Check if mass factor type is "EQUAL".
//         if ( !massDistributionType.compare( "EQUAL" ) )
//         {
//             // Fill map of equal mass factors.
//             for ( unsigned int i = 0; i < selectedSimulationNumbers.size( ); i++ )
//             {
//                 massFactors[ selectedSimulationNumbers.at( i ) ]
//                         = massDistributionParameters.at( 0 );
//             }
//         }

//         // Else, check if mass factor type is "POWERLAW".
//         else if ( !massDistributionType.compare( "POWERLAW" ) )
//         {
//             // Fill map of power-law distributed mass factors.
//             for ( unsigned int i = 0; i < selectedSimulationNumbers.size( ); i++ )
//             {
//                 double powerLawMassFactor = generatePowerLawMassFactors( );

//                 while ( powerLawMassFactor > 1.0 )
//                 {
//                     powerLawMassFactor = generatePowerLawMassFactors( );
//                 }

//                 massFactors[ selectedSimulationNumbers.at( i ) ]
//                         = powerLawMassFactor / dohnanyiPowerLawIndex;
//             }
//         }

//         // Else, the mass distribution type is unknown, throw a runtime error.
//         else
//         {
//             throw_exception( enable_error_info(
//                                  runtime_error( "Mass distribution type unknown!" ) ) );
//         }

// //        // DEBUG.
// //        for ( MassFactors::const_iterator iteratorMassFactors = massFactors.begin( );
// //              iteratorMassFactors != massFactors.end( ); iteratorMassFactors++ )
// //            std::cout << iteratorMassFactors->first
// //                      << ", " << iteratorMassFactors->second << std::endl;

//         // Compile aggregate kick table to database. To avoid locking of the database, this
//         // section is thread-critical, so will be executed one-by-one by multiple threads.
//         KickTable aggregateKickTable;

// #pragma omp critical( accessDatabase )
//         {
//             aggregateKickTable = getAggregateKickTable(
//                         databasePath, caseData.randomWalkSimulationDuration,
//                         selectedSimulationNumbers, massFactors );
//         }

//         ///////////////////////////////////////////////////////////////////////////

//         ///////////////////////////////////////////////////////////////////////////

//         // Execute random walk simulation.

//         // Declare Mab propagation history. This stores the propagation history of the
//         // action variables only (semi-major axis, eccentricity, inclination).
//         ActionPropagationHistory keplerianActionElementsHistory;

//         // Set Mab initial state in Keplerian elements.
//         keplerianActionElementsHistory[ aggregateKickTable.begin( )->first ]
//                 = caseData.mabStateInKeplerianElementsAtT0.segment( 0, 3 );

//         // Declare iterator to previous state in Keplerian elements.
//         ActionPropagationHistory::iterator iteratorPreviousKeplerianActionElements
//                 = keplerianActionElementsHistory.begin( );

//         // Loop through aggregate kick table and execute kicks on Mab. Store results in propagation
//         // history.
//         for ( KickTable::iterator iteratorKickTable = aggregateKickTable.begin( );
//               iteratorKickTable != aggregateKickTable.end( ); iteratorKickTable++ )
//         {
//             // Execute kick and store results in propagation history.
//             keplerianActionElementsHistory[ iteratorKickTable->first ]
//                     = executeKick( iteratorPreviousKeplerianActionElements->second,
//                                    iteratorKickTable->second );

//             iteratorPreviousKeplerianActionElements
//                     = keplerianActionElementsHistory.find( iteratorKickTable->first );
//         }

// //        // DEBUG.
// //        // Declare stringstream for Keplerian action elements history filename.
// //        const double outputDataPrecision = std::numeric_limits< double >::digits10;
// //        std::stringstream keplerianActionElementsHistoryFilenameStream;
// //        keplerianActionElementsHistoryFilenameStream
// //                << "keplerianActionElementsHistory_case" << caseData.caseNumber
// //                << "_perturbers" << perturberPopulation << ".dat";

//         ///////////////////////////////////////////////////////////////////////////

//         ///////////////////////////////////////////////////////////////////////////

//         // Select a random observation period epoch.
//         const double observationPeriodEpoch = generatorObservationPeriodEpoch( );

// //        // DEBUG.
// //        std::cout << std::setprecision( outputDataPrecision )
// //                  << observationPeriodEpoch << std::endl;

//         // Compute maximum eccentricity change in propagation history.

//         // Declare map of epoch-window average eccentricities.
//         DoubleKeyDoubleValueMap epochWindowAverageEccentricities;

//         // Loop over observation period and compute average eccentricities per epoch window.
//         for ( unsigned int windowNumber = 0; windowNumber < numberOfEpochWindows; windowNumber++ )
//         {
//             const double epochWindowCenter = observationPeriodEpoch
//                     + windowNumber * epochWindowSpacing;

// //            // DEBUG.
// //            std::cout << std::setprecision( outputDataPrecision )
// //                      << epochWindowCenter << std::endl;

//             epochWindowAverageEccentricities[ epochWindowCenter ]
//                     = computeAverageEccentricityForEpochWindow(
//                         keplerianActionElementsHistory,
//                         epochWindowCenter - epochWindowSize / 2.0,
//                         epochWindowCenter + epochWindowSize / 2.0 );
//         }

// //        // DEBUG.
// //        for ( DoubleKeyDoubleValueMap::iterator iteratorWindow
// //              = epochWindowAverageEccentricities.begin( );
// //              iteratorWindow != epochWindowAverageEccentricities.end( );
// //              iteratorWindow++ )
// //        {
// //            std::cout << std::setprecision( outputDataPrecision )
// //                      << iteratorWindow->first << ", "
// //                      << iteratorWindow->second << std::endl;
// //        }

//         // Compute map of maximum eccentricity change during propagation history.
//         double maximumEccentricityChange = ( std::max_element(
//                     epochWindowAverageEccentricities.begin( ),
//                     epochWindowAverageEccentricities.end( ),
//                     compareDoubleKeyDoubleValueElements ) )->second
//                 - ( std::min_element(
//                         epochWindowAverageEccentricities.begin( ),
//                         epochWindowAverageEccentricities.end( ),
//                         compareDoubleKeyDoubleValueElements ) )->second;

// //        // DEBUG.
// //        std::cout << std::setprecision( outputDataPrecision )
// //                  << maximumEccentricityChange << std::endl;

//         ///////////////////////////////////////////////////////////////////////////

//         ///////////////////////////////////////////////////////////////////////////

//         // Compute maximum inclination change in propagation history.

//         // Declare map of epoch-window average inclinations.
//         DoubleKeyDoubleValueMap epochWindowAverageInclinations;

//         // Loop over observation period and compute average inclinations per epoch window.
//         for ( unsigned int windowNumber = 0; windowNumber < numberOfEpochWindows; windowNumber++ )
//         {
//             const double epochWindowCenter = observationPeriodEpoch
//                     + windowNumber * epochWindowSpacing;

//             epochWindowAverageInclinations[ epochWindowCenter ]
//                     = computeAverageInclinationForEpochWindow(
//                         keplerianActionElementsHistory,
//                         epochWindowCenter - epochWindowSize / 2.0,
//                         epochWindowCenter + epochWindowSize / 2.0 );
//         }

// //        // DEBUG.
// //        for ( DoubleKeyDoubleValueMap::iterator iteratorWindow
// //              = epochWindowAverageInclinations.begin( );
// //              iteratorWindow != epochWindowAverageInclinations.end( );
// //              iteratorWindow++ )
// //        {
// //            std::cout << std::setprecision( outputDataPrecision )
// //                      << iteratorWindow->first << ", "
// //                      << iteratorWindow->second << std::endl;
// //        }

//         // Compute map of maximum inclination change during propagation history.
//         double maximumInclinationChange = ( std::max_element(
//                     epochWindowAverageInclinations.begin( ),
//                     epochWindowAverageInclinations.end( ),
//                     compareDoubleKeyDoubleValueElements ) )->second
//                 - ( std::min_element(
//                         epochWindowAverageInclinations.begin( ),
//                         epochWindowAverageInclinations.end( ),
//                         compareDoubleKeyDoubleValueElements ) )->second;

// //        // DEBUG.
// //        std::cout << std::setprecision( outputDataPrecision )
// //                  << maximumInclinationChange << std::endl;

//           ///////////////////////////////////////////////////////////////////////////

//         ///////////////////////////////////////////////////////////////////////////

//         // Compute maximum longitude residual change in propagation history.

//         // Compute longitude history.
//         DoubleKeyDoubleValueMap longitudeHistory
//                 = computeLongitudeHistory( keplerianActionElementsHistory,
//                                            caseData.uranusGravitationalParameter );

//         // Compute reduced longitude history.
//         DoubleKeyDoubleValueMap reducedLongitudeHistory
//                 = reduceLongitudeHistory( longitudeHistory, observationPeriodEpoch,
//                                           epochWindowSpacing, epochWindowSize,
//                                           numberOfEpochWindows );

// //        // DEBUG.
// //        // Write test particle (reduced) longitude history to file.
// //        writeDataMapToTextFile( longitudeHistory,
// //                                "longitudeHistory.dat",
// //                                "/Users/kartikkumar/Desktop", "",
// //                                outputDataPrecision,
// //                                outputDataPrecision, "," );
// //        writeDataMapToTextFile( reducedLongitudeHistory,
// //                                "reducedLongitudeHistory.dat",
// //                                "/Users/kartikkumar/Desktop", "",
// //                                outputDataPrecision,
// //                                outputDataPrecision, "," );

//         // Set input data for simple linear regression.
//         SimpleLinearRegression longitudeHistoryRegression( reducedLongitudeHistory );

//         // Compute linear fit.
//         longitudeHistoryRegression.computeFit( );

//         // Clear longitude residuals history.
//         DoubleKeyDoubleValueMap longitudeResidualsHistory;

//         // Generate longitude history residuals by subtracting linear fit from data.
//         for ( DoubleKeyDoubleValueMap::iterator iteratorReducedLongitudeHistory
//               = reducedLongitudeHistory.begin( );
//               iteratorReducedLongitudeHistory != reducedLongitudeHistory.end( );
//               iteratorReducedLongitudeHistory++ )
//         {
//             longitudeResidualsHistory[ iteratorReducedLongitudeHistory->first ]
//                     = iteratorReducedLongitudeHistory->second
//                     - longitudeHistoryRegression.getCoefficientOfConstantTerm( )
//                     - longitudeHistoryRegression.getCoefficientOfLinearTerm( )
//                     * iteratorReducedLongitudeHistory->first;
//         }

// //        // DEBUG.
// //        writeDataMapToTextFile( longitudeResidualsHistory,
// //                                "longitudeResidualsHistory.dat",
// //                                "/Users/kartikkumar/Desktop", "",
// //                                outputDataPrecision,
// //                                outputDataPrecision, "," );

//         // Declare map of epoch-window average longitude residuals.
//         DoubleKeyDoubleValueMap epochWindowAverageLongitudeResiduals;

//         for ( unsigned int windowNumber = 0; windowNumber < numberOfEpochWindows; windowNumber++ )
//         {
//             const double epochWindowCenter = epochWindowSize / 2.0
//                     + windowNumber * epochWindowSpacing;

//             epochWindowAverageLongitudeResiduals[ epochWindowCenter ]
//                     = computeAverageLongitudeResidualForEpochWindow(
//                         longitudeResidualsHistory,
//                         epochWindowCenter - epochWindowSize / 2.0,
//                         epochWindowCenter + epochWindowSize / 2.0 );
//         }

// //        // DEBUG.
// //        for ( DoubleKeyDoubleValueMap::iterator iteratorWindow
// //              = epochWindowAverageLongitudeResiduals.begin( );
// //              iteratorWindow != epochWindowAverageLongitudeResiduals.end( );
// //              iteratorWindow++ )
// //        {
// //            std::cout << std::setprecision( outputDataPrecision )
// //                      << iteratorWindow->first << ", "
// //                      << iteratorWindow->second << std::endl;
// //        }

//         // Compute map of maximum longitude residual change during propagation history.
//         double maximumLongitudeResidualChange = ( std::max_element(
//                     epochWindowAverageLongitudeResiduals.begin( ),
//                     epochWindowAverageLongitudeResiduals.end( ),
//                     compareDoubleKeyDoubleValueElements ) )->second
//                 - ( std::min_element(
//                         epochWindowAverageLongitudeResiduals.begin( ),
//                         epochWindowAverageLongitudeResiduals.end( ),
//                         compareDoubleKeyDoubleValueElements ) )->second;

// //        // DEBUG.
// //        std::cout << std::setprecision( outputDataPrecision )
// //                  << maximumLongitudeResidualChange << std::endl;

//         ///////////////////////////////////////////////////////////////////////////

//         ///////////////////////////////////////////////////////////////////////////

//         // Write random walk table to database. To avoid locking of the database, this section is
//         // thread-critical, so will be executed one-by-one by multiple threads.
// #pragma omp critical( accessDatabase )
//         {
//             // Populate kick table in database.
//             populateRandomWalkTable( databasePath, selectedSimulationNumbers, massFactors,
//                                      massDistributionType, massDistributionParameters,
//                                      observationPeriod, epochWindowSize, numberOfEpochWindows,
//                                      maximumEccentricityChange, maximumLongitudeResidualChange,
//                                      maximumInclinationChange );
//         }

//         ///////////////////////////////////////////////////////////////////
//     }

//     ///////////////////////////////////////////////////////////////////////////

    return 0;
}
