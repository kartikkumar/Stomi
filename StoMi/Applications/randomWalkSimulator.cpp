/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <sstream>

#include <omp.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/make_shared.hpp>

#include <Assist/Astrodynamics/astrodynamicsBasics.h>
#include <Assist/Astrodynamics/hillSphere.h>
#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/commonTypedefs.h>
#include <Assist/Basics/comparisonFunctions.h>
#include <Assist/InputOutput/basicInputOutput.h>
#include <Assist/Mathematics/statistics.h>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/InputOutput/dictionaryTools.h>
#include <Tudat/InputOutput/fieldType.h>
#include <Tudat/InputOutput/separatedParser.h>
#include <Tudat/InputOutput/parsedDataVectorUtilities.h>
#include <Tudat/Mathematics/Statistics/simpleLinearRegression.h>

#include "StoMi/Astrodynamics/randomWalkFunctions.h"
#include "StoMi/Database/databaseReadFunctions.h" 
#include "StoMi/Database/databaseWriteFunctions.h"  
#include "StoMi/Database/testParticleCase.h"
#include "StoMi/Database/testParticleKick.h"
#include "StoMi/Database/randomWalkCase.h"
#include "StoMi/Database/randomWalkInput.h"
#include "StoMi/InputOutput/dictionaries.h"


#include <SQLiteCpp/SQLiteCpp.h> 

//! Execute random walk simulations.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    using std::advance;
    using std::cout;
    using std::endl;
    using std::ostringstream;
    using std::string;

    using boost::iequals;
    using namespace boost::filesystem; 
    using boost::make_shared;   

    using namespace SQLite;

    using namespace assist::astrodynamics;
    using namespace assist::basics;
    using namespace assist::input_output;
    using namespace assist::mathematics;

    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_mathematics::mathematical_constants;
    using namespace tudat::input_output;
    using namespace tudat::input_output::dictionary;
    using namespace tudat::input_output::field_types::general;
    using namespace tudat::input_output::parsed_data_vector_utilities;
    using namespace tudat::statistics;

    using namespace stomi::astrodynamics;
    using namespace stomi::database;
    using namespace stomi::input_output;

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
    const string databasePath = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "DATABASE" ) );
    cout << "Database                                                  " 
         << databasePath << endl;

    const string randomWalkCaseName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "CASE" ) );
    cout << "Random walk case                                          "  
         << randomWalkCaseName << endl;  

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
                findEntry( dictionary, "FILEOUTPUTDIRECTORY" ), "" ) + "/";
    cout << "File output directory                                     "
         << fileOutputDirectory << endl;

    const string monteCarloRunsToExecute = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "MONTECARLORUNSTOEXECUTE" ), "ALL" );
    cout << "Monte Carlo runs to execute                               "
         << monteCarloRunsToExecute << endl;         

    const string randomWalkCaseTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKCASETABLENAME" ), "random_walk_case" );
    cout << "Random walk case table                                    " 
         << randomWalkCaseTableName << endl;

    const string randomWalkInputTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKINPUTTABLENAME" ), "random_walk_input" );
    cout << "Random walk input table                                   " 
         << randomWalkInputTableName << endl;         

    const string randomWalkPerturberTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKPERTURBERTABLENAME" ),
                "random_walk_perturbers" );
    cout << "Random walk perturber table                               "
         << randomWalkPerturberTableName << endl;

    const string randomWalkOutputTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKOUTPUTTABLENAME" ), "random_walk_output" );
    cout << "Random walk output table                                  "
              << randomWalkOutputTableName << endl; 

    const string testParticleCaseTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLECASETABLENAME" ), "test_particle_case" );
    cout << "Test particle case table                                  "
         << testParticleCaseTableName << endl;                

    const string testParticleKickTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLEKICKTABLENAME" ), "test_particle_kicks" );
    cout << "Test particle kick table                                  "
         << testParticleKickTableName << endl;          

    // Retrieve and store random walk case data from database.
    RandomWalkCasePointer randomWalkCase;

    // Case data is extracted in local scope from the database, with overwritten parameters 
    // extracted from the input file, to ensure that none of these parameters are used globally
    // elsewhere in this file.
    {
        const RandomWalkCasePointer caseDataFromDatabase = getRandomWalkCase( 
            databasePath, randomWalkCaseName, randomWalkCaseTableName ); 

        const double perturberDensity = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "PERTURBERDENSITY" ),
                    caseDataFromDatabase->perturberDensity );    
        cout << "Perturber density                                         "
             << perturberDensity << " perturbers per R_Hill" << endl;                         

        const double perturberRingMass = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "PERTURBERRINGMASS" ),
                    caseDataFromDatabase->perturberRingMass );    
        cout << "Perturber ring mass                                       "
             << perturberRingMass << " M_PerturbedBody" << endl;              

        const double observationPeriod = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "OBSERVATIONPERIOD" ),
                    caseDataFromDatabase->observationPeriod, &convertJulianYearsToSeconds );
        cout << "Observation period                                        "
             << convertSecondsToJulianYears( observationPeriod ) << " yrs" << endl; 

        const unsigned int numberOfEpochWindows = extractParameterValue< unsigned int >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "NUMBEROFEPOCHWINDOWS" ),
                    caseDataFromDatabase->numberOfEpochWindows );    
        cout << "Number of epoch windows                                   "
             << numberOfEpochWindows << endl;                                                 

        const double epochWindowSize = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "EPOCHWINDOWSIZE" ),
                    caseDataFromDatabase->epochWindowSize, &convertJulianDaysToSeconds  );
        cout << "Epoch window size                                         "
             << convertSecondsToJulianDays( epochWindowSize ) << " days" << endl;     

        // Store case data with possible overwritten data.
        randomWalkCase = make_shared< RandomWalkCase >( 
            caseDataFromDatabase->caseId, randomWalkCaseName, 
            caseDataFromDatabase->testParticleCaseId, perturberDensity, perturberRingMass,
            observationPeriod, numberOfEpochWindows, epochWindowSize );
    }

    // Check that all required parameters have been set.
    checkRequiredParameters( dictionary );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Fetch case and input data from database.

    cout << endl;
    cout << "****************************************************************************" << endl;
    cout << "Database operations" << endl;
    cout << "****************************************************************************" << endl;
    cout << endl;

    // Generate output message.
    cout << "Fetching test particle case data from database ..." << endl;      

    // Retrieve and store case data.
    const TestParticleCasePointer testParticleCaseData = getTestParticleCase( 
        databasePath, randomWalkCase->testParticleCaseId, testParticleCaseTableName ); 

    // Generate output message to indicate that case data was fetched successfully.
    cout << "Test particle case data fetched successfully from database!" << endl;   

    // Generate output message.
    cout << "Fetching random walk input data from database ..." << endl;   

    // Check if all incomplete Monte Carlo runs are to be run and fetch the input table, else only 
    // fetch the requested Monte Carlo run IDs.
    RandomWalkInputTable randomWalkInputTable;

    if ( iequals( monteCarloRunsToExecute, "ALL" )  )
    {
        cout << "Fetching all incomplete random walk Monte Carlo runs ..." << endl;    

        // Get entire random walk input table from database.
        randomWalkInputTable = getCompleteRandomWalkInputTable(
                    databasePath, randomWalkCase->caseId, 
                    randomWalkInputTableName, randomWalkPerturberTableName );
    }

    else
    {
        cout << "Fetching all requested random walk Monte Carlo runs ..." << endl;    

        // Get selected random walk input table from database.
        randomWalkInputTable = getSelectedRandomWalkInputTable(
                    databasePath, randomWalkCase->caseId, monteCarloRunsToExecute, 
                    randomWalkInputTableName, randomWalkPerturberTableName );
    }

    // Generate output message to indicate that the input table was fetched successfully.
    cout << "Random walk input data (" << randomWalkInputTable.size( )
         << " rows) fetched successfully from database!" << endl;    

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Compute derived parameters.

    cout << endl;
    cout << "****************************************************************************" << endl;
    cout << "Derived parameters" << endl;
    cout << "****************************************************************************" << endl;
    cout << endl;

    // Compute epoch window spacing [s].
    const double epochWindowSpacing 
        = randomWalkCase->observationPeriod / ( randomWalkCase->numberOfEpochWindows - 1 );
    cout << "Epoch window spacing                                      " 
         << convertSecondsToJulianDays( epochWindowSpacing ) << " days" << endl;

    // Compute mass of perturbed body [kg].
    const double perturbedBodyMass = computeMassOfSphere(
                testParticleCaseData->perturbedBodyRadius, 
                testParticleCaseData->perturbedBodyBulkDensity );
    cout << "Perturbed body mass                                       " 
         << perturbedBodyMass << " kg" << endl;

    // Compute perturbed body's gravitational parameter [m^3 s^-2].
    const double perturbedBodyGravitationalParameter
            = computeGravitationalParameter( perturbedBodyMass );
    cout << "Perturbed body gravitational parameter                    " 
         << perturbedBodyGravitationalParameter << " m^3 s^-2" << endl;
         
    // Compute perturber population using density and semi-major axis distribution limits.
    // Note, in the floor() function, adding 0.5 is a workaround for the fact that there is no
    // round() function in C++03 (it is available in C++11).
    ConvertHillRadiiToMeters convertHillRadiiToMeters( 
        testParticleCaseData->centralBodyGravitationalParameter, 
        perturbedBodyGravitationalParameter, 
        testParticleCaseData->perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) );
    const double perturberDensityInMeters 
        = randomWalkCase->perturberDensity / convertHillRadiiToMeters( 1.0 );
    const unsigned int perturberPopulation = std::floor( 
        2.0 * testParticleCaseData->semiMajorAxisDistributionLimit 
        * perturberDensityInMeters + 0.5 );
    cout << "Perturber population                                      " 
         << perturberPopulation << endl;

    // Compute perturber mass ratio.
    // Note, for the random walk simulations, the mass ratio is equal for all perturbers.
    const double perturberMassRatio = randomWalkCase->perturberRingMass / perturberPopulation;
    cout << "Perturber mass ratio                                      "
         << perturberMassRatio << endl;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Execute Monte Carlo simulation.
    cout << endl;
    cout << "****************************************************************************" << endl;
    cout << "Simulation loop" << endl;
    cout << "****************************************************************************" << endl;
    cout << endl;

    // Execute simulation loop.
    cout << "Starting simulation loop ... " << endl;
    cout << randomWalkInputTable.size( ) << " Monte Carlo simulations queued for execution ..." 
         << endl;
    cout << endl;    

    // Loop over input table.
#pragma omp parallel for num_threads( numberOfThreads )
    for ( unsigned int i = 0; i < randomWalkInputTable.size( ); i++ )
    {
        ///////////////////////////////////////////////////////////////////////////

        // Set input table iterator and emit output message.

        // Set input table iterator for current simulation wrt to start of input table and counter.
        RandomWalkInputTable::iterator iteratorRandomWalkInputTable 
            = randomWalkInputTable.begin( );
        advance( iteratorRandomWalkInputTable, i );

        // Emit output message.
#pragma omp critical( outputToConsole )
        {
            cout << "Monte Carlo run " << iteratorRandomWalkInputTable->monteCarloRunId 
                 << " on thread " << omp_get_thread_num( ) + 1 << " / " 
                 << omp_get_num_threads( ) << endl;
        }

        ///////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////

        // Fetch kick table based on test particle simulation IDs for Monte Carlo run.
        TestParticleKickTable kickTable;
#pragma omp critical( accessKickTable )
        {        
            kickTable = getTestParticleKickTable( 
                databasePath, testParticleCaseData->randomWalkSimulationPeriod,
                iteratorRandomWalkInputTable->testParticleSimulationIds, 
                testParticleKickTableName );
        }

        // Check if output mode is set to "FILE".
        // If so, open output file and write kick table data.
        // Check if the output directory exists: if not, create it.
        if ( iequals( outputMode, "FILE" ) )
        {
            // Check if output directory exists.
            if ( !exists( fileOutputDirectory ) )
            {
               std::cerr << "Directory does not exist. Will be created." << std::endl;
               create_directories( fileOutputDirectory );
            }        

            // Declare file handler.
            std::ofstream kickTableFile;

            // Set up and write file header to file.
            std::ostringstream kickTableFilename;
            kickTableFilename << fileOutputDirectory << "monteCarloRun" 
                              << iteratorRandomWalkInputTable->monteCarloRunId 
                              << "_kickTable.csv";    

            kickTableFile.open( kickTableFilename.str( ).c_str( ) );

            kickTableFile << "kickId,simulationId,conjunctionEpoch,conjunctionDistance,"
                          << "preConjunctionEpoch,preConjunctionDistance,"
                          << "preConjunctionSemiMajorAxis,preConjunctionEccentricity,"
                          << "preConjunctionInclination,preConjunctionArgumentOfPeriapsis,"
                          << "preConjunctionLongitudeOfAscendingNode,preConjunctionTrueAnomaly"
                          << "postConjunctionEpoch,postConjunctionDistance,"
                          << "postConjunctionSemiMajorAxis,postConjunctionEccentricity,"
                          << "postConjunctionInclination,postConjunctionArgumentOfPeriapsis,"
                          << "postConjunctionLongitudeOfAscendingNode,postConjunctionTrueAnomaly"                          
                          << endl;
            kickTableFile << "# [-],[-],[s],[m],[s],[m],[m],[-],[rad],[rad],[rad],[rad]," 
                          << "[s],[m],[m],[-],[rad],[rad],[rad],[rad]" << endl;       

            // Write kick table to file.
            for ( TestParticleKickTable::iterator iteratorTestParticleKicks = kickTable.begin( );
                  iteratorTestParticleKicks != kickTable.end( );
                  iteratorTestParticleKicks++ )
            {
                kickTableFile 
                    << iteratorTestParticleKicks->kickId << "," 
                    << iteratorTestParticleKicks->testParticleSimulationId << ",";
                kickTableFile
                    << std::setprecision( std::numeric_limits< double >::digits10 )
                    << iteratorTestParticleKicks->conjunctionEpoch << ","
                    << iteratorTestParticleKicks->conjunctionDistance << ","
                    << iteratorTestParticleKicks->preConjunctionEpoch << ","
                    << iteratorTestParticleKicks->preConjunctionDistance << ","
                    << iteratorTestParticleKicks->preConjunctionStateInKeplerianElements( 
                        semiMajorAxisIndex ) << ","
                    << iteratorTestParticleKicks->preConjunctionStateInKeplerianElements( 
                        eccentricityIndex ) << ","
                    << iteratorTestParticleKicks->preConjunctionStateInKeplerianElements( 
                        inclinationIndex ) << ","
                    << iteratorTestParticleKicks->preConjunctionStateInKeplerianElements( 
                        argumentOfPeriapsisIndex ) << ","
                    << iteratorTestParticleKicks->preConjunctionStateInKeplerianElements( 
                        longitudeOfAscendingNodeIndex ) << ","
                    << iteratorTestParticleKicks->preConjunctionStateInKeplerianElements( 
                        trueAnomalyIndex ) << ","
                    << iteratorTestParticleKicks->postConjunctionEpoch << ","
                    << iteratorTestParticleKicks->postConjunctionDistance << ","
                    << iteratorTestParticleKicks->postConjunctionStateInKeplerianElements( 
                        semiMajorAxisIndex ) << ","
                    << iteratorTestParticleKicks->postConjunctionStateInKeplerianElements( 
                        eccentricityIndex ) << ","
                    << iteratorTestParticleKicks->postConjunctionStateInKeplerianElements( 
                        inclinationIndex ) << ","
                    << iteratorTestParticleKicks->postConjunctionStateInKeplerianElements( 
                        argumentOfPeriapsisIndex ) << ","
                    << iteratorTestParticleKicks->postConjunctionStateInKeplerianElements( 
                        longitudeOfAscendingNodeIndex ) << ","
                    << iteratorTestParticleKicks->postConjunctionStateInKeplerianElements( 
                        trueAnomalyIndex ) << endl;
            }                 

            // Close file handler.
            kickTableFile.close( );
        }

        ///////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////

        // Execute random walk simulation.

        // Declare perturbed body propagation history. This stores the propagation history of the
        // action variables only (semi-major axis, eccentricity, inclination).
        DoubleKeyVector3dValueMap keplerianActionElementsHistory;

        // Set perturbed body initial state (actions) in Keplerian elements.
        keplerianActionElementsHistory[ 0.0 ]
                = testParticleCaseData->perturbedBodyStateInKeplerianElementsAtT0.segment( 0, 3 );

        // Declare iterator to previous state in Keplerian elements.
        DoubleKeyVector3dValueMap::iterator iteratorPreviousKeplerianElements
                = keplerianActionElementsHistory.begin( );

        // Loop through aggregate kick table and execute kicks on perturbed body. 
        for ( TestParticleKickTable::iterator iteratorKickTable = kickTable.begin( );
              iteratorKickTable != kickTable.end( ); iteratorKickTable++ )
        {
            // Execute kick and store results in propagation history.
            keplerianActionElementsHistory[ iteratorKickTable->conjunctionEpoch ]
                    = executeKick( iteratorPreviousKeplerianElements->second,
                                   iteratorKickTable, perturberMassRatio );

            advance( iteratorPreviousKeplerianElements, 1 );
        }

        // Check if output mode is set to "FILE".
        // If so, open output file and write header content.
        if ( iequals( outputMode, "FILE" ) )
        {
            ostringstream keplerianActionElementsFilename;
            keplerianActionElementsFilename << "monteCarloRun" 
                                            << iteratorRandomWalkInputTable->monteCarloRunId
                                            << "_keplerianActionElements.csv"; 
            
            ostringstream keplerianActionElementsFileHeader;
            keplerianActionElementsFileHeader << "epoch,semiMajorAxis,eccentricity,inclination" 
                                              << endl;
            keplerianActionElementsFileHeader << "# [s],[m],[-],[rad]" << endl;                                            

            writeDataMapToTextFile( keplerianActionElementsHistory,
                                    keplerianActionElementsFilename.str( ),
                                    fileOutputDirectory, 
                                    keplerianActionElementsFileHeader.str( ),
                                    std::numeric_limits< double >::digits10, 
                                    std::numeric_limits< double >::digits10,
                                    "," ); 
        }        

        ///////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////

        // Compute average longitude residual and maximum longitude residual change in observation
        // period.

        // Declare average longitude residual in observation period [-].
        double averageLongitudeResidual = TUDAT_NAN;

        // Declare maximum longitude residual change in observation period [-].
        double maximumLongitudeResidualChange = TUDAT_NAN;

        // Declare map of average longitude residuals per window [rad].
        DoubleKeyDoubleValueMap averageLongitudeResiduals;

        {
            // Populate temporary map with epochs and semi-major axes.
            DoubleKeyDoubleValueMap semiMajorAxisHistory;

            for ( DoubleKeyVector3dValueMap::iterator iteratorKeplerianActionElements 
                  = keplerianActionElementsHistory.begin( );
                  iteratorKeplerianActionElements != keplerianActionElementsHistory.end( );
                  iteratorKeplerianActionElements++ )
            {
                semiMajorAxisHistory[ iteratorKeplerianActionElements->first ]
                    = iteratorKeplerianActionElements->second( semiMajorAxisIndex );
            }  

            // Compute longitude history.
            DoubleKeyDoubleValueMap longitudeHistory
                    = computeLongitudeHistory( 
                        semiMajorAxisHistory,
                        testParticleCaseData->centralBodyGravitationalParameter );      

            // Compute reduced longitude history (data is reduced to only the epoch windows).
            DoubleKeyDoubleValueMap reducedLongitudeHistory
                    = reduceLongitudeHistory( 
                        longitudeHistory, 
                        iteratorRandomWalkInputTable->observationPeriodStartEpoch,
                        epochWindowSpacing, 
                        randomWalkCase->epochWindowSize,
                        randomWalkCase->numberOfEpochWindows ); 

            // Set input data for simple linear regression.
            SimpleLinearRegression longitudeHistoryRegression( reducedLongitudeHistory );
         
            // Compute linear fit.
            longitudeHistoryRegression.computeFit( );        
                
            // Store longitude residuals history.
            DoubleKeyDoubleValueMap longitudeResidualsHistory;

            // Generate longitude history residuals by subtracting linear fit from data.
            for ( DoubleKeyDoubleValueMap::iterator iteratorReducedLongitudeHistory
                    = reducedLongitudeHistory.begin( );
                  iteratorReducedLongitudeHistory != reducedLongitudeHistory.end( );
                  iteratorReducedLongitudeHistory++ )
            {
                longitudeResidualsHistory[ iteratorReducedLongitudeHistory->first ]
                        = iteratorReducedLongitudeHistory->second
                          - longitudeHistoryRegression.getCoefficientOfConstantTerm( )
                          - longitudeHistoryRegression.getCoefficientOfLinearTerm( )
                          * iteratorReducedLongitudeHistory->first;
            }           

            // Loop over observation period and compute average longitude residuals per epoch window.
            for ( int windowNumber = 0; 
                  windowNumber < randomWalkCase->numberOfEpochWindows; 
                  windowNumber++ )
            {
                const double epochWindowCenter 
                    = iteratorRandomWalkInputTable->observationPeriodStartEpoch 
                      + windowNumber * epochWindowSpacing;

                averageLongitudeResiduals[ epochWindowCenter ] 
                    = computeStepFunctionWindowAverage( 
                        longitudeResidualsHistory, 
                        epochWindowCenter - 0.5 * randomWalkCase->epochWindowSize, 
                        epochWindowCenter + 0.5 * randomWalkCase->epochWindowSize );
            }

             // Compute average longitude residual during propagation history.
            double sumLongitudeResiduals = 0.0;

            for ( DoubleKeyDoubleValueMap::iterator iteratorAverageLongitudeResiduals 
                    = averageLongitudeResiduals.begin( );
                  iteratorAverageLongitudeResiduals != averageLongitudeResiduals.end( );
                  iteratorAverageLongitudeResiduals++ )
            {
                sumLongitudeResiduals += iteratorAverageLongitudeResiduals->second;
            }

            averageLongitudeResidual 
                = sumLongitudeResiduals / randomWalkCase->numberOfEpochWindows;

            // Compute maximum longitude residual change during propagation history.
            maximumLongitudeResidualChange 
                = ( std::max_element( averageLongitudeResiduals.begin( ), 
                                      averageLongitudeResiduals.end( ),
                                      CompareDoubleKeyDoubleValueMapValues( ) ) )->second
                  - ( std::min_element( averageLongitudeResiduals.begin( ), 
                                        averageLongitudeResiduals.end( ),
                                        CompareDoubleKeyDoubleValueMapValues( ) ) )->second;            

            // Check if output mode is set to "FILE".
            // If so, open output file and write header content.
            if ( iequals( outputMode, "FILE" ) )
            {
                ostringstream longitudeHistoryFilename;
                longitudeHistoryFilename << "monteCarloRun"
                                         << iteratorRandomWalkInputTable->monteCarloRunId
                                         << "_longitudeHistory.csv"; 
                
                ostringstream longitudeHistoryFileHeader;
                longitudeHistoryFileHeader << "epoch,longitude" << endl;
                longitudeHistoryFileHeader << "# [s],[rad]" << endl;                                            

                writeDataMapToTextFile( longitudeHistory,
                                        longitudeHistoryFilename.str( ),
                                        fileOutputDirectory, 
                                        longitudeHistoryFileHeader.str( ),
                                        std::numeric_limits< double >::digits10, 
                                        std::numeric_limits< double >::digits10,
                                        "," );             

                ostringstream reducedLongitudeHistoryFilename;
                reducedLongitudeHistoryFilename << "monteCarloRun"
                                                << iteratorRandomWalkInputTable->monteCarloRunId
                                                << "_reducedLongitudeHistory.csv"; 
                
                ostringstream reducedLongitudeHistoryFileHeader;
                reducedLongitudeHistoryFileHeader << "epoch,longitude" << endl;
                reducedLongitudeHistoryFileHeader << "# [s],[rad]" << endl;                                            

                writeDataMapToTextFile( reducedLongitudeHistory,
                                        reducedLongitudeHistoryFilename.str( ),
                                        fileOutputDirectory, 
                                        reducedLongitudeHistoryFileHeader.str( ),
                                        std::numeric_limits< double >::digits10, 
                                        std::numeric_limits< double >::digits10,
                                        "," );            

                ostringstream longitudeResidualsFilename;
                longitudeResidualsFilename << "monteCarloRun"
                                           << iteratorRandomWalkInputTable->monteCarloRunId
                                           << "_longitudeResiduals.csv"; 
                
                ostringstream longitudeResidualsFileHeader;
                longitudeResidualsFileHeader << "epoch,longitudeResidual" << endl;
                longitudeResidualsFileHeader << "# [s],[rad]" << endl;                                            

                writeDataMapToTextFile( longitudeResidualsHistory,
                                        longitudeResidualsFilename.str( ),
                                        fileOutputDirectory, 
                                        longitudeResidualsFileHeader.str( ),
                                        std::numeric_limits< double >::digits10, 
                                        std::numeric_limits< double >::digits10,
                                        "," );                                                 
            }                       
        }

        ///////////////////////////////////////////////////////////////////////////        

        ///////////////////////////////////////////////////////////////////////////

        // Compute average eccentricity and maximum eccentricity change in observation window.

        // Declare average eccentricity in observation period [-].
        double averageEccentricity = TUDAT_NAN;

        // Declare maximum eccentricity change in observation period [-].
        double maximumEccentricityChange = TUDAT_NAN;

        // Declare map of average eccentricities per window [-].
        DoubleKeyDoubleValueMap averageEccentricities;

        {
            // Populate temporary map with epochs and eccentricities.
            DoubleKeyDoubleValueMap eccentricityHistory;

            for ( DoubleKeyVector3dValueMap::iterator iteratorKeplerianActionElements 
                  = keplerianActionElementsHistory.begin( );
                  iteratorKeplerianActionElements != keplerianActionElementsHistory.end( );
                  iteratorKeplerianActionElements++ )
            {
                eccentricityHistory[ iteratorKeplerianActionElements->first ]
                    = iteratorKeplerianActionElements->second( eccentricityIndex );
            }            

            // Loop over observation period and compute average eccentricities per epoch window.
            for ( int windowNumber = 0; 
                  windowNumber < randomWalkCase->numberOfEpochWindows; 
                  windowNumber++ )
            {
                const double epochWindowCenter 
                    = iteratorRandomWalkInputTable->observationPeriodStartEpoch
                      + windowNumber * epochWindowSpacing;

                averageEccentricities[ epochWindowCenter ] 
                    = computeStepFunctionWindowAverage( 
                        eccentricityHistory, 
                        epochWindowCenter - 0.5 * randomWalkCase->epochWindowSize, 
                        epochWindowCenter + 0.5 * randomWalkCase->epochWindowSize );
            }

            // Compute average eccentricity during propagation history.
            double sumEccentricities = 0.0;

            for ( DoubleKeyDoubleValueMap::iterator iteratorAverageEccentricities 
                    = averageEccentricities.begin( );
                  iteratorAverageEccentricities != averageEccentricities.end( );
                  iteratorAverageEccentricities++ )
            {
                sumEccentricities += iteratorAverageEccentricities->second;
            }

            averageEccentricity = sumEccentricities / randomWalkCase->numberOfEpochWindows;

            // Compute maximum eccentricity change during propagation history.
            maximumEccentricityChange 
                = ( std::max_element( averageEccentricities.begin( ), 
                                      averageEccentricities.end( ),
                                      CompareDoubleKeyDoubleValueMapValues( ) ) )->second
                  - ( std::min_element( averageEccentricities.begin( ), 
                                        averageEccentricities.end( ),
                                        CompareDoubleKeyDoubleValueMapValues( ) ) )->second;                                                
        }

        ///////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////

        // Compute average inclination and maximum inclination change in observation window.

        // Declare average inclination in observation period [rad].
        double averageInclination = TUDAT_NAN;

        // Declare maximum inclination change in observation period [rad].
        double maximumInclinationChange = TUDAT_NAN;

        // Declare map of average inclinations per window [rad].
        DoubleKeyDoubleValueMap averageInclinations;        

        {
            // Populate temporary map with epochs and inclinations.
            DoubleKeyDoubleValueMap inclinationHistory;

            for ( DoubleKeyVector3dValueMap::iterator iteratorKeplerianActionElements 
                  = keplerianActionElementsHistory.begin( );
                  iteratorKeplerianActionElements != keplerianActionElementsHistory.end( );
                  iteratorKeplerianActionElements++ )
            {
                inclinationHistory[ iteratorKeplerianActionElements->first ]
                    = iteratorKeplerianActionElements->second( inclinationIndex );
            }            

            // Loop over observation period and compute average inclinations per epoch window.
            for ( int windowNumber = 0; 
                  windowNumber < randomWalkCase->numberOfEpochWindows; 
                  windowNumber++ )
            {
                const double epochWindowCenter 
                    = iteratorRandomWalkInputTable->observationPeriodStartEpoch
                      + windowNumber * epochWindowSpacing;

                averageInclinations[ epochWindowCenter ] 
                    = computeStepFunctionWindowAverage( 
                        inclinationHistory, 
                        epochWindowCenter - randomWalkCase->epochWindowSize * 0.5, 
                        epochWindowCenter + randomWalkCase->epochWindowSize * 0.5 );
            }

            // Compute average inclination during propagation history.
            double sumInclinations = 0.0;

            for ( DoubleKeyDoubleValueMap::iterator iteratorAverageInclinations 
                    = averageInclinations.begin( );
                  iteratorAverageInclinations != averageInclinations.end( );
                  iteratorAverageInclinations++ )
            {
                sumInclinations += iteratorAverageInclinations->second;
            }

            averageInclination = sumInclinations / randomWalkCase->numberOfEpochWindows;

            // Compute maximum inclination change during propagation history.
            maximumInclinationChange 
                = ( std::max_element( averageInclinations.begin( ), 
                                      averageInclinations.end( ),
                                      CompareDoubleKeyDoubleValueMapValues( ) ) )->second
                  - ( std::min_element( averageInclinations.begin( ), 
                                        averageInclinations.end( ),
                                        CompareDoubleKeyDoubleValueMapValues( ) ) )->second;        
        }

        ///////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////

        // Write epoch-windoow average values to file.

        // Check if output mode is set to "FILE".
        // If so, open output file and write header content.
        if ( iequals( outputMode, "FILE" ) )
        {
            // Declare file handler.
            std::ofstream epochWindowAveragesFile;   

            ostringstream epochWindowAveragesFilename;
            epochWindowAveragesFilename << fileOutputDirectory << "monteCarloRun" 
                                        << iteratorRandomWalkInputTable->monteCarloRunId
                                        << "_epochWindowAverages.csv";
            
            epochWindowAveragesFile.open( epochWindowAveragesFilename.str( ).c_str( ) );                                        

            epochWindowAveragesFile << "epoch,longitudeResidual,eccentricity,inclination" << endl;
            epochWindowAveragesFile << "# [s],[rad],[-],[rad]" << endl; 

            // Loop through epoch-window averages and write data to file.
            DoubleKeyDoubleValueMap::iterator iteratorEpochWindowAverageLongitudeResiduals
                = averageLongitudeResiduals.begin( );
            DoubleKeyDoubleValueMap::iterator iteratorEpochWindowAverageEccentricities
                = averageEccentricities.begin( );
            DoubleKeyDoubleValueMap::iterator iteratorEpochWindowAverageInclinations
                = averageInclinations.begin( );
                                                
            for ( int i = 0; i < randomWalkCase->numberOfEpochWindows; i++ )
            {
                epochWindowAveragesFile 
                    << std::setprecision( std::numeric_limits< double >::digits10 )
                    << iteratorEpochWindowAverageLongitudeResiduals->first << ","
                    << iteratorEpochWindowAverageLongitudeResiduals->second << ","
                    << iteratorEpochWindowAverageEccentricities->second << ","
                    << iteratorEpochWindowAverageInclinations->second << endl;

                iteratorEpochWindowAverageLongitudeResiduals++;    
                iteratorEpochWindowAverageEccentricities++;
                iteratorEpochWindowAverageInclinations++;
            }

            // Close file handler.
            epochWindowAveragesFile.close( );
        }        

        ///////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////

        // Write random walk output to database. To avoid locking of the database, this section is
        // thread-critical, so will be executed one-by-one by multiple threads.

        // Check if output mode is set to "DATABASE".
        if ( iequals( outputMode, "DATABASE" ) )
        {
    #pragma omp critical( writeToDatabase )
            {
                // Populate output table in database.
                populateRandomWalkOutputTable( 
                    databasePath, iteratorRandomWalkInputTable->monteCarloRunId,
                    averageLongitudeResidual, maximumLongitudeResidualChange,
                    averageEccentricity, maximumEccentricityChange,
                    averageInclination, maximumInclinationChange,
                    randomWalkOutputTableName, randomWalkInputTableName );
            }
        }

        ///////////////////////////////////////////////////////////////////
    }

    ///////////////////////////////////////////////////////////////////////////

    return EXIT_SUCCESS;
}
