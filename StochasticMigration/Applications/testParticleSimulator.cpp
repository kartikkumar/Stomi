/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120330    K. Kumar          File created from mabSystemIntegrator.cpp.
 *      120331    K. Kumar          Updated file to work with openmp and serialize write
 *                                  operations to database.
 *      120403    K. Kumar          Updated code to write kick tables to temoporary swap files;
 *                                  the database is updated at the end, so there is no concurrent
 *                                  access, which leads to problems.
 *      120502    K. Kumar          Modified filename to mabSystemSimulator.cpp.
 *      120503    K. Kumar          Updated code to not save Mab's state; to reduce memory
 *                                  footprint.
 *      120522    K. Kumar          Changed filename to testParticleSimulator.cpp.
 *      120629    K. Kumar          Updated propagation loop to used fixed output interval.
 *      120808    K. Kumar          Updated to new dictionary-based input file system.
 *      120830    K. Kumar          Cleaned up code; implemented composite state derivative model.
 *
 *    References
 *      Kumar, K., de Pater, I., Showalter, M.R. In prep, 2013.
 *
 *    Notes
 *
 */

//#include <iomanip>
#include <iostream>
//#include <limits>
#include <string>

//#include <Eigen/Core>

// #include <cmath>
// #include <iostream>
// #include <limits>
// #include <stdexcept>
// #include <string>
// #include <utility>

// #include <boost/assign/list_of.hpp>
// #include <boost/bind.hpp>
// #include <boost/exception/all.hpp>
#include <boost/make_shared.hpp>
// #include <boost/shared_ptr.hpp>

// #include <omp.h>

#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/basics.h>
#include <Assist/InputOutput/basicInputOutput.h>

// #include <TudatCore/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
// #include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

// #include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
// #include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>
// #include <Tudat/Astrodynamics/Gravitation/centralJ2GravityModel.h>
// #include <Tudat/Astrodynamics/StateDerivativeModels/cartesianStateDerivativeModel.h>
// #include <Tudat/Astrodynamics/StateDerivativeModels/compositeStateDerivativeModel.h>
#include <Tudat/InputOutput/dictionaryTools.h>
#include <Tudat/InputOutput/fieldType.h>
#include <Tudat/InputOutput/parsedDataVectorUtilities.h>
#include <Tudat/InputOutput/separatedParser.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h>
// #include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include "StochasticMigration/Basics/basics.h"
#include "StochasticMigration/Database/databaseReadFunctions.h"
#include "StochasticMigration/Database/testParticleCase.h"
// #include "StochasticMigration/Database/testParticleInput.h"
#include "StochasticMigration/InputOutput/dictionaries.h"

//! Execute test particle simulations.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{

    ///////////////////////////////////////////////////////////////////////////

    // Declare using-statements and type definitions.

//     using std::fabs;
//     using std::make_pair;
//     using std::numeric_limits;
//     using std::runtime_error;

//     using boost::assign::list_of;
//     using boost::bind;
//     using boost::enable_error_info;
    using boost::make_shared;
//     using boost::shared_ptr;
//     using boost::throw_exception;

    using namespace assist::astrodynamics;
    using namespace assist::basics;
    using namespace assist::input_output;

//     using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_astrodynamics::physical_constants;
    using namespace tudat::basic_astrodynamics::unit_conversions;
    using namespace tudat::basic_mathematics;
//     using namespace tudat::basic_mathematics::mathematical_constants;
//     using namespace tudat::gravitation;
    using namespace tudat::input_output;
    using namespace tudat::input_output::dictionary;
    using namespace tudat::input_output::field_types::general;
    using namespace tudat::input_output::parsed_data_vector_utilities;
    using namespace tudat::numerical_integrators;
//     using namespace tudat::state_derivative_models;

//     using namespace general_tools::astrodynamics;
//     using namespace general_tools::basics;
//     using namespace general_tools::input_output;

    using namespace stochastic_migration::basics;
    using namespace stochastic_migration::database;
    using namespace stochastic_migration::input_output;

//     // Typedefs.
//     typedef CompositeStateDerivativeModel< double, Vector12d, Vector6d >
//             CompositeStateDerivativeModel12d;
//     typedef shared_ptr< CompositeStateDerivativeModel12d > CompositeStateDerivativeModel12dPointer;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up input deck.

    // Check number of input parameters is correct (the numberOfInputs variable includes the
    // application itself, so one is subtracted from this number).
    checkNumberOfInputArguments( numberOfInputs - 1 );

    // Get input parameter dictionary.
    const DictionaryPointer dictionary = getTestParticleSimulatorDictionary( );

    // Read and filter input stream.
    std::string filteredInput = readAndFilterInputFile( inputArguments[ 1 ] );

    // Declare a separated parser.
    SeparatedParser parser( std::string( ": " ), 2, parameterName, parameterValue );

    // Parse filtered data.
    const ParsedDataVectorPtr parsedData = parser.parse( filteredInput );

    std::cout << std::endl;
    std::cout << "****************************************************************************" 
              << std::endl;
    std::cout << "Input parameters" << std::endl;
    std::cout << "****************************************************************************" 
              << std::endl;
    std::cout << std::endl;

    // Extract input parameters.
    const std::string applicationMode = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "APPLICATIONMODE" ), "BULKSIMULATION" );
    std::cout << "Application mode                                          " 
              << applicationMode << std::endl;

    const std::string databasePath = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "DATABASE" ) );
    std::cout << "Database                                                  "
              << databasePath << std::endl;

    const std::string caseName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "CASE" ) );
    std::cout << "Case                                                      " 
              << caseName << std::endl;              

    const int numberOfThreads = extractParameterValue< int >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "NUMBEROFTHREADS" ), 1 );
    std::cout << "Number of threads                                         "
              << numberOfThreads << std::endl;

    const std::string fileOutputDirectory = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "FILEOUTPUTDIRECTORY" ), "NO FILE OUTPUT" );
    std::cout << "File output directory                                     "
              << fileOutputDirectory << std::endl;

    const std::string simulationsToExecute = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "SIMULATIONSTOEXECUTE" ), "ALL" );
    std::cout << "Simulations to execute                                    "
              << simulationsToExecute << std::endl;

    const std::string testParticleCaseTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLECASETABLENAME" ), "test_particle_case" );
    std::cout << "Test particle case table                                  "
              << testParticleCaseTableName << std::endl;

    const std::string testParticleInputTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLEINPUTTABLENAME" ), "test_particle_input" );
    std::cout << "Test particle input table                                 "
              << testParticleInputTableName << std::endl;

    const std::string testParticleKickTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLEKICKTABLENAME" ), "test_particle_kicks" );
    std::cout << "Test particle kick table                                  "
              << testParticleKickTableName << std::endl;

    // Retrieve and store test particle case data from database.
    const TestParticleCasePointer caseDataFromDatabase = getTestParticleCase(
                databasePath, caseName, testParticleCaseTableName );

    // Check if any case parameters are overwritten by user input.
    const double randomWalkSimulationDuration = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKSIMULATIONDURATION" ),
                caseDataFromDatabase->randomWalkSimulationDuration, &convertJulianYearsToSeconds );
    std::cout << "Random walk simulation duration                           " 
              << randomWalkSimulationDuration / JULIAN_YEAR << " yrs" << std::endl;

    const double synodicPeriodLimit = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "SYNODICPERIODLIMIT" ),
                caseDataFromDatabase->synodicPeriodLimit, &convertJulianYearsToSeconds );
    std::cout << "Synodic period limit                                      " 
              << synodicPeriodLimit / JULIAN_YEAR << " yrs" << std::endl;

    const double outputInterval = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "OUTPUTINTERVAL" ),
                caseDataFromDatabase->outputInterval, &convertHoursToSeconds< double > );
    std::cout << "Output interval                                           " 
              << convertSecondsToHours( outputInterval ) << " hrs" << std::endl;

    const double startUpIntegrationDuration = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "STARTUPINTEGRATIONDURATION" ),
                caseDataFromDatabase->startUpIntegrationDuration, &convertJulianYearsToSeconds );
    std::cout << "Start-up integration duration                             " 
              << startUpIntegrationDuration / JULIAN_YEAR << " yrs" << std::endl;

    const double conjunctionEventDetectionDistance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CONJUNCTIONEVENTDETECTIONDISTANCE" ), 
                caseDataFromDatabase->conjunctionEventDetectionDistance );
    std::cout << "Conjunction event detection distance                      " 
              << convertMetersToKilometers( conjunctionEventDetectionDistance ) 
              << " km" << std::endl;

    const double oppositionEventDetectionDistance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "OPPOSITIONEVENTDETECTIONDISTANCE" ),
                caseDataFromDatabase->oppositionEventDetectionDistance );
    std::cout << "Opposition event detection distance                       " 
              << convertMetersToKilometers( oppositionEventDetectionDistance ) 
              << " km" << std::endl;

    const double centralBodyGravitationalParameter = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYGRAVITATIONALPARAMETER" ),
                caseDataFromDatabase->centralBodyGravitationalParameter );
    std::cout << "Central body gravitational parameter                      " 
              << centralBodyGravitationalParameter << " m^3 s^-2" << std::endl;

    const double centralBodyJ2GravityCoefficient = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYJ2GRAVITYCOEFFICIENT" ), 
                caseDataFromDatabase->centralBodyJ2GravityCoefficient );
    std::cout << "Central body J2 gravity coefficient                       "
              << centralBodyJ2GravityCoefficient << std::endl;

    const double centralBodyEquatorialRadius = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYEQUATORIALRADIUS" ),
                caseDataFromDatabase->centralBodyEquatorialRadius );
    std::cout << "Central body equatorial radius                            "
              << convertMetersToKilometers( centralBodyEquatorialRadius ) << " km" << std::endl;

    const double perturbedBodyRadius = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYRADIUS" ),
                caseDataFromDatabase->perturbedBodyRadius, &convertKilometersToMeters< double > );
    std::cout << "Perturbed body radius                                     " 
              << convertMetersToKilometers( perturbedBodyRadius ) << " km" << std::endl;

    const double perturbedBodyBulkDensity = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYBULKDENSITY" ),
                caseDataFromDatabase->perturbedBodyBulkDensity );
    std::cout << "Perturbed body bulk density                               " 
              << perturbedBodyBulkDensity << " kg m^-3" << std::endl;

    Vector6d perturbedBodyStateInKeplerianElementsAtT0( 6 );

    perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYSEMIMAJORAXISATT0" ),
                caseDataFromDatabase->perturbedBodyStateInKeplerianElementsAtT0( 
                    semiMajorAxisIndex ),
                &convertKilometersToMeters< double > );
    std::cout << "Perturbed body semi-major axis at TO                      "
              << convertMetersToKilometers( 
                    perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) ) 
              << " km" << std::endl;

    perturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYECCENTRICITYATT0" ), 
                caseDataFromDatabase->perturbedBodyStateInKeplerianElementsAtT0( 
                    eccentricityIndex ) );
    std::cout << "Perturbed body eccentricity at TO                         "
              << perturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex ) << std::endl;

    perturbedBodyStateInKeplerianElementsAtT0( inclinationIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYINCLINATIONATT0" ),
                caseDataFromDatabase->perturbedBodyStateInKeplerianElementsAtT0( 
                    inclinationIndex ),
                 &convertDegreesToRadians< double > );
    std::cout << "Perturbed body inclination at TO                          "
              << convertRadiansToDegrees( 
                    perturbedBodyStateInKeplerianElementsAtT0( inclinationIndex ) ) 
              << " deg" << std::endl;

    perturbedBodyStateInKeplerianElementsAtT0( argumentOfPeriapsisIndex )
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYARGUMENTOFPERIAPSISATT0" ),
                caseDataFromDatabase->perturbedBodyStateInKeplerianElementsAtT0( 
                    argumentOfPeriapsisIndex ),
                &convertDegreesToRadians< double > );
    std::cout << "Perturbed body argument of periapsis at TO                "
              << convertRadiansToDegrees( 
                    perturbedBodyStateInKeplerianElementsAtT0( argumentOfPeriapsisIndex ) ) 
              << " deg" << std::endl;

    perturbedBodyStateInKeplerianElementsAtT0( longitudeOfAscendingNodeIndex )
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYLONGITUDEOFASCENDINGNODEATT0" ),
                caseDataFromDatabase->perturbedBodyStateInKeplerianElementsAtT0( 
                    longitudeOfAscendingNodeIndex ),
                &convertDegreesToRadians< double > );
    std::cout << "Perturbed body longitude of ascending node at TO          "
              << convertRadiansToDegrees( 
                    perturbedBodyStateInKeplerianElementsAtT0( longitudeOfAscendingNodeIndex ) ) 
              << " deg" << std::endl;

    perturbedBodyStateInKeplerianElementsAtT0( trueAnomalyIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYTRUEANOMALYATT0" ),
                caseDataFromDatabase->perturbedBodyStateInKeplerianElementsAtT0(
                    trueAnomalyIndex ),
                 &convertDegreesToRadians< double > );
    std::cout << "Perturbed body true anomaly at TO                         "
              << convertRadiansToDegrees( 
                    perturbedBodyStateInKeplerianElementsAtT0( trueAnomalyIndex ) ) 
              << " deg" << std::endl;

    const std::string numericalIntegratorType = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "NUMERICALINTEGRATORTYPE" ),
                caseDataFromDatabase->numericalIntegratorType );
    std::cout << "Numerical integrator type                                 "
              << numericalIntegratorType << std::endl;

    const double initialStepSize = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INITIALSTEPSIZE" ), 
                caseDataFromDatabase->initialStepSize );
    std::cout << "Initial step size                                         "
              << initialStepSize << " s" << std::endl;

    const double numericalIntegratorRelativeTolerance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RUNGEKUTTARELATIVEERRORTOLERANCE" ),
                caseDataFromDatabase->numericalIntegratorRelativeTolerance );
    std::cout << "Numerical integrator relative tolerance                   " 
              << numericalIntegratorRelativeTolerance << std::endl;

    const double numericalIntegratorAbsoluteTolerance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RUNGEKUTTAABSOLUTEERRORTOLERANCE" ), 
                caseDataFromDatabase->numericalIntegratorAbsoluteTolerance );
    std::cout << "Numerical integrator absolute tolerance                   " 
              << numericalIntegratorAbsoluteTolerance << std::endl;

    // Store case data with possible overwritten data.
    TestParticleCasePointer testParticleCase = make_shared< TestParticleCase >(
        TestParticleCase( 
            caseDataFromDatabase->caseId, caseDataFromDatabase->caseName, 
            randomWalkSimulationDuration, synodicPeriodLimit, outputInterval, 
            startUpIntegrationDuration, conjunctionEventDetectionDistance, 
            oppositionEventDetectionDistance, centralBodyGravitationalParameter,
            centralBodyJ2GravityCoefficient, centralBodyEquatorialRadius, 
            caseDataFromDatabase->semiMajorAxisDistributionLimit,
            caseDataFromDatabase->eccentricityDistributionMean, 
            caseDataFromDatabase->eccentricityDistributionAngle,
            caseDataFromDatabase->eccentricityDistributionFullWidthHalfMaximum,
            caseDataFromDatabase->inclinationDistributionMean, 
            caseDataFromDatabase->inclinationDistributionAngle,
            caseDataFromDatabase->inclinationDistributionFullWidthHalfMaximum, 
            perturbedBodyRadius, perturbedBodyBulkDensity, 
            perturbedBodyStateInKeplerianElementsAtT0, numericalIntegratorType, initialStepSize,
            numericalIntegratorRelativeTolerance, numericalIntegratorRelativeTolerance ) );

    // Check that all required parameters have been set.
    checkRequiredParameters( dictionary );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Compute derived parameters.

    std::cout << std::endl;
    std::cout << "****************************************************************************" 
              << std::endl;
    std::cout << "Derived parameters" << std::endl;
    std::cout << "****************************************************************************" 
              << std::endl;
    std::cout << std::endl;

    // Compute mass of perturbed body [kg].
    const double perturbedBodyMass = computeMassOfSphere(
                testParticleCase->perturbedBodyRadius, 
                testParticleCase->perturbedBodyBulkDensity );
    std::cout << "Perturbed body mass                                       " 
              << perturbedBodyMass << " kg" << std::endl;

    // Compute perturbed body's gravitational parameter [m^3 s^-2].
    const double perturbedBodyGravitationalParameter
            = computeGravitationalParameter( perturbedBodyMass );
    std::cout << "Perturbed body gravitational parameter                    " 
              << perturbedBodyGravitationalParameter << " m^3 s^-2" << std::endl;

    // Set coefficient set selected for numerical integrator.
    RungeKuttaCoefficients rungeKuttaCoefficients 
        = getRungeKuttaCoefficients( testParticleCase->numericalIntegratorType );

    ///////////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////////

    // Fetch test particle input table.    

    // Check if all incomplete simulations are to be run and fetch the input table, else only fetch
    // the requested test particle simulation numbers.
    // TestParticleInputTablePointer testParticleInputTable;

    // if ( iequals( applicationMode, "BULKSIMULATION" )  )
    // {
    //     // Get entire test particle input table from database.
    //     testParticleInputTable = getTestParticleInputTable(
    //                 databasePath, false, testParticleInputTableName );
    // }

//     else if ( iequals( applicationMode, "SINGLESIMULATION" ) )
//     {
//         // Get selected test particle input table from database.
//         testParticleInputTable = getTestParticleInputTable(
//                     databasePath, testParticleSimulations, testParticleInputTableName );
//     }

//     // DEBUG.
//     if ( iequals( debugMode, "ON" ) )
//     {
//         std::cout << "test particle simulation number, completed, semi-major axis [m], "
//                   << "eccentricity [-], inclination [rad], argument of periapsis [rad], "
//                   << "longitude of ascending node [rad], true anomaly [rad]" << std::endl;

//         for ( unsigned int i = 0; i < testParticleInputTable.size( ); i++ )
//         {
//             std::cout << testParticleInputTable.at( i ).simulationNumber << ", "
//                       << testParticleInputTable.at( i ).isCompleted << ", "
//                       << testParticleInputTable.at( i ).initialStateInKeplerianElements(
//                              semiMajorAxisIndex ) << ", "
//                       << testParticleInputTable.at( i ).initialStateInKeplerianElements(
//                              eccentricityIndex ) << ", "
//                       << testParticleInputTable.at( i ).initialStateInKeplerianElements(
//                              inclinationIndex ) << ", "
//                       << testParticleInputTable.at( i ).initialStateInKeplerianElements(
//                              argumentOfPeriapsisIndex ) << ", "
//                       << testParticleInputTable.at( i ).initialStateInKeplerianElements(
//                              longitudeOfAscendingNodeIndex ) << ", "
//                       << testParticleInputTable.at( i ).initialStateInKeplerianElements(
//                              trueAnomalyIndex ) << std::endl;
//         }

//         std::cout << std::endl;
//     }

    ///////////////////////////////////////////////////////////////////////////  





//     ///////////////////////////////////////////////////////////////////////////

//     // Execute simulation loop.
// #pragma omp parallel for num_threads( numberOfThreads ) schedule( static, 1 )
//     for ( unsigned int i = 0; i < testParticleInputTable.size( ); i++ )
//     {

// #pragma omp critical( outputToConsole )
//         {
//             std::cout << "Run " << i + 1 << " / " << testParticleInputTable.size( )
//                       << " (simulation " << testParticleInputTable.at( i ).simulationNumber << ")"
//                       << " on thread " << omp_get_thread_num( ) + 1
//                       << " / " <<  omp_get_num_threads( ) << std::endl;
//         }

//         ///////////////////////////////////////////////////////////////////////////

//         // Create perturbed body and test particle.

//         // Convert test particle initial state in Keplerian elements to Cartesian elements.
//         const Vector6d testParticleInitialState
//                 = convertKeplerianToCartesianElements(
//                     testParticleInputTable.at( i ).initialStateInKeplerianElements,
//                     testParticleCase->centralBodyGravitationalParameter );

//         // Convert perturbed body initial state in Keplerian elements to Cartesian elements.
//         const Vector6d perturbedBodyInitialState
//                 = convertKeplerianToCartesianElements(
//                     testParticleCase->perturbedBodyStateInKeplerianElementsAtT0,
//                     testParticleCase->centralBodyGravitationalParameter );

//         // Create perturbed body and test particle.
//         BodyPointer perturbedBody = make_shared< Body >(
//                     "Perturbed body", perturbedBodyInitialState );
//         BodyPointer testParticle = make_shared< Body >(
//                     "Test particle", testParticleInitialState );

//         ///////////////////////////////////////////////////////////////////////////

//         ///////////////////////////////////////////////////////////////////////////

//         // Set up dynamics.

//         // Set acceleration models for perturbed body and test particle.
//         AccelerationModel3dPointer perturbedBodyGravityOnTestParticle
//                 = make_shared< CentralGravitationalAccelerationModel3d >(
//                     bind( &Body::getCurrentPosition, testParticle ),
//                     testParticleCase->perturbedBodyGravitationalParameter,
//                     bind( &Body::getCurrentPosition, perturbedBody ) );

//         AccelerationModel3dPointer centralBodyGravityOnPerturbedBody;
//         AccelerationModel3dPointer centralBodyGravityOnTestParticle;

//         if ( fabs( testParticleCase->centralBodyJ2GravityCoefficient )
//              < numeric_limits< double >::min( ) )
//         {
//             centralBodyGravityOnPerturbedBody
//                     = make_shared< CentralGravitationalAccelerationModel3d >(
//                         bind( &Body::getCurrentPosition, perturbedBody ),
//                         testParticleCase->centralBodyGravitationalParameter );

//             centralBodyGravityOnTestParticle
//                     = make_shared< CentralGravitationalAccelerationModel3d >(
//                         bind( &Body::getCurrentPosition, testParticle ),
//                         testParticleCase->centralBodyGravitationalParameter );
//         }

//         else
//         {
//             centralBodyGravityOnPerturbedBody
//                     = make_shared< CentralJ2GravitationalAccelerationModel3d >(
//                         bind( &Body::getCurrentPosition, perturbedBody ),
//                         testParticleCase->centralBodyGravitationalParameter,
//                         testParticleCase->centralBodyJ2GravityCoefficient,
//                         testParticleCase->centralBodyEquatorialRadius );

//             centralBodyGravityOnTestParticle
//                     = make_shared< CentralJ2GravitationalAccelerationModel3d >(
//                         bind( &Body::getCurrentPosition, testParticle ),
//                         testParticleCase->centralBodyGravitationalParameter,
//                         testParticleCase->centralBodyJ2GravityCoefficient,
//                         testParticleCase->centralBodyEquatorialRadius );
//         }

//         // Create lists of acceleration models to provide to state derivative models.
//         CartesianStateDerivativeModel6d::AccelerationModelPointerVector
//                 perturbedBodyAccelerationList = list_of( centralBodyGravityOnPerturbedBody );

//         CartesianStateDerivativeModel6d::AccelerationModelPointerVector
//                 testParticleAccelerationList = list_of( centralBodyGravityOnTestParticle )(
//                     perturbedBodyGravityOnTestParticle );

//         // Set Cartesian state derivative models for perturbed body and test particle.
//         CartesianStateDerivativeModel6dPointer perturbedBodyStateDerivative
//                 = make_shared< CartesianStateDerivativeModel6d >(
//                     perturbedBodyAccelerationList, &updateNothing< double, Vector6d > );

//         CartesianStateDerivativeModel6dPointer testParticleStateDerivative
//                 = make_shared< CartesianStateDerivativeModel6d >(
//                     testParticleAccelerationList, &updateNothing< double, Vector6d > );

//         // Construct and return composite state derivative model.
//         CompositeStateDerivativeModel12d::VectorStateDerivativeModelMap stateDerivativeModelMap;
//         stateDerivativeModelMap[ make_pair( 0, 6 ) ]
//                 = bind( &CartesianStateDerivativeModel6d::computeStateDerivative,
//                         perturbedBodyStateDerivative, _1, _2 );
//         stateDerivativeModelMap[ make_pair( 6, 6 ) ]
//                 = bind( &CartesianStateDerivativeModel6d::computeStateDerivative,
//                         testParticleStateDerivative, _1, _2 );

//         // Construct data updater that disassembles the composite state and updates the states of
//         // the bodies.
//         DataUpdaterPointer dataUpdater = make_shared< DataUpdater >( perturbedBody, testParticle );

//         // Construct composite state derivative model.
//         CompositeStateDerivativeModel12dPointer stateDerivativeModel
//                 = make_shared< CompositeStateDerivativeModel12d >(
//                     stateDerivativeModelMap,
//                     bind( &DataUpdater::updateTimeAndCompositeState, dataUpdater, _1, _2 ) );

//         ///////////////////////////////////////////////////////////////////////////

//         ///////////////////////////////////////////////////////////////////////////

//         // Set up integrator.

//         // Declare Runge-Kutta, variable-stepsize, integrator.
//         RungeKuttaVariableStepSizeIntegratorXdPointer integrator
//                 = make_shared< RungeKuttaVariableStepSizeIntegratorXd >(
//                     rungeKuttaCoefficients,
//                     bind( &CompositeStateDerivativeModel12d::computeStateDerivative,
//                           stateDerivativeModel, _1, _2 ),
//                     0.0,
//                     ( Eigen::VectorXd( 12 ) << perturbedBody->getCurrentState( ),
//                       testParticle->getCurrentState( ) ).finished( ),
//                     numeric_limits< double >::epsilon( ),
//                     numeric_limits< double >::max( ),
//                     testParticleCase->numericalnumericalIntegratorRelativeTolerance,
//                     testParticleCase->numericalnumericalIntegratorAbsoluteTolerance );

//         // Numerically integrate motion of test particle up to end of start-up period.
//         integrator->integrateTo( testParticleCase->startUpIntegrationDuration,
//                                  testParticleCase->initialStepSize );

//         ///////////////////////////////////////////////////////////////////////////

//         ///////////////////////////////////////////////////////////////////////////

//         // Compute system parameters at end of start-up period (T0).

//         // Convert test particle's state to Keplerian elements.
//         const Vector6d testParticleStateAfterStartUp = convertCartesianToKeplerianElements(
//                     testParticle->getCurrentState( ),
//                     testParticleCase->centralBodyGravitationalParameter );

//         // Compute orbital period of perturbed body [s].
//         const double orbitalPeriodOfPerturbedBody = computeKeplerOrbitalPeriod(
//                     testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
//                         semiMajorAxisIndex ),
//                     testParticleCase->centralBodyGravitationalParameter );

//         // Compute orbital period of test particle [s].
//         const double orbitalPeriodOfTestParticle = computeKeplerOrbitalPeriod(
//                     testParticleStateAfterStartUp( semiMajorAxisIndex ),
//                     testParticleCase->centralBodyGravitationalParameter );

//         // Compute synodic period of test particle's motion with respect to perturbed body [s].
//         const double synodicPeriod = computeSynodicPeriod(
//                     orbitalPeriodOfPerturbedBody, orbitalPeriodOfTestParticle );

//         // Compute perturbed body's initial energy.
//         const double perturbedBodyInitialEnergy = computeKeplerEnergy(
//                     testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
//                         semiMajorAxisIndex ),
//                     testParticleCase->centralBodyGravitationalParameter );

//         // Compute perturbed body's initial angular momentum.
//         const double perturbedBodyInitialAngularMomentum = computeKeplerAngularMomentum(
//                     testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
//                         semiMajorAxisIndex ),
//                     testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
//                         eccentricityIndex ),
//                     testParticleCase->centralBodyGravitationalParameter );

//         // DEBUG.
//         if ( iequals( debugMode, "ON" ) )
//         {
//             std::cout << "Test particle orbital period [hr]: "
//                       << convertSecondsToHours( orbitalPeriodOfTestParticle ) << std::endl;
//             std::cout << "Perturbed body orbital period [hr]: "
//                       << convertSecondsToHours( orbitalPeriodOfPerturbedBody ) << std::endl;
//             std::cout << "Test particle synodic period [yr]: "
//                       << synodicPeriod / JULIAN_YEAR << std::endl;
//             std::cout << "Perturbed body's initial energy [m^2 s^-2]: "
//                       << perturbedBodyInitialEnergy << std::endl;
//             std::cout << "Perturbed body's initial angular momentum [m^2 s^-1]: "
//                       << perturbedBodyInitialAngularMomentum << std::endl;
//         }

//         // Numerically integrate motion of test particle back by one synodic period
//         // (TMinusSynodicPeriod).
//         integrator->integrateTo( testParticleCase->startUpIntegrationDuration - synodicPeriod,
//                                  -testParticleCase->initialStepSize );

//         // Copy integrator in case the full history needs to be written to file.
//         RungeKuttaVariableStepSizeIntegratorXdPointer integratorCopy =
//                 make_shared< RungeKuttaVariableStepSizeIntegratorXd >( *integrator );

//         ///////////////////////////////////////////////////////////////////////////

//         ///////////////////////////////////////////////////////////////////////////

//         // Numerically integrate system from TMinusSynodicPeriod to
//         // simulationDuration + synodicPeriod.

//         // Declare error in perturbed body's energy.
//         double perturbedBodyEnergyError = TUDAT_NAN;

//         // Declare error in perturbed body's angular momentum.
//         double perturbedBodyAngularMomentumError = TUDAT_NAN;

// //        // Propagate system and generate test particle kick table.
// //        KickTable kickTable = propagateSystemAndGenerateKickTable(
// //                    mab, testParticle, initialStepSize, simulationDuration, synodicPeriod,
// //                    startUpIntegrationDuration, minimumDistanceCrossing, maximumDistanceCrossing,
// //                    integrator, testParticleCase->uranusGravitationalParameter, mabInitialEnergy,
// //                    mabInitialAngularMomentum, perturbedBodyEnergyError, perturbedBodyAngularMomentumError );

// ////        // DEBUG.
// ////        std::cout << "Mab energy error: " << perturbedBodyEnergyError << std::endl;
// ////        std::cout << "Mab angular momentum error: " << perturbedBodyAngularMomentumError << std::endl;

// //        ///////////////////////////////////////////////////////////////////////////

// //        // Save output. To avoid locking of the database, this section is thread-critical, so will
// //        // be executed one-by-one by multiple threads.

// //#pragma omp critical( writeToDatabase )
// //        {
// //            // Check if output should be written to database or to files.
// //            if ( outputDirectory.empty( ) )
// //            {
// //                // Populate kick table in database.
// //                mab_simulations::database_functions::populateKickTable(
// //                            databasePath, testParticleInputTable.at( i ).simulationNumber,
// //                            kickTable, perturbedBodyEnergyError, perturbedBodyAngularMomentumError );
// //            }

// //            // If file output is required, re-integrate the system and generate file output at
// //            // fixed output intervals defined by the user.
// //            else
// //            {
// //                // Update time and state stored in Mab and test particle bodies to
// //                // TMinusSynodicPeriod.
// //                dataUpdater->updateTimeAndCompositeState(
// //                            integratorCopy->getCurrentIndependentVariable( ),
// //                            integratorCopy->getCurrentState( ) );

// //                // Propagate system and write output generated to files.
// //                propagateSystemAndGenerateFileOutput(
// //                            simulationDuration, synodicPeriod, outputInterval, initialStepSize,
// //                            startUpIntegrationDuration, integratorCopy, mab, testParticle,
// //                            caseNumber, testParticleInputTable.at( i ).simulationNumber, outputDirectory,
// //                            std::numeric_limits< double >::digits10,
// //                            kickTable, testParticleCase->uranusGravitationalParameter );
// //            }
// //        }

//         std::cout << std::endl;

//         /////////////////////////////////////////////////////////////////////////
//     }

    return 0;
}
