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

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
 #include <utility>

#include <Eigen/Core>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <omp.h>

#include <Assist/Astrodynamics/body.h>
#include <Assist/Astrodynamics/dataUpdater.h>
#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/basics.h>
#include <Assist/Basics/commonTypedefs.h>
#include <Assist/InputOutput/basicInputOutput.h>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>
#include <Tudat/Astrodynamics/Gravitation/centralJ2GravityModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/cartesianStateDerivativeModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/compositeStateDerivativeModel.h>
#include <Tudat/InputOutput/dictionaryTools.h>
#include <Tudat/InputOutput/fieldType.h>
#include <Tudat/InputOutput/parsedDataVectorUtilities.h>
#include <Tudat/InputOutput/separatedParser.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include "StochasticMigration/Basics/basics.h"
#include "StochasticMigration/Database/databaseReadFunctions.h"
#include "StochasticMigration/Database/testParticleCase.h"
#include "StochasticMigration/Database/testParticleInput.h"
#include "StochasticMigration/InputOutput/dictionaries.h"

//! Execute test particle simulations.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{

    ///////////////////////////////////////////////////////////////////////////

    // Declare using-statements and type definitions.
    using std::cout;
    using std::endl;
    using std::fabs;
    using std::numeric_limits;
    using std::string;
    using std::make_pair;

    using boost::assign::list_of;
    using boost::bind;
    using boost::iequals;
    using boost::make_shared;
    using boost::shared_ptr;

    using namespace Eigen;

    using namespace assist::astrodynamics;
    using namespace assist::basics;
    using namespace assist::input_output;

    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_astrodynamics::physical_constants;
    using namespace tudat::basic_astrodynamics::unit_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::input_output;
    using namespace tudat::input_output::dictionary;
    using namespace tudat::input_output::field_types::general;
    using namespace tudat::input_output::parsed_data_vector_utilities;
    using namespace tudat::numerical_integrators;
    using namespace tudat::state_derivative_models;

    using namespace stochastic_migration::basics;
    using namespace stochastic_migration::database;
    using namespace stochastic_migration::input_output;

    // Typedefs.
    typedef CompositeStateDerivativeModel< double, Vector12d, Vector6d > 
            CompositeStateDerivativeModel12d;
    typedef shared_ptr< CompositeStateDerivativeModel12d > CompositeStateDerivativeModel12dPointer;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up input deck.

    // Check number of input parameters is correct (the numberOfInputs variable includes the
    // application itself, so one is subtracted from this number).
    checkNumberOfInputArguments( numberOfInputs - 1 );

    // Get input parameter dictionary.
    const DictionaryPointer dictionary = getTestParticleSimulatorDictionary( );

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
    const string applicationMode = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "APPLICATIONMODE" ), "BULKSIMULATION" );
    cout << "Application mode                                          " 
         << applicationMode << endl;

    const string databasePath = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "DATABASE" ) );
    cout << "Database                                                  "
         << databasePath << endl;

    const string caseName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "CASE" ) );
    cout << "Case                                                      " 
         << caseName << endl;              

    const int numberOfThreads = extractParameterValue< int >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "NUMBEROFTHREADS" ), 1 );
    cout << "Number of threads                                         "
         << numberOfThreads << endl;

    const string fileOutputDirectory = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "FILEOUTPUTDIRECTORY" ), "NO FILE OUTPUT" );
    cout << "File output directory                                     "
         << fileOutputDirectory << endl;

    const string simulationsToExecute = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "SIMULATIONSTOEXECUTE" ), "ALL" );
    cout << "Simulations to execute                                    "
         << simulationsToExecute << endl;

    const string testParticleCaseTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLECASETABLENAME" ), "test_particle_case" );
    cout << "Test particle case table                                  "
         << testParticleCaseTableName << endl;

    const string testParticleInputTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLEINPUTTABLENAME" ), "test_particle_input" );
    cout << "Test particle input table                                 "
         << testParticleInputTableName << endl;

    const string testParticleKickTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLEKICKTABLENAME" ), "test_particle_kicks" );
    cout << "Test particle kick table                                  "
         << testParticleKickTableName << endl;

    // Retrieve and store test particle case data from database.
    const TestParticleCasePointer caseDataFromDatabase = getTestParticleCase(
                databasePath, caseName, testParticleCaseTableName );

    // Check if any case parameters are overwritten by user input.
    const double randomWalkSimulationDuration = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKSIMULATIONDURATION" ),
                caseDataFromDatabase->randomWalkSimulationDuration, &convertJulianYearsToSeconds );
    cout << "Random walk simulation duration                           " 
         << randomWalkSimulationDuration / JULIAN_YEAR << " yrs" << endl;

    const double synodicPeriodLimit = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "SYNODICPERIODLIMIT" ),
                caseDataFromDatabase->synodicPeriodLimit, &convertJulianYearsToSeconds );
    cout << "Synodic period limit                                      " 
         << synodicPeriodLimit / JULIAN_YEAR << " yrs" << endl;

    const double outputInterval = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "OUTPUTINTERVAL" ),
                caseDataFromDatabase->outputInterval, &convertHoursToSeconds< double > );
    cout << "Output interval                                           " 
         << convertSecondsToHours( outputInterval ) << " hrs" << endl;

    const double startUpIntegrationDuration = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "STARTUPINTEGRATIONDURATION" ),
                caseDataFromDatabase->startUpIntegrationDuration, &convertJulianYearsToSeconds );
    cout << "Start-up integration duration                             " 
         << startUpIntegrationDuration / JULIAN_YEAR << " yrs" << endl;

    const double conjunctionEventDetectionDistance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CONJUNCTIONEVENTDETECTIONDISTANCE" ), 
                caseDataFromDatabase->conjunctionEventDetectionDistance );
    cout << "Conjunction event detection distance                      " 
         << convertMetersToKilometers( conjunctionEventDetectionDistance ) << " km" << endl;

    const double oppositionEventDetectionDistance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "OPPOSITIONEVENTDETECTIONDISTANCE" ),
                caseDataFromDatabase->oppositionEventDetectionDistance );
    cout << "Opposition event detection distance                       " 
         << convertMetersToKilometers( oppositionEventDetectionDistance ) << " km" << endl;

    const double centralBodyGravitationalParameter = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYGRAVITATIONALPARAMETER" ),
                caseDataFromDatabase->centralBodyGravitationalParameter );
    cout << "Central body gravitational parameter                      " 
         << centralBodyGravitationalParameter << " m^3 s^-2" << endl;

    const double centralBodyJ2GravityCoefficient = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYJ2GRAVITYCOEFFICIENT" ), 
                caseDataFromDatabase->centralBodyJ2GravityCoefficient );
    cout << "Central body J2 gravity coefficient                       "
         << centralBodyJ2GravityCoefficient << endl;

    const double centralBodyEquatorialRadius = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYEQUATORIALRADIUS" ),
                caseDataFromDatabase->centralBodyEquatorialRadius );
    cout << "Central body equatorial radius                            "
         << convertMetersToKilometers( centralBodyEquatorialRadius ) << " km" << endl;

    const double perturbedBodyRadius = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYRADIUS" ),
                caseDataFromDatabase->perturbedBodyRadius, &convertKilometersToMeters< double > );
    cout << "Perturbed body radius                                     " 
         << convertMetersToKilometers( perturbedBodyRadius ) << " km" << endl;

    const double perturbedBodyBulkDensity = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYBULKDENSITY" ),
                caseDataFromDatabase->perturbedBodyBulkDensity );
    cout << "Perturbed body bulk density                               " 
         << perturbedBodyBulkDensity << " kg m^-3" << endl;

    Vector6d perturbedBodyStateInKeplerianElementsAtT0( 6 );

    perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYSEMIMAJORAXISATT0" ),
                caseDataFromDatabase->perturbedBodyStateInKeplerianElementsAtT0( 
                    semiMajorAxisIndex ),
                &convertKilometersToMeters< double > );
    cout << "Perturbed body semi-major axis at TO                      "
         << convertMetersToKilometers( 
                perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) ) 
         << " km" << endl;

    perturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYECCENTRICITYATT0" ), 
                caseDataFromDatabase->perturbedBodyStateInKeplerianElementsAtT0( 
                    eccentricityIndex ) );
    cout << "Perturbed body eccentricity at TO                         "
         << perturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex ) << endl;

    perturbedBodyStateInKeplerianElementsAtT0( inclinationIndex ) 
            = extractParameterValue< double >( 
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYINCLINATIONATT0" ),
                caseDataFromDatabase->perturbedBodyStateInKeplerianElementsAtT0( 
                    inclinationIndex ),
                 &convertDegreesToRadians< double > );
    cout << "Perturbed body inclination at TO                          "
         << convertRadiansToDegrees( 
                perturbedBodyStateInKeplerianElementsAtT0( inclinationIndex ) ) 
         << " deg" << endl;

    perturbedBodyStateInKeplerianElementsAtT0( argumentOfPeriapsisIndex )
            = extractParameterValue< double >( 
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYARGUMENTOFPERIAPSISATT0" ),
                caseDataFromDatabase->perturbedBodyStateInKeplerianElementsAtT0( 
                    argumentOfPeriapsisIndex ),
                &convertDegreesToRadians< double > );
    cout << "Perturbed body argument of periapsis at TO                "
         << convertRadiansToDegrees( 
                perturbedBodyStateInKeplerianElementsAtT0( argumentOfPeriapsisIndex ) ) 
         << " deg" << endl;

    perturbedBodyStateInKeplerianElementsAtT0( longitudeOfAscendingNodeIndex )
            = extractParameterValue< double >( 
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYLONGITUDEOFASCENDINGNODEATT0" ),
                caseDataFromDatabase->perturbedBodyStateInKeplerianElementsAtT0( 
                    longitudeOfAscendingNodeIndex ),
                &convertDegreesToRadians< double > );
    cout << "Perturbed body longitude of ascending node at TO          "
         << convertRadiansToDegrees( 
               perturbedBodyStateInKeplerianElementsAtT0( longitudeOfAscendingNodeIndex ) ) 
         << " deg" << endl;

    perturbedBodyStateInKeplerianElementsAtT0( trueAnomalyIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYTRUEANOMALYATT0" ),
                caseDataFromDatabase->perturbedBodyStateInKeplerianElementsAtT0(
                    trueAnomalyIndex ),
                 &convertDegreesToRadians< double > );
    cout << "Perturbed body true anomaly at TO                         "
         << convertRadiansToDegrees( 
               perturbedBodyStateInKeplerianElementsAtT0( trueAnomalyIndex ) ) 
         << " deg" << endl;

    const string numericalIntegratorType = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "NUMERICALINTEGRATORTYPE" ),
                caseDataFromDatabase->numericalIntegratorType );
    cout << "Numerical integrator type                                 "
         << numericalIntegratorType << endl;

    const double initialStepSize = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INITIALSTEPSIZE" ), 
                caseDataFromDatabase->initialStepSize );
    cout << "Initial step size                                         "
         << initialStepSize << " s" << endl;

    const double numericalIntegratorRelativeTolerance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RUNGEKUTTARELATIVEERRORTOLERANCE" ),
                caseDataFromDatabase->numericalIntegratorRelativeTolerance );
    cout << "Numerical integrator relative tolerance                   " 
         << numericalIntegratorRelativeTolerance << endl;

    const double numericalIntegratorAbsoluteTolerance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RUNGEKUTTAABSOLUTEERRORTOLERANCE" ), 
                caseDataFromDatabase->numericalIntegratorAbsoluteTolerance );
    cout << "Numerical integrator absolute tolerance                   " 
         << numericalIntegratorAbsoluteTolerance << endl;

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

    cout << endl;
    cout << "****************************************************************************" << endl;
    cout << "Derived parameters" << endl;
    cout << "****************************************************************************" << endl;
    cout << endl;

    // Compute mass of perturbed body [kg].
    const double perturbedBodyMass = computeMassOfSphere(
                testParticleCase->perturbedBodyRadius, 
                testParticleCase->perturbedBodyBulkDensity );
    cout << "Perturbed body mass                                       " 
         << perturbedBodyMass << " kg" << endl;

    // Compute perturbed body's gravitational parameter [m^3 s^-2].
    const double perturbedBodyGravitationalParameter
            = computeGravitationalParameter( perturbedBodyMass );
    cout << "Perturbed body gravitational parameter                    " 
         << perturbedBodyGravitationalParameter << " m^3 s^-2" << endl;

    // Set coefficient set selected for numerical integrator.
    RungeKuttaCoefficients rungeKuttaCoefficients 
        = getRungeKuttaCoefficients( testParticleCase->numericalIntegratorType );

    // Set start epoch for numerical integration.
    const double startEpoch = 0.0;            

    ///////////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////////

    // Fetch test particle input table.

    cout << endl;
    cout << "****************************************************************************" << endl;
    cout << "Database operations" << endl;
    cout << "****************************************************************************" << endl;
    cout << endl;

    // Generate output message.
    cout << "Fetching test particle input data from database ..." << endl;    

    // Check if all incomplete simulations are to be run and fetch the input table, else only fetch
    // the requested test particle simulation numbers.
    TestParticleInputTable testParticleInputTable;

    if ( iequals( applicationMode, "BULKSIMULATION" )  )
    {
        cout << "Fetching all incomplete test particle simulations ..." << endl;    

        // Get entire test particle input table from database.
        testParticleInputTable = getCompleteTestParticleInputTable(
                    databasePath, testParticleCase->caseId, testParticleInputTableName );
    }

    else if ( iequals( applicationMode, "SINGLESIMULATION" ) )
    {

        cout << "Fetching all requested test particle simulations ..." << endl;    

        // Get selected test particle input table from database.
        testParticleInputTable = getSelectedTestParticleInputTable(
                    databasePath, testParticleCase->caseId, 
                    simulationsToExecute, testParticleInputTableName );
    }

    cout << "Test particle input data (" << testParticleInputTable.size( )
         << " rows) fetched successfully from database!" << endl;    

    ///////////////////////////////////////////////////////////////////////////  

    ///////////////////////////////////////////////////////////////////////////

    cout << endl;
    cout << "****************************************************************************" << endl;
    cout << "Simulation loop" << endl;
    cout << "****************************************************************************" << endl;
    cout << endl;

    // Execute simulation loop.
    cout << "Starting simulation loop ... " << endl;
    cout << endl;

    TestParticleInputTable::iterator iteratorInputTable;
    iteratorInputTable = testParticleInputTable.begin( );

#pragma omp parallel for num_threads( numberOfThreads ) schedule( static, 1 )
    for ( unsigned int i = 0; i < testParticleInputTable.size( ); i++ )
    {

#pragma omp critical( outputToConsole )
        {
            cout << "Run " << i + 1 << " / " << testParticleInputTable.size( )
                 << " (simulation " << iteratorInputTable->simulationNumber << ")"
                 << " on thread " << omp_get_thread_num( ) + 1
                 << " / " << omp_get_num_threads( ) << endl;
        }

        ///////////////////////////////////////////////////////////////////////////

        // Create perturbed body and test particle.

        // Convert test particle initial state in Keplerian elements to Cartesian elements.
        const Vector6d testParticleInitialState
                = convertKeplerianToCartesianElements(
                    iteratorInputTable->initialStateInKeplerianElements,
                    testParticleCase->centralBodyGravitationalParameter );

        // Convert perturbed body initial state in Keplerian elements to Cartesian elements.
        const Vector6d perturbedBodyInitialState
                = convertKeplerianToCartesianElements(
                    testParticleCase->perturbedBodyStateInKeplerianElementsAtT0,
                    testParticleCase->centralBodyGravitationalParameter );

        // Create perturbed body and test particle.
        BodyPointer perturbedBody = make_shared< Body >(
                    "Perturbed body", perturbedBodyInitialState );
        BodyPointer testParticle = make_shared< Body >(
                    "Test particle", testParticleInitialState );

        ///////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////

        // Set up dynamics.

        // Set acceleration models for perturbed body and test particle.
        AccelerationModel3dPointer perturbedBodyGravityOnTestParticle
                = make_shared< CentralGravitationalAccelerationModel3d >(
                    bind( &Body::getCurrentPosition, testParticle ),
                    perturbedBodyGravitationalParameter,
                    bind( &Body::getCurrentPosition, perturbedBody ) );

        AccelerationModel3dPointer centralBodyGravityOnPerturbedBody;
        AccelerationModel3dPointer centralBodyGravityOnTestParticle;

        if ( fabs( testParticleCase->centralBodyJ2GravityCoefficient )
             < numeric_limits< double >::min( ) )
        {
            centralBodyGravityOnPerturbedBody
                    = make_shared< CentralGravitationalAccelerationModel3d >(
                        bind( &Body::getCurrentPosition, perturbedBody ),
                        testParticleCase->centralBodyGravitationalParameter );

            centralBodyGravityOnTestParticle
                    = make_shared< CentralGravitationalAccelerationModel3d >(
                        bind( &Body::getCurrentPosition, testParticle ),
                        testParticleCase->centralBodyGravitationalParameter );
        }

        else
        {
            centralBodyGravityOnPerturbedBody
                    = make_shared< CentralJ2GravitationalAccelerationModel >(
                        bind( &Body::getCurrentPosition, perturbedBody ),
                        testParticleCase->centralBodyGravitationalParameter,
                        testParticleCase->centralBodyJ2GravityCoefficient,
                        testParticleCase->centralBodyEquatorialRadius );

            centralBodyGravityOnTestParticle
                    = make_shared< CentralJ2GravitationalAccelerationModel >(
                        bind( &Body::getCurrentPosition, testParticle ),
                        testParticleCase->centralBodyGravitationalParameter,
                        testParticleCase->centralBodyJ2GravityCoefficient,
                        testParticleCase->centralBodyEquatorialRadius );
        }

        // Create lists of acceleration models to provide to state derivative models.
        CartesianStateDerivativeModel6d::AccelerationModelPointerVector
                perturbedBodyAccelerationList = list_of( centralBodyGravityOnPerturbedBody );

        CartesianStateDerivativeModel6d::AccelerationModelPointerVector
                testParticleAccelerationList = list_of( centralBodyGravityOnTestParticle )(
                    perturbedBodyGravityOnTestParticle );

        // Set Cartesian state derivative models for perturbed body and test particle.
        CartesianStateDerivativeModel6dPointer perturbedBodyStateDerivative
                = make_shared< CartesianStateDerivativeModel6d >(
                    perturbedBodyAccelerationList, &updateNothing< double, Vector6d > );

        CartesianStateDerivativeModel6dPointer testParticleStateDerivative
                = make_shared< CartesianStateDerivativeModel6d >(
                    testParticleAccelerationList, &updateNothing< double, Vector6d > );

        // Construct and return composite state derivative model.
        CompositeStateDerivativeModel12d::VectorStateDerivativeModelMap stateDerivativeModelMap;
        stateDerivativeModelMap[ make_pair( 0, 6 ) ]
                = bind( &CartesianStateDerivativeModel6d::computeStateDerivative,
                        perturbedBodyStateDerivative, _1, _2 );
        stateDerivativeModelMap[ make_pair( 6, 6 ) ]
                = bind( &CartesianStateDerivativeModel6d::computeStateDerivative,
                        testParticleStateDerivative, _1, _2 );

        // Construct data updater that disassembles the composite state and updates the states of
        // the bodies.
        DataUpdaterPointer dataUpdater = make_shared< DataUpdater >( perturbedBody, testParticle );

        // Construct composite state derivative model.
        CompositeStateDerivativeModel12dPointer stateDerivativeModel
                = make_shared< CompositeStateDerivativeModel12d >(
                    stateDerivativeModelMap,
                    bind( &DataUpdater::updateTimeAndCompositeState, dataUpdater, _1, _2 ) );

        ///////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////

        // Set up integrator.

        // Declare Runge-Kutta, variable-stepsize, integrator.
        RungeKuttaVariableStepSizeIntegratorXdPointer integrator
                = make_shared< RungeKuttaVariableStepSizeIntegratorXd >(
                    rungeKuttaCoefficients,
                    bind( &CompositeStateDerivativeModel12d::computeStateDerivative,
                          stateDerivativeModel, _1, _2 ),
                    startEpoch,
                    ( VectorXd( 12 ) << perturbedBody->getCurrentState( ),
                      testParticle->getCurrentState( ) ).finished( ),
                    numeric_limits< double >::epsilon( ),
                    numeric_limits< double >::max( ),
                    testParticleCase->numericalIntegratorRelativeTolerance,
                    testParticleCase->numericalIntegratorAbsoluteTolerance );

        // Numerically integrate motion of test particle up to end of start-up period.
        integrator->integrateTo( testParticleCase->startUpIntegrationDuration,
                                 testParticleCase->initialStepSize );

        ///////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////

        // Compute system parameters at end of start-up period (T0).

        // Convert test particle's state to Keplerian elements.
        const Vector6d testParticleStateAfterStartUp = convertCartesianToKeplerianElements(
                    testParticle->getCurrentState( ),
                    testParticleCase->centralBodyGravitationalParameter );

        // Compute orbital period of perturbed body [s].
        const double orbitalPeriodOfPerturbedBody = computeKeplerOrbitalPeriod(
                    testParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                        semiMajorAxisIndex ),
                    testParticleCase->centralBodyGravitationalParameter );

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
//             cout << "Test particle orbital period [hr]: "
//                       << convertSecondsToHours( orbitalPeriodOfTestParticle ) << endl;
//             cout << "Perturbed body orbital period [hr]: "
//                       << convertSecondsToHours( orbitalPeriodOfPerturbedBody ) << endl;
//             cout << "Test particle synodic period [yr]: "
//                       << synodicPeriod / JULIAN_YEAR << endl;
//             cout << "Perturbed body's initial energy [m^2 s^-2]: "
//                       << perturbedBodyInitialEnergy << endl;
//             cout << "Perturbed body's initial angular momentum [m^2 s^-1]: "
//                       << perturbedBodyInitialAngularMomentum << endl;
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
// ////        cout << "Mab energy error: " << perturbedBodyEnergyError << endl;
// ////        cout << "Mab angular momentum error: " << perturbedBodyAngularMomentumError << endl;

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

//         cout << endl;

//         /////////////////////////////////////////////////////////////////////////

        iteratorInputTable++;
    }

    return 0;
}
