/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <algorithm>
#include <cmath>
#include <fstream> 
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <utility>

#include <Eigen/Core>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <omp.h>

#include <Assist/Astrodynamics/astrodynamicsBasics.h>
#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/commonTypedefs.h>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>
#include <Tudat/Astrodynamics/Gravitation/centralJ2GravityModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/cartesianStateDerivativeModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/compositeStateDerivativeModel.h>
#include <Tudat/InputOutput/dictionaryTools.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include "Stomi/ApplicationModes/testParticleSimulator.h"
#include "Stomi/Astrodynamics/body.h"
#include "Stomi/Astrodynamics/dataUpdater.h"
#include "Stomi/Astrodynamics/propagationDataPoint.h"
#include "Stomi/Database/databaseReadFunctions.h"
#include "Stomi/Database/databaseWriteFunctions.h"
#include "Stomi/Database/testParticleCase.h"
#include "Stomi/Database/testParticleInput.h"
#include "Stomi/Database/testParticleKick.h"
#include "Stomi/InputOutput/dictionaries.h"
#include "Stomi/Mathematics/basicMathematics.h"

namespace stomi
{
namespace application_modes
{

//! Execute test particle simulator application mode.
void executeTestParticleSimulator( 
    const std::string databasePath, 
    const tudat::input_output::parsed_data_vector_utilities::ParsedDataVectorPtr parsedData )
{

    ///////////////////////////////////////////////////////////////////////////

    // Declare using-statements and type definitions.
    
    using std::advance;
    using std::cerr;
    using std::cout;
    using std::endl;
    using std::fabs;
    using std::max_element;
    using std::min_element;
    using std::numeric_limits;
    using std::ofstream;
    using std::ostringstream;
    using std::setprecision;
    using std::string;
    using std::make_pair;

    using boost::assign::list_of;
    using boost::bind;
    using boost::iequals;
    using namespace boost::filesystem;
    using boost::make_shared;
    using boost::shared_ptr;

    using namespace Eigen;

    using namespace assist::astrodynamics;
    using namespace assist::basics;

    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_astrodynamics::unit_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::input_output::dictionary;
    using namespace tudat::numerical_integrators;
    using namespace tudat::state_derivative_models;

    using namespace stomi::astrodynamics;
    using namespace stomi::database;
    using namespace stomi::input_output;
    using namespace stomi::mathematics;    

    // Typedefs.
    typedef CompositeStateDerivativeModel< double, Vector12d, Vector6d > 
        CompositeStateDerivativeModel12d;
    typedef shared_ptr< CompositeStateDerivativeModel12d > CompositeStateDerivativeModel12dPointer;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Extract input parameters.

    // Get dictionary.
    const DictionaryPointer dictionary = getTestParticleSimulatorDictionary( );

    // Print database path to console.
    cout << "Database                                                  " 
         << databasePath << endl;

    // Extract required parameters.
    const string testParticleCaseName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), 
                findEntry( dictionary, "TESTPARTICLECASE" ) );
    cout << "Test particle case                                        " 
         << testParticleCaseName << endl; 

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

    const string testParticleSimulations = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLESIMULATIONS" ), "ALL" );
    cout << "Simulations to execute                                    "
         << testParticleSimulations << endl;
         
    const double outputInterval = extractParameterValue< double >(
            parsedData->begin( ), parsedData->end( ),
            findEntry( dictionary, "OUTPUTINTERVAL" ), 3600.0 );
    cout << "Output interval                                           " 
         << outputInterval << " s" << endl;     

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
    TestParticleCasePointer testParticleCase;

    // Test particle case data is extracted in local scope from the database, with overwritten 
    // parameters extracted from the configuration file, to ensure that none of these parameters 
    // are used globally elsewhere in this file.
    {
        const TestParticleCasePointer retrievedTestParticleCase = getTestParticleCase(
                databasePath, testParticleCaseName, testParticleCaseTableName );

        // Check if any test particle case parameters are overwritten by user input.
        const double randomWalkSimulationPeriod = extractParameterValue< double >(
            parsedData->begin( ), parsedData->end( ),
            findEntry( dictionary, "RANDOMWALKSIMULATIONPERIOD" ),
            retrievedTestParticleCase->randomWalkSimulationPeriod, &convertJulianYearsToSeconds );
        cout << "Random walk simulation period                             " 
             << convertSecondsToJulianYears( randomWalkSimulationPeriod ) << " yrs" << endl;

        const double centralBodyGravitationalParameter = extractParameterValue< double >(
            parsedData->begin( ), parsedData->end( ),
            findEntry( dictionary, "CENTRALBODYGRAVITATIONALPARAMETER" ),
            retrievedTestParticleCase->centralBodyGravitationalParameter );
        cout << "Central body gravitational parameter                      " 
             << centralBodyGravitationalParameter << " m^3 s^-2" << endl;  

        const double perturbedBodyRadius = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "PERTURBEDBODYRADIUS" ),
                    retrievedTestParticleCase->perturbedBodyRadius, 
                    &convertKilometersToMeters< double > );
        cout << "Perturbed body radius                                     " 
             << convertMetersToKilometers( perturbedBodyRadius ) << " km" << endl;

        const double perturbedBodyBulkDensity = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "PERTURBEDBODYBULKDENSITY" ),
                    retrievedTestParticleCase->perturbedBodyBulkDensity );
        cout << "Perturbed body bulk density                               " 
             << perturbedBodyBulkDensity << " kg m^-3" << endl;                

        Vector6d perturbedBodyStateInKeplerianElementsAtT0( 6 );

        perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) 
                = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "PERTURBEDBODYSEMIMAJORAXISATT0" ),
                    retrievedTestParticleCase->perturbedBodyStateInKeplerianElementsAtT0( 
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
                    retrievedTestParticleCase->perturbedBodyStateInKeplerianElementsAtT0( 
                        eccentricityIndex ) );
        cout << "Perturbed body eccentricity at TO                         "
             << perturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex ) << endl;

        perturbedBodyStateInKeplerianElementsAtT0( inclinationIndex ) 
                = extractParameterValue< double >( 
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "PERTURBEDBODYINCLINATIONATT0" ),
                    retrievedTestParticleCase->perturbedBodyStateInKeplerianElementsAtT0( 
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
                    retrievedTestParticleCase->perturbedBodyStateInKeplerianElementsAtT0( 
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
                    retrievedTestParticleCase->perturbedBodyStateInKeplerianElementsAtT0( 
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
                    retrievedTestParticleCase->perturbedBodyStateInKeplerianElementsAtT0(
                        trueAnomalyIndex ),
                     &convertDegreesToRadians< double > );
        cout << "Perturbed body true anomaly at TO                         "
             << convertRadiansToDegrees( 
                   perturbedBodyStateInKeplerianElementsAtT0( trueAnomalyIndex ) ) 
             << " deg" << endl;

        const double synodicPeriodMaximum = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "SYNODICPERIODMAXIMUM" ),
                    retrievedTestParticleCase->synodicPeriodMaximum, &convertJulianYearsToSeconds );
        cout << "Synodic period limit                                      " 
             << convertSecondsToJulianYears( synodicPeriodMaximum ) << " yrs" << endl;

        const double startUpIntegrationPeriod = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "STARTUPINTEGRATIONPERIOD" ),
                    retrievedTestParticleCase->startUpIntegrationPeriod, 
                    &convertJulianYearsToSeconds );
        cout << "Start-up integration duration                             " 
             << convertSecondsToJulianYears( startUpIntegrationPeriod ) << " yrs" << endl;

        const double centralBodyJ2GravityCoefficient = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "CENTRALBODYJ2GRAVITYCOEFFICIENT" ), 
                    retrievedTestParticleCase->centralBodyJ2GravityCoefficient );
        cout << "Central body J2 gravity coefficient                       "
             << centralBodyJ2GravityCoefficient << endl;

        const double centralBodyEquatorialRadius = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "CENTRALBODYEQUATORIALRADIUS" ),
                    retrievedTestParticleCase->centralBodyEquatorialRadius );
        cout << "Central body equatorial radius                            "
             << convertMetersToKilometers( centralBodyEquatorialRadius ) << " km" << endl;    

        const double conjunctionEventDetectionDistance = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "CONJUNCTIONEVENTDETECTIONDISTANCE" ), 
                    retrievedTestParticleCase->conjunctionEventDetectionDistance );
        cout << "Conjunction event detection distance                      " 
             << convertMetersToKilometers( conjunctionEventDetectionDistance ) << " km" << endl;

        const double oppositionEventDetectionDistance = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "OPPOSITIONEVENTDETECTIONDISTANCE" ),
                    retrievedTestParticleCase->oppositionEventDetectionDistance );
        cout << "Opposition event detection distance                       " 
             << convertMetersToKilometers( oppositionEventDetectionDistance ) << " km" << endl;

        string defaultNumericalIntegratorType;

        if ( retrievedTestParticleCase->numericalIntegratorType == DOPRI853 )
        {
            defaultNumericalIntegratorType = "DOPRI853";
        }

        else if ( retrievedTestParticleCase->numericalIntegratorType == RKF78 )
        {
            defaultNumericalIntegratorType = "RKF78";
        }

        string numericalIntegratorType = extractParameterValue< string >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "NUMERICALINTEGRATORTYPE" ), 
                    defaultNumericalIntegratorType );
        cout << "Numerical integrator type                                 "
             << numericalIntegratorType << endl;

        const double numericalIntegratorInitialStepSize = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "INITIALSTEPSIZE" ), 
                    retrievedTestParticleCase->numericalIntegratorInitialStepSize );
        cout << "Initial step size                                         "
             << numericalIntegratorInitialStepSize << " s" << endl;

        const double numericalIntegratorRelativeTolerance = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "RUNGEKUTTARELATIVEERRORTOLERANCE" ),
                    retrievedTestParticleCase->numericalIntegratorRelativeTolerance );
        cout << "Numerical integrator relative tolerance                   " 
             << numericalIntegratorRelativeTolerance << endl;

        const double numericalIntegratorAbsoluteTolerance = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "RUNGEKUTTAABSOLUTEERRORTOLERANCE" ), 
                    retrievedTestParticleCase->numericalIntegratorAbsoluteTolerance );
        cout << "Numerical integrator absolute tolerance                   " 
             << numericalIntegratorAbsoluteTolerance << endl;             

    // Store case data with possible overwritten data.
    testParticleCase = make_shared< TestParticleCase >(
        TestParticleCase( 
            retrievedTestParticleCase->testParticleCaseId, testParticleCaseName, 
            randomWalkSimulationPeriod, centralBodyGravitationalParameter, perturbedBodyRadius, 
            perturbedBodyBulkDensity, perturbedBodyStateInKeplerianElementsAtT0, 
            retrievedTestParticleCase->semiMajorAxisDistributionLimit,
            synodicPeriodMaximum, startUpIntegrationPeriod, centralBodyJ2GravityCoefficient, 
            centralBodyEquatorialRadius, conjunctionEventDetectionDistance, 
            oppositionEventDetectionDistance, 
            retrievedTestParticleCase->eccentricityDistributionMean, 
            retrievedTestParticleCase->eccentricityDistributionStandardDeviation,
            retrievedTestParticleCase->inclinationDistributionMean, 
            retrievedTestParticleCase->inclinationDistributionStandardDeviation, 
            numericalIntegratorType, numericalIntegratorInitialStepSize, 
            numericalIntegratorRelativeTolerance, numericalIntegratorRelativeTolerance ) );
    }

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
    cout << "Start epoch for numerical integration (T0)                " 
         << startEpoch << " yrs" << endl;    

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
    // the requested test particle simulation IDs.
    TestParticleInputTable testParticleInputTable;

    if ( iequals( testParticleSimulations, "ALL" )  )
    {
        cout << "Fetching all incomplete test particle simulations ..." << endl;    

        // Get entire test particle input table from database.
        testParticleInputTable = getCompleteTestParticleInputTable(
             databasePath, testParticleCase->testParticleCaseId, testParticleInputTableName );
    }

    else
    {
        cout << "Fetching all requested test particle simulations ..." << endl;    

        // Get selected test particle input table from database.
        testParticleInputTable = getSelectedTestParticleInputTable(
            databasePath, testParticleCase->testParticleCaseId, 
            testParticleSimulations, testParticleInputTableName );
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
    cout << testParticleInputTable.size( ) << " test particle simulations queued for execution ..."
         << endl;
    cout << endl;

#pragma omp parallel for num_threads( numberOfThreads )
    for ( unsigned int i = 0; i < testParticleInputTable.size( ); i++ )
    {
        /////////////////////////////////////////////////////////////////////////

        // Set input table iterator and emit output message.

        // Set input table iterator for current simulation wrt to start of input table and counter.
        TestParticleInputTable::iterator iteratorInputTable = testParticleInputTable.begin( );
        advance( iteratorInputTable, i );

        // Emit output message.
#pragma omp critical( outputToConsole )
        {
            cout << "Test particle simulation ID " << iteratorInputTable->testParticleSimulationId
                 << " on thread " << omp_get_thread_num( ) + 1 << " / " << omp_get_num_threads( ) 
                 << endl;
        }

        ///////////////////////////////////////////////////////////////////////////

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

        // Set up integrator and perform start-up integration.

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
        // Set a flag that indicates whether the start-up period is non-zero, i.e., numerical
        // integration was performed.
        bool isStartup = false;

        if ( integrator->getCurrentIndependentVariable( ) 
                < testParticleCase->startUpIntegrationPeriod )
        {
            isStartup = true;
        }

        while ( integrator->getCurrentIndependentVariable( ) 
                < testParticleCase->startUpIntegrationPeriod )
        {
            integrator->performIntegrationStep( 
                testParticleCase->numericalIntegratorInitialStepSize );
        }   

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

        // Compute orbital period of test particle [s].
        const double orbitalPeriodOfTestParticle = computeKeplerOrbitalPeriod(
                    testParticleStateAfterStartUp( semiMajorAxisIndex ),
                    testParticleCase->centralBodyGravitationalParameter );

        // Compute synodic period of test particle's motion with respect to perturbed body [s].
        const double synodicPeriod = computeSynodicPeriod(
                    orbitalPeriodOfPerturbedBody, orbitalPeriodOfTestParticle );

        ///////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////

        // Numerically integrate system to TPlusSimulationAndSynodicPeriod.
        // Write output to file if the OUTPUTMODE variable is set to "FILE".

        // Set the next step size, depending on whether the start-up integration was performed.
        double nextStepSize = 0.0;

        // Check if the start-up integration period was performed. If so, set the next step for the
        // numerical integrator based on the internally computed value.
        if ( isStartup ) { nextStepSize = integrator->getNextStepSize( ); }

        // Else, set it to the initial step size specified in the case data.
        else { nextStepSize = testParticleCase->numericalIntegratorInitialStepSize; }

        // Propagate test-particle-perturbed-body system and retrieve table of kicks experienced by
        // the test particle.

        // Set time shift (this is the amount of time by which all the epochs have to be shifted
        // so that the kicks lie in the range [0.0, randomWalkSimulationPeriod]).
        const double timeShift = -( integrator->getCurrentIndependentVariable( ) + synodicPeriod );

        // Compute mutual distance between test particle and perturbed body [m].
        double mutualDistance = ( perturbedBody->getCurrentPosition( ) 
                - testParticle->getCurrentPosition( ) ).norm( );

        // Integrate until the first opposition event is detected.
        while ( mutualDistance < testParticleCase->oppositionEventDetectionDistance )
        {
            // Integrate one step forwards.
            integrator->performIntegrationStep( nextStepSize );

            // Compute mutual distance based on new test particle and perturbed body states.
            mutualDistance = ( perturbedBody->getCurrentPosition( ) 
                - testParticle->getCurrentPosition( ) ).norm( );  

            // Update next step size.
            nextStepSize = integrator->getNextStepSize( );  
        }        

        // Create a table of propagation data points.
        PropagationDataPointTable dataPoints;

        // Set data point table iterator.
        PropagationDataPointTable::iterator iteratorDataPoint;                   

        // Create tables of opposition and conjunction events.
        PropagationDataPointTable oppositionEvents;
        PropagationDataPointTable conjunctionEvents; 

        // Set up output files.
        ofstream mutualDistanceFile;
        ofstream keplerianElementsFile;
        ofstream cartesianElementsFile;
        ofstream keplerianElementsPerturbedBodyFile;
        ofstream cartesianElementsPerturbedBodyFile;            

        // Check if output mode is set to "FILE".
        // If so, open output files and write header content.
        // Check if the output directory exists: if not, create it.
        if ( iequals( outputMode, "FILE" ) )
        {
            // Check if output directory exists.
            if ( !exists( fileOutputDirectory ) )
            {
               cerr << "Directory does not exist. Will be created." << endl;
               create_directories( fileOutputDirectory );
            }

            ostringstream mutualDistanceFilename;
            mutualDistanceFilename << fileOutputDirectory << "testParticleSimulation" 
                                   << iteratorInputTable->testParticleSimulationId
                                   << "_mutualDistance.csv";            
            mutualDistanceFile.open( mutualDistanceFilename.str( ).c_str( ) );
            mutualDistanceFile << "epoch,mutualDistance" << endl;
            mutualDistanceFile << "# [s],[m]" << endl;

            ostringstream keplerianElementsFilename;
            keplerianElementsFilename << fileOutputDirectory << "testParticleSimulation" 
                                      << iteratorInputTable->testParticleSimulationId 
                                      << "_keplerianElements.csv";
            keplerianElementsFile.open( keplerianElementsFilename.str( ).c_str( ) );
            keplerianElementsFile << "epoch,semiMajorAxis,eccentricity,inclination,"
                                  << "argumentofPeriapsis,longitudeOfAscendingNode,trueAnomaly"
                                  << endl;
            keplerianElementsFile << "# [s],[m],[-],[deg],[deg],[deg],[deg]" << endl;

            ostringstream cartesianElementsFilename;
            cartesianElementsFilename << fileOutputDirectory << "testParticleSimulation" 
                                      << iteratorInputTable->testParticleSimulationId 
                                      << "_cartesianElements.csv";
            cartesianElementsFile.open( cartesianElementsFilename.str( ).c_str( ) );
            cartesianElementsFile << "epoch,xPosition,yPosition,zPosition,"
                                  << "xVelocity,yVelocity,zVelocity"
                                  << endl;
            cartesianElementsFile << "# [s],[m],[m],[m],[m s^-1],[m s^-1],[m s^-1]" << endl;     

            ostringstream keplerianElementsPerturbedBodyFilename;
            keplerianElementsPerturbedBodyFilename << fileOutputDirectory 
                                                   << "testParticleSimulation" 
                                                   << iteratorInputTable->testParticleSimulationId 
                                                   << "_keplerianElementsPerturbedBody.csv";
            keplerianElementsPerturbedBodyFile.open( 
                keplerianElementsPerturbedBodyFilename.str( ).c_str( ) );
            keplerianElementsPerturbedBodyFile 
                << "epoch,semiMajorAxis,eccentricity,inclination,"
                << "argumentofPeriapsis,longitudeOfAscendingNode,trueAnomaly" << endl;
            keplerianElementsPerturbedBodyFile << "# [s],[m],[-],[deg],[deg],[deg],[deg]" 
                                               << endl;

            ostringstream cartesianElementsPerturbedBodyFilename;
            cartesianElementsPerturbedBodyFilename << fileOutputDirectory 
                                                   << "testParticleSimulation" 
                                                   << iteratorInputTable->testParticleSimulationId 
                                                   << "_cartesianElementsPerturbedBody.csv";
            cartesianElementsPerturbedBodyFile.open( 
                cartesianElementsPerturbedBodyFilename.str( ).c_str( ) );
            cartesianElementsPerturbedBodyFilename << "epoch,xPosition,yPosition,zPosition,"
                                                   << "xVelocity,yVelocity,zVelocity" << endl;
            cartesianElementsPerturbedBodyFilename 
                << "# [s],[m],[m],[m],[m s^-1],[m s^-1],[m s^-1]" << endl;
        }

        // Set output counter.
        unsigned int outputCounter = 0;        

        // Set flag indicating if an opposition event has been detected to true.
        bool isOppositionEventDetected = true;

        // Integrate until the end of the simulation period and detect local maxima and minima.
        while ( integrator->getCurrentIndependentVariable( )
                < testParticleCase->startUpIntegrationPeriod
                + testParticleCase->randomWalkSimulationPeriod + 2.0 * synodicPeriod )
        {            
            // Check if the current event detected is an opposition event and the end of opposition
            // event is detected (i.e., start of conjunction event).
            if ( isOppositionEventDetected 
                 && mutualDistance < testParticleCase->conjunctionEventDetectionDistance )
            {
                // Set flag to false.
                isOppositionEventDetected = false;

                // Find data point at opposition event.
                iteratorDataPoint = max_element( 
                    dataPoints.begin( ), dataPoints.end( ), compareMutualDistances );

                // Add opposition event to table.
                oppositionEvents.insert( PropagationDataPoint( *iteratorDataPoint ) );

                // Clear all data points from table.
                dataPoints.clear( );    
            }

            // Else, the current event must be a conjunction event. Check if end of conjunction 
            // event is detected (i.e., start of opposition event).
            else if ( !isOppositionEventDetected
                      && mutualDistance > testParticleCase->oppositionEventDetectionDistance )
            {
                // Set flag to true.
                isOppositionEventDetected = true;

                // Find data point at conjunction event.
                iteratorDataPoint = min_element( 
                    dataPoints.begin( ), dataPoints.end( ), compareMutualDistances );

                // Add conjunction event to table.
                conjunctionEvents.insert( PropagationDataPoint( *iteratorDataPoint ) );            

                // Clear all data points from table.
                dataPoints.clear( );     
            }

            // Integrate one step forward and store data in table.
            integrator->performIntegrationStep( nextStepSize );

            // Update next step size.
            nextStepSize = integrator->getNextStepSize( );

            // Compute mutual distance between test particle and perturbed body.
            mutualDistance = ( perturbedBody->getCurrentPosition( ) 
                - testParticle->getCurrentPosition( ) ).norm( );

            // Add current data point to table.
            dataPoints.insert( PropagationDataPoint( 
                testParticle->getCurrentTime( ) + timeShift,
                mutualDistance,
                convertCartesianToKeplerianElements( 
                    testParticle->getCurrentState( ),
                    testParticleCase->centralBodyGravitationalParameter ),
                convertCartesianToKeplerianElements( 
                    perturbedBody->getCurrentState( ),
                    testParticleCase->centralBodyGravitationalParameter ) ) );

            // Advance data point iterator.
            if ( dataPoints.size( ) == 1 )
            {
                iteratorDataPoint = dataPoints.begin( );
            }

            else
            {
                iteratorDataPoint++;
            }

            // Check if output mode is set to "FILE".
            // If so, write mutual distance and Keplerian elements to file.
            if ( iequals( outputMode, "FILE" ) )
            {
                // Only write data points that are spaced by more than specified by the 
                // outputInterval variable.
                if ( iteratorDataPoint->epoch > outputInterval * outputCounter - synodicPeriod )
                {
                    mutualDistanceFile 
                         << setprecision( numeric_limits< double >::digits10 )
                         << iteratorDataPoint->epoch << "," 
                         << iteratorDataPoint->mutualDistance << endl;

                    keplerianElementsFile 
                         << setprecision( numeric_limits< double >::digits10 )
                         << iteratorDataPoint->epoch << "," 
                         << iteratorDataPoint->testParticleStateInKeplerianElements( 
                                semiMajorAxisIndex ) << ","
                         << iteratorDataPoint->testParticleStateInKeplerianElements( 
                                eccentricityIndex ) << ","
                         << iteratorDataPoint->testParticleStateInKeplerianElements( 
                                inclinationIndex ) << ","
                         << iteratorDataPoint->testParticleStateInKeplerianElements( 
                                argumentOfPeriapsisIndex ) << ","
                         << iteratorDataPoint->testParticleStateInKeplerianElements( 
                                longitudeOfAscendingNodeIndex ) << ","
                         << iteratorDataPoint->testParticleStateInKeplerianElements( 
                                trueAnomalyIndex ) << endl; 

                    cartesianElementsFile 
                         << setprecision( numeric_limits< double >::digits10 )
                         << iteratorDataPoint->epoch << "," 
                         << testParticle->getCurrentState( )( xPositionIndex ) << ","
                         << testParticle->getCurrentState( )( yPositionIndex ) << ","
                         << testParticle->getCurrentState( )( zPositionIndex ) << ","
                         << testParticle->getCurrentState( )( xVelocityIndex ) << ","
                         << testParticle->getCurrentState( )( yVelocityIndex ) << ","
                         << testParticle->getCurrentState( )( zVelocityIndex ) << endl;                                            

                    keplerianElementsPerturbedBodyFile 
                         << setprecision( numeric_limits< double >::digits10 )
                         << iteratorDataPoint->epoch << "," 
                         << iteratorDataPoint->perturbedBodyStateInKeplerianElements( 
                                semiMajorAxisIndex ) << ","
                         << iteratorDataPoint->perturbedBodyStateInKeplerianElements( 
                                eccentricityIndex ) << ","
                         << iteratorDataPoint->perturbedBodyStateInKeplerianElements( 
                                inclinationIndex ) << ","
                         << iteratorDataPoint->perturbedBodyStateInKeplerianElements( 
                                argumentOfPeriapsisIndex ) << ","
                         << iteratorDataPoint->perturbedBodyStateInKeplerianElements( 
                                longitudeOfAscendingNodeIndex ) << ","
                         << iteratorDataPoint->perturbedBodyStateInKeplerianElements( 
                                trueAnomalyIndex ) << endl;        

                    cartesianElementsPerturbedBodyFile 
                         << setprecision( numeric_limits< double >::digits10 )
                         << iteratorDataPoint->epoch << "," 
                         << testParticle->getCurrentState( )( xPositionIndex ) << ","
                         << testParticle->getCurrentState( )( yPositionIndex ) << ","
                         << testParticle->getCurrentState( )( zPositionIndex ) << ","
                         << testParticle->getCurrentState( )( xVelocityIndex ) << ","
                         << testParticle->getCurrentState( )( yVelocityIndex ) << ","
                         << testParticle->getCurrentState( )( zVelocityIndex ) << endl; 

                    outputCounter++;
                }               
            }
        }         

        // Close output files if output mode is set to "FILE".
        if ( iequals( outputMode, "FILE" ) )
        {
            mutualDistanceFile.close( ); 
            keplerianElementsFile.close( );
            cartesianElementsFile.close( );
            keplerianElementsPerturbedBodyFile.close( );
            cartesianElementsPerturbedBodyFile.close( );               
        }
      
        // Check if there is a problem with the detection algorithm by checking if either the list
        // of conjunction of opposition events is empty.
        // In case this is true, emit a warning message and skip the current simulation.
        if ( conjunctionEvents.size( ) == 0 || oppositionEvents.size( ) == 0 ) 
        {
             cerr << "WARNING: Zero conjunction or opposition events were detected!" << endl;
        }      

        // Check if epoch of last conjunction event is greater than that of the last opposition 
        // event.
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
        // window [0, randomWalkSimulationPeriod].
        // If so, delete the last conjunction event and the last opposition event.
        // Repeat until the last conjunction epoch is within the window.
        while ( conjunctionEvents.rbegin( )->epoch > testParticleCase->randomWalkSimulationPeriod )
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
        // window [0.0, randomWalkSimulationPeriod].
        // If so, delete the first conjunction event and the first opposition event.
        // Repeat until the first epoch is within the window.
        while ( conjunctionEvents.begin( )->epoch < 0.0 )
        {
            // Erase first element.
            conjunctionEvents.erase( conjunctionEvents.begin( ) );

            // Erase last element.
            oppositionEvents.erase( oppositionEvents.begin( ) );
        }    

        // Check if output mode is set to "FILE".
        // If so, generate output files containing opposition and conjunction event data.
        if ( iequals( outputMode, "FILE" ) )
        {
            // Set up and populate opposition events output file.
            ostringstream oppositionEventsFilename;
            oppositionEventsFilename << fileOutputDirectory << "testParticleSimulation" 
                                     << iteratorInputTable->testParticleSimulationId 
                                     << "_oppositionEvents.csv";
            ofstream oppositionEventsFile( oppositionEventsFilename.str( ).c_str( ) );
            oppositionEventsFile << "epoch,mutualDistance" << endl;
            oppositionEventsFile << "# [s],[m]" << endl;        

            for ( PropagationDataPointTable::iterator iteratorOppositionEvents 
                  = oppositionEvents.begin( );
                  iteratorOppositionEvents != oppositionEvents.end( );
                  iteratorOppositionEvents++ )
            {
                oppositionEventsFile 
                    << setprecision( numeric_limits< double >::digits10 )
                    << iteratorOppositionEvents->epoch << "," 
                    << iteratorOppositionEvents->mutualDistance << endl;
            }

            oppositionEventsFile.close( );

            // Set up and populate conjunction events output file.
            ostringstream conjunctionEventsFilename;
            conjunctionEventsFilename << fileOutputDirectory << "testParticleSimulation" 
                                      << iteratorInputTable->testParticleSimulationId 
                                      << "_conjunctionEvents.csv";
            ofstream conjunctionEventsFile( conjunctionEventsFilename.str( ).c_str( ) );
            conjunctionEventsFile << "epoch,mutualDistance" << endl;
            conjunctionEventsFile << "# [s],[m]" << endl;        

            for ( PropagationDataPointTable::iterator iteratorConjunctionEvents 
                  = conjunctionEvents.begin( );
                  iteratorConjunctionEvents != conjunctionEvents.end( );
                  iteratorConjunctionEvents++ )
            {
                conjunctionEventsFile 
                    << setprecision( numeric_limits< double >::digits10 )
                    << iteratorConjunctionEvents->epoch << "," 
                    << iteratorConjunctionEvents->mutualDistance << endl;
            }

            conjunctionEventsFile.close( );
        }

        ///////////////////////////////////////////////////////////////////////////

        // Generate test particle kick table from conjunction and opposition events.

        // Declare test particle kick table.
        TestParticleKickTable testParticleKickTable;

        // Loop through tables of conjunction and opposition events to generate entries in test 
        // particle kick table.
        PropagationDataPointTable::iterator iteratorOppositionEventBefore 
            = oppositionEvents.begin( );
        PropagationDataPointTable::iterator iteratorOppositionEventAfter 
            = oppositionEvents.begin( );
        iteratorOppositionEventAfter++;    

        for ( PropagationDataPointTable::iterator iteratorConjunctionEvents 
              = conjunctionEvents.begin( );
              iteratorConjunctionEvents != conjunctionEvents.end( );
              iteratorConjunctionEvents++ )
        {
            // Add new test particle kick to table.
            testParticleKickTable.insert( new TestParticleKick( 
                0, iteratorInputTable->testParticleSimulationId, 
                iteratorConjunctionEvents->epoch, iteratorConjunctionEvents->mutualDistance,  
                iteratorOppositionEventBefore->epoch,
                iteratorOppositionEventBefore->mutualDistance,
                iteratorOppositionEventBefore->testParticleStateInKeplerianElements,
                iteratorOppositionEventAfter->epoch,
                iteratorOppositionEventAfter->mutualDistance,
                iteratorOppositionEventAfter->testParticleStateInKeplerianElements ) );

            // Advance iterators.
            iteratorOppositionEventBefore = iteratorOppositionEventAfter;
            iteratorOppositionEventAfter++;
        }

        // Write test particle kick table to database or file based on OUTPUTMODE parameter value.

        // Check if output mode is set to "FILE".
        if ( iequals( outputMode, "FILE" ) )
        {
            // Set up test particle kick table output file.
            ostringstream testParticleKickTableFilename;
            testParticleKickTableFilename << fileOutputDirectory << "testParticleSimulation" 
                              << iteratorInputTable->testParticleSimulationId << "_kickTable.csv";
            ofstream testParticleKickTableFile( testParticleKickTableFilename.str( ).c_str( ) );
            testParticleKickTableFile 
                << "conjunctionEpoch,conjunctionDistance,preConjunctionEpoch,"
                << "preConjunctionDistance,preConjunctionSemiMajorAxis,"
                << "preConjunctionEccentricity,preConjunctionInclination," 
                << "preConjunctionArgumentOfPeriapsis, "
                << "preConjunctionLongitudeOfAscendingNode,preConjunctionTrueAnomaly, "
                << "postConjunctionEpoch,postConjunctionDistance,"
                << "postConjunctionSemiMajorAxis,postConjunctionEccentricity,"
                << "postConjunctionInclination,postConjunctionArgumentOfPeriapsis, "
                << "postConjunctionLongitudeOfAscendingNode,postConjunctionTrueAnomaly, "
                << endl;
            testParticleKickTableFile << "# [s],[m],[s],[m],[m],[-],[rad],[rad],[rad],[rad], "
                                      << "[s],[m],[m],[-],[rad],[rad],[rad],[rad]" << endl;        

            // Loop through test particle kick table and write data to file.
            for ( TestParticleKickTable::iterator iteratorKickTable 
                    = testParticleKickTable.begin( );
                  iteratorKickTable != testParticleKickTable.end( );
                  iteratorKickTable++ )
            {
                testParticleKickTableFile 
                    << setprecision( numeric_limits< double >::digits10 )
                    << iteratorKickTable->conjunctionEpoch << "," 
                    << iteratorKickTable->conjunctionDistance << ","
                    << iteratorKickTable->preConjunctionEpoch << ","
                    << iteratorKickTable->preConjunctionDistance << ","
                    << iteratorKickTable->preConjunctionStateInKeplerianElements( 
                                            semiMajorAxisIndex ) << ","
                    << iteratorKickTable->preConjunctionStateInKeplerianElements( 
                                            eccentricityIndex ) << ","
                    << iteratorKickTable->preConjunctionStateInKeplerianElements( 
                                            inclinationIndex ) << ","
                    << iteratorKickTable->preConjunctionStateInKeplerianElements( 
                                            argumentOfPeriapsisIndex ) << ","
                    << iteratorKickTable->preConjunctionStateInKeplerianElements( 
                                            longitudeOfAscendingNodeIndex ) << ","
                    << iteratorKickTable->preConjunctionStateInKeplerianElements( 
                                            trueAnomalyIndex ) << ","
                    << iteratorKickTable->postConjunctionEpoch << ","
                    << iteratorKickTable->postConjunctionDistance << ","
                    << iteratorKickTable->postConjunctionStateInKeplerianElements( 
                                            semiMajorAxisIndex ) << ","
                    << iteratorKickTable->postConjunctionStateInKeplerianElements( 
                                            eccentricityIndex ) << ","
                    << iteratorKickTable->postConjunctionStateInKeplerianElements( 
                                            inclinationIndex ) << ","
                    << iteratorKickTable->postConjunctionStateInKeplerianElements( 
                                            argumentOfPeriapsisIndex ) << ","
                    << iteratorKickTable->postConjunctionStateInKeplerianElements( 
                                            longitudeOfAscendingNodeIndex ) << ","
                    << iteratorKickTable->postConjunctionStateInKeplerianElements( 
                                            trueAnomalyIndex ) << endl;
            }

            // Close output file.
            testParticleKickTableFile.close( );
        }

        // Else, check if output mode is set to "DATABASE".
        else if ( iequals( outputMode, "DATABASE" ) )
        {
            // To avoid locking of the database, this section is thread-critical, so will be 
            // executed one-by-one by multiple threads.
#pragma omp critical( writeKickTableToDatabase )
            {
                // Write kick table to database.
                populateTestParticleKickTable( databasePath, 
                    iteratorInputTable->testParticleSimulationId, testParticleKickTable, 
                    testParticleKickTableName, testParticleInputTableName );
            }
        }

        ///////////////////////////////////////////////////////////////////////////
    } // simulation for-loop

    ///////////////////////////////////////////////////////////////////////////
}

} // namespace application_modes
} // namespace stomi
