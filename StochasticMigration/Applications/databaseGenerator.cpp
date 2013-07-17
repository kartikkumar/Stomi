/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120502    K. Kumar          File imported from old code (mabSystemCaseGenerator.cpp and
 *                                  mabSystemCaseSimulationGenerator.cpp).
 *      120515    K. Kumar          Added support for text input file.
 *      120520    K. Kumar          Replaced sqlite3pp interface with native C API for SQLite;
 *                                  tweaked database operations for performance.
 *      120529    K. Kumar          Added use of custom database API for SQLite.
 *      120607    K. Kumar          Changed semi-major axis distribution from normal to uniform.
 *      120704    K. Kumar          Updated to new input file functionality.
 *      120808    K. Kumar          Updated to new dictionary-based input file system.
 *      130217    K. Kumar          Renamed file, updated code to work with refactored code for
 *                                  StochasticMigration project.
 *      130218    K. Kumar          Updated "encounter" to "conjunction".
 *      130702    K. Kumar          Updated code to use SQLiteCpp library; added output messages.
 *      130702    K. Kumar          Completed update of code; updated table schemas for test 
 *                                  particle input and output.
 *      130704    K. Kumar          Updated case table schema; updated variable-naming.
 *      130705    K. Kumar          Updated database structure to store multiple cases in a single
 *                                  file.
 *      130715    K. Kumar          Updated test particle case table schema, split input parameters
 *                                  into required and optional categories, changed input for 
 *                                  SMALIMIT to Hill radii.
 *
 *    References
 *      Kumar, K., de Pater, I., Showalter, M.R. In prep, 2013.
 *
 *    Notes
 *      TODO:
 *          - encapsulate code in try-catch blocks to capture exceptions.
 *          - execute verification of existing case data against input parameters provided to
 *            ensure consistency of inputs and possibly warn user.
 *          - expand code to enable interactive interface for user to provide inputs and select
 *            options (capture command line user input).          
 *
 */

#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <Eigen/Core>

#include <SQLiteC++.h>

#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/Basics/basics.h>
#include <Assist/InputOutput/basicInputOutput.h>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h>
#include <Tudat/InputOutput/dictionaryTools.h>
#include <Tudat/InputOutput/fieldType.h>
#include <Tudat/InputOutput/separatedParser.h>
#include <Tudat/InputOutput/parsedDataVectorUtilities.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include "StochasticMigration/Astrodynamics/hillSphere.h"
#include "StochasticMigration/InputOutput/dictionaries.h"

//! Execute stochastic migration database generator.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    ///////////////////////////////////////////////////////////////////////////

    // Declare using-statements.
    
    using namespace boost::random;

    using namespace assist::astrodynamics;
    using namespace assist::basics;
    using namespace assist::input_output;

    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_astrodynamics::physical_constants;
    using namespace tudat::basic_astrodynamics::unit_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::basic_mathematics::mathematical_constants;
    using namespace tudat::input_output::dictionary;
    using namespace tudat::input_output;
    using namespace tudat::input_output::field_types::general;
    using namespace tudat::input_output::parsed_data_vector_utilities;

    using namespace stochastic_migration::astrodynamics;
    using namespace stochastic_migration::input_output;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up input deck.

    // Check number of input parameters is correct (the numberOfInputs variable includes the
    // application itself, so one is subtracted from this number).
    checkNumberOfInputArguments( numberOfInputs - 1 );

    // Get input parameter dictionary.
    DictionaryPointer dictionary = getDatabaseGeneratorDictionary( );

    // Read and filter input stream (this can't be declared const because the parser's parse
    // function is not const-correct at the moment).
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

    // Extract required parameters.
    const std::string caseName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "CASE" ) );
    std::cout << "Case                                                      " 
              << caseName << std::endl;

    const std::string databasePath = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ), 
                findEntry( dictionary, "DATABASEPATH" ) );
    std::cout << "Database                                                  "
              << databasePath << std::endl;

    const double numberOfSimulations = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "NUMBEROFSIMULATIONS" ) );
    std::cout << "Number of simulations                                     " 
              << numberOfSimulations << std::endl;

    const double randomWalkSimulationPeriod = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKSIMULATIONPERIOD" ),
                TUDAT_NAN, &convertJulianYearsToSeconds );
    std::cout << "Random walk simulation period                             " 
              << randomWalkSimulationPeriod / JULIAN_YEAR << " yrs" << std::endl;

    const double centralBodyGravitationalParameter = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYGRAVITATIONALPARAMETER" ) );
    std::cout << "Central body gravitational parameter                      " 
              << centralBodyGravitationalParameter << " m^3 s^-2" << std::endl;

    const double perturbedBodyRadius = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYRADIUS" ) );
    std::cout << "Perturbed body radius                                     " 
              << convertMetersToKilometers( perturbedBodyRadius ) << " km" << std::endl;

    const double perturbedBodyBulkDensity = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYBULKDENSITY" ) );
    std::cout << "Perturbed body bulk density                               " 
              << perturbedBodyBulkDensity << " kg m^-3" << std::endl;

    // Compute mass of perturbed body [kg].
    const double perturbedBodyMass = computeMassOfSphere(
                perturbedBodyRadius, perturbedBodyBulkDensity );

    // Compute perturbed body's gravitational parameter [m^3 s^-2].
    const double perturbedBodyGravitationalParameter
            = computeGravitationalParameter( perturbedBodyMass );       

    Vector6d perturbedBodyStateInKeplerianElementsAtT0( 6 );

    perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYSEMIMAJORAXISATT0" ) );
    std::cout << "Perturbed body semi-major axis at TO                      "
              << perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex )
              << " m" << std::endl;

    perturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYECCENTRICITYATT0" ) );
    std::cout << "Perturbed body eccentricity at TO                         "
              << perturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex ) << std::endl;

    perturbedBodyStateInKeplerianElementsAtT0( inclinationIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYINCLINATIONATT0" ), TUDAT_NAN,
                &convertDegreesToRadians< double > );
    std::cout << "Perturbed body inclination at TO                          "
              << convertRadiansToDegrees(
                    perturbedBodyStateInKeplerianElementsAtT0( inclinationIndex ) ) 
              << " deg" << std::endl;

    perturbedBodyStateInKeplerianElementsAtT0( argumentOfPeriapsisIndex )
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYARGUMENTOFPERIAPSISATT0" ), TUDAT_NAN,
                 &convertDegreesToRadians< double > );
    std::cout << "Perturbed body argument of periapsis at TO                "
              << convertRadiansToDegrees( 
                    perturbedBodyStateInKeplerianElementsAtT0( argumentOfPeriapsisIndex ) ) 
              << " deg" << std::endl;

    perturbedBodyStateInKeplerianElementsAtT0( longitudeOfAscendingNodeIndex )
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYLONGITUDEOFASCENDINGNODEATT0" ), TUDAT_NAN,
                &convertDegreesToRadians< double > );
    std::cout << "Perturbed body longitude of ascending node at TO          "
              << convertRadiansToDegrees( 
                    perturbedBodyStateInKeplerianElementsAtT0( longitudeOfAscendingNodeIndex ) ) 
              << " deg" << std::endl;

    perturbedBodyStateInKeplerianElementsAtT0( trueAnomalyIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYTRUEANOMALYATT0" ), TUDAT_NAN,
                &convertDegreesToRadians< double > );
    std::cout << "Perturbed body true anomaly at TO                         "
              << convertRadiansToDegrees( 
                    perturbedBodyStateInKeplerianElementsAtT0( trueAnomalyIndex ) ) 
              << " deg" << std::endl;

    const double semiMajorAxisDistributionLimit = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "SEMIMAJORAXISDISTRIBUTIONLIMIT" ), TUDAT_NAN,
                ConvertHillRadiiToMeters( 
                    centralBodyGravitationalParameter, 
                    perturbedBodyGravitationalParameter,
                    perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) ) );
    std::cout << "Semi-major axis limit                                     " 
              << semiMajorAxisDistributionLimit << " m" << std::endl;

    // Extract optional parameters (parameters that take on default values if they are not  
    // specified in the input file).
    const double synodicPeriodMaximum = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "SYNODICPERIODMAXIMUM" ),
                randomWalkSimulationPeriod, &convertJulianYearsToSeconds );
    std::cout << "Synodic period maximum                                    " 
              << synodicPeriodMaximum / JULIAN_YEAR << " yrs" << std::endl;

    const double startUpIntegrationPeriod = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "STARTUPINTEGRATIONPERIOD" ), 0.0,
                &convertJulianYearsToSeconds );
    std::cout << "Start-up integration period                               " 
              << startUpIntegrationPeriod / JULIAN_YEAR << " yrs" << std::endl;

    const double centralBodyJ2GravityCoefficient = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYJ2GRAVITYCOEFFICIENT" ), 0.0 );
    std::cout << "Central body J2 gravity coefficient                       "
              << centralBodyJ2GravityCoefficient << std::endl;

    const double centralBodyEquatorialRadius = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYEQUATORIALRADIUS" ), 0.0 );
    std::cout << "Central body equatorial radius                            "
              << convertMetersToKilometers( centralBodyEquatorialRadius ) << " km" << std::endl;

    const double conjunctionEventDetectionDistance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CONJUNCTIONEVENTDETECTIONDISTANCE" ), 
                1.5 * perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) );
    std::cout << "Conjunction event detection distance                      " 
              << convertMetersToKilometers( conjunctionEventDetectionDistance ) 
              << " km" << std::endl;

    const double oppositionEventDetectionDistance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "OPPOSITIONEVENTDETECTIONDISTANCE" ),
                0.5 * perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) );
    std::cout << "Opposition event detection distance                       " 
              << convertMetersToKilometers( oppositionEventDetectionDistance ) 
              << " km" << std::endl;

    const double eccentricityDistributionMean = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "ECCENTRICITYDISTRIBUTIONMEAN" ), 0.0 );
    std::cout << "Eccentricity distribution mean                            " 
              << eccentricityDistributionMean << std::endl;

    const double eccentricityDistributionAngle = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "ECCENTRICITYDISTRIBUTIONANGLE" ), PI / 4.0,
                &convertDegreesToRadians< double > );
    std::cout << "Eccentricity distribution angle                           " 
              << convertRadiansToDegrees( eccentricityDistributionAngle ) << " deg" << std::endl;

    const double eccentricityDistributionFullWidthHalfMaxmimum = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "ECCENTRICITYDISTRIBUTIONFULLWIDTHHALFMAXIMUM" ), 0.0 );
    std::cout << "Eccentricity distribution Full-Width Half-Maximum         " 
              << eccentricityDistributionFullWidthHalfMaxmimum << std::endl;

    const double inclinationDistributionMean = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INCLINATIONDISTRIBUTIONMEAN" ), 0.0,
                &convertDegreesToRadians< double > );
    std::cout << "Inclination distribution mean                             "
              << convertRadiansToDegrees( inclinationDistributionMean ) << " deg" << std::endl;

    const double inclinationDistributionAngle = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INCLINATIONDISTRIBUTIONANGLE" ), PI / 4.0,
                &convertDegreesToRadians< double > );
    std::cout << "Inclination distribution angle                            " 
              << convertRadiansToDegrees( inclinationDistributionAngle ) << " deg" << std::endl;

    const double inclinationDistributionFullWidthHalfMaxmimum = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INCLINATIONDISTRIBUTIONFULLWIDTHHALFMAXIMUM" ), 0.0,
                &convertDegreesToRadians< double > );
    std::cout << "Inclination distribution Full-Width Half-Maximum          " 
              << convertRadiansToDegrees( inclinationDistributionFullWidthHalfMaxmimum ) 
              << " deg" << std::endl;

    const std::string numericalIntegratorType = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "NUMERICALINTEGRATORTYPE" ), "DOPRI853" );
    std::cout << "Numerical integrator type                                 "
              << numericalIntegratorType << std::endl;

    const double initialStepSize = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INITIALSTEPSIZE" ), 60.0 );
    std::cout << "Initial step size                                         "
              << initialStepSize << " s" << std::endl;

    const double numericalIntegratorRelativeTolerance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RUNGEKUTTARELATIVEERRORTOLERANCE" ), 1.0e-12 );
    std::cout << "Numerical integrator relative tolerance                   " 
              << numericalIntegratorRelativeTolerance << std::endl;

    const double numericalIntegratorAbsoluteTolerance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RUNGEKUTTAABSOLUTEERRORTOLERANCE" ), 1.0e-15 );
    std::cout << "Numerical integrator absolute tolerance                   " 
              << numericalIntegratorAbsoluteTolerance << std::endl;

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

    const std::string randomWalkMonteCarloRunTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKMONTECARLORUNTABLENAME" ),
                "random_walk_monte_carlo_runs" );
    std::cout << "Random walk Monte Carlo run table                         "
              << randomWalkMonteCarloRunTableName << std::endl;

    const std::string randomWalkPerturberTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKPERTURBERTABLENAME" ),
                "random_walk_perturbers" );
    std::cout << "Random walk perturber table                               "
              << randomWalkPerturberTableName << std::endl;

    const std::string randomWalkOutputTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKOUTPUTTABLENAME" ), "random_walk_output" );
    std::cout << "Random walk output table                                  "
              << randomWalkOutputTableName << std::endl;

    // Check that all required parameters have been set.
    checkRequiredParameters( dictionary );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Define random number generators.

    // Define a random number generator and initialize it with a reproducible seed (current cpu
    // time).
    GlobalRandomNumberGeneratorType randomNumberGenerator = getGlobalRandomNumberGenerator( );

    // Define an uniform random number distribution for test particle mean anomaly [rad,rad].
    uniform_real_distribution< > meanAnomalyDistribution( 0.0, 2.0 * PI );

    // Define variate generator for mean anomaly values using the random number
    // generator and uniform distribution of mean anomaly.
    variate_generator< GlobalRandomNumberGeneratorType&, uniform_real_distribution< > >
            generateMeanAnomaly( randomNumberGenerator, meanAnomalyDistribution );

    // Define an uniform random number distribution for test particle semi-major axis, centered at
    // perturbed body's semi-major axis [m].
    uniform_real_distribution< > semiMajorAxisDistribution(
                -semiMajorAxisDistributionLimit 
                + perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ),
                semiMajorAxisDistributionLimit 
                + perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) );

    // Define variate generator for semi-major axis values using the random number generator
    // and uniform distribution of semi-major axis.
    variate_generator< GlobalRandomNumberGeneratorType&, uniform_real_distribution< > >
            generateSemiMajorAxis( randomNumberGenerator, semiMajorAxisDistribution );

    // Define normal random number distributions for test particle components of eccentricity
    // vector (h_e = e*cos( AoP ), k_e = e*sin( AoP ) ).
    normal_distribution< > distributionOfXComponentOfEccentricityVector(
                eccentricityDistributionMean * std::cos( eccentricityDistributionAngle ),
                convertFullWidthHalfMaximumToStandardDeviation( 
                    eccentricityDistributionFullWidthHalfMaxmimum ) );

    normal_distribution< > distributionOfYComponentOfEccentricityVector(
                eccentricityDistributionMean * std::sin( eccentricityDistributionAngle ),
                convertFullWidthHalfMaximumToStandardDeviation( 
                    eccentricityDistributionFullWidthHalfMaxmimum ) );

    // Define variate generators for h_e- and k_e-values using the random number generator
    // and normal distribution of h_e- and k_e-values.
    variate_generator< GlobalRandomNumberGeneratorType&, normal_distribution< > >
            generateXComponentOfEccentricityVector(
                randomNumberGenerator, distributionOfXComponentOfEccentricityVector );

    variate_generator< GlobalRandomNumberGeneratorType&, normal_distribution< > >
            generateYComponentOfEccentricityVector(
                randomNumberGenerator, distributionOfYComponentOfEccentricityVector );

    // Define normal random number distributions for test particle components of inclination
    // vector (h_i = i*cos( RAAN ), k_i = i*sin( RAAN ) ).
    normal_distribution< > distributionOfXComponentOfInclinationVector(
                inclinationDistributionMean * std::cos( inclinationDistributionAngle ),
                convertFullWidthHalfMaximumToStandardDeviation( 
                    inclinationDistributionFullWidthHalfMaxmimum ) );

    normal_distribution< > distributionOfYComponentOfInclinationVector(
                inclinationDistributionMean * std::sin( inclinationDistributionAngle ),
                convertFullWidthHalfMaximumToStandardDeviation( 
                    inclinationDistributionFullWidthHalfMaxmimum ) );

    // Define variate generators for h_i- and k_i-values using the random number generator
    // and normal distribution of h_i- and k_i-values.
    variate_generator< GlobalRandomNumberGeneratorType&, normal_distribution< > >
            generateXComponentOfInclinationVector(
                randomNumberGenerator, distributionOfXComponentOfInclinationVector );

    variate_generator< GlobalRandomNumberGeneratorType&, normal_distribution< > >
            generateYComponentOfInclinationVector(
                randomNumberGenerator, distributionOfYComponentOfInclinationVector );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Open (and if necessary create) database.

    std::cout << std::endl;
    std::cout << "****************************************************************************" 
              << std::endl;
    std::cout << "Database operations" << std::endl;
    std::cout << "****************************************************************************" 
              << std::endl;
    std::cout << std::endl;

    // Check if database file already exists.
    bool isDatabaseCreated = false;
    if ( boost::filesystem::exists( databasePath ) )
    {
        std::cout << "Opening existing database at " << databasePath << " ..." << std::endl;
    }

    else
    {
        isDatabaseCreated = true;
        std::cout << "WARNING: Database does not exist!" << std::endl;
        std::cout << "Creating database at " << databasePath << " ..." << std::endl;
    }

    // Create/open database.          
    SQLite::Database database( databasePath.c_str( ), SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE );

    std::cout << "SQLite database file '" << database.getFilename().c_str();
    if ( isDatabaseCreated ) { std::cout << "' created &"; }
    std::cout << " opened successfully ..." << std::endl;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up test particle case table.

    // Create table if it doesn't exist.
    if ( !database.tableExists( testParticleCaseTableName.c_str( ) ) )
    {
        std::cout << "Table '" << testParticleCaseTableName << "' does not exist ..." << std::endl;
        std::cout << "Creating table ... " << std::endl;

        std::ostringstream testParticleCaseTableCreate;
        testParticleCaseTableCreate
            << "CREATE TABLE IF NOT EXISTS " << testParticleCaseTableName << " ("
            << "\"caseId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"                
            << "\"caseName\" TEXT NOT NULL,"
            << "\"randomWalkSimulationPeriod\" REAL NOT NULL,"
            << "\"centralBodyGravitationalParameter\" REAL NOT NULL,"
            << "\"perturbedBodyRadius\" REAL NOT NULL,"
            << "\"perturbedBodyBulkDensity\" REAL NOT NULL,"
            << "\"perturbedBodySemiMajorAxisAtT0\" REAL NOT NULL,"
            << "\"perturbedBodyEccentricityAtT0\" REAL NOT NULL,"
            << "\"perturbedBodyInclinationAtT0\" REAL NOT NULL,"
            << "\"perturbedBodyArgumentOfPeriapsisAtT0\" REAL NOT NULL,"
            << "\"perturbedBodyLongitudeOfAscendingNodeAtT0\" REAL NOT NULL,"
            << "\"perturbedBodyTrueAnomalyAtT0\" REAL NOT NULL,"
            << "\"semiMajorAxisDistributionLimit\" REAL NOT NULL,"
            << "\"synodicPeriodMaximum\" REAL NOT NULL,"
            << "\"startUpIntegrationPeriod\" REAL NOT NULL,"
            << "\"centralBodyJ2GravityCoefficient\" REAL NOT NULL,"
            << "\"centralBodyEquatorialRadius\" REAL NOT NULL,"
            << "\"conjunctionEventDetectionDistance\" REAL NOT NULL,"
            << "\"oppositionEventDetectionDistance\" REAL NOT NULL,"
            << "\"eccentricityDistributionMean\" REAL NOT NULL,"
            << "\"eccentricityDistributionAngle\" REAL NOT NULL,"
            << "\"eccentricityDistributionFullWidthHalfMaxmimum\" REAL NOT NULL,"
            << "\"inclinationDistributionMean\" REAL NOT NULL,"
            << "\"inclinationDistributionAngle\" REAL NOT NULL,"
            << "\"inclinationDistributionFullWidthHalfMaxmimum\" REAL NOT NULL,"
            << "\"numericalIntegratorType\" TEXT NOT NULL,"
            << "\"initialStepSize\" REAL NOT NULL,"
            << "\"relativeTolerance\" REAL NOT NULL,"
            << "\"absoluteTolerance\" REAL NOT NULL);";

        // Execute command to create table.
        database.exec( testParticleCaseTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( testParticleCaseTableName.c_str( ) ) )
        {
            std::cout << "Table '" << testParticleCaseTableName << "' successfully created!"
                      << std::endl; 
        }

        else
        {
            std::ostringstream tableCreateError;
            tableCreateError << "Error: Creating table '" << testParticleCaseTableName 
            << "'' failed!";
            throw std::runtime_error( tableCreateError.str( ).c_str( ) );

        }
    }

    else
    {
        std::cout << "Table '" << testParticleCaseTableName 
                  << "' already exists ... skipping creating table ..." << std::endl;

    }

    // Check data present in table.

    // Check if case is already present in table.
    std::ostringstream testParticleCaseCheck;
    testParticleCaseCheck << "SELECT COUNT( * ) FROM " << testParticleCaseTableName 
                          << " WHERE \"caseName\" = \"" << caseName << "\"";
    int numberOfCaseRows = database.execAndGet( testParticleCaseCheck.str( ).c_str( ) );
    int caseId = 0;

    if ( numberOfCaseRows > 1 )
    {
        std::ostringstream numberOfCaseRowsError;
        numberOfCaseRowsError << "Error: Table '" << testParticleCaseTableName << "' contains " 
                              << numberOfCaseRows << " rows for case '" << caseName << "'!";
        throw std::runtime_error( numberOfCaseRowsError.str( ).c_str( ) );
    }

    else if ( numberOfCaseRows == 1 )
    {
        std::cout << "Table '" << testParticleCaseTableName 
                  << "' contains 1 row of data for case '" << caseName << "' ... "
                  << "skipping populating table ... " << std::endl;
    }

    // Write test particle case data to table.
    else if ( numberOfCaseRows == 0 )
    {
        std::cout << "No data present in table '" << testParticleCaseTableName 
                  << "' for case '" << caseName << "' ... " << std::endl;
        std::cout << "Populating table ... " << std::endl;

        // Create stringstream with test particle case data insert command.
        // For floating-point values, ensure the data is written to the stream at full precision.
        std::ostringstream testParticleCaseDataInsert;
        testParticleCaseDataInsert
            << "INSERT INTO " << testParticleCaseTableName << " VALUES ("
            << "NULL,"
            << "\"" << caseName << "\",";
        testParticleCaseDataInsert 
            << std::setprecision( std::numeric_limits< double >::digits10 )
            << randomWalkSimulationPeriod << ","
            << centralBodyGravitationalParameter << ","
            << perturbedBodyRadius << ","
            << perturbedBodyBulkDensity << ","
            << perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) << ","
            << perturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex ) << ","
            << perturbedBodyStateInKeplerianElementsAtT0( inclinationIndex ) << ","
            << perturbedBodyStateInKeplerianElementsAtT0( argumentOfPeriapsisIndex ) << ","
            << perturbedBodyStateInKeplerianElementsAtT0( longitudeOfAscendingNodeIndex ) << ","
            << perturbedBodyStateInKeplerianElementsAtT0( trueAnomalyIndex ) << ","            
            << semiMajorAxisDistributionLimit << ","
            << synodicPeriodMaximum << ","
            << startUpIntegrationPeriod << ","
            << centralBodyJ2GravityCoefficient << ","
            << centralBodyEquatorialRadius << ","            
            << conjunctionEventDetectionDistance << ","
            << oppositionEventDetectionDistance << ","            
            << eccentricityDistributionMean << ","
            << eccentricityDistributionAngle << ","
            << eccentricityDistributionFullWidthHalfMaxmimum << ","
            << inclinationDistributionMean << ","
            << inclinationDistributionAngle << ","
            << inclinationDistributionFullWidthHalfMaxmimum << ",";
        testParticleCaseDataInsert    
            << "\"" << numericalIntegratorType << "\",";
        testParticleCaseDataInsert
            << std::setprecision( std::numeric_limits< double >::digits10 )
            << initialStepSize << ","
            << numericalIntegratorRelativeTolerance << ","
            << numericalIntegratorAbsoluteTolerance
            << ");";

        // Insert test particle case data.
        database.exec( testParticleCaseDataInsert.str( ).c_str( ) );

        // Check that there is only one row present in the table.
        numberOfCaseRows = database.execAndGet( testParticleCaseCheck.str( ).c_str( ) );
        if ( numberOfCaseRows == 1 )
        {
            std::cout << "Table '" << testParticleCaseTableName << "' populated successfully!" 
                      << std::endl; 
        }

        else
        {
            std::ostringstream numberOfCaseRowsError;
            numberOfCaseRowsError << "Error: Table '" << testParticleCaseTableName 
                                  << "' contains " << numberOfCaseRows << " rows for case '"
                                  << caseName << "'!";
            std::runtime_error( numberOfCaseRowsError.str( ).c_str( ) );
        }
    }

    // Retrieve and output case id.
    std::ostringstream testParticleCaseId;          
    testParticleCaseId << "SELECT \"caseId\" FROM " << testParticleCaseTableName
                       << " WHERE \"caseName\" = \"" << caseName << "\"";
    caseId = database.execAndGet( testParticleCaseId.str( ).c_str( ) );
    std::cout << "Case ID is " << caseId << " for case '" << caseName << "'" << std::endl;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up test particle input table.
    if ( !database.tableExists( testParticleInputTableName.c_str( ) ) )
    {
        std::cout << "Table '" << testParticleInputTableName << "' does not exist ..." 
                  << std::endl;
        std::cout << "Creating table ... " << std::endl;

        // Create table.
        std::ostringstream testParticleInputTableCreate;
        testParticleInputTableCreate
            << "CREATE TABLE IF NOT EXISTS " << testParticleInputTableName << " ("
            << "\"simulationId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"caseId\" INTEGER NOT NULL,"
            << "\"completed\" INTEGER NOT NULL,"
            << "\"semiMajorAxis\" REAL NOT NULL,"
            << "\"eccentricity\" REAL NOT NULL,"
            << "\"inclination\" REAL NOT NULL,"
            << "\"argumentOfPeriapsis\" REAL NOT NULL,"
            << "\"longitudeOfAscendingNode\" REAL NOT NULL,"
            << "\"trueAnomaly\" REAL NOT NULL);";

        // Execute command to create table.
        database.exec( testParticleInputTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( testParticleInputTableName.c_str( ) ) )
        {
            std::cout << "Table '" << testParticleInputTableName << "' successfully created!"
                      << std::endl; 
        }

        else
        {
            std::ostringstream tableCreateError;
            tableCreateError << "Error: Creating table '" << testParticleInputTableName 
            << "'' failed!";
            throw std::runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        std::cout << "Table '" << testParticleInputTableName 
                  << "' already exists ... skipping creating table ..." << std::endl;
    }

    // Check data present in table.

    // Check how many rows are present in table.
    std::ostringstream testParticleInputRowCount;
    testParticleInputRowCount << "SELECT COUNT( * ) FROM " << testParticleInputTableName 
                              << " WHERE \"caseId\" = " << caseId;
    int inputTableRows = database.execAndGet( testParticleInputRowCount.str( ).c_str( ) );

    if ( inputTableRows > 0 )
    {
        std::cout << "Table '" << testParticleInputTableName << "' contains "
                  << inputTableRows << " rows for case '" << caseName << "' ... " << std::endl;
    }

    // Populate table.
    std::cout << "Populating input table with data for " << numberOfSimulations 
              << " new simulations ... " << std::endl;

    // Set up database transaction.
    SQLite::Transaction testParticleInputTableTransaction( database );

    // Set up test particle input table insert statement.
    std::ostringstream testParticleInputTableInsert;
    testParticleInputTableInsert << "INSERT INTO " << testParticleInputTableName
                                 << " VALUES (NULL, " << caseId << ", 0, :semiMajorAxis, "
                                 << ":eccentricity, :inclination, :argumentOfPeriapsis, "
                                 << ":longitudeOfAscendingNode, :trueAnomaly);";

    // Compile a SQL query.
    SQLite::Statement testParticleInputTableInsertQuery( 
        database, testParticleInputTableInsert.str( ).c_str( ) );

    // Generate test particle input data and populate table.
    for ( int simulationNumber = 0; simulationNumber < numberOfSimulations; simulationNumber++ )
    {
        // Set random eccentricity vector.
        const Eigen::Vector2d eccentricityVector( generateXComponentOfEccentricityVector( ),
                                                  generateYComponentOfEccentricityVector( ) );

        // Compute eccentricity [-].
        const double eccentricity = eccentricityVector.norm( );

        // Compute argument of periapsis [rad].
        const double argumentOfPeriapsis = std::atan2( eccentricityVector.y( ),
                                                       eccentricityVector.x( ) );

        // Set random inclination vector.
        const Eigen::Vector2d inclinationVector( generateXComponentOfInclinationVector( ),
                                                 generateYComponentOfInclinationVector( ) );

        // Compute inclination [rad].
        const double inclination = inclinationVector.norm( );

        // Compute longitude ascension of ascending node [rad].
        const double longitudeOfAscendingNode = std::atan2( inclinationVector.y( ),
                                                            inclinationVector.x( ) );

        // Set random mean anomaly [rad].
        const double meanAnomaly = generateMeanAnomaly( );

        // Convert mean anomaly to eccentric anomaly [rad].
        ConvertMeanAnomalyToEccentricAnomaly convertMeanAnomalyToEccentricAnomaly(
                    eccentricity, meanAnomaly );
        const double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly.convert( );

        // Convert eccentric anomaly to true anomaly.
        const double trueAnomaly = convertEccentricAnomalyToTrueAnomaly(
                    eccentricAnomaly, eccentricity );

        // Compute orbital period of perturbed body.
        const double orbitalPeriodOfPerturbedBody = computeKeplerOrbitalPeriod(
                    perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ),
                    centralBodyGravitationalParameter );

        // Set random semi-major axis [m].
        // If synodic period is too long, i.e., test particle semi-major axis is too close to
        // perturbed body, regenerate value.
        double synodicPeriodOfTestParticle = TUDAT_NAN;
        double semiMajorAxis = TUDAT_NAN;
        do
        {
            // Set random semi-major axis [m].
            semiMajorAxis = generateSemiMajorAxis( );

            // Compute orbital period of test particle [s].
            const double orbitalPeriodOfTestParticle = computeKeplerOrbitalPeriod(
                        semiMajorAxis, centralBodyGravitationalParameter );

            // Compute synodic period of test particle's motion with respect to perturbed body [s].
            synodicPeriodOfTestParticle = computeSynodicPeriod(
                        orbitalPeriodOfPerturbedBody, orbitalPeriodOfTestParticle );
        }
        while ( synodicPeriodOfTestParticle > synodicPeriodMaximum );

        // Bind values to prepared SQLite statement.
        testParticleInputTableInsertQuery.bind( ":semiMajorAxis", semiMajorAxis );
        testParticleInputTableInsertQuery.bind( ":eccentricity", eccentricity );
        testParticleInputTableInsertQuery.bind( ":inclination", inclination );
        testParticleInputTableInsertQuery.bind( ":argumentOfPeriapsis", argumentOfPeriapsis );
        testParticleInputTableInsertQuery.bind( 
            ":longitudeOfAscendingNode", longitudeOfAscendingNode );
        testParticleInputTableInsertQuery.bind( ":trueAnomaly", trueAnomaly );

        // Execute insert query.
        testParticleInputTableInsertQuery.exec( );

        // Reset query.
        testParticleInputTableInsertQuery.reset( );
    }

    // Commit transaction.
    testParticleInputTableTransaction.commit();

    ///////////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////////

    // Set up test particle kick table.
    if ( !database.tableExists( testParticleKickTableName.c_str( ) ) )
    {
        std::cout << "Table '" << testParticleKickTableName << "' does not exist ..." 
                  << std::endl;
        std::cout << "Creating table ... " << std::endl;

        // Create table.
        std::ostringstream testParticleKickTableCreate;
        testParticleKickTableCreate
            << "CREATE TABLE IF NOT EXISTS " << testParticleKickTableName << " ("
            << "\"kickId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"simulationId\" INTEGER NOT NULL,"
            << "\"conjunctionEpoch\" REAL NOT NULL,"
            << "\"conjunctionDistance\" REAL NOT NULL,"
            << "\"conjunctionDuration\" REAL NOT NULL,"
            << "\"preConjunctionEpoch\" REAL NOT NULL,"
            << "\"preconjunctionEventDetectionDistance\" REAL NOT NULL,"
            << "\"preConjunctionSemiMajorAxis\" REAL NOT NULL,"
            << "\"preConjunctionEccentricity\" REAL NOT NULL,"
            << "\"preConjunctionInclination\" REAL NOT NULL,"
            << "\"postConjunctionEpoch\" REAL NOT NULL,"
            << "\"postconjunctionEventDetectionDistance\" REAL NOT NULL,"
            << "\"postConjunctionSemiMajorAxis\" REAL NOT NULL,"
            << "\"postConjunctionEccentricity\" REAL NOT NULL,"
            << "\"postConjunctionInclination\" REAL NOT NULL);";

        // Execute command to create table.
        database.exec( testParticleKickTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( testParticleKickTableName.c_str( ) ) )
        {
            std::cout << "Table '" << testParticleKickTableName << "' successfully created!"
                      << std::endl; 
        }

        else
        {
            std::ostringstream tableCreateError;
            tableCreateError << "Error: Creating table '" << testParticleKickTableName 
            << "'' failed!";
            throw std::runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        std::cout << "Table '" << testParticleKickTableName 
                  << "' already exists ... skipping creating table ..." << std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up random walk Monte Carlo run table.
    if ( !database.tableExists( randomWalkMonteCarloRunTableName.c_str( ) ) )
    {
        std::cout << "Table '" << randomWalkMonteCarloRunTableName << "' does not exist ..." 
                  << std::endl;
        std::cout << "Creating table ... " << std::endl;

        // Create table.
        std::ostringstream randomWalkMonteCarloRunTableCreate;
        randomWalkMonteCarloRunTableCreate
            << "CREATE TABLE IF NOT EXISTS " << randomWalkMonteCarloRunTableName << " ("
            << "\"monteCarloRunId\" INTEGER PRIMARY KEY NOT NULL,"
            << "\"perturberPopulation\" INTEGER NOT NULL,"
            << "\"massDistributionType\" TEXT NOT NULL,"
            << "\"massDistributionParameter1\" REAL NOT NULL,"
            << "\"massDistributionParameter2\" REAL NOT NULL,"
            << "\"observationPeriod\" REAL NOT NULL,"
            << "\"epochWindowSize\" REAL NOT NULL,"
            << "\"numberOfEpochWindows\" REAL NOT NULL);";

        // Execute command to create table.
        database.exec( randomWalkMonteCarloRunTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( randomWalkMonteCarloRunTableName.c_str( ) ) )
        {
            std::cout << "Table '" << randomWalkMonteCarloRunTableName 
                      << "' successfully created!" << std::endl; 
        }

        else
        {
            std::ostringstream tableCreateError;
            tableCreateError << "Error: Creating table '" << randomWalkMonteCarloRunTableName 
            << "'' failed!";
            throw std::runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        std::cout << "Table '" << randomWalkMonteCarloRunTableName 
                  << "' already exists ... skipping creating table ..." << std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up random walk perturber table.
    if ( !database.tableExists( randomWalkPerturberTableName.c_str( ) ) )
    {
        std::cout << "Table '" << randomWalkPerturberTableName << "' does not exist ..." 
                  << std::endl;
        std::cout << "Creating table ... " << std::endl;

        // Create table.
        std::ostringstream randomWalkPerturberTableCreate;
        randomWalkPerturberTableCreate
            << "CREATE TABLE IF NOT EXISTS " << randomWalkPerturberTableName << " ("
            << "\"perturberId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"monteCarloRunId\" INTEGER NOT NULL,"
            << "\"testParticleSimulationId\" INTEGER NOT NULL,"
            << "\"massFactor\" REAL NOT NULL);";

        // Execute command to create table.
        database.exec( randomWalkPerturberTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( randomWalkPerturberTableName.c_str( ) ) )
        {
            std::cout << "Table '" << randomWalkPerturberTableName 
                      << "' successfully created!" << std::endl; 
        }

        else
        {
            std::ostringstream tableCreateError;
            tableCreateError << "Error: Creating table '" << randomWalkPerturberTableName 
            << "'' failed!";
            throw std::runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        std::cout << "Table '" << randomWalkPerturberTableName 
                  << "' already exists ... skipping creating table ..." << std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up random walk output table.
    if ( !database.tableExists( randomWalkOutputTableName.c_str( ) ) )
    {
        std::cout << "Table '" << randomWalkOutputTableName << "' does not exist ..." << std::endl;
        std::cout << "Creating table ... " << std::endl;

        // Create table.
        std::ostringstream randomWalkOutputTableCreate;
        randomWalkOutputTableCreate
            << "CREATE TABLE IF NOT EXISTS " << randomWalkOutputTableName << " ("
            << "\"outputId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"monteCarloRunId\" INTEGER NOT NULL,"
            << "\"maximumEccentricityChange\" REAL NOT NULL,"
            << "\"maximumLongitudeResidualChange\" REAL NOT NULL,"
            << "\"maximumInclinationChange\" REAL NOT NULL);";

        // Execute command to create table.
        database.exec( randomWalkOutputTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( randomWalkOutputTableName.c_str( ) ) )
        {
            std::cout << "Table '" << randomWalkOutputTableName 
                      << "' successfully created!" << std::endl; 
        }

        else
        {
            std::ostringstream tableCreateError;
            tableCreateError << "Error: Creating table '" << randomWalkOutputTableName 
            << "'' failed!";
            throw std::runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        std::cout << "Table '" << randomWalkOutputTableName 
                  << "' already exists ... skipping creating table ..." << std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Close database.

    // Database will be automatically closed after this statement 
    // (when object goes out of scope its destructor will be called).
    std::cout << "SQLite database file '" << database.getFilename().c_str() 
              << "' closed successfully ..." << std::endl;
    std::cout << std::endl;

    ///////////////////////////////////////////////////////////////////////////

    return 0;
}
