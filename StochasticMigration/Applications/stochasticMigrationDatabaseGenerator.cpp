/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
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
 *      130702    K. Kumar          Completed update of code.
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

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <Eigen/Core>

#include <sqlite3.h>

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

#include "StochasticMigration/InputOutput/dictionaries.h"

//! Execute stochastic migration database generator.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    // Using statements.
    using boost::iequals;
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

    using namespace stochastic_migration::input_output;

    ///////////////////////////////////////////////////////////////////////////

    // Input deck.

    // Check number of input parameters is correct (the numberOfInputs variable includes the
    // application itself, so one is subtracted from this number).
    checkNumberOfInputArguments( numberOfInputs - 1 );

    // Get input parameter dictionary.
    DictionaryPointer dictionary = getStochasticMigrationDatabaseGeneratorDictionary( );

    // Read and filter input stream (this can't be declared const because the parser's parse
    // function is not const-correct at the moment).
    std::string filteredInput = readAndFilterInputFile( inputArguments[ 1 ] );

    // Declare a separated parser.
    SeparatedParser parser( std::string( ": " ), 2, parameterName, parameterValue );

    // Parse filtered data.
    const ParsedDataVectorPtr parsedData = parser.parse( filteredInput );

    // Extract input parameters.
    const int caseNumber = extractParameterValue< int >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "CASE" ) );

    const std::string databasePath = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "DATABASE" ) );

    const double numberOfSimulations = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "NUMBEROFSIMULATIONS" ) );

    const double randomWalkDuration = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKDURATION" ),
                50.0 * JULIAN_YEAR, &convertJulianYearsToSeconds );

    const double synodicPeriodLimit = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "SYNODICPERIODLIMIT" ),
                50.0 * JULIAN_YEAR, &convertJulianYearsToSeconds );

    const double outputInterval = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "OUTPUTINTERVAL" ),
                convertHoursToSeconds( 4.0 ), &convertHoursToSeconds< double > );

    const double startUpIntegrationDuration = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "STARTUPINTEGRATIONDURATION" ), 0.0,
                &convertJulianYearsToSeconds );

    const double conjunctionEventDetectionDistance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CONJUNCTIONEVENTDETECTIONDISTANCE" ), 3.0e7 );

    const double oppositionEventDetectionDistance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "OPPOSITIONEVENTDETECTIONDISTANCE" ), 1.5e8 );

    const double centralBodyGravitationalParameter = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYGRAVITATIONALPARAMETER" ), 5.793966e15 );

    const double centralBodyJ2GravityCoefficient = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYJ2GRAVITYCOEFFICIENT" ), 0.0 );

    const double centralBodyEquatorialRadius = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYEQUATORIALRADIUS" ), 0.0 );

    const double semiMajorAxisLimit = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "SEMIMAJORAXISLIMIT" ), TUDAT_NAN,
                &convertKilometersToMeters< double > );

    const double eccentricityMean = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "ECCENTRICITYMEAN" ), 0.0 );

    const double eccentricityAngle = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "ECCENTRICITYANGLE" ), PI / 4.0,
                &convertDegreesToRadians< double > );

    const double eccentricityFWHM = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "ECCENTRICITYFWHM" ) );

    const double inclinationMean = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INCLINATIONMEAN" ), 0.0,
                &convertDegreesToRadians< double > );

    const double inclinationAngle = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INCLINATIONANGLE" ), PI / 4.0,
                &convertDegreesToRadians< double > );

    const double inclinationFWHM = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INCLINATIONFWHM" ), 0.0,
                &convertDegreesToRadians< double > );

    const double perturbedBodyRadius = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYRADIUS" ),
                1.2e4, &convertKilometersToMeters< double > );

    const double perturbedBodyBulkDensity = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYBULKDENSITY" ), 1500.0 );

    Vector6d perturbedBodyKeplerianElementsAtT0( 6 );

    perturbedBodyKeplerianElementsAtT0( semiMajorAxisIndex ) = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYSEMIMAJORAXISATT0" ), 9.7736e7,
                &convertKilometersToMeters< double > );

    perturbedBodyKeplerianElementsAtT0( eccentricityIndex ) = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYECCENTRICITYATT0" ), 0.00254 );

    perturbedBodyKeplerianElementsAtT0( inclinationIndex ) = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYINCLINATIONATT0" ),
                convertDegreesToRadians( 0.14 ), &convertDegreesToRadians< double > );

    perturbedBodyKeplerianElementsAtT0( argumentOfPeriapsisIndex )
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYARGUMENTOFPERIAPSISATT0" ),
                convertDegreesToRadians( 18.9594 ), &convertDegreesToRadians< double > );

    perturbedBodyKeplerianElementsAtT0( longitudeOfAscendingNodeIndex )
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYLONGITUDEOFASCENDINGNODEATT0" ),
                convertDegreesToRadians( 251.932 ), &convertDegreesToRadians< double > );

    perturbedBodyKeplerianElementsAtT0( trueAnomalyIndex ) = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYTRUEANOMALYATT0" ),
                convertDegreesToRadians( 354.516 ), &convertDegreesToRadians< double > );

    const std::string numericalIntegratorType = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "NUMERICALINTEGRATORTYPE" ), "DOPRI853" );

    const double initialStepSize = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INITIALSTEPSIZE" ), 60.0 );

    const double integratorRelativeTolerance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RUNGEKUTTARELATIVEERRORTOLERANCE" ), 1.0e-12 );

    const double integratorAbsoluteTolerance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RUNGEKUTTAABSOLUTEERRORTOLERANCE" ), 1.0e-15 );

    const std::string testParticleCaseTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLECASETABLENAME" ), "test_particle_case" );

    const std::string testParticleInputTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLEINPUTTABLENAME" ), "test_particle_input" );

    const std::string testParticleKickTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLEKICKTABLENAME" ), "test_particle_kicks" );

    const std::string randomWalkMonteCarloRunTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKMONTECARLORUNTABLENAME" ),
                "random_walk_monte_carlo_runs" );

    const std::string randomWalkPerturberTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKPERTURBERTABLENAME" ),
                "random_walk_perturbers" );

    const std::string randomWalkOutputTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKOUTPUTTABLENAME" ), "random_walk_output" );

    // Compute mass of perturbed body [kg].
    const double perturbedBodyMass = computeMassOfSphere(
                perturbedBodyRadius, perturbedBodyBulkDensity );

    // Compute perturbed body's gravitational parameter [m^3 s^-2].
    const double perturbedBodyGravitationalParameter
            = computeGravitationalParameter( perturbedBodyMass );

    std::cout << std::endl;
    std::cout << "****************************************************************************" 
              << std::endl;
    std::cout << "Input parameters provided" << std::endl;
    std::cout << "****************************************************************************" 
              << std::endl;
    std::cout << std::endl;
    std::cout << "Case                              " << caseNumber << std::endl;
    std::cout << "Database                          " << databasePath << std::endl;
    std::cout << "Number of simulations             " << numberOfSimulations << std::endl;
    std::cout << "Random walk duration              " << randomWalkDuration / JULIAN_YEAR << " yrs" 
              << std::endl;
    std::cout << "Synodic period limit              " << synodicPeriodLimit / JULIAN_YEAR << " yrs"
              << std::endl;
    std::cout << "Output interval                   " << convertSecondsToHours( outputInterval ) 
              << " hrs" << std::endl;
    std::cout << "Start-up integration duration     " << startUpIntegrationDuration / JULIAN_YEAR
              << " yrs" << std::endl;
    std::cout << conjunctionEventDetectionDistance << std::endl;
    std::cout << oppositionEventDetectionDistance << std::endl;
    std::cout << centralBodyGravitationalParameter << std::endl;
    std::cout << centralBodyJ2GravityCoefficient << std::endl;
    std::cout << centralBodyEquatorialRadius << std::endl;
    std::cout << semiMajorAxisLimit << std::endl;
    std::cout << eccentricityFWHM << std::endl;
    std::cout << eccentricityMean << std::endl;
    std::cout << eccentricityAngle << std::endl;
    std::cout << inclinationFWHM << std::endl;
    std::cout << inclinationMean << std::endl;
    std::cout << inclinationAngle << std::endl;
    std::cout << perturbedBodyRadius << std::endl;
    std::cout << perturbedBodyBulkDensity << std::endl;
    std::cout << perturbedBodyKeplerianElementsAtT0 << std::endl;
    std::cout << numericalIntegratorType << std::endl;
    std::cout << initialStepSize << std::endl;
    std::cout << integratorRelativeTolerance << ", "
              << integratorAbsoluteTolerance << std::endl;
    std::cout << perturbedBodyMass << std::endl;
    std::cout << perturbedBodyGravitationalParameter << std::endl;
    std::cout << testParticleCaseTableName << std::endl;
    std::cout << testParticleInputTableName << std::endl;
    std::cout << testParticleKickTableName << std::endl;
    std::cout << randomWalkMonteCarloRunTableName << std::endl;
    std::cout << randomWalkPerturberTableName << std::endl;
    std::cout << randomWalkOutputTableName << std::endl;
    std::cout << std::endl;

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
                -semiMajorAxisLimit + perturbedBodyKeplerianElementsAtT0( semiMajorAxisIndex ),
                semiMajorAxisLimit + perturbedBodyKeplerianElementsAtT0( semiMajorAxisIndex ) );

    // Define variate generator for semi-major axis values using the random number generator
    // and uniform distribution of semi-major axis.
    variate_generator< GlobalRandomNumberGeneratorType&, uniform_real_distribution< > >
            generateSemiMajorAxis( randomNumberGenerator, semiMajorAxisDistribution );

    // Define normal random number distributions for test particle components of eccentricity
    // vector (h_e = e*cos( AoP ), k_e = e*sin( AoP ) ).
    normal_distribution< > distributionOfXComponentOfEccentricityVector(
                eccentricityMean * std::cos( eccentricityAngle ),
                convertFullWidthHalfMaximumToStandardDeviation( eccentricityFWHM ) );

    normal_distribution< > distributionOfYComponentOfEccentricityVector(
                eccentricityMean * std::sin( eccentricityAngle ),
                convertFullWidthHalfMaximumToStandardDeviation( eccentricityFWHM ) );

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
                inclinationMean * std::cos( inclinationAngle ),
                convertFullWidthHalfMaximumToStandardDeviation( inclinationFWHM ) );

    normal_distribution< > distributionOfYComponentOfInclinationVector(
                inclinationMean * std::sin( inclinationAngle ),
                convertFullWidthHalfMaximumToStandardDeviation( inclinationFWHM ) );

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
            << "\"case\" INTEGER PRIMARY KEY NOT NULL,"
            << "\"randomWalkSimulationDuration\" REAL NOT NULL,"
            << "\"synodicPeriodLimit\" REAL NOT NULL,"
            << "\"outputInterval\" REAL NOT NULL,"
            << "\"startUpIntegrationDuration\" REAL NOT NULL,"
            << "\"conjunctionEventDetectionDistance\" REAL NOT NULL,"
            << "\"oppositionEventDetectionDistance\" REAL NOT NULL,"
            << "\"centralBodyGravitationalParameter\" REAL NOT NULL,"
            << "\"centralBodyJ2GravityCoefficient\" REAL NOT NULL,"
            << "\"centralBodyEquatorialRadius\" REAL NOT NULL,"
            << "\"limitSemiMajorAxisDistribution\" REAL NOT NULL,"
            << "\"meanEccentricityDistribution\" REAL NOT NULL,"
            << "\"fullWidthHalfMaxmimumEccentricityDistribution\" REAL NOT NULL,"
            << "\"meanInclinationDistribution\" REAL NOT NULL,"
            << "\"fullWidthHalfMaxmimumInclinationDistribution\" REAL NOT NULL,"
            << "\"perturbedBodyGravitationalParameter\" REAL NOT NULL,"
            << "\"perturbedBodySemiMajorAxisAtT0\" REAL NOT NULL,"
            << "\"perturbedBodyEccentricityAtT0\" REAL NOT NULL,"
            << "\"perturbedBodyInclinationAtT0\" REAL NOT NULL,"
            << "\"perturbedBodyArgumentOfPeriapsisAtT0\" REAL NOT NULL,"
            << "\"perturbedBodyLongitudeOfAscendingNodeAtT0\" REAL NOT NULL,"
            << "\"perturbedBodyTrueAnomalyAtT0\" REAL NOT NULL,"
            << "\"numericalIntegratorType\" TEXT NOT NULL,"
            << "\"relativeTolerance\" REAL NOT NULL,"
            << "\"absoluteTolerance\" REAL NOT NULL,"
            << "\"initialStepSize\" REAL NOT NULL);";

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

    // Check how many rows are present in table.
    std::ostringstream testParticleCaseRowCount;
    testParticleCaseRowCount << "SELECT COUNT( * ) FROM " << testParticleCaseTableName;
    int caseTableRows = database.execAndGet( testParticleCaseRowCount.str( ).c_str( ) );

    if ( caseTableRows > 1 )
    {
        std::ostringstream caseTableRowsError;
        caseTableRowsError << "Error: Table '" << testParticleCaseTableName 
                           << "'' contains " << caseTableRows << " rows!";
        throw std::runtime_error( caseTableRowsError.str( ).c_str( ) );
    }

    else if ( caseTableRows == 1 )
    {
        std::cout << "Table '" << testParticleCaseTableName << "' contains 1 row of data ... "
                  << "skipping populating table ... " << std::endl;;
    }

    // Write test particle case data to table.
    else if ( caseTableRows == 0 )
    {
        std::cout << "No data present in table '" << testParticleCaseTableName << "' ... " 
                  << std::endl;
        std::cout << "Populating table ... " << std::endl;

        // Create stringstream with test particle case data insert command.
        std::stringstream testParticleCaseDataInsert;
        testParticleCaseDataInsert
            << "INSERT INTO " << testParticleCaseTableName << " VALUES ("
            << caseNumber << ","
            << randomWalkDuration << ","
            << synodicPeriodLimit << ","
            << outputInterval << ","
            << startUpIntegrationDuration << ","
            << conjunctionEventDetectionDistance << ","
            << oppositionEventDetectionDistance << ","
            << centralBodyGravitationalParameter << ","
            << centralBodyJ2GravityCoefficient << ","
            << centralBodyEquatorialRadius << ","
            << semiMajorAxisLimit << ","
            << eccentricityMean << ","
            << eccentricityFWHM << ","
            << inclinationMean << ","
            << inclinationFWHM << ","
            << perturbedBodyGravitationalParameter << ","
            << perturbedBodyKeplerianElementsAtT0( semiMajorAxisIndex ) << ","
            << perturbedBodyKeplerianElementsAtT0( eccentricityIndex ) << ","
            << perturbedBodyKeplerianElementsAtT0( inclinationIndex ) << ","
            << perturbedBodyKeplerianElementsAtT0( argumentOfPeriapsisIndex ) << ","
            << perturbedBodyKeplerianElementsAtT0( longitudeOfAscendingNodeIndex ) << ","
            << perturbedBodyKeplerianElementsAtT0( trueAnomalyIndex ) << ","
            << "\"" << numericalIntegratorType << "\", "
            << integratorRelativeTolerance << ", "
            << integratorAbsoluteTolerance << ", "
            << initialStepSize << ");";

        // Insert test particle case data.
        database.exec( testParticleCaseDataInsert.str( ).c_str( ) );

        // Check that there is only one row present in the table.
        caseTableRows = database.execAndGet( testParticleCaseRowCount.str( ).c_str( ) );
        if ( caseTableRows == 1 )
        {
            std::cout << "Table '" << testParticleCaseTableName << "' populated successfully!" 
                      << std::endl; 
        }

        else
        {
            std::ostringstream caseTableRowsError;
            caseTableRowsError << "Error: Table '" << testParticleCaseTableName << "' contains" 
                               << caseTableRows << "rows!";
            std::runtime_error( caseTableRowsError.str( ).c_str( ) );
        }
    }

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
            << "\"testParticleSimulation\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"completed\" INTEGER NOT NULL,"
            << "\"semiMajorAxis\" REAL NOT NULL,"
            << "\"eccentricity\" REAL NOT NULL,"
            << "\"inclination\" REAL NOT NULL,"
            << "\"argumentOfPeriapsis\" REAL NOT NULL,"
            << "\"longitudeOfAscendingNode\" REAL NOT NULL,"
            << "\"trueAnomaly\" REAL NOT NULL,"
            << "\"perturbedBodyEnergyError\" REAL NULL,"
            << "\"perturbedBodyAngularMomentumError\" REAL NULL);";

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
    testParticleInputRowCount << "SELECT COUNT( * ) FROM " << testParticleInputTableName;
    int inputTableRows = database.execAndGet( testParticleInputRowCount.str( ).c_str( ) );

    if ( inputTableRows > 0 )
    {
        std::cout << "Table '" << testParticleInputTableName << "' contains "
                  << inputTableRows << " rows ... " << std::endl;
    }

    // Populate table.
    std::cout << "Populating table with input data for " << numberOfSimulations 
              << " new simulations ... " << std::endl;

    // Set up database transaction.
    SQLite::Transaction testParticleInputTableTransaction( database );

    // Set up test particle input table insert statement.
    std::ostringstream testParticleInputTableInsert;
    testParticleInputTableInsert << "INSERT INTO " << testParticleInputTableName
                                 << " VALUES (NULL, 0, :semiMajorAxis, :eccentricity, "
                                 << ":inclination, :argumentOfPeriapsis, "
                                 << ":longitudeOfAscendingNode, :trueAnomaly, NULL, NULL);";

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
                    perturbedBodyKeplerianElementsAtT0( semiMajorAxisIndex ),
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
        while ( synodicPeriodOfTestParticle > synodicPeriodLimit );

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
            << "\"key\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"testParticleSimulation\" INTEGER NOT NULL,"
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
            << "\"monteCarloRun\" INTEGER PRIMARY KEY NOT NULL,"
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
            << "\"key\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"monteCarloRun\" INTEGER NOT NULL,"
            << "\"testParticleSimulation\" INTEGER NOT NULL,"
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
            << "\"key\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"monteCarloRun\" INTEGER NOT NULL,"
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
