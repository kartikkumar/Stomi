/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <SQLiteCpp/SQLiteCpp.h>

#include <Assist/Astrodynamics/astrodynamicsBasics.h>
#include <Assist/Astrodynamics/hillSphere.h> 
#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/InputOutput/basicInputOutput.h>
#include <Assist/Mathematics/statistics.h>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h>
#include <Tudat/InputOutput/dictionaryTools.h>
#include <Tudat/InputOutput/fieldType.h>
#include <Tudat/InputOutput/separatedParser.h>
#include <Tudat/InputOutput/parsedDataVectorUtilities.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include "Stomi/InputOutput/dictionaries.h"

//! Execute test particle database generator.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    ///////////////////////////////////////////////////////////////////////////

    // Declare using-statements.

    using std::cout;
    using std::endl;
    using std::ostringstream;
    using std::runtime_error;
    using std::scientific;
    using std::setprecision;
    using std::string;
    using std::stringstream;
    
    using namespace boost::filesystem;
    using namespace boost::random;

    using namespace SQLite;

    using namespace assist::astrodynamics;
    using namespace assist::input_output;
    using namespace assist::mathematics;

    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_astrodynamics::unit_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::basic_mathematics::mathematical_constants;
    using namespace tudat::input_output::dictionary;
    using namespace tudat::input_output;
    using namespace tudat::input_output::field_types::general;
    using namespace tudat::input_output::parsed_data_vector_utilities;

    using namespace stomi::input_output;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up input deck.

    // Check number of input parameters is correct (the numberOfInputs variable includes the
    // application itself, so one is subtracted from this number).
    checkNumberOfInputArguments( numberOfInputs - 1 );

    // Get input parameter dictionary.
    DictionaryPointer dictionary = getTestParticleDatabaseGeneratorDictionary( );

    // Read and filter input stream (this can't be declared const because the parser's parse
    // function is not const-correct at the moment).
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
                parsedData->begin( ), parsedData->end( ), 
                findEntry( dictionary, "DATABASEPATH" ) );
    cout << "Database                                                  " << databasePath << endl;

    const string caseName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "CASE" ) );
    cout << "Case                                                      " << caseName << endl;

    const double numberOfSimulations = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "NUMBEROFSIMULATIONS" ) );
    cout << "Number of simulations                                     "
         << numberOfSimulations << endl;

    const double randomWalkSimulationPeriod = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKSIMULATIONPERIOD" ),
                TUDAT_NAN, &convertJulianYearsToSeconds );
    cout << "Random walk simulation period                             " 
         << convertSecondsToJulianYears( randomWalkSimulationPeriod ) << " yrs" << endl;

    const double centralBodyGravitationalParameter = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYGRAVITATIONALPARAMETER" ) );
    cout << "Central body gravitational parameter                      " 
         << centralBodyGravitationalParameter << " m^3 s^-2" << endl;

    const double perturbedBodyRadius = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYRADIUS" ) );
    cout << "Perturbed body radius                                     " 
         << convertMetersToKilometers( perturbedBodyRadius ) << " km" << endl;

    const double perturbedBodyBulkDensity = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYBULKDENSITY" ) );
    cout << "Perturbed body bulk density                               " 
         << perturbedBodyBulkDensity << " kg m^-3" << endl;

    // Compute mass of perturbed body [kg].
    const double perturbedBodyMass = computeMassOfSphere(
                perturbedBodyRadius, perturbedBodyBulkDensity );
    cout << "Perturbed body mass                                       " 
         << perturbedBodyMass << " kg" << endl;    

    // Compute perturbed body's gravitational parameter [m^3 s^-2].
    const double perturbedBodyGravitationalParameter
            = computeGravitationalParameter( perturbedBodyMass );  
    cout << "Perturbed body gravitational parameter                    " 
         << perturbedBodyGravitationalParameter << " m^3 s^-2" << endl;    

    Vector6d perturbedBodyStateInKeplerianElementsAtT0( 6 );

    perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYSEMIMAJORAXISATT0" ) );
    cout << "Perturbed body semi-major axis at TO                      "
         << perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) << " m" << endl;

    perturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYECCENTRICITYATT0" ) );
    cout << "Perturbed body eccentricity at TO                         "
         << perturbedBodyStateInKeplerianElementsAtT0( eccentricityIndex ) << endl;

    perturbedBodyStateInKeplerianElementsAtT0( inclinationIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYINCLINATIONATT0" ), TUDAT_NAN,
                &convertDegreesToRadians< double > );
    cout << "Perturbed body inclination at TO                          "
         << convertRadiansToDegrees( 
                perturbedBodyStateInKeplerianElementsAtT0( inclinationIndex ) ) << " deg" << endl;

    perturbedBodyStateInKeplerianElementsAtT0( argumentOfPeriapsisIndex )
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYARGUMENTOFPERIAPSISATT0" ), TUDAT_NAN,
                 &convertDegreesToRadians< double > );
    cout << "Perturbed body argument of periapsis at TO                "
         << convertRadiansToDegrees( 
                perturbedBodyStateInKeplerianElementsAtT0( argumentOfPeriapsisIndex ) ) 
         << " deg" << endl;

    perturbedBodyStateInKeplerianElementsAtT0( longitudeOfAscendingNodeIndex )
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYLONGITUDEOFASCENDINGNODEATT0" ), TUDAT_NAN,
                &convertDegreesToRadians< double > );
    cout << "Perturbed body longitude of ascending node at TO          "
         << convertRadiansToDegrees( 
                perturbedBodyStateInKeplerianElementsAtT0( longitudeOfAscendingNodeIndex ) ) 
         << " deg" << endl;

    perturbedBodyStateInKeplerianElementsAtT0( trueAnomalyIndex ) 
            = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBEDBODYTRUEANOMALYATT0" ), TUDAT_NAN,
                &convertDegreesToRadians< double > );
    cout << "Perturbed body true anomaly at TO                         "
         << convertRadiansToDegrees( 
                perturbedBodyStateInKeplerianElementsAtT0( trueAnomalyIndex ) ) << " deg" << endl;

    const double semiMajorAxisDistributionLimit = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "SEMIMAJORAXISDISTRIBUTIONLIMIT" ), TUDAT_NAN,
                ConvertHillRadiiToMeters( 
                    centralBodyGravitationalParameter, 
                    perturbedBodyGravitationalParameter,
                    perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) ) );
    cout << "Semi-major axis limit                                     " 
         << semiMajorAxisDistributionLimit << " m" << endl;

    // Extract optional parameters (parameters that take on default values if they are not  
    // specified in the input file).

    const double synodicPeriodMaximum = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "SYNODICPERIODMAXIMUM" ),
                randomWalkSimulationPeriod, &convertJulianYearsToSeconds );
    cout << "Synodic period maximum                                    " 
         << convertSecondsToJulianYears( synodicPeriodMaximum ) << " yrs" << endl;

    const double startUpIntegrationPeriod = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "STARTUPINTEGRATIONPERIOD" ), 0.0,
                &convertJulianYearsToSeconds );
    cout << "Start-up integration period                               " 
         << convertSecondsToJulianYears( startUpIntegrationPeriod ) << " yrs" << endl;

    const double centralBodyJ2GravityCoefficient = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYJ2GRAVITYCOEFFICIENT" ), 0.0 );
    cout << "Central body J2 gravity coefficient                       "
         << centralBodyJ2GravityCoefficient << endl;

    const double centralBodyEquatorialRadius = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CENTRALBODYEQUATORIALRADIUS" ), 0.0 );
    cout << "Central body equatorial radius                            "
         << convertMetersToKilometers( centralBodyEquatorialRadius ) << " km" << endl;

    const double conjunctionEventDetectionDistance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "CONJUNCTIONEVENTDETECTIONDISTANCE" ), 
                0.5 * perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) );
    cout << "Conjunction event detection distance                      " 
         << convertMetersToKilometers( conjunctionEventDetectionDistance ) << " km" << endl;

    const double oppositionEventDetectionDistance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "OPPOSITIONEVENTDETECTIONDISTANCE" ),
                1.5 * perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) );
    cout << "Opposition event detection distance                       " 
         << convertMetersToKilometers( oppositionEventDetectionDistance ) << " km" << endl;

    const double eccentricityDistributionMean = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "ECCENTRICITYDISTRIBUTIONMEAN" ), 0.0 );
    cout << "Eccentricity distribution mean                            " 
         << eccentricityDistributionMean << endl;

    const double eccentricityDistributionFullWidthHalfMaxmimum = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "ECCENTRICITYDISTRIBUTIONFULLWIDTHHALFMAXIMUM" ), 0.0 );
    cout << "Eccentricity distribution Full-Width Half-Maximum         " 
         << eccentricityDistributionFullWidthHalfMaxmimum << endl;

    const double inclinationDistributionMean = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INCLINATIONDISTRIBUTIONMEAN" ), 0.0,
                &convertDegreesToRadians< double > );
    cout << "Inclination distribution mean                             "
         << convertRadiansToDegrees( inclinationDistributionMean ) << " deg" << endl;

    const double inclinationDistributionFullWidthHalfMaxmimum = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INCLINATIONDISTRIBUTIONFULLWIDTHHALFMAXIMUM" ), 0.0,
                &convertDegreesToRadians< double > );
    cout << "Inclination distribution Full-Width Half-Maximum          " 
         << convertRadiansToDegrees( inclinationDistributionFullWidthHalfMaxmimum ) 
         << " deg" << endl;

    const string numericalIntegratorType = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "NUMERICALINTEGRATORTYPE" ), "DOPRI853" );
    cout << "Numerical integrator type                                 "
         << numericalIntegratorType << endl;

    const double numericalIntegratorInitialStepSize = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "INITIALSTEPSIZE" ), 60.0 );
    cout << "Initial step size                                         " 
         << numericalIntegratorInitialStepSize << " s" << endl;

    const double numericalIntegratorRelativeTolerance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RUNGEKUTTARELATIVEERRORTOLERANCE" ), 1.0e-12 );
    cout << "Numerical integrator relative tolerance                   " 
         << numericalIntegratorRelativeTolerance << endl;

    // Compute absolute tolerance with margin and round off.
    const double absoluteTolerance 
        = perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) * 1.0e-15;
    stringstream absoluteToleranceRoundedBuffer;
    absoluteToleranceRoundedBuffer << setprecision( 0 ) << scientific << absoluteTolerance;
    double absoluteToleranceRounded = TUDAT_NAN;
    absoluteToleranceRoundedBuffer >> absoluteToleranceRounded;

    const double numericalIntegratorAbsoluteTolerance = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RUNGEKUTTAABSOLUTEERRORTOLERANCE" ), 
                absoluteToleranceRounded );
    cout << "Numerical integrator absolute tolerance                   " 
         << numericalIntegratorAbsoluteTolerance << endl;

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

    // Define normal random number distribution for test particle components of eccentricity
    // vector (h_e = e*cos( AoP ), k_e = e*sin( AoP ) ).
    normal_distribution< > distributionOfEccentricityVectorComponent(
                eccentricityDistributionMean * std::sqrt( 2 ) * 0.5,
                convertFullWidthHalfMaximumToStandardDeviation( 
                    eccentricityDistributionFullWidthHalfMaxmimum ) );

    // Define variate generator for h_e- and k_e-values using the random number generator
    // and normal distribution of h_e- and k_e-values.
    variate_generator< GlobalRandomNumberGeneratorType&, normal_distribution< > >
            generateEccentricityVectorComponent(
                randomNumberGenerator, distributionOfEccentricityVectorComponent );

    // Define normal random number distribution for test particle components of inclination
    // vector (h_i = i*cos( RAAN ), k_i = i*sin( RAAN ) ).
    normal_distribution< > distributionOfInclinationVectorComponent(
                inclinationDistributionMean * std::sqrt( 2 ) * 0.5,
                convertFullWidthHalfMaximumToStandardDeviation( 
                    inclinationDistributionFullWidthHalfMaxmimum ) );

    // Define variate generator for h_i- and k_i-values using the random number generator
    // and normal distribution of h_i- and k_i-values.
    variate_generator< GlobalRandomNumberGeneratorType&, normal_distribution< > >
            generateInclinationVectorComponent(
                randomNumberGenerator, distributionOfInclinationVectorComponent );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Open (and if necessary create) database.

    cout << endl;
    cout << "****************************************************************************" << endl;
    cout << "Database operations" << endl;
    cout << "****************************************************************************" << endl;
    cout << endl;

    // Check if database file already exists.
    bool isDatabaseCreated = false;
    if ( exists( databasePath ) )
    {
        cout << "Opening existing database at " << databasePath << " ..." << endl;
    }

    else
    {
        isDatabaseCreated = true;
        cout << "WARNING: Database does not exist!" << endl;
        cout << "Creating database at " << databasePath << " ..." << endl;
    }

    // Create/open database.          
    Database database( databasePath.c_str( ), SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE );

    cout << "SQLite database file '" << database.getFilename( ).c_str( );
    if ( isDatabaseCreated ) { cout << "' created &"; }
    cout << " opened successfully ..." << endl;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up test particle case table.

    // Create table if it doesn't exist.
    if ( !database.tableExists( testParticleCaseTableName.c_str( ) ) )
    {
        cout << "Table '" << testParticleCaseTableName << "' does not exist ..." << endl;
        cout << "Creating table ... " << endl;

        ostringstream testParticleCaseTableCreate;
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
            << "\"eccentricityDistributionFullWidthHalfMaxmimum\" REAL NOT NULL,"
            << "\"inclinationDistributionMean\" REAL NOT NULL,"
            << "\"inclinationDistributionFullWidthHalfMaxmimum\" REAL NOT NULL,"
            << "\"numericalIntegratorType\" TEXT NOT NULL,"
            << "\"numericalIntegratorInitialStepSize\" REAL NOT NULL,"
            << "\"numericalIntegratorRelativeTolerance\" REAL NOT NULL,"
            << "\"numericalIntegratorAbsoluteTolerance\" REAL NOT NULL);";

        // Execute command to create table.
        database.exec( testParticleCaseTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( testParticleCaseTableName.c_str( ) ) )
        {
            cout << "Table '" << testParticleCaseTableName << "' successfully created!" << endl; 
        }

        else
        {
            ostringstream tableCreateError;
            tableCreateError << "Error: Creating table '" << testParticleCaseTableName 
                             << "'' failed!";
            throw runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        cout << "Table '" << testParticleCaseTableName 
             << "' already exists ... skipping creating table ..." << endl;
    }

    // Check data present in table.

    // Check if case is already present in table.
    ostringstream testParticleCaseCheck;
    testParticleCaseCheck << "SELECT COUNT(*) FROM " << testParticleCaseTableName 
                          << " WHERE \"caseName\" = \"" << caseName << "\"";
    int numberOfCaseRows = database.execAndGet( testParticleCaseCheck.str( ).c_str( ) );

    if ( numberOfCaseRows > 1 )
    {
        ostringstream numberOfCaseRowsError;
        numberOfCaseRowsError << "Error: Table '" << testParticleCaseTableName << "' contains " 
                              << numberOfCaseRows << " rows for case '" << caseName << "'!";
        throw runtime_error( numberOfCaseRowsError.str( ).c_str( ) );
    }

    else if ( numberOfCaseRows == 1 )
    {
        cout << "Table '" << testParticleCaseTableName << "' contains 1 row of data for case '" 
             << caseName << "' ... " << "skipping populating table ... " << endl;
    }

    // Write test particle case data to table.
    else if ( numberOfCaseRows == 0 )
    {
        cout << "No data present in table '" << testParticleCaseTableName << "' for case '" 
             << caseName << "' ... " << endl;
        cout << "Populating table ... " << endl;

        // Create stringstream with test particle case data insert command.
        // For floating-point values, ensure the data is written to the stream at full precision.
        ostringstream testParticleCaseDataInsert;
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
            << eccentricityDistributionFullWidthHalfMaxmimum << ","
            << inclinationDistributionMean << ","
            << inclinationDistributionFullWidthHalfMaxmimum << ",";
        testParticleCaseDataInsert    
            << "\"" << numericalIntegratorType << "\",";
        testParticleCaseDataInsert
            << std::setprecision( std::numeric_limits< double >::digits10 )
            << numericalIntegratorInitialStepSize << ","
            << numericalIntegratorRelativeTolerance << ","
            << numericalIntegratorAbsoluteTolerance
            << ");";

        // Insert test particle case data.
        database.exec( testParticleCaseDataInsert.str( ).c_str( ) );

        // Check that there is only one row present in the table.
        numberOfCaseRows = database.execAndGet( testParticleCaseCheck.str( ).c_str( ) );
        if ( numberOfCaseRows == 1 )
        {
            cout << "Table '" << testParticleCaseTableName << "' populated successfully!" << endl; 
        }

        else
        {
            ostringstream numberOfCaseRowsError;
            numberOfCaseRowsError << "Error: Table '" << testParticleCaseTableName << "' contains "
                                  << numberOfCaseRows << " rows for case '" << caseName << "'!";
            throw runtime_error( numberOfCaseRowsError.str( ).c_str( ) );
        }
    }

    // Retrieve and output case id.
    ostringstream testParticleCaseId;          
    testParticleCaseId << "SELECT \"caseId\" FROM " << testParticleCaseTableName 
                       << " WHERE \"caseName\" = \"" << caseName << "\"";
    const int caseId = database.execAndGet( testParticleCaseId.str( ).c_str( ) );
    cout << "Case ID is " << caseId << " for case '" << caseName << "'" << endl;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up test particle input table.
    if ( !database.tableExists( testParticleInputTableName.c_str( ) ) )
    {
        cout << "Table '" << testParticleInputTableName << "' does not exist ..." << endl;
        cout << "Creating table ... " << endl;

        // Create table.
        ostringstream testParticleInputTableCreate;
        testParticleInputTableCreate
            << "CREATE TABLE IF NOT EXISTS " << testParticleInputTableName << " ("
            << "\"simulationId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"testParticleCaseId\" INTEGER NOT NULL,"
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
            cout << "Table '" << testParticleInputTableName << "' successfully created!" << endl; 
        }

        else
        {
            ostringstream tableCreateError;
            tableCreateError << "Error: Creating table '" << testParticleInputTableName
                             << "'' failed!";
            throw runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        cout << "Table '" << testParticleInputTableName 
             << "' already exists ... skipping creating table ..." << endl;
    }

    // Check data present in table.

    // Check how many rows are present in table.
    ostringstream testParticleInputRowCount;
    testParticleInputRowCount << "SELECT COUNT(*) FROM " << testParticleInputTableName 
                              << " WHERE \"caseId\" = " << caseId;
    const int inputTableRows = database.execAndGet( testParticleInputRowCount.str( ).c_str( ) );

    if ( inputTableRows > 0 )
    {
        cout << "Table '" << testParticleInputTableName << "' contains " << inputTableRows 
             << " rows for case '" << caseName << "' ... " << endl;
    }

    // Populate table.
    cout << "Populating input table with " << numberOfSimulations << " new simulations ... " 
         << endl;

    // Set up database transaction.
    Transaction testParticleInputTableTransaction( database );

    // Set up test particle input table insert statement.
    ostringstream testParticleInputTableInsert;
    testParticleInputTableInsert << "INSERT INTO " << testParticleInputTableName
                                 << " VALUES (NULL, " << caseId << ", 0, :semiMajorAxis, "
                                 << ":eccentricity, :inclination, :argumentOfPeriapsis, "
                                 << ":longitudeOfAscendingNode, :trueAnomaly);";

    // Compile a SQL query.
    Statement testParticleInputTableInsertQuery( 
        database, testParticleInputTableInsert.str( ).c_str( ) );

    // Generate test particle input data and populate table.
    for ( int simulationNumber = 0; simulationNumber < numberOfSimulations; simulationNumber++ )
    {
        // Set random eccentricity vector.
        const Eigen::Vector2d eccentricityVector( generateEccentricityVectorComponent( ),
                                                  generateEccentricityVectorComponent( ) );

        // Compute eccentricity [-].
        const double eccentricity = eccentricityVector.norm( );

        // Compute argument of periapsis [rad].
        const double argumentOfPeriapsis = std::atan2( eccentricityVector.y( ),
                                                       eccentricityVector.x( ) );

        // Set random inclination vector.
        const Eigen::Vector2d inclinationVector( generateInclinationVectorComponent( ),
                                                 generateInclinationVectorComponent( ) );

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
    testParticleInputTableTransaction.commit( );

    ///////////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////////

    // Set up test particle kick table.
    if ( !database.tableExists( testParticleKickTableName.c_str( ) ) )
    {
        cout << "Table '" << testParticleKickTableName << "' does not exist ..." << endl;
        cout << "Creating table ... " << endl;

        // Create table.
        ostringstream testParticleKickTableCreate;
        testParticleKickTableCreate
            << "CREATE TABLE IF NOT EXISTS " << testParticleKickTableName << " ("
            << "\"kickId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"testParticleSimulationId\" INTEGER NOT NULL,"
            << "\"conjunctionEpoch\" REAL NOT NULL,"
            << "\"conjunctionDistance\" REAL NOT NULL,"
            << "\"preConjunctionEpoch\" REAL NOT NULL,"
            << "\"preConjunctionDistance\" REAL NOT NULL,"            
            << "\"preConjunctionSemiMajorAxis\" REAL NOT NULL,"
            << "\"preConjunctionEccentricity\" REAL NOT NULL,"
            << "\"preConjunctionInclination\" REAL NOT NULL,"
            << "\"preConjunctionArgumentOfPeriapsis\" REAL NOT NULL,"
            << "\"preConjunctionLongitudeOfAscendingNode\" REAL NOT NULL,"
            << "\"preConjunctionTrueAnomaly\" REAL NOT NULL,"
            << "\"postConjunctionEpoch\" REAL NOT NULL,"
            << "\"postConjunctionDistance\" REAL NOT NULL,"            
            << "\"postConjunctionSemiMajorAxis\" REAL NOT NULL,"
            << "\"postConjunctionEccentricity\" REAL NOT NULL,"
            << "\"postConjunctionInclination\" REAL NOT NULL,"
            << "\"postConjunctionArgumentOfPeriapsis\" REAL NOT NULL,"
            << "\"postConjunctionLongitudeOfAscendingNode\" REAL NOT NULL,"
            << "\"postConjunctionTrueAnomaly\" REAL NOT NULL)";

        // Execute command to create table.
        database.exec( testParticleKickTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( testParticleKickTableName.c_str( ) ) )
        {
            cout << "Table '" << testParticleKickTableName << "' successfully created!" << endl; 
        }

        else
        {
            ostringstream tableCreateError;
            tableCreateError << "Error: Creating table '" << testParticleKickTableName 
                             << "'' failed!";
            throw runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        cout << "Table '" << testParticleKickTableName 
             << "' already exists ... skipping creating table ..." << endl;
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Close database.

    // Database will be automatically closed after the application is terminated. 
    // (when object goes out of scope its destructor will be called).
    cout << "SQLite database file '" << database.getFilename( ).c_str( ) 
         << "' closed successfully ..." << endl;
    cout << endl;

    ///////////////////////////////////////////////////////////////////////////

    // If program is successfully completed, return 0.
    return EXIT_SUCCESS;
}

/*    
 *    TODO:
 *      - encapsulate code in try-catch blocks to capture exceptions.
 *      - execute verification of existing case data against input parameters provided to
 *        ensure consistency of inputs and possibly warn user.
 *      - expand code to enable interactive interface for user to provide inputs and select
 *        options (capture command line user input). 
 */ 
