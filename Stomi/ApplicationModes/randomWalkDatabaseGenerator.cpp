/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <boost/filesystem.hpp> 
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <SQLiteCpp/SQLiteCpp.h> 

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>

#include <Tudat/InputOutput/dictionaryTools.h>

#include <Assist/Astrodynamics/astrodynamicsBasics.h>
#include <Assist/Astrodynamics/hillSphere.h>
#include <Assist/Astrodynamics/unitConversions.h>

#include "Stomi/ApplicationModes/randomWalkDatabaseGenerator.h"
#include "Stomi/Database/databaseReadFunctions.h" 
#include "Stomi/Database/testParticleCase.h" 
#include "Stomi/InputOutput/dictionaries.h"

namespace stomi
{
namespace application_modes
{

//! Execute random walk database generator.
void executeRandomWalkDatabaseGenerator( 
    const std::string databasePath, 
    const tudat::input_output::parsed_data_vector_utilities::ParsedDataVectorPtr parsedData )
{
    ///////////////////////////////////////////////////////////////////////////

    // Declare using-statements.

    using std::advance;
    using std::cout;
    using std::endl;
    using std::generate;
    using std::ostringstream;
    using std::numeric_limits;
    using std::runtime_error;
    using std::setprecision;
    using std::string;
    using std::vector;

    using namespace boost::filesystem;
    using namespace boost::random;

    using namespace SQLite;

    using namespace assist::astrodynamics;
    // using namespace assist::input_output;

    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::input_output::dictionary;

    using namespace stomi::database;
    using namespace stomi::input_output;    

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Extract input parameters.

    // Get dictionary.
    const DictionaryPointer dictionary = getRandomWalkDatabaseGeneratorDictionary( );

    // Print database path to console.
    cout << "Database                                                  " 
         << databasePath << endl;

    // Extract required parameters.
    const string randomWalkRunName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), 
                findEntry( dictionary, "RANDOMWALKRUN" ) );
    cout << "Run                                                       " 
         << randomWalkRunName << endl;

    const string testParticleCaseName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), 
                findEntry( dictionary, "TESTPARTICLECASE" ) );
    cout << "Test particle case                                        " 
         << testParticleCaseName << endl;

    const unsigned int monteCarloPopulation = extractParameterValue< unsigned int >(
                parsedData->begin( ), parsedData->end( ), 
                findEntry( dictionary, "MONTECARLOPOPULATION" ) );
    cout << "Monte Carlo population                                    "
         << monteCarloPopulation << endl;    

    const double perturberRingNumberDensity = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBERRINGNUMBERDENSITY" ) );    
    cout << "Perturber ring number density                             " 
         << perturberRingNumberDensity << " perturbers per R_Hill" << endl;

    const double perturberRingMass = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBERRINGMASS" ) );    
    cout << "Perturber ring mass                                       "
         << perturberRingMass << " M_PerturbedBody" << endl;      

    const double observationPeriod = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "OBSERVATIONPERIOD" ),
                TUDAT_NAN, &convertJulianYearsToSeconds );
    cout << "Observation period                                        "
         << convertSecondsToJulianYears( observationPeriod ) << " yrs" << endl; 

    const unsigned int numberOfEpochWindows = extractParameterValue< unsigned int >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "NUMBEROFEPOCHWINDOWS" ) );    
    cout << "Number of epoch windows                                   "
         << numberOfEpochWindows << endl;                                                 

    const double epochWindowSize = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "EPOCHWINDOWSIZE" ),
                TUDAT_NAN, &convertJulianDaysToSeconds  );
    cout << "Epoch window size                                         "
         << convertSecondsToJulianDays( epochWindowSize ) << " days" << endl;                    

    // Extract optional parameters (parameters that take on default values if they are not  
    // specified in the configuration file).

    const string randomWalkRunTable = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKRUNTABLENAME" ),
                "random_walk_run" );
    cout << "Random walk run table                                     " 
         << randomWalkRunTable << endl;

    const string randomWalkInputTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKINPUTTABLENAME" ),
                "random_walk_input" );
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

    const string testParticleInputTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLEINPUTTABLENAME" ), "test_particle_input" );
    cout << "Test particle input table                                 "
         << testParticleInputTableName << endl;                

    // Check that all required parameters have been set.
    checkRequiredParameters( dictionary );

    ///////////////////////////////////////////////////////////////////////////  

    ///////////////////////////////////////////////////////////////////////////

    cout << endl;
    cout << "****************************************************************************" << endl;
    cout << "Database operations" << endl;
    cout << "****************************************************************************" << endl;
    cout << endl;

    // Open (and if necessary create) database.

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

    // Fetch test particle case and input data from database.

    // Generate output message.
    cout << "Fetching test particle case data from database ..." << endl;   

    // Retrieve and store test particle case data.
    const TestParticleCasePointer testParticleCaseData = getTestParticleCase( 
        databasePath, testParticleCaseName, testParticleCaseTableName ); 

    // Generate output message to indicate that test particle case data was fetched successfully.
    cout << "Test particle case data fetched successfully from database!" << endl;   

    // Generate output message.
    cout << "Fetching test particle input data from database ..." << endl;   

    // Get entire test particle input table from database. Only test particle simulations that 
    // are complete are fetched for the given test particle case ID.
    const TestParticleInputTable testParticleInputTable 
        = getCompleteTestParticleInputTable( 
            databasePath, testParticleCaseData->testParticleCaseId, 
            testParticleInputTableName, true );    

    // Generate output message to indicate that test particle input table was fetched successfully.
    cout << "Test particle input data (" << testParticleInputTable.size( )
         << " rows) fetched successfully from database!" << endl; 

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Compute derived parameters.

    // Compute mass of perturbed body [kg].
    const double perturbedBodyMass = computeMassOfSphere(
                testParticleCaseData->perturbedBodyRadius, 
                testParticleCaseData->perturbedBodyBulkDensity );

    // Compute perturbed body's gravitational parameter [m^3 s^-2].
    const double perturbedBodyGravitationalParameter
            = computeGravitationalParameter( perturbedBodyMass );

    // Compute perturber population using density and semi-major axis distribution limits.
    // Note, in the floor() function, adding 0.5 is a workaround for the fact that there is no
    // round() function in C++03 (it is available in C++11).
    ConvertHillRadiiToMeters convertHillRadiiToMeters( 
        testParticleCaseData->centralBodyGravitationalParameter, 
        perturbedBodyGravitationalParameter, 
        testParticleCaseData->perturbedBodyStateInKeplerianElementsAtT0( semiMajorAxisIndex ) );
    const double perturberRingNumberDensityInMeters = perturberRingNumberDensity / convertHillRadiiToMeters( 1.0 );
    const unsigned int perturberPopulation = std::floor( 
        2.0 * testParticleCaseData->semiMajorAxisDistributionLimit 
        * perturberRingNumberDensityInMeters + 0.5 );

    // Check if desired perturber population is greater than number of completed simulations 
    // fetched, from the database input table, throw a run-time error.
    if ( perturberPopulation > testParticleInputTable.size( ) )
    {
        throw runtime_error(
            "ERROR: Perturber population > Number of completed test particle simulations" );
    }

    ///////////////////////////////////////////////////////////////////////////  

    ///////////////////////////////////////////////////////////////////////////

    // Define random number generators.

    // Define a random number generator and initialize it with a reproducible seed (current cpu
    // time).
    GlobalRandomNumberGeneratorType randomNumberGenerator = getGlobalRandomNumberGenerator( );

    // Define an uniform random number distribution to randomly select the start epoch of the 
    // observation period during random walk simulation.
    uniform_real_distribution< > observationPeriodStartDistribution(
                epochWindowSize / 2.0, 
                testParticleCaseData->randomWalkSimulationPeriod 
                - observationPeriod - epochWindowSize / 2.0 );

    // Define variate generator for the first epoch in the observation period.
    variate_generator< GlobalRandomNumberGeneratorType&, uniform_real_distribution< > > 
    generateObservationPeriodStartEpoch( 
        randomNumberGenerator, observationPeriodStartDistribution );

    // Define an uniform random number distribution for test particle simulation ID indices in 
    // input table.
    uniform_int_distribution< > testParticleSimulationIdIndexDistribution( 
        0, testParticleInputTable.size( ) - 1 );

    // Define variate generator for test particle simulation ID indices using the random number
    // generator and uniform distribution of test particle simulation ID indices in input table.
    variate_generator< GlobalRandomNumberGeneratorType&, uniform_int_distribution< > > 
    generateTestParticleSimulationId( randomNumberGenerator, 
        testParticleSimulationIdIndexDistribution );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up random walk run table.
    if ( !database.tableExists( randomWalkRunTable.c_str( ) ) )
    {
        cout << "Table '" << randomWalkRunTable << "' does not exist ..." << endl;
        cout << "Creating table ... " << endl;

        // Create table.
        ostringstream randomWalkRunTableCreate;
        randomWalkRunTableCreate
            << "CREATE TABLE IF NOT EXISTS " << randomWalkRunTable << " ("
            << "\"randomWalkRunId\" INTEGER PRIMARY KEY NOT NULL,"
            << "\"randomWalkRunName\" TEXT NOT NULL,"
            << "\"testParticleCaseId\" INTEGER NOT NULL,"
            << "\"perturberRingNumberDensity\" INTEGER NOT NULL,"
            << "\"perturberRingMass\" INTEGER NOT NULL,"            
            << "\"observationPeriod\" REAL NOT NULL,"
            << "\"numberOfEpochWindows\" INTEGER NOT NULL,"
            << "\"epochWindowSize\" REAL NOT NULL);";

        // Execute command to create table.
        database.exec( randomWalkRunTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( randomWalkRunTable.c_str( ) ) )
        {
            cout << "Table '" << randomWalkRunTable << "' successfully created!" << endl; 
        }

        else
        {
            ostringstream tableCreateError;
            tableCreateError << "ERROR: Creating table '" << randomWalkRunTable 
                             << "'' failed!";
            throw runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        cout << "Table '" << randomWalkRunTable 
             << "' already exists ... skipping creating table ..." << endl;
    }

    // Check data present in table.

    // Check if run is already present in table.
    ostringstream randomWalkRunCheck;
    randomWalkRunCheck << "SELECT COUNT(*) FROM " << randomWalkRunTable 
                        << " WHERE \"randomWalkRunName\" = \"" << randomWalkRunName << "\"";
    int numberOfRunRows = database.execAndGet( randomWalkRunCheck.str( ).c_str( ) );

    if ( numberOfRunRows > 1 )
    {
        ostringstream numberOfRunRowsError;
        numberOfRunRowsError << "ERROR: Table '" << randomWalkRunTable << "' contains " 
                              << numberOfRunRows << " rows for random walk run '" 
                              << randomWalkRunName << "'!";
        throw runtime_error( numberOfRunRowsError.str( ).c_str( ) );
    }

    else if ( numberOfRunRows == 1 )
    {
        cout << "Table '" << randomWalkRunTable << "' contains 1 row of data for run '" 
             << randomWalkRunName << "' ... " << "skipping populating table ... " << endl;
    }

    // Write random walk run data to table.
    else if ( numberOfRunRows == 0 )
    {
        cout << "No data present in table '" << randomWalkRunTable << "' for random walk run '" 
             << randomWalkRunName << "' ... " << endl;
        cout << "Populating table ... " << endl;

        // Create stringstream with random walk run data insert command.
        // For floating-point values, ensure the data is written to the stream at full precision.
        ostringstream randomWalkRunDataInsert;
        randomWalkRunDataInsert
            << "INSERT INTO " << randomWalkRunTable << " VALUES ("
            << "NULL,"
            << "\"" << randomWalkRunName << "\","
            << testParticleCaseData->testParticleCaseId << ",";
        randomWalkRunDataInsert 
            << setprecision( numeric_limits< double >::digits10 )
            << perturberRingNumberDensity << ","
            << perturberRingMass << ","
            << observationPeriod << ",";
        randomWalkRunDataInsert
            << numberOfEpochWindows << ",";
        randomWalkRunDataInsert
            << setprecision( numeric_limits< double >::digits10 )
            << epochWindowSize
            << ");";

        // Insert random walk run data.
        database.exec( randomWalkRunDataInsert.str( ).c_str( ) );

        // Check that there is only one row present in the table.
        numberOfRunRows = database.execAndGet( randomWalkRunCheck.str( ).c_str( ) );
        if ( numberOfRunRows == 1 )
        {
            cout << "Table '" << randomWalkRunTable << "' populated successfully!" << endl; 
        }

        else
        {
            ostringstream numberOfRunRowsError;
            numberOfRunRowsError << "ERROR: Table '" << randomWalkRunTable << "' contains "
                                 << numberOfRunRows << " rows for random walk run '" 
                                 << randomWalkRunName << "'!";
            throw runtime_error( numberOfRunRowsError.str( ).c_str( ) );
        }
    }    

    // Retrieve and output run id.
    ostringstream randomWalkRunIdQuery;          
    randomWalkRunIdQuery << "SELECT \"randomWalkRunId\" FROM " << randomWalkRunTable 
                         << " WHERE \"randomWalkRunName\" = \"" << randomWalkRunName << "\"";
    const int randomWalkRunId = database.execAndGet( randomWalkRunIdQuery.str( ).c_str( ) );
    cout << "Run ID is " << randomWalkRunId << " for random walk run '" 
         << randomWalkRunName << "'" << endl;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up random walk perturber table.
    if ( !database.tableExists( randomWalkPerturberTableName.c_str( ) ) )
    {
        cout << "Table '" << randomWalkPerturberTableName << "' does not exist ..." << endl;
        cout << "Creating table ... " << endl;

        // Create table.
        ostringstream randomWalkPerturberTableCreate;
        randomWalkPerturberTableCreate
            << "CREATE TABLE IF NOT EXISTS " << randomWalkPerturberTableName << " ("
            << "\"randomWalkPerturberId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"randomWalkSimulationId\" INTEGER NOT NULL,"
            << "\"testParticleSimulationId\" INTEGER NOT NULL);";

        // Execute command to create table.
        database.exec( randomWalkPerturberTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( randomWalkPerturberTableName.c_str( ) ) )
        {
            cout << "Table '" << randomWalkPerturberTableName << "' successfully created!" 
                 << endl; 
        }

        else
        {
            ostringstream tableCreateError;
            tableCreateError << "ERROR: Creating table '" << randomWalkPerturberTableName 
                             << "'' failed!";
            throw runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        cout << "Table '" << randomWalkPerturberTableName 
             << "' already exists ... skipping creating table ..." << endl;
    }

    ///////////////////////////////////////////////////////////////////////////    

    ///////////////////////////////////////////////////////////////////////////

    // Set up random walk input table.
    if ( !database.tableExists( randomWalkInputTableName.c_str( ) ) )
    {
        cout << "Table '" << randomWalkInputTableName << "' does not exist ..." << endl;
        cout << "Creating table ... " << endl;

        // Create table.
        ostringstream randomWalkInputTableCreate;
        randomWalkInputTableCreate
            << "CREATE TABLE IF NOT EXISTS " << randomWalkInputTableName << " ("
            << "\"randomWalkSimulationId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"randomWalkRunId\" INTEGER NOT NULL,"
            << "\"completed\" INTEGER NOT NULL,"
            << "\"observationPeriodStartEpoch\" REAL NOT NULL);";

        // Execute command to create table.
        database.exec( randomWalkInputTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( randomWalkInputTableName.c_str( ) ) )
        {
            cout << "Table '" << randomWalkInputTableName << "' successfully created!" << endl; 
        }

        else
        {
            ostringstream tableCreateError;
            tableCreateError << "ERROR: Creating table '" << randomWalkInputTableName 
                             << "'' failed!";
            throw runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        cout << "Table '" << randomWalkInputTableName 
             << "' already exists ... skipping creating table ..." << endl;
    }

    // Check data present in table.

    // Check how many rows are present in table.
    ostringstream randomWalkInputRowCount;
    randomWalkInputRowCount << "SELECT COUNT(*) FROM " << randomWalkInputTableName 
                            << " WHERE \"randomWalkRunId\" == " << randomWalkRunId;
    const int inputTableRows = database.execAndGet( randomWalkInputRowCount.str( ).c_str( ) );

    if ( inputTableRows > 0 )
    {
        cout << "Table '" << randomWalkInputTableName << "' contains " << inputTableRows 
             << " rows for random walk run '" << randomWalkRunName 
             << "' ... skipping populating table ..." << endl;
    }

    else
    {
        // Populate table.
        // Table containing list of random walk perturbers is populated simultaneously.
        cout << "Populating input table with " << monteCarloPopulation 
             << " new random walk simulations ... " << endl;
        cout << "Populating perturber table with " << monteCarloPopulation 
             << " new data for random walk simulations ... " << endl;         

        // Set up database transaction.
        Transaction randomWalkInputTableTransaction( database );

        // Set up random walk input table insert statement.
        ostringstream randomWalkInputTableInsert;
        randomWalkInputTableInsert << "INSERT INTO " << randomWalkInputTableName
                                   << " VALUES (NULL, " << randomWalkRunId 
                                   << ", 0, :observationPeriodStartEpoch);";

        // Compile a SQL query.
        Statement randomWalkInputTableInsertQuery( 
            database, randomWalkInputTableInsert.str( ).c_str( ) );

        // Set up random walk perturber table insert statement.
        ostringstream randomWalkPerturberTableInsert;
        randomWalkPerturberTableInsert << "INSERT INTO " << randomWalkPerturberTableName
                                       << " VALUES (NULL, :randomWalkSimulationId" 
                                       << ", :testParticleSimulationId);";

        // Compile a SQL query.
        Statement randomWalkPerturberTableInsertQuery( 
            database, randomWalkPerturberTableInsert.str( ).c_str( ) );

        // Generate random walk input data and populate table.
        // Also populate perturber table for each random walk simulation.
        for ( unsigned int randomWalkSimulation = 0; 
              randomWalkSimulation < monteCarloPopulation; randomWalkSimulation++ )
        { 
            // Bind values to prepared SQLite statement.
            randomWalkInputTableInsertQuery.bind( ":observationPeriodStartEpoch", 
                                                  generateObservationPeriodStartEpoch( ) );

            // Execute insert query.
            randomWalkInputTableInsertQuery.exec( );

            // Reset query.
            randomWalkInputTableInsertQuery.reset( );

            // Get row ID for last input row inserted in table.
            const int lastInputTableRowId = database.getLastInsertRowid( );

            // Get random walk simulation ID corresponding to row ID.
            ostringstream randomWalkSimulationRowIdQuery;
            randomWalkSimulationRowIdQuery << "SELECT \"randomWalkSimulationId\" FROM " 
                                           << randomWalkInputTableName 
                                           << " WHERE rowid == " << lastInputTableRowId;
            const int randomWalkSimulationId = database.execAndGet( 
                randomWalkSimulationRowIdQuery.str( ).c_str( ) );

            // Select test particle simulation ID indices to generate list of random walk 
            // perturbers.

            // Declare vector containing test particle simulation IDs.
            vector< unsigned int > testParticleSimulationIdIndices( perturberPopulation );

            // Generate vector of randomly selected test particle simulation IDs.
            generate( testParticleSimulationIdIndices.begin( ), 
                      testParticleSimulationIdIndices.end( ),
                      generateTestParticleSimulationId );

            // Check if the test particle simulation IDs are unique, and if not, generate new 
            // random number.
            for ( unsigned int i = 0; i < testParticleSimulationIdIndices.size( ); i++ )
            {
                for ( unsigned int j = 0; j < testParticleSimulationIdIndices.size( ); j++ )
                {
                    // If inner and outer loop point to the same element, skip.
                    if ( i == j )
                    {
                        continue;
                    }

                    // Else, check if the elements are equal, and if they are generate a new 
                    // test particle simulation ID index and restart looping.
                    else if ( testParticleSimulationIdIndices.at( j ) 
                              == testParticleSimulationIdIndices.at( i ) )
                    {
                        testParticleSimulationIdIndices.at( j ) 
                            = generateTestParticleSimulationId( );
                        i = 0;
                        break;
                    }
                }
            }       

            // Loop through list of test particle simulation ID indices and write data to 
            // random walk perturber table.
            for ( unsigned int k = 0; k < testParticleSimulationIdIndices.size( ); k++ ) 
            {
                // Bind values to prepared SQLite statement.
                randomWalkPerturberTableInsertQuery.bind( 
                    ":randomWalkSimulationId", randomWalkSimulationId );
                
                TestParticleInputTable::iterator iteratorTestParticleInputTable 
                    = testParticleInputTable.begin( );
                advance( iteratorTestParticleInputTable, testParticleSimulationIdIndices.at( k ) );

                randomWalkPerturberTableInsertQuery.bind( 
                    ":testParticleSimulationId",
                     iteratorTestParticleInputTable->testParticleSimulationId );  

                // Execute insert query.
                randomWalkPerturberTableInsertQuery.exec( );

                // Reset query.
                randomWalkPerturberTableInsertQuery.reset( );            
            }
        }

        // Commit transaction.
        randomWalkInputTableTransaction.commit( );
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up random walk output table.
    if ( !database.tableExists( randomWalkOutputTableName.c_str( ) ) )
    {
        cout << "Table '" << randomWalkOutputTableName << "' does not exist ..." << endl;
        cout << "Creating table ... " << endl;

        // Create table.
        ostringstream randomWalkOutputTableCreate;
        randomWalkOutputTableCreate
            << "CREATE TABLE IF NOT EXISTS " << randomWalkOutputTableName << " ("
            << "\"randomWalkOutputId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"randomWalkSimulationId\" INTEGER NOT NULL,"
            << "\"averageLongitudeResidual\" REAL NOT NULL,"
            << "\"maximumLongitudeResidualChange\" REAL NOT NULL,"            
            << "\"averageEccentricity\" REAL NOT NULL,"            
            << "\"maximumEccentricityChange\" REAL NOT NULL,"
            << "\"averageInclination\" REAL NOT NULL,"
            << "\"maximumInclinationChange\" REAL NOT NULL);";

        // Execute command to create table.
        database.exec( randomWalkOutputTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( randomWalkOutputTableName.c_str( ) ) )
        {
            cout << "Table '" << randomWalkOutputTableName << "' successfully created!" << endl; 
        }

        else
        {
            ostringstream tableCreateError;
            tableCreateError << "ERROR: Creating table '" << randomWalkOutputTableName
                             << "'' failed!";
            throw runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        cout << "Table '" << randomWalkOutputTableName 
             << "' already exists ... skipping creating table ..." << endl;
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Close database.

    // Database will be automatically closed after the application is terminated. 
    // (when object goes out of scope its destructor will be called).
    cout << "SQLite database file '" << database.getFilename( ).c_str( ) 
         << "' closed successfully ..." << endl;

    ///////////////////////////////////////////////////////////////////////////
}

} // namespace application_modes
} // namespace stomi 
