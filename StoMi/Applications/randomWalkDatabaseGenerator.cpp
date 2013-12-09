/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string> 
#include <vector>

#include <boost/filesystem.hpp> 
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <SQLiteC++.h> 

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>

#include <Tudat/InputOutput/dictionaryTools.h>
#include <Tudat/InputOutput/fieldType.h>
#include <Tudat/InputOutput/separatedParser.h>
#include <Tudat/InputOutput/parsedDataVectorUtilities.h>

#include <Assist/Astrodynamics/astrodynamicsbasics.h>
#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/InputOutput/basicInputOutput.h> 

#include "StochasticMigration/Astrodynamics/hillSphere.h"
#include "StochasticMigration/Database/databaseReadFunctions.h" 
#include "StochasticMigration/Database/testParticleCase.h" 
#include "StochasticMigration/InputOutput/dictionaries.h"

//! Execute random walk database generator.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{

    ///////////////////////////////////////////////////////////////////////////

    // Declare using-statements.

    using std::cout;
    using std::endl;
    using std::generate;
    using std::ostringstream;
    using std::runtime_error;
    using std::string;
    using std::vector;

    using namespace boost::filesystem;
    using namespace boost::random;

    using namespace SQLite;

    using namespace assist::astrodynamics;
    using namespace assist::input_output;

    using namespace tudat::basic_astrodynamics::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::input_output::dictionary;
    using namespace tudat::input_output;
    using namespace tudat::input_output::field_types::general;
    using namespace tudat::input_output::parsed_data_vector_utilities;

    using namespace stochastic_migration::astrodynamics;
    using namespace stochastic_migration::database;
    using namespace stochastic_migration::input_output;    

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up input deck.

    // Check number of input parameters is correct (the numberOfInputs variable includes the
    // application itself, so one is subtracted from this number).
    checkNumberOfInputArguments( numberOfInputs - 1 );

    // Get input parameter dictionary.
    DictionaryPointer dictionary = getRandomWalkDatabaseGeneratorDictionary( );

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

    // Extract required parameters.
    const string databasePath = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), 
                findEntry( dictionary, "DATABASEPATH" ) );
    cout << "Database                                                  " << databasePath << endl;

    const string caseName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "CASE" ) );
    cout << "Case                                                      " << caseName << endl;

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

    const double perturberDensity = extractParameterValue< double >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "PERTURBERDENSITY" ) );    
    cout << "Perturber density                                         "
         << perturberDensity << " perturbers per R_Hill" << endl;                         

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
    // specified in the input file).

    const string randomWalkCaseTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKCASETABLENAME" ),
                "random_walk_case" );
    cout << "Random walk case table                                    " 
         << randomWalkCaseTableName << endl;

    const string randomWalkInputTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKINPUTTABLENAME" ),
                "random_walk_input" );
    cout << "Random walk input table                                   " 
         << randomWalkInputTableName << endl;          

    const string randomWalkWindowsTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKWINDOWSTABLENAME" ),
                "random_walk_windows" );
    cout << "Random walk windows table                                 " 
         << randomWalkWindowsTableName << endl;      

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

    // Fetch case and input data from database.

    // Generate output message.
    cout << "Fetching test particle case data from database ..." << endl;   

    // Retrieve and store case data.
    const TestParticleCasePointer testParticleCaseData = getTestParticleCase( 
        databasePath, testParticleCaseName, testParticleCaseTableName ); 

    // Generate output message to indicate that case data was fetched successfully.
    cout << "Test particle case data fetched successfully from database!" << endl;   

    // Generate output message.
    cout << "Fetching test particle input data from database ..." << endl;   

    // Get entire test particle input table from database. Only test particle simulations that 
    // are complete are fetched for the given case ID.
    const TestParticleInputTable testParticleInputTable = getCompleteTestParticleInputTable(
                databasePath, testParticleCaseData->caseId, testParticleInputTableName, true );    

    // Generate output message to indicate that input table was fetched successfully.
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
    const double perturberDensityInMeters = perturberDensity / convertHillRadiiToMeters( 1.0 );
    const unsigned int perturberPopulation = std::floor( 
        2.0 * testParticleCaseData->semiMajorAxisDistributionLimit 
        * perturberDensityInMeters + 0.5 );

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

    // Define variate generator for simulation ID indices using the random number
    // generator and uniform distribution of simulation ID indices in input table.
    variate_generator< GlobalRandomNumberGeneratorType&, uniform_int_distribution< > > 
    generateTestParticleSimulationId( randomNumberGenerator, 
        testParticleSimulationIdIndexDistribution );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up random walk case table.
    if ( !database.tableExists( randomWalkCaseTableName.c_str( ) ) )
    {
        cout << "Table '" << randomWalkCaseTableName << "' does not exist ..." << endl;
        cout << "Creating table ... " << endl;

        // Create table.
        ostringstream randomWalkCaseTableCreate;
        randomWalkCaseTableCreate
            << "CREATE TABLE IF NOT EXISTS " << randomWalkCaseTableName << " ("
            << "\"caseId\" INTEGER PRIMARY KEY NOT NULL,"
            << "\"caseName\" TEXT NOT NULL,"
            << "\"testParticleCaseId\" INTEGER NOT NULL,"
            << "\"perturberDensity\" INTEGER NOT NULL,"
            << "\"perturberRingMass\" INTEGER NOT NULL,"            
            << "\"observationPeriod\" REAL NOT NULL,"
            << "\"numberOfEpochWindows\" INTEGER NOT NULL,"
            << "\"epochWindowSize\" REAL NOT NULL);";

        // Execute command to create table.
        database.exec( randomWalkCaseTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( randomWalkCaseTableName.c_str( ) ) )
        {
            cout << "Table '" << randomWalkCaseTableName << "' successfully created!" << endl; 
        }

        else
        {
            ostringstream tableCreateError;
            tableCreateError << "Error: Creating table '" << randomWalkCaseTableName 
                             << "'' failed!";
            throw runtime_error( tableCreateError.str( ).c_str( ) );
        }
    }

    else
    {
        cout << "Table '" << randomWalkCaseTableName 
             << "' already exists ... skipping creating table ..." << endl;
    }

    // Check data present in table.

    // Check if case is already present in table.
    ostringstream randomWalkCaseCheck;
    randomWalkCaseCheck << "SELECT COUNT(*) FROM " << randomWalkCaseTableName 
                        << " WHERE \"caseName\" = \"" << caseName << "\"";
    int numberOfCaseRows = database.execAndGet( randomWalkCaseCheck.str( ).c_str( ) );

    if ( numberOfCaseRows > 1 )
    {
        ostringstream numberOfCaseRowsError;
        numberOfCaseRowsError << "Error: Table '" << randomWalkCaseTableName << "' contains " 
                              << numberOfCaseRows << " rows for case '" << caseName << "'!";
        throw runtime_error( numberOfCaseRowsError.str( ).c_str( ) );
    }

    else if ( numberOfCaseRows == 1 )
    {
        cout << "Table '" << randomWalkCaseTableName << "' contains 1 row of data for case '" 
             << caseName << "' ... " << "skipping populating table ... " << endl;
    }

    // Write test particle case data to table.
    else if ( numberOfCaseRows == 0 )
    {
        cout << "No data present in table '" << randomWalkCaseTableName << "' for case '" 
             << caseName << "' ... " << endl;
        cout << "Populating table ... " << endl;

        // Create stringstream with random walk case data insert command.
        // For floating-point values, ensure the data is written to the stream at full precision.
        ostringstream randomWalkCaseDataInsert;
        randomWalkCaseDataInsert
            << "INSERT INTO " << randomWalkCaseTableName << " VALUES ("
            << "NULL,"
            << "\"" << caseName << "\","
            << testParticleCaseData->caseId << ",";
        randomWalkCaseDataInsert 
            << std::setprecision( std::numeric_limits< double >::digits10 )
            << perturberDensity << ","
            << perturberRingMass << ","
            << observationPeriod << ",";
        randomWalkCaseDataInsert
            << numberOfEpochWindows << ",";
        randomWalkCaseDataInsert
            << std::setprecision( std::numeric_limits< double >::digits10 )
            << epochWindowSize
            << ");";

        // Insert random walk case data.
        database.exec( randomWalkCaseDataInsert.str( ).c_str( ) );

        // Check that there is only one row present in the table.
        numberOfCaseRows = database.execAndGet( randomWalkCaseCheck.str( ).c_str( ) );
        if ( numberOfCaseRows == 1 )
        {
            cout << "Table '" << randomWalkCaseTableName << "' populated successfully!" << endl; 
        }

        else
        {
            ostringstream numberOfCaseRowsError;
            numberOfCaseRowsError << "Error: Table '" << randomWalkCaseTableName << "' contains "
                                  << numberOfCaseRows << " rows for case '" << caseName << "'!";
            throw runtime_error( numberOfCaseRowsError.str( ).c_str( ) );
        }
    }    

    // Retrieve and output case id.
    ostringstream randomWalkCaseId;          
    randomWalkCaseId << "SELECT \"caseId\" FROM " << randomWalkCaseTableName 
                     << " WHERE \"caseName\" = \"" << caseName << "\"";
    const int caseId = database.execAndGet( randomWalkCaseId.str( ).c_str( ) );
    cout << "Case ID is " << caseId << " for case '" << caseName << "'" << endl;

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
            << "\"perturberId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"monteCarloRunId\" INTEGER NOT NULL,"
            << "\"testParticleSimulationId\" INTEGER NOT NULL);";

        // Execute command to create table.
        database.exec( randomWalkPerturberTableCreate.str( ).c_str( ) );

        // Check that the table was created successfully.
        if ( database.tableExists( randomWalkPerturberTableName.c_str( ) ) )
        {
            cout << "Table '" << randomWalkPerturberTableName << "' successfully created!" << endl; 
        }

        else
        {
            ostringstream tableCreateError;
            tableCreateError << "Error: Creating table '" << randomWalkPerturberTableName 
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
        << "\"monteCarloRunId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
        << "\"randomWalkCaseId\" INTEGER NOT NULL,"
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
            tableCreateError << "Error: Creating table '" << randomWalkInputTableName 
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
                            << " WHERE \"randomWalkCaseId\" == " << caseId;
    const int inputTableRows = database.execAndGet( randomWalkInputRowCount.str( ).c_str( ) );

    if ( inputTableRows > 0 )
    {
        cout << "Table '" << randomWalkInputTableName << "' contains " << inputTableRows 
             << " rows for case '" << caseName << "' ... " << endl;
    }

    // Populate table.
    // Table containing list of perturbers is populated simultaneously.
    cout << "Populating input table with " << monteCarloPopulation 
         << " new Monte Carlo simulations ... " << endl;
    cout << "Populating perturber table with " << monteCarloPopulation 
         << " new data for Monte Carlo simulations ... " << endl;         

    // Set up database transaction.
    Transaction randomWalkInputTableTransaction( database );

    // Set up random walk input table insert statement.
    ostringstream randomWalkInputTableInsert;
    randomWalkInputTableInsert << "INSERT INTO " << randomWalkInputTableName
                               << " VALUES (NULL, " << caseId 
                               << ", 0, :observationPeriodStartEpoch);";

    // Compile a SQL query.
    Statement randomWalkInputTableInsertQuery( 
        database, randomWalkInputTableInsert.str( ).c_str( ) );

    // Set up random walk perturber table insert statement.
    ostringstream randomWalkPerturberTableInsert;
    randomWalkPerturberTableInsert << "INSERT INTO " << randomWalkPerturberTableName
                                   << " VALUES (NULL, :monteCarloRunId" 
                                   << ", :testParticleSimulationId);";

    // Compile a SQL query.
    Statement randomWalkPerturberTableInsertQuery( 
        database, randomWalkPerturberTableInsert.str( ).c_str( ) );

    // Generate random walk input data and populate table.
    // Also populate perturber table for each Monte Carlo run.
    for ( unsigned int monteCarloRun = 0; monteCarloRun < monteCarloPopulation; monteCarloRun++ )
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

        // Get Monte Carlo run ID corresponding to row ID.
        ostringstream randomWalkMonteCarloRowIdQuery;
        randomWalkMonteCarloRowIdQuery << "SELECT \"monteCarloRunId\" FROM " 
                                       << randomWalkInputTableName 
                                       << " WHERE rowid == " << lastInputTableRowId;
        const int monteCarloRunId = database.execAndGet( 
            randomWalkMonteCarloRowIdQuery.str( ).c_str( ) );

        // Select simulation ID indices (test particle simulation indices in input 
        // table retrieved from database) to generate list of perturbers.

        // Declare vector containing test particle simulation IDs.
        vector< unsigned int > testParticleSimulationIdIndices( perturberPopulation );

        // Generate vector of randomly selected test particle simulation IDs.
        generate( testParticleSimulationIdIndices.begin( ), testParticleSimulationIdIndices.end( ),
                  generateTestParticleSimulationId );

        // Check if the test particle simulation IDs are unique, and if not, generate new random
        // number.
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
                // simulation ID index and restart looping.
                else if ( testParticleSimulationIdIndices.at( j ) 
                          == testParticleSimulationIdIndices.at( i ) )
                {
                    testParticleSimulationIdIndices.at( j ) = generateTestParticleSimulationId( );
                    i = 0;
                    break;
                }
            }
        }       

        // Loop through list of test particle simulation ID indices and write data to perturber 
        // table.
        for ( unsigned int k = 0; k < testParticleSimulationIdIndices.size( ); k++ ) 
        {
            // Bind values to prepared SQLite statement.
            randomWalkPerturberTableInsertQuery.bind( ":monteCarloRunId", monteCarloRunId );
            
            TestParticleInputTable::iterator iteratorTestParticleInputTable 
                = testParticleInputTable.begin( );
            std::advance( iteratorTestParticleInputTable, 
                testParticleSimulationIdIndices.at( k ) );

            randomWalkPerturberTableInsertQuery.bind( 
                ":testParticleSimulationId",
                 iteratorTestParticleInputTable->simulationId );  

            // Execute insert query.
            randomWalkPerturberTableInsertQuery.exec( );

            // Reset query.
            randomWalkPerturberTableInsertQuery.reset( );            
        }
    }

    // Commit transaction.
    randomWalkInputTableTransaction.commit( );

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
            << "\"key\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
            << "\"monteCarloRunId\" INTEGER NOT NULL,"
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
            tableCreateError << "Error: Creating table '" << randomWalkOutputTableName
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
    cout << endl;

    ///////////////////////////////////////////////////////////////////////////

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