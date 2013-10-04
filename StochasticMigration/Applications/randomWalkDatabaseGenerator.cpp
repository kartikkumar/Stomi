/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string> 

#include <boost/filesystem.hpp> 

#include <SQLiteC++.h> 

#include <Tudat/InputOutput/dictionaryTools.h>
#include <Tudat/InputOutput/fieldType.h>
#include <Tudat/InputOutput/separatedParser.h>
#include <Tudat/InputOutput/parsedDataVectorUtilities.h>

#include <Assist/Astrodynamics/unitConversions.h>
#include <Assist/InputOutput/basicInputOutput.h> 

#include "StochasticMigration/InputOutput/dictionaries.h" 

//! Execute random walk database generator.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{

    ///////////////////////////////////////////////////////////////////////////

    // Declare using-statements.

    using std::cout;
    using std::endl;
    using std::ostringstream;
    using std::runtime_error;
    using std::string;

    using namespace boost::filesystem;

    using namespace SQLite;

    using namespace assist::astrodynamics;
    using namespace assist::input_output;

    using namespace tudat::input_output::dictionary;
    using namespace tudat::input_output;
    using namespace tudat::input_output::field_types::general;
    using namespace tudat::input_output::parsed_data_vector_utilities;

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

    const std::string randomWalkCaseTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKCASETABLENAME" ),
                "random_walk_case" );
    cout << "Random walk case table                                    " 
         << randomWalkCaseTableName << endl;

    const std::string randomWalkPerturberTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKPERTURBERTABLENAME" ),
                "random_walk_perturbers" );
    cout << "Random walk perturber table                               "
         << randomWalkPerturberTableName << endl;

    const std::string randomWalkOutputTableName = extractParameterValue< std::string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "RANDOMWALKOUTPUTTABLENAME" ), "random_walk_output" );
    cout << "Random walk output table                                  "
         << randomWalkOutputTableName << endl;

    const string testParticleCaseTableName = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ),
                findEntry( dictionary, "TESTPARTICLECASETABLENAME" ), "test_particle_case" );
    cout << "Test particle case table                                  "
         << testParticleCaseTableName << endl;         

    // Check that all required parameters have been set.
    checkRequiredParameters( dictionary );

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
            << "\"monteCarloPopulation\" INTEGER NOT NULL,"
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

        // Get ID for test particle case.
        ostringstream testParticleCaseIdQuery;
        testParticleCaseIdQuery << "SELECT caseId FROM " << testParticleCaseTableName 
                                << " WHERE \"caseName\" = \"" << testParticleCaseName << "\"";
        const int testParticleCaseId = database.execAndGet( 
            testParticleCaseIdQuery.str( ).c_str( ) );

        // Create stringstream with random walk case data insert command.
        // For floating-point values, ensure the data is written to the stream at full precision.
        ostringstream randomWalkCaseDataInsert;
        randomWalkCaseDataInsert
            << "INSERT INTO " << randomWalkCaseTableName << " VALUES ("
            << "NULL,"
            << "\"" << caseName << "\","
            << testParticleCaseId << ","
            << monteCarloPopulation << ",";
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

//     ///////////////////////////////////////////////////////////////////////////

//     // Set up random walk perturber table.
//     if ( !database.tableExists( randomWalkPerturberTableName.c_str( ) ) )
//     {
//         cout << "Table '" << randomWalkPerturberTableName << "' does not exist ..." 
//                   << endl;
//         cout << "Creating table ... " << endl;

//         // Create table.
//         ostringstream randomWalkPerturberTableCreate;
//         randomWalkPerturberTableCreate
//             << "CREATE TABLE IF NOT EXISTS " << randomWalkPerturberTableName << " ("
//             << "\"perturberId\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
//             << "\"monteCarloRunId\" INTEGER NOT NULL,"
//             << "\"testParticleSimulationId\" INTEGER NOT NULL,"
//             << "\"massFactor\" REAL NOT NULL);";

// // << "\"observationPeriodStartEpoch\" REAL NOT NULL

//         // Execute command to create table.
//         database.exec( randomWalkPerturberTableCreate.str( ).c_str( ) );

//         // Check that the table was created successfully.
//         if ( database.tableExists( randomWalkPerturberTableName.c_str( ) ) )
//         {
//             cout << "Table '" << randomWalkPerturberTableName 
//                       << "' successfully created!" << endl; 
//         }

//         else
//         {
//             ostringstream tableCreateError;
//             tableCreateError << "Error: Creating table '" << randomWalkPerturberTableName 
//             << "'' failed!";
//             throw runtime_error( tableCreateError.str( ).c_str( ) );
//         }
//     }

//     else
//     {
//         cout << "Table '" << randomWalkPerturberTableName 
//                   << "' already exists ... skipping creating table ..." << endl;
//     }

//     ///////////////////////////////////////////////////////////////////////////

//     ///////////////////////////////////////////////////////////////////////////

//     // Set up random walk output table.
//     if ( !database.tableExists( randomWalkOutputTableName.c_str( ) ) )
//     {
//         cout << "Table '" << randomWalkOutputTableName << "' does not exist ..." << endl;
//         cout << "Creating table ... " << endl;

//         // Create table.
//         ostringstream randomWalkOutputTableCreate;
//         randomWalkOutputTableCreate
//             << "CREATE TABLE IF NOT EXISTS " << randomWalkOutputTableName << " ("
//             << "\"key\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
//             << "\"monteCarloRunId\" INTEGER NOT NULL,"
//             << "\"maximumEccentricityChange\" REAL NOT NULL,"
//             << "\"maximumLongitudeResidualChange\" REAL NOT NULL,"
//             << "\"maximumInclinationChange\" REAL NOT NULL);";

//         // Execute command to create table.
//         database.exec( randomWalkOutputTableCreate.str( ).c_str( ) );

//         // Check that the table was created successfully.
//         if ( database.tableExists( randomWalkOutputTableName.c_str( ) ) )
//         {
//             cout << "Table '" << randomWalkOutputTableName 
//                       << "' successfully created!" << endl; 
//         }

//         else
//         {
//             ostringstream tableCreateError;
//             tableCreateError << "Error: Creating table '" << randomWalkOutputTableName 
//             << "'' failed!";
//             throw runtime_error( tableCreateError.str( ).c_str( ) );
//         }
//     }

//     else
//     {
//         cout << "Table '" << randomWalkOutputTableName 
//                   << "' already exists ... skipping creating table ..." << endl;
//     }

//     ///////////////////////////////////////////////////////////////////////////


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