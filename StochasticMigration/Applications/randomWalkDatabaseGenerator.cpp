/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130924    K. Kumar          File created based on code in old databaseGenerator.cpp file.
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


//! Execute random walk database generator.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{

    //     const std::string randomWalkMonteCarloRunTableName = extractParameterValue< std::string >(
//                 parsedData->begin( ), parsedData->end( ),
//                 findEntry( dictionary, "RANDOMWALKMONTECARLORUNTABLENAME" ),
//                 "random_walk_monte_carlo_runs" );
//     cout << "Random walk Monte Carlo run table                         "
//               << randomWalkMonteCarloRunTableName << endl;

//     const std::string randomWalkPerturberTableName = extractParameterValue< std::string >(
//                 parsedData->begin( ), parsedData->end( ),
//                 findEntry( dictionary, "RANDOMWALKPERTURBERTABLENAME" ),
//                 "random_walk_perturbers" );
//     cout << "Random walk perturber table                               "
//               << randomWalkPerturberTableName << endl;

//     const std::string randomWalkOutputTableName = extractParameterValue< std::string >(
//                 parsedData->begin( ), parsedData->end( ),
//                 findEntry( dictionary, "RANDOMWALKOUTPUTTABLENAME" ), "random_walk_output" );
//     cout << "Random walk output table                                  "
//               << randomWalkOutputTableName << endl;




    //     ///////////////////////////////////////////////////////////////////////////

//     // Set up random walk Monte Carlo run table.
//     if ( !database.tableExists( randomWalkMonteCarloRunTableName.c_str( ) ) )
//     {
//         cout << "Table '" << randomWalkMonteCarloRunTableName << "' does not exist ..." 
//                   << endl;
//         cout << "Creating table ... " << endl;

//         // Create table.
//         ostringstream randomWalkMonteCarloRunTableCreate;
//         randomWalkMonteCarloRunTableCreate
//             << "CREATE TABLE IF NOT EXISTS " << randomWalkMonteCarloRunTableName << " ("
//             << "\"monteCarloRunId\" INTEGER PRIMARY KEY NOT NULL,"
//             << "\"caseId\" INTEGER NOT NULL,"
//             << "\"perturberDensity\" INTEGER NOT NULL,"
//             << "\"perturberRingMass\" INTEGER NOT NULL,"            
//             << "\"observationPeriod\" REAL NOT NULL,"
//             << "\"numberOfEpochs\" REAL NOT NULL,"
//             << "\"epochWindowSize\" REAL NOT NULL);";

//         // Execute command to create table.
//         database.exec( randomWalkMonteCarloRunTableCreate.str( ).c_str( ) );

//         // Check that the table was created successfully.
//         if ( database.tableExists( randomWalkMonteCarloRunTableName.c_str( ) ) )
//         {
//             cout << "Table '" << randomWalkMonteCarloRunTableName 
//                       << "' successfully created!" << endl; 
//         }

//         else
//         {
//             ostringstream tableCreateError;
//             tableCreateError << "Error: Creating table '" << randomWalkMonteCarloRunTableName 
//             << "'' failed!";
//             throw runtime_error( tableCreateError.str( ).c_str( ) );
//         }
//     }

//     else
//     {
//         cout << "Table '" << randomWalkMonteCarloRunTableName 
//                   << "' already exists ... skipping creating table ..." << endl;
//     }

//     ///////////////////////////////////////////////////////////////////////////

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
}