/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/timer/timer.hpp>

#include <Assist/InputOutput/basicInputOutput.h>

#include <Tudat/InputOutput/dictionaryTools.h>
#include <Tudat/InputOutput/fieldType.h>
#include <Tudat/InputOutput/separatedParser.h>
#include <Tudat/InputOutput/parsedDataVectorUtilities.h>

#include "Stomi/ApplicationModes/randomWalkDatabaseGenerator.h"
#include "Stomi/ApplicationModes/randomWalkSimulator.h" 
#include "Stomi/ApplicationModes/testParticleDatabaseGenerator.h"
#include "Stomi/ApplicationModes/testParticleSimulator.h"
#include "Stomi/InputOutput/dictionaries.h"

//! Execute Stomi.
int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    ///////////////////////////////////////////////////////////////////////////

    // Start timer. Timer automatically ends when this object goes out of scope.
    boost::timer::auto_cpu_timer timer;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Declare using-statements.
    using std::cout;
    using std::endl;
    using std::runtime_error;
    using std::string;

    using boost::iequals;

    using namespace assist::input_output;

    using namespace tudat::input_output;
    using namespace tudat::input_output::dictionary;
    using namespace tudat::input_output::field_types::general;
    using namespace tudat::input_output::parsed_data_vector_utilities;

    using namespace stomi::application_modes;
    using namespace stomi::input_output;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up input deck.

    // Check number of input parameters is correct (the numberOfInputs variable includes the
    // application itself, so one is subtracted from this number).
    checkNumberOfInputArguments( numberOfInputs - 1 );

    // Read and filter input stream (this can't be declared const because the parser's parse
    // function is not const-correct at the moment).
    string filteredInput = readAndFilterInputFile( inputArguments[ 1 ] );

    // Declare a separated parser.
    SeparatedParser parser( string( ": " ), 2, parameterName, parameterValue );

    // Parse filtered data.
    const ParsedDataVectorPtr parsedData = parser.parse( filteredInput );

    // Get general parameters dictionary.
    DictionaryPointer dictionary = getGeneralParametersDictionary( );

    // Extract application mode.
    const string applicationMode = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), findEntry( dictionary, "MODE" ) );

    const string databasePath = extractParameterValue< string >(
                parsedData->begin( ), parsedData->end( ), 
                findEntry( dictionary, "DATABASEPATH" ) );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Call selected application mode.
    cout << endl;
    cout << "----------------------------------------------------------------------------" << endl;
    cout << endl;
    cout << "                                    STOMI                                   " << endl;
    cout << "                                    2.0.0                                   " << endl;
    cout << endl;
    cout << "          Copyright (c) 2010-2014, Delft University of Technology           " << endl;
    cout << "           Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)           " << endl;
    cout << endl;
    cout << "----------------------------------------------------------------------------" << endl;

    cout << endl;
    cout << "****************************************************************************" << endl;
    cout << "Input parameters" << endl;
    cout << "****************************************************************************" << endl;
    cout << endl;

    cout << endl;
    cout << "Application mode                                          ";

    if ( iequals( applicationMode, "TPDB" ) )
    {
        cout << "TEST PARTICLE DATABASE GENERATOR" << endl;
        executeTestParticleDatabaseGenerator( databasePath, parsedData );
    }    

    else if ( iequals( applicationMode, "TPSIM" ) )
    {
        cout << "TEST PARTICLE SIMULATOR" << endl;
        executeTestParticleSimulator( databasePath, parsedData );             
    } 

    else if ( iequals( applicationMode, "RWDB" ) )
    {
        cout << "RANDOM WALK DATABASE GENERATOR" << endl;
        executeRandomWalkDatabaseGenerator( databasePath, parsedData );        
    }    

    else if ( iequals( applicationMode, "RWSIM" ) )
    {
        cout << "RANDOM WALK SIMULATOR" << endl;
        executeRandomWalkSimulator( databasePath, parsedData );             
    }    

    else
    {
        throw runtime_error( "ERROR: Application mode not recognized!" );
    }

    cout << endl;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Print timing information.
    cout << "Timing information: ";
    
    // If program is successfully completed, return 0.
    return EXIT_SUCCESS;

    ///////////////////////////////////////////////////////////////////////////
}

/*    
 *    TODO:
 *      - encapsulate code in try-catch blocks to capture exceptions.
 *      - execute verification of existing data against input parameters provided to ensure
 *        consistency of inputs and possibly warn user.
 *      - expand code to enable interactive interface for user to provide inputs and select
 *        options (capture command line user input). 
 */ 
