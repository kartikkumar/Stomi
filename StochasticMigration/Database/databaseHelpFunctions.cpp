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
 *      130217    K. Kumar          File created.
 *      130329    K. Kumar          Updated references to GeneralTools to Assist.
 *
 *    References
 *
 *    Notes
 *
 */

#include <sstream>
#include <stdexcept>

#include <boost/make_shared.hpp>

#include "StochasticMigration/Database/databaseHelpFunctions.h"

namespace stochastic_migration
{
namespace database
{

//! Using statements.
using std::runtime_error;

using boost::enable_error_info;
using boost::make_shared;
using boost::throw_exception;

using namespace assist::database;

//! Initiate SQLite3 database connector.
Sqlite3DatabaseConnectorPointer initiateDatabaseConnector(
        const std::string& databaseAbsolutePath )
{
    // Open database connection.
    Sqlite3DatabaseConnectorPointer databaseConnector
            = make_shared< Sqlite3DatabaseConnector >( databaseAbsolutePath );

    // Begin database transaction.
    databaseConnector->beginTransaction( );

    // Return database connector.
    return databaseConnector;
}

//! Terminate SQLite3 database connector.
void terminateDatabaseConnector( Sqlite3DatabaseConnectorPointer databaseConnector )
{
    // Finalize fetched data to ensure memory is not leaked.
    databaseConnector->finalizeStatement( );

    // End database transaction.
    databaseConnector->endTransaction( );

    // Disconnect from SQLite database.
    databaseConnector->closeDatabase( );
}

//! Throw database error, including status returned by handler.
void throwDatabaseError( Sqlite3DatabaseConnectorPointer databaseConnector,
                         const unsigned int databaseStatus )
{
    // Terminate database connector cleanly.
    terminateDatabaseConnector( databaseConnector );

    // Set error message.
    std::ostringstream errorMessage;
    errorMessage << "SQLite3 Error: " << databaseStatus << "!";

    // Throw run-time error.
    throw runtime_error( errorMessage.str( ) );
}

} // namespace database
} // namespace stochastic_migration
