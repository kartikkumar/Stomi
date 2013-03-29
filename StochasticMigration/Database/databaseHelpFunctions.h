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
 *      The functions given here are used by both read and write database functions
 *      (StochasticMigration/Database/databaseReadFunctions.h &
 *      StochasticMigration/Database/databaseWriteFunctions.h ).
 *
 *      WARNING: There are no unit tests for these functions as yet.
 *
 */

#ifndef STOCHASTIC_MIGRATION_DATABASE_HELP_FUNCTIONS_H
#define STOCHASTIC_MIGRATION_DATABASE_HELP_FUNCTIONS_H

#include <string>

#include <Assist/Database/sqlite3DatabaseConnector.h>

namespace stochastic_migration
{
namespace database
{

//! Initiate SQLite3 database connector.
/*!
 * Initiates a new database connector, with the databased opened based on the input path. An SQLite
 * transaction is also started. Database connector should be terminated after use by calling
 * terminateDatabaseConnector().
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \return Shared-pointer to new database connector.
 * \sa initiateDatabaseConnector().
 */
assist::database::Sqlite3DatabaseConnectorPointer initiateDatabaseConnector(
        const std::string& databaseAbsolutePath );

//! Terminate SQLite3 database connector.
/*!
 * Terminates an existing database connector safely. Also ends an SQLITE transaction that is
 * assumed to have been started when initiateDatabaseConnector() is called.
 * \param databaseConnector Shared-pointer to active database connector.
 * \sa initiateDatabaseConnector().
 */
void terminateDatabaseConnector(
        assist::database::Sqlite3DatabaseConnectorPointer databaseConnector );

//! Throw database error, including status returned by handler.
/*!
 * Throws a run-time error that includes the status returned by the SQLite handler. The error code
 * returned corresponds to the standard SQLite error codes (SQLite C Interface, 2013). Before the
 * run-time error is thrown, the database is closed cleanly by calling
 * terminateDatabaseConnector().
 * \param databaseConnector Shared-pointer to active database connector.
 * \param databaseStatus Error code status of SQLite database handler.
 * \sa terminateDatabaseConnector().
 */
void throwDatabaseError(
        assist::database::Sqlite3DatabaseConnectorPointer databaseConnector,
        const unsigned int databaseStatus );

} // namespace database
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_DATABASE_HELP_FUNCTIONS_H
