/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOMI_RANDOM_WALK_DATABASE_GENERATOR_H
#define STOMI_RANDOM_WALK_DATABASE_GENERATOR_H

#include <string>

#include <Tudat/InputOutput/parsedDataVectorUtilities.h>

namespace stomi
{
namespace application_modes
{

//! Execute random walk database generator application mode.
/*!
 * Executes random walk database generator application mode.
 * \param databasePath Absolute path to SQLite database.
 * \param parsedData Pointer to vector of data parsed from configuration file provided by user.
 */
void executeRandomWalkDatabaseGenerator( 
    const std::string databasePath, 
    const tudat::input_output::parsed_data_vector_utilities::ParsedDataVectorPtr parsedData );

} // namespace application_modes
} // namespace stomi

#endif // STOMI_RANDOM_WALK_DATABASE_GENERATOR_H
