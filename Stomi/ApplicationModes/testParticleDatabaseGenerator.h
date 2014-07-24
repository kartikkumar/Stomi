/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOMI_TEST_PARTICLE_DATABASE_GENERATOR_H
#define STOMI_TEST_PARTICLE_DATABASE_GENERATOR_H

#include <string>

#include <Tudat/InputOutput/parsedDataVectorUtilities.h>

namespace stomi
{
namespace application_modes
{

//! Execute test particle database generator application mode.
/*!
 * Executes test particle database generator application mode.
 * \param databasePath Absolute path to SQLite database.
 * \param parsedData Pointer to vector of data parsed from configuration file provided by user.
 */
void executeTestParticleDatabaseGenerator( 
    const std::string databasePath, 
    const tudat::input_output::parsed_data_vector_utilities::ParsedDataVectorPtr parsedData );

} // namespace application_modes
} // namespace stomi

#endif // STOMI_TEST_PARTICLE_DATABASE_GENERATOR_H
