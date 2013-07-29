/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120808    K. Kumar          File created.
 *      130217    K. Kumar          Updated "mab simulations" references to "stochastic migration".
 *      130704    K. Kumar          Updated function-naming.
 *
 *    References
 *
 *    Notes
 *      Unit tests need to be added to test the lexical dictionaries defined.
 *
 */

#ifndef STOCHASTIC_MIGRATION_DICTIONARIES_H
#define STOCHASTIC_MIGRATION_DICTIONARIES_H

#include <Tudat/InputOutput/dictionaryTools.h>

namespace stochastic_migration
{
namespace input_output
{

//! Get dictionary for databaseGenerator application.
/*!
 * Returns standard dictionary for StochasticMigrationDatabaseGenerator application.
 * \return Shared-pointer to new dictionary for StochasticMigrationDatabaseGenerator application.
 */
tudat::input_output::dictionary::DictionaryPointer getDatabaseGeneratorDictionary( );

//! Get dictionary for testParticleSimulator application.
/*!
 * Returns standard dictionary for TestParticleSimulator application.
 * \return Shared-pointer to new dictionary for TestParticleSimulator application.
 */
tudat::input_output::dictionary::DictionaryPointer getTestParticleSimulatorDictionary( );

//! Get dictionary for randomWalkSimulator application.
/*!
 * Returns standard dictionary for RandomWalkSimulator application.
 * \return Shared-pointer to new dictionary for RandomWalkSimulator application.
 */
tudat::input_output::dictionary::DictionaryPointer getRandomWalkSimulatorDictionary( );

} // namespace input_output
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_DICTIONARIES_H
