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
 *      120808    K. Kumar          File created.
 *      130217    K. Kumar          Updated "mab simulations" references to "stochastic migration".
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

//! Get dictionary for StochasticMigrationDatabaseGenerator application.
/*!
 * Returns standard dictionary for StochasticMigrationDatabaseGenerator application.
 * \return Shared-pointer to new dictionary for StochasticMigrationDatabaseGenerator application.
 */
tudat::input_output::dictionary::DictionaryPointer
getStochasticMigrationDatabaseGeneratorDictionary( );

//! Get dictionary for TestParticleSimulator application.
/*!
 * Returns standard dictionary for TestParticleSimulator application.
 * \return Shared-pointer to new dictionary for TestParticleSimulator application.
 */
tudat::input_output::dictionary::DictionaryPointer getTestParticleSimulatorDictionary( );

//! Get dictionary for RandomWalkSimulator application.
/*!
 * Returns standard dictionary for RandomWalkSimulator application.
 * \return Shared-pointer to new dictionary for RandomWalkSimulator application.
 */
tudat::input_output::dictionary::DictionaryPointer getRandomWalkSimulatorDictionary( );

} // namespace input_output
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_DICTIONARIES_H
