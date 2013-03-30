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
 *      120402    K. Kumar          File created from old mabSimulationDatabaseFunctions.h.
 *      130212    K. Kumar          Added planetary_rings namespace; updated Doxygen documentation.
 *      130214    K. Kumar          Split file to contain only read-functions (write-functions
 *                                  ported to new file).
 *      130217    K. Kumar          Updated "mab simulations" references to "stochastic migration".
 *
 *    References
 *      SQLite C Interface. Result codes, http://www.sqlite.org/c3ref/c_abort.html, last accessed:
 *          17th Feb, 2013.
 *
 *    Notes
 *
 */

#ifndef STOCHASTIC_MIGRATION_DATABASE_READ_FUNCTIONS_H
#define STOCHASTIC_MIGRATION_DATABASE_READ_FUNCTIONS_H

#include <string>

// #include "StochasticMigration/Database/databaseHelpFunctions.h"
// #include "StochasticMigration/Database/randomWalkMonteCarloRun.h"
// #include "StochasticMigration/Database/randomWalkPerturber.h"
#include "StochasticMigration/Database/testParticleCase.h"
#include "StochasticMigration/Database/testParticleInput.h"
// #include "StochasticMigration/Database/testParticleKick.h"

namespace stochastic_migration
{
namespace database
{

//! Get test particle case data.
/*!
 * Returns case data, used as metadata for a set of test particle simulations, retrieved from
 * simulation database.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param testParticleCaseTableName String name of test particle case table in database (default is
 *          set to "test_particle_case").
 * \return Test particle case data, stored in a shared-point to a TestParticleCase object.
 */
TestParticleCasePointer getTestParticleCase(
        const std::string& databaseAbsolutePath,
        const std::string& testParticleCaseTableName = "test_particle_case" );

//! Get test particle input table.
/*!
 * Returns table of input data for test particle simulations, retrieved from simulation database.
 * Completed (true) or incomplete (false) simulations can be retrieved.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param isCompleted Flag indicating whether to retrieve completed or incomplete simulations
 *          (default is set to retrieve input data for simulations that haven't been completed).
 * \param testParticleInputTableName String name of test particle input table in database (default
 *          is set to "test_particle_input").
 * \return Test particle input table, as a set of TestParticleInput pointers.
 */
// TestParticleInputTable getTestParticleInputTable(
//         const std::string& databaseAbsolutePath, bool isCompleted = false,
//         const std::string& testParticleInputTableName = "test_particle_input" );

// //! Get test particle input table.
// /*!
//  * Returns table of input data for test particle simulations, retrieved from simulation database.
//  * A list of specific test particle simulation numbers can be requested to be retrieved from the
//  * SQLite database.
//  * \param databaseAbsolutePath Absolute path to simulation database.
//  * \param testParticleSimulationNumbers List of test particle simulation numbers to select from
//  *          database.
//  * \param testParticleInputTableName String name of test particle input table in database (default
//  *          is set to "test_particle_input").
//  * \return Test particle input table, as a vector of TestParticleInput objects.
//  */
// TestParticleInputTable getTestParticleInputTable(
//         const std::string& databaseAbsolutePath, const std::string& testParticleSimulationNumbers,
//         const std::string& testParticleInputTableName = "test_particle_input" );

// //! Get test particle kick table.
// !
//  * Returns table of test particle kick data, aggregated based on specified test particle simulation
//  * numbers.
//  * \param databaseAbsolutePath Absolute path to simulation database.
//  * \param randomWalkDuration Duration of random walk [s].
//  * \param testParticleSimulationNumbersAndMassFactors Map of specified test particle simulation
//  *          numbers to include in kick table and mass factors associated with each.
//  * \param testParticleKickTableName String name of test particle kick table in database (default
//  *          is set to "test_particle_kicks").
//  * \return Table of test particle kick data, aggregated from selected test particle simulation
//  *          numbers, as a vector of TestParticleKick objects.
 
// TestParticleKickTable getTestParticleKickTable(
//         const std::string& databaseAbsolutePath, const double randomWalkDuration,
//         const TestParticleSimulationNumbersAndMassFactors&
//         testParticleSimulationNumbersAndMassFactors,
//         const std::string& testParticleKickTableName = "test_particle_kicks" );

// //! Get table of random walk Monte Carlo runs.
// /*!
//  * Returns random_walk_runs table retrieved from SQLite simulation database, aggregated from the
//  * Monte Carlo runs requested.
//  * \param databaseAbsolutePath Absolute path to simulation database.
//  * \param monteCarloRun Monte Carlo runs to retrieve data for.
//  * \param randomWalkMonteCarloRunTableName String name of random walk Monte Carlo run table in
//  *          database (default is set to "random_walk_monte_carlo_runs").
//  * \return Monte Carlo run table, stored in a vector of RandomWalkMonteCarloRun objects.
//  */
// RandomWalkMonteCarloRunTable getRandomWalkMonteCarloRunsTable(
//         const std::string& databaseAbsolutePath, const std::vector< unsigned int >& monteCarloRuns,
//         const std::string& randomWalkMonteCarloRunTableName = "random_walk_monte_carlo_runs" );

// //! Get table of selected perturbers for random walk Monte Carlo run.
// /*!
//  * Returns table of selected perturbers retrieved from SQLite simulation database, aggregated from
//  * the test particle simulation numbers and associated mass factors for the requested Monte Carlo
//  * run.
//  * \param databaseAbsolutePath Absolute path to simulation database.
//  * \param monteCarloRun Monte Carlo run to retrieve list of selected perturbers for.
//  * \param randomWalkPerturberTableName String name of random walk perturber table in database
//  *          (default is set to "random_walk_perturber_selection").
//  * \return Random walk perturber selection table, stored in a vector of RandomWalkPerturber
//  *          objects.
//  */
// RandomWalkPerturberTable getRandomWalkPerturberTable(
//         const std::string& databaseAbsolutePath, const unsigned int monteCarloRun,
//         const std::string& randomWalkPerturberTableName = "random_walk_perturbers" );

} // namespace database
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_DATABASE_READ_FUNCTIONS_H
