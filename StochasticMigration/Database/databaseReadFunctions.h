/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120402    K. Kumar          File created from old mabSimulationDatabaseFunctions.h.
 *      130212    K. Kumar          Added planetary_rings namespace; updated Doxygen documentation.
 *      130214    K. Kumar          Split file to contain only read-functions (write-functions
 *                                  ported to new file).
 *      130217    K. Kumar          Updated "mab simulations" references to "stochastic migration".
 *      130918    K. Kumar          Uncommented and update function to get kick table.
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
#include <vector>

// #include "StochasticMigration/Database/randomWalkMonteCarloRun.h"
// #include "StochasticMigration/Database/randomWalkPerturber.h"
#include "StochasticMigration/Database/testParticleCase.h"
#include "StochasticMigration/Database/testParticleInput.h"
#include "StochasticMigration/Database/testParticleKick.h"

namespace stochastic_migration
{
namespace database
{

//! Get test particle case data.
/*!
 * Returns case data, used as metadata for a set of test particle simulations, retrieved from
 * simulation database.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param caseName Name stored in table for test particle case.
 * \param testParticleCaseTableName String name of test particle case table in database.
 * \return Test particle case data, stored in a shared-point to a TestParticleCase object.
 */
TestParticleCasePointer getTestParticleCase( const std::string& databaseAbsolutePath, 
                                             const std::string& caseName,
                                             const std::string& testParticleCaseTableName );

//! Get complete test particle input table.
/*!
 * Returns table of input data for test particle simulations, retrieved from simulation database,
 * for a given case, defined by a string-name.
 * Completed (true) or incomplete (false) simulations can be retrieved.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param caseId ID stored in table for requested test particle case.
 * \param testParticleInputTableName String name of test particle input table in database.
 * \param isCompleted Flag indicating whether to retrieve completed or incomplete simulations
 *          (default is set to retrieve input data for simulations that haven't been completed).
 * \return Test particle input table, as a set of TestParticleInput pointers.
 */
TestParticleInputTable getCompleteTestParticleInputTable(
        const std::string& databaseAbsolutePath, const int caseId,
        const std::string& testParticleInputTableName, bool isCompleted = false );

//! Get selected test particle input table.
/*!
 * Returns table of input data for test particle simulations, retrieved from simulation database,
 * for a given case, defined by a string-name.
 * A list of specific test particle simulation numbers can be requested to be retrieved from the
 * SQLite database.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param caseId ID stored in table for requested test particle case.
 * \param testParticleSimulationNumbers List of test particle simulation numbers to select from
 *          database.
 * \param testParticleInputTableName String name of test particle input table in database.
 * \return Test particle input table, as a vector of TestParticleInput objects.
 */
TestParticleInputTable getSelectedTestParticleInputTable(
        const std::string& databaseAbsolutePath, const int caseId,
        const std::string& testParticleSimulationNumbers,
        const std::string& testParticleInputTableName );

//! Get test particle kick table.
/*!
 * Returns table of test particle kick data, aggregated based on specified test particle simulation
 * IDs.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param randomWalkSimulationPeriod Duration of random walk [s].
 * \param testParticleSimulationIds Vector of specified test particle simulation IDs to include in 
 *          kick table.
 * \param testParticleKickTableName String name of test particle kick table in database.
 * \return Table of test particle kick data, aggregated from selected test particle simulation IDs.
 */
TestParticleKickTable getTestParticleKickTable(
        const std::string& databaseAbsolutePath, const double randomWalkSimulationPeriod, 
        const std::vector< int >& selectedSimulationIds, 
        const std::string& testParticleKickTableName );

// //! Get table of random walk Monte Carlo runs.
// /*!
//  * Returns random_walk_runs table retrieved from SQLite simulation database, aggregated from the
//  * Monte Carlo runs requested.
//  * \param databaseAbsolutePath Absolute path to simulation database.
//  * \param monteCarloRuns Monte Carlo runs to retrieve data for.
//  * \param randomWalkMonteCarloRunTableName String name of random walk Monte Carlo run table in
//  *          database (default is set to "random_walk_monte_carlo_runs").
//  * \return Monte Carlo run table, stored in a set of RandomWalkMonteCarloRun objects.
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
//  * \return Random walk perturber selection table, stored in a set of RandomWalkPerturber
//  *          objects.
//  */
// RandomWalkPerturberTable getRandomWalkPerturberTable(
//         const std::string& databaseAbsolutePath, const unsigned int monteCarloRun,
//         const std::string& randomWalkPerturberTableName = "random_walk_perturbers" );

} // namespace database
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_DATABASE_READ_FUNCTIONS_H
