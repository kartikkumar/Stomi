/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOMI_DATABASE_READ_FUNCTIONS_H
#define STOMI_DATABASE_READ_FUNCTIONS_H

#include <string>
#include <vector>

#include "StoMi/Database/randomWalkCase.h"
#include "StoMi/Database/randomWalkInput.h"
#include "StoMi/Database/testParticleCase.h"
#include "StoMi/Database/testParticleInput.h"
#include "StoMi/Database/testParticleKick.h"

namespace stomi
{
namespace database
{

//! Get test particle case data.
/*!
 * Returns case data, used as metadata for a set of test particle simulations, retrieved from
 * simulation database.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param caseId ID stored in table for test particle case.
 * \param testParticleCaseTableName String name of test particle case table in database.
 * \return Test particle case data, stored in a shared-point to a TestParticleCase object.
 */
TestParticleCasePointer getTestParticleCase( const std::string& databaseAbsolutePath, 
                                             const int caseId,
                                             const std::string& testParticleCaseTableName );

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
 * A list of specific test particle simulation IDs can be requested to be retrieved from the SQLite
 * database.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param caseId ID stored in table for requested test particle case.
 * \param testParticleSimulationIds List of test particle simulation IDs to select from database.
 * \param testParticleInputTableName String name of test particle input table in database.
 * \return Test particle input table, as a vector of TestParticleInput objects.
 */
TestParticleInputTable getSelectedTestParticleInputTable(
        const std::string& databaseAbsolutePath, const int caseId,
        const std::string& testParticleSimulationIds,
        const std::string& testParticleInputTableName );

//! Get test particle kick table.
/*!
 * Returns table of test particle kick data, aggregated based on specified test particle simulation
 * IDs.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param randomWalkSimulationPeriod Duration of random walk [s].
 * \param selectedSimulationIds Vector of specified test particle simulation IDs to include in kick
 *          table.
 * \param testParticleKickTableName String name of test particle kick table in database.
 * \return Table of test particle kick data, aggregated from selected test particle simulation IDs.
 */
TestParticleKickTable getTestParticleKickTable(
        const std::string& databaseAbsolutePath, const double randomWalkSimulationPeriod, 
        const std::vector< int >& selectedSimulationIds, 
        const std::string& testParticleKickTableName );

//! Get random walk case data.
/*!
 * Returns case data, used as metadata for a set of random walk simulations, retrieved from
 * simulation database.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param caseName Name stored in table for random walk case.
 * \param randomWalkCaseTableName String name of random walk case table in database.
 * \return Random Walk case data, stored in a shared-point to a RandomWalkCase object.
 */
RandomWalkCasePointer getRandomWalkCase( const std::string& databaseAbsolutePath, 
                                         const std::string& caseName,
                                         const std::string& randomWalkCaseTableName );

//! Get complete random walk input table.
/*!
 * Returns table of input data for random walk simulations, retrieved from simulation database,
 * for a given case, defined by a string-name.
 * Completed (true) or incomplete (false) simulations can be retrieved.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param caseId ID stored in table for requested random walk case.
 * \param randomWalkInputTableName String name of random walk input table in database.
 * \param randomWalkPerturberTableName String name of random walk perturber table in database.
 * \param isCompleted Flag indicating whether to retrieve completed or incomplete simulations
 *          (default is set to retrieve input data for simulations that haven't been completed).
 * \return Random walk input table, as a set of RandomWalkInput pointers.
 */
RandomWalkInputTable getCompleteRandomWalkInputTable(
        const std::string& databaseAbsolutePath, const int caseId,
        const std::string& randomWalkInputTableName, 
        const std::string& randomWalkPerturberTableName, bool isCompleted = false );

//! Get selected random walk input table.
/*!
 * Returns table of input data for random walk simulations, retrieved from simulation database,
 * for a given case, defined by a string-name.
 * A list of specific random walk Monte Carlo run IDs can be requested to be retrieved from the
 * SQLite database.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param caseId ID stored in table for requested random walk case.
 * \param monteCarloRunIds List of Monte Carlo run IDs to select from database.
 * \param randomWalkInputTableName String name of random walk input table in database.
 * \param randomWalkPerturberTableName String name of random walk perturber table in database. 
 * \return Random walk input table, as a vector of RandomWalkInput objects.
 */
RandomWalkInputTable getSelectedRandomWalkInputTable(
        const std::string& databaseAbsolutePath, const int caseId,
        const std::string& monteCarloRunIds,
        const std::string& randomWalkInputTableName,
        const std::string& randomWalkPerturberTableName );

//! Get list of selected perturbers for random walk Monte Carlo run.
/*!
 * Returns selected perturbers, given as test particle simulation IDs retrieved from SQLite
 * simulation database, for the specified Monte Carlo run ID.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param monteCarloRunId Monte Carlo run ID to retrieve list of selected perturbers for.
 * \param randomWalkPerturberTableName String name of random walk perturber table in database.
 * \return Vector containing list of test particle simulation IDs.
 */
std::vector< int > getRandomWalkPerturberList(
        const std::string& databaseAbsolutePath, const unsigned int monteCarloRunId,
        const std::string& randomWalkPerturberTableName );

} // namespace database
} // namespace stomi

#endif // STOMI_DATABASE_READ_FUNCTIONS_H

/*
 *    References
 *      SQLite C Interface. Result codes, http://www.sqlite.org/c3ref/c_abort.html, last accessed:
 *          17th Feb, 2013.
 */
