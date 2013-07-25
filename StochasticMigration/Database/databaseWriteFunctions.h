/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130214    K. Kumar          File created from old databaseFunctions.h.
 *      130217    K. Kumar          updated "mab simulations" references to "stochastic migration".
 *
 *    References
 *
 *    Notes
 *
 */

#include <string>

#include "StochasticMigration/Database/testParticleKick.h"

#ifndef STOCHASTIC_MIGRATION_DATABASE_WRITE_FUNCTIONS_H
#define STOCHASTIC_MIGRATION_DATABASE_WRITE_FUNCTIONS_H

namespace stochastic_migration
{
namespace database
{

//! Populate test particle kick table.
/*!
 * Populates table of test particle kicks in SQLite database.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param kickTable Table containing kick data as set of pointers to TestParticleKick objects.
 * \param testParticleKickTableName String name of test particle kick table in database.
 * \param testParticleInputTableName String name of test particle input table in database.
 */
void populateTestParticleKickTable( const std::string& databaseAbsolutePath,
                                    const TestParticleKickTable& kickTable,
                                    const std::string& testParticleKickTableName, 
                                    const std::string& testParticleInputTableName );

//! Populate random walk Monte Carlo run and output tables.
/*!
 * Populates SQLite database with random walk Monte Carlo run and output data. The random walk
 * Monte Carlo run serves as metadata for the random walk simulations. The random walk Monte Carlo
 * output data is the data generated through execution of the Monte Carlo run.
 * \param databaseAbsolutePath Absolute path to simulation database.
 * \param testParticleSimulationNumbersAndMassFactors Map of specified test particle simulation
 *          numbers and mass factors associated with each, used during execution of random walk
 *          Monte Carlo simulation.
 * \param massDistributionType Mass distribution type (EQUAL, UNIFORM, LINEAR, or POWERLAW). [ONLY
 *          EQUAL IMPLEMENTED AT THE MOMENT].
 * \param massDistributionParameters Mass distribution parameters (type EQUAL has only one
 *          parameter).
 * \param maximumEccentricityChange Maximum eccentricity change that perturber body underwent
 *          during random walk Monte Carlo simulation. The eccentricity change is computed based on
 *          the method that the observation metrics are specified.
 * \param maximumLongitudeResidualChange Maximum longitude residual changes that perturber body
 *          underwent during random walk Monte Carlo simulation [rad]. The longitude residual
 *          change is computed based on the method that the observation metrics are specified.
 * \param maximumInclinationChange Maximum inclination changes that perturber body underwent during
 *          random walk Monte Carlo simulation [rad]. The inclination change is computed based on
 *          the method that the observation metrics are specified.
 */
// void populateRandomWalkRunAndOutputTables(
//         const std::string& databaseAbsolutePath,
//         const TestParticleSimulationNumbersAndMassFactors&
//         testParticleSimulationNumbersAndMassFactors,
//         const std::string& massDistributionType,
//         const std::vector< double > massDistributionParameters, const double observationPeriod,
//         const double epochWindowSize, const int numberOfEpochWindows,
//         const double maximumEccentricityChange, const double maximumLongitudeResidualChange,
//         const double maximumInclinationChange );

} // namespace database
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_DATABASE_WRITE_FUNCTIONS_H
