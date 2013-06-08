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
 * \param testParticleSimulationNumber Test particle simulation number that kick data corresponds
 *          to.
 * \param testParticleKickTable Table in SQLite database containing kick data for specified test
 *          particle simulation number.
 * \param perturbedBodyEnergyError Error in perturbed body's energy due to numerical integration
 *          (as a percentage).
 * \param perturbedBodyAngularMomentumError Error in perturbed body's angular momentum due to
 *          numerical integration (as a percentage).
 * \param testParticleKickTableName String name of test particle kick table in database (default
 *          is set to "test_particle_kicks").
 * \param testParticleInputTableName String name of test particle input table in database (default
 *          is set to "test_particle_input").
 */
void populateTestParticleKickTable(
        const std::string& databaseAbsolutePath, const int testParticleSimulationNumber,
        const TestParticleKickTable& testParticleKickTable, const double perturbedBodyEnergyError,
        const double perturbedBodyAngularMomentumError,
        const std::string& testParticleKickTableName = "test_particle_kicks",
        const std::string& testParticleInputTableName = "test_particle_input" );

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
