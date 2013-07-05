/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120808    K. Kumar          File created.
 *      130212    K. Kumar          Migrated code to GeneralTools project and new basics.h file;
 *                                  updated writeKickTableToFile() function; updated documentation.
 *      130217    K. Kumar          Updated "mab simulations" references to "stochastic migration".
 *
 *    References
 *
 *    Notes
 *      The implementation of the writeTestParticleKickPointerTableToFile() function can be improved by
 *      incorporating a string file header for the output file as an argument of the function.
 *
 */

#ifndef STOCHASTIC_MIGRATION_BASIC_INPUT_OUTPUT_H
#define STOCHASTIC_MIGRATION_BASIC_INPUT_OUTPUT_H

#include <limits>
#include <string>

#include <boost/filesystem.hpp>

// #include "StochasticMigration/Database/testParticleKick.h"

namespace stochastic_migration
{
namespace input_output
{

// //! Write test particle kick table to file.
// !
//  * Writes test particle kick table data stored in a Boost associative container to file.
//  * \param testParticleKickTable Table of kicks, extracted from test particle simulations, stored as
//  *          TestParticleKickPointer objects in a Boost associative container set.
//  * \param outputFilename Output filename.
//  * \param outputDirectory Output directory, where output file is saved to (absolute path). Default
//  *          value is an empty string, meaning the output file will be saved somewhere in the build
//  *          directory of the project (specified by the user)
//  * \param delimiter Delimiter string used to delimit column data in file (default=",").
//  * \param outputPrecision Number of digits of precision of floating-point data written to output
//  *          file (default=maximum number of significant digits for a double-precision number).
 
// void writeTestParticleKickTableToFile(
//        const database::TestParticleKickTable& testParticleKickTable,
//        const std::string& outputFilename, const boost::filesystem::path& outputDirectory = "",
//        const std::string& delimiter = ",",
//        const int outputPrecision = std::numeric_limits< double >::digits10 );

} // namespace input_output
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_BASIC_INPUT_OUTPUT_H
