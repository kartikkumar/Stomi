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

#include "StochasticMigration/Database/testParticleKick.h"

namespace stochastic_migration
{
namespace input_output
{

//! Write test particle kick table to file.
/*!
 * Writes test particle kick table data stored in a STL set to file.
 * \param testParticleKickTable Table of kicks, extracted from test particle simulations, stored as
 *          TestParticleKickPointer objects in a STL set.
 * \param outputFilename Output filename.
 * \param outputDirectory Output directory, where output file is saved to (absolute path). Default
 *          value is an empty string, meaning the output file will be saved somewhere in the build
 *          directory of the project (specified by the user)
 * \param delimiter Delimiter string used to delimit column data in file (default=",").
 * \param outputPrecision Number of digits of precision of floating-point data written to output
 *          file (default=maximum number of significant digits for a double-precision number).
 */
void writeTestParticleKickPointerTableToFile(
       const database::TestParticleKickPointerTable& testParticleKickTable,
       const std::string& outputFilename, const boost::filesystem::path& outputDirectory = "",
       const std::string& delimiter = ",",
       const int outputPrecision = std::numeric_limits< double >::digits10 );

} // namespace input_output
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_BASIC_INPUT_OUTPUT_H
