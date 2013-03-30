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
 *      120515    K. Kumar          File created.
 *      120520    K. Kumar          Added function to parse random walk input file.
 *      120523    K. Kumar          Added function to parse random walk verification input file.
 *      130212    K. Kumar          Migrated code to GeneralTools project and new basics.h file;
 *                                  updated writeKickTableToFile() function; updated documentation;
 *                                  modified algorithm to use updated vector-definition of kick
 *                                  table.
 *      130217    K. Kumar          Updated "mab simulations" references to "stochastic migration".
 *      130218    K. Kumar          Updated "encounter" to "conjunction".
 *      130329    K. Kumar          Updated looping through kick table to work with ptr_set.
 *
 *    References
 *
 *    Notes
 *      The implementation of the writeTestParticleKickPointerTableToFile() function can be improved by
 *      incorporating a string file header for the output file as an argument of the function.
 *
 *      The test particle kick table includes mass factor as a member variable of TestParticleKickPointer,
 *      but this is not written to the output file.
 *
 */

#include <fstream>
#include <iostream>
#include <string>

#include "StochasticMigration/InputOutput/basicInputOutput.h"

namespace stochastic_migration
{
namespace input_output
{

//! Write test particle kick table to file.
void writeTestParticleKickTableToFile(
       const database::TestParticleKickTable& testParticleKickTable,
       const std::string& outputFilename, const boost::filesystem::path& outputDirectory,
       const std::string& delimiter, const int outputPrecision )
{
    // Check if output directory exists.
    if ( !boost::filesystem::exists( outputDirectory ) )
    {
       std::cerr << "Directory does not exist. Will be created." << std::endl;
       boost::filesystem::create_directories( outputDirectory );
    }

    // Open output file.
    std::ofstream outputFile;
    std::string outputDirectoryAndFilename = outputDirectory.string( ) + "/" + outputFilename;
    outputFile.open( outputDirectoryAndFilename.c_str( ) );

    // Set precision in output file.
    outputFile.precision( outputPrecision );

    // Loop over kick table.
    for ( database::TestParticleKickTable::iterator iteratorKickTable 
          = testParticleKickTable.begin( );
          iteratorKickTable != testParticleKickTable.end( ); iteratorKickTable++ )
    {        
       // Print map key and value to output file.
       outputFile << iteratorKickTable->conjunctionEpoch << delimiter
                  << iteratorKickTable->conjunctionDistance << delimiter
                  << iteratorKickTable->conjunctionDuration << delimiter
                  << iteratorKickTable->preConjunctionEpoch << delimiter
                  << iteratorKickTable->preConjunctionDistance << delimiter
                  << iteratorKickTable->preConjunctionSemiMajorAxis << delimiter
                  << iteratorKickTable->preConjunctionEccentricity << delimiter
                  << iteratorKickTable->preConjunctionInclination << delimiter
                  << iteratorKickTable->postConjunctionEpoch << delimiter
                  << iteratorKickTable->postConjunctionDistance << delimiter
                  << iteratorKickTable->postConjunctionSemiMajorAxis << delimiter
                  << iteratorKickTable->postConjunctionEccentricity << delimiter
                  << iteratorKickTable->postConjunctionInclination << std::endl;
    }

   // Close output file.
   outputFile.close( );
}

} // namespace input_output
} // namespace stochastic_migration
