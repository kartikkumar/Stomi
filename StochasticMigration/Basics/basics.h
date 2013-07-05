/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130212    K. Kumar          File created from code in basicInputOutput.h.
 *      130217    K. Kumar          Updated "mab simulations" references to "stochastic migration".
 *
 *    References
 *
 *    Notes
 *      A unit test is needed for the getStochasticMigrationRootPath() function.
 *
 */

#ifndef STOCHASTIC_MIGRATION_BASICS_H
#define STOCHASTIC_MIGRATION_BASICS_H

#include <string>

namespace stochastic_migration
{
namespace basics
{

//! Get root-path for StochasticMigration directory.
/*!
 * Returns root-path corresponding with root-directory of StochasticMigration as a string with
 * trailing slash included.
 * \return StochasticMigration root-path.
 */
static inline std::string getStochasticMigrationRootPath( )
{
#ifdef STOCHASTIC_MIGRATION_CUSTOM_ROOT_PATH
    return std::string( STOCHASTIC_MIGRATION_CUSTOM_ROOT_PATH );
#else
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path in the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return filePath_.substr( 0, filePath_.length( ) -
                                std::string( "Basics/basics.h" ).length( ) );
#endif
}

} // namespace basics
} // namespace stochastic_migration

#endif // STOCHASTIC_MIGRATION_BASICS_H
