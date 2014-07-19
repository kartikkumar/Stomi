/*    
 *    Copyright (c) 2010-2014, Delft University of Technology
 *    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#ifndef STOMI_ROOT_PATH_H
#define STOMI_ROOT_PATH_H

#include <string>

namespace stomi
{
namespace input_output
{

//! Get root-path for Stomi directory.
/*!
 * Returns root-path corresponding with root-directory of Stomi as a string with
 * trailing slash included.
 * \return Stomi root-path.
 */
static inline std::string getStomiRootPath( )
{
#ifdef STOMI_CUSTOM_ROOT_PATH
    return std::string( STOMI_CUSTOM_ROOT_PATH );
#else
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path in the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return filePath_.substr( 0, filePath_.length( ) -
                                std::string( "InputOutput/rootPath.h" ).length( ) );
#endif
}

} // namespace basics
} // namespace stomi

#endif // STOMI_ROOT_PATH_H
