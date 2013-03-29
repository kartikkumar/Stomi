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
 *      130217    K. Kumar          Updated "mab simulations" references to "stochastic migration".
 *
 *    References
 *
 *    Notes
 *
 */

#include <boost/assign/list_of.hpp>
#include <boost/make_shared.hpp>

#include "StochasticMigration/InputOutput/dictionaries.h"

namespace stochastic_migration
{
namespace input_output
{

using namespace tudat::input_output::dictionary;
using boost::assign::list_of;
using boost::make_shared;

//! Get dictionary for StochasticMigrationDatabaseGenerator application.
tudat::input_output::dictionary::DictionaryPointer
getStochasticMigrationDatabaseGeneratorDictionary( )
{
    DictionaryPointer dictionary = make_shared< Dictionary >( );

    addEntry( dictionary, "DEBUGMODE",                             0, 0, list_of( "DEBUG" ) );
    addEntry( dictionary, "CASE",                                  1, 0 );
    addEntry( dictionary, "DATABASE",                              1, 0, list_of( "DB" ) );
    addEntry( dictionary, "NUMBEROFSIMULATIONS",                   1, 0, list_of( "POPULATION" ) );
    addEntry( dictionary, "RANDOMWALKDURATION",                    0, 0 );
    addEntry( dictionary, "SYNODICPERIODLIMIT",                    0, 0 );
    addEntry( dictionary, "OUTPUTINTERVAL",                        0, 0 );
    addEntry( dictionary, "STARTUPINTEGRATIONDURATION",            0, 0, list_of( "STARTUP" ) );
    addEntry( dictionary, "CONJUNCTIONEVENTDETECTIONDISTANCE",     0, 0,
                 list_of( "CONJUNCTIONDISTANCE ") );
    addEntry( dictionary, "OPPOSITIONEVENTDETECTIONDISTANCE",      0, 0,
                 list_of( "OPPOSITIONDISTANCE ") );
    addEntry( dictionary, "CENTRALBODYGRAVITATIONALPARAMETER",     0, 0, list_of( "GRAVPARAM" ) );
    addEntry( dictionary, "CENTRALBODYJ2GRAVITYCOEFFICIENT",       0, 0, list_of( "J2" ) );
    addEntry( dictionary, "CENTRALBODYEQUATORIALRADIUS",           0, 0, list_of( "RADIUS" ) );
    addEntry( dictionary, "SEMIMAJORAXISLIMIT",                    1, 0, list_of( "SMALIMIT" ) );
    addEntry( dictionary, "ECCENTRICITYMEAN",                      0, 0, list_of( "ECCMEAN" ) );
    addEntry( dictionary, "ECCENTRICITYANGLE",                     0, 0, list_of( "ECCANGLE" ) );
    addEntry( dictionary, "ECCENTRICITYFWHM",                      1, 0, list_of( "ECCFWHM" ) );
    addEntry( dictionary, "INCLINATIONMEAN",                       0, 0, list_of( "INCMEAN" ) );
    addEntry( dictionary, "INCLINATIONANGLE",                      0, 0, list_of( "INCANGLE" ) );
    addEntry( dictionary, "INCLINATIONFWHM",                       1, 0, list_of( "INCFWHM" ) );
    addEntry( dictionary, "PERTURBEDBODYRADIUS",                   0, 0 );
    addEntry( dictionary, "PERTURBEDBODYBULKDENSITY",              0, 0 );
    addEntry( dictionary, "PERTURBEDBODYSEMIMAJORAXISATT0",        0, 0, list_of( "SMA0" ) );
    addEntry( dictionary, "PERTURBEDBODYECCENTRICITYATT0",         0, 0, list_of( "ECC0" ) );
    addEntry( dictionary, "PERTURBEDBODYINCLINATIONATT0",          0, 0, list_of( "INC0" ) );
    addEntry( dictionary, "PERTURBEDBODYARGUMENTOFPERIAPSISATT0",  0, 0, list_of( "AOP0" ) );
    addEntry( dictionary, "PERTURBEDBODYLONGITUDEOFASCENDINGNODEATT0",
           0, 0, list_of( "RAAN0" ) );
    addEntry( dictionary, "PERTURBEDBODYTRUEANOMALYATT0",          0, 0, list_of( "TRAN0" ) );
    addEntry( dictionary, "NUMERICALINTEGRATORTYPE",               0, 0, list_of( "INTEGRATOR" ) );
    addEntry( dictionary, "INITIALSTEPSIZE",                       0, 0, list_of( "STEPSIZE0" ) );
    addEntry( dictionary, "RUNGEKUTTARELATIVEERRORTOLERANCE",      0, 0, list_of( "RELTOL" ) );
    addEntry( dictionary, "RUNGEKUTTAABSOLUTEERRORTOLERANCE",      0, 0, list_of( "ABSTOL" ) );
    addEntry( dictionary, "TESTPARTICLECASETABLENAME",             0, 0 );
    addEntry( dictionary, "TESTPARTICLEINPUTTABLENAME",            0, 0 );
    addEntry( dictionary, "TESTPARTICLEKICKTABLENAME",             0, 0 );
    addEntry( dictionary, "RANDOMWALKMONTECARLORUNTABLENAME",      0, 0 );
    addEntry( dictionary, "RANDOMWALKPERTURBERTABLENAME",          0, 0 );
    addEntry( dictionary, "RANDOMWALKOUTPUTTABLENAME",             0, 0 );

    return dictionary;
}

//! Get dictionary for TestParticleSimulator application.
DictionaryPointer getTestParticleSimulatorDictionary( )
{
    DictionaryPointer dictionary = make_shared< Dictionary >( );

    addEntry( dictionary, "APPLICATIONMODE",                   0, 0, list_of( "MODE" ) );
    addEntry( dictionary, "DEBUGMODE",                         0, 0, list_of( "DEBUG" ) );
    addEntry( dictionary, "DATABASE",                          1, 0, list_of( "DB" ) );
    addEntry( dictionary, "NUMBEROFTHREADS",                   0, 0, list_of( "THREADS" ) );
    addEntry( dictionary, "OUTPUTDIRECTORY",                   0, 0 );
    addEntry( dictionary, "TESTPARTICLESIMULATIONS",           0, 0, list_of( "SIMULATIONS" ) );
    addEntry( dictionary, "TESTPARTICLECASETABLENAME",         0, 0 );
    addEntry( dictionary, "TESTPARTICLEINPUTTABLENAME",        0, 0 );
    addEntry( dictionary, "TESTPARTICLEKICKTABLENAME",         0, 0 );
    addEntry( dictionary, "CONJUNCTIONEVENTDETECTIONDISTANCE", 0, 0,
                 list_of( "CONJUNCTIONDISTANCE ") );
    addEntry( dictionary, "OPPOSITIONEVENTDETECTIONDISTANCE",  0, 0,
                 list_of( "OPPOSITIONDISTANCE ") );
    addEntry( dictionary, "DURATION",                          0, 0 );
    addEntry( dictionary, "OUTPUTINTERVAL",                    0, 0 );
    addEntry( dictionary, "STARTUPINTEGRATIONDURATION",        0, 0, list_of( "STARTUP" ) );
    addEntry( dictionary, "CASE",                              0, 0 );
    addEntry( dictionary, "NUMERICALINTEGRATORTYPE",           0, 0, list_of( "INTEGRATOR" ) );
    addEntry( dictionary, "INITIALSTEPSIZE",                   0, 0, list_of( "STEPSIZE0" )  );
    addEntry( dictionary, "RUNGEKUTTARELATIVEERRORTOLERANCE",  0, 0, list_of( "RELTOL" ) );
    addEntry( dictionary, "RUNGEKUTTAABSOLUTEERRORTOLERANCE",  0, 0, list_of( "ABSTOL" ) );

    return dictionary;
}

//! Get dictionary for RandomWalkSimulator application.
DictionaryPointer getRandomWalkSimulatorDictionary( )
{
    DictionaryPointer dictionary = make_shared< Dictionary >( );

    addEntry( dictionary, "DATABASE",                         1, 0, list_of( "DB" ) );
    addEntry( dictionary, "NUMBEROFTHREADS",                  0, 0, list_of( "THREADS" ) );
    addEntry( dictionary, "OUTPUTDIRECTORY",                  0, 0 );
    addEntry( dictionary, "MONTECARLOPOPULATION",             1, 0, list_of( "MCPOP" ) );
    addEntry( dictionary, "PERTURBERPOPULATION",              0, 0, list_of( "PERTURBERPOP" ) );
    addEntry( dictionary, "MASSDISTRIBUTIONTYPE",             0, 0, list_of( "MASSDIST" ) );
    addEntry( dictionary, "MASSDISTRIBUTIONPARAMETERS",       0, 0, list_of( "MASSDISTPARAM" ) );
    addEntry( dictionary, "OBSERVATIONPERIOD",                0, 0, list_of( "POBS" ) );
    addEntry( dictionary, "EPOCHWINDOWSIZE",                  0, 0  );
    addEntry( dictionary, "NUMBEROFEPOCHWINDOWS",             0, 0, list_of( "NUMEPOCHS" ) );

    return dictionary;
}

} // namespace input_output
} // namespace stochastic_migration
