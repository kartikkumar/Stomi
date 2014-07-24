/*    
 * Copyright (c) 2010-2014, Delft University of Technology
 * Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 * See http://bit.ly/12SHPLR for license details.
 */

#include <boost/assign/list_of.hpp>
#include <boost/make_shared.hpp>

#include "Stomi/InputOutput/dictionaries.h"

namespace stomi
{
namespace input_output
{

using namespace tudat::input_output::dictionary;
using boost::assign::list_of;
using boost::make_shared;

//! Get general parameters dictionary.
DictionaryPointer getGeneralParametersDictionary( )
{
    DictionaryPointer dictionary = make_shared< Dictionary >( );

    // Add required parameters.
    addEntry( dictionary, "MODE",                                  1, 0 );  
    addEntry( dictionary, "DATABASEPATH",                          1, 0, list_of( "DB" ) ); 

    return dictionary;       
}

//! Get dictionary for test particle database generator application mode.
DictionaryPointer getTestParticleDatabaseGeneratorDictionary( )
{
    DictionaryPointer dictionary = make_shared< Dictionary >( );

    // Add required parameters.
    addEntry( dictionary, "TESTPARTICLECASE",                      1, 0, list_of( "TPCASE" ) );    
    addEntry( dictionary, "MONTECARLOPOPULATION",                  1, 0, list_of( "MCPOP" ) );
    addEntry( dictionary, "RANDOMWALKSIMULATIONPERIOD",            1, 0, list_of( "TRANDOM" ) );
    addEntry( dictionary, "CENTRALBODYGRAVITATIONALPARAMETER",     1, 0, list_of( "GRAVPARAM" ) );
    addEntry( dictionary, "PERTURBEDBODYRADIUS",                   1, 0, list_of( "RPERTURBED") );
    addEntry( dictionary, "PERTURBEDBODYBULKDENSITY",              1, 0, 
                list_of( "RHOPERTURBED" ) );
    addEntry( dictionary, "PERTURBEDBODYSEMIMAJORAXISATT0",        1, 0, list_of( "SMA0" ) );
    addEntry( dictionary, "PERTURBEDBODYECCENTRICITYATT0",         1, 0, list_of( "ECC0" ) );
    addEntry( dictionary, "PERTURBEDBODYINCLINATIONATT0",          1, 0, list_of( "INC0" ) );
    addEntry( dictionary, "PERTURBEDBODYARGUMENTOFPERIAPSISATT0",  1, 0, list_of( "AOP0" ) );
    addEntry( dictionary, "PERTURBEDBODYLONGITUDEOFASCENDINGNODEATT0",
                0, 0, list_of( "LAN0" ) );
    addEntry( dictionary, "PERTURBEDBODYTRUEANOMALYATT0",          1, 0, list_of( "TRAN0" ) );
    addEntry( dictionary, "SEMIMAJORAXISDISTRIBUTIONLIMIT",        1, 0, list_of( "SMALIMIT" ) );

    // Add optional parameters.
    addEntry( dictionary, "SYNODICPERIODMAXIMUM",                  0, 0, 
        list_of( "TSYNODICMAX" ) );
    addEntry( dictionary, "STARTUPINTEGRATIONPERIOD",              0, 0, list_of( "TSTARTUP" ) );
    addEntry( dictionary, "CENTRALBODYJ2GRAVITYCOEFFICIENT",       0, 0, list_of( "J2" ) );
    addEntry( dictionary, "CENTRALBODYEQUATORIALRADIUS",           0, 0, list_of( "RCENTRAL" ) );
    addEntry( dictionary, "CONJUNCTIONEVENTDETECTIONDISTANCE",     0, 0, 
                list_of( "DCONJUNCTION ") );
    addEntry( dictionary, "OPPOSITIONEVENTDETECTIONDISTANCE",      0, 0,
                list_of( "DOPPOSITION ") );
    addEntry( dictionary, "ECCENTRICITYDISTRIBUTIONMEAN",          0, 0, list_of( "ECCMEAN" ) );
    addEntry( dictionary, "ECCENTRICITYDISTRIBUTIONSTANDARDDEVIATION", 
                0, 0, list_of( "ECCSIGMA" ) );
    addEntry( dictionary, "INCLINATIONDISTRIBUTIONMEAN",           0, 0, list_of( "INCMEAN" ) );
    addEntry( dictionary, "INCLINATIONDISTRIBUTIONSTANDARDDEVIATION",                       
                0, 0, list_of( "INCSIGMA" ) );
    addEntry( dictionary, "NUMERICALINTEGRATORTYPE",               0, 0, list_of( "INTEGRATOR" ) );
    addEntry( dictionary, "INITIALSTEPSIZE",                       0, 0, list_of( "STEPSIZE" ) );
    addEntry( dictionary, "RUNGEKUTTARELATIVEERRORTOLERANCE",      0, 0, list_of( "RELTOL" ) );
    addEntry( dictionary, "RUNGEKUTTAABSOLUTEERRORTOLERANCE",      0, 0, list_of( "ABSTOL" ) );
    addEntry( dictionary, "TESTPARTICLECASETABLENAME",             
                0, 0, list_of( "TPCASETABLE" ) );
    addEntry( dictionary, "TESTPARTICLEINPUTTABLENAME",            
                0, 0, list_of( "TPINPUTTABLE" ) );
    addEntry( dictionary, "TESTPARTICLEKICKTABLENAME",             
                0, 0, list_of( "TPKICKTABLE" ) );
    addEntry( dictionary, "RANDOMWALKRUNTABLENAME",                0, 0, 
                list_of( "RWRUNTABLE" ) );
    addEntry( dictionary, "RANDOMWALKINPUTTABLENAME",              0, 0, 
                list_of( "RWINPUTTABLE" ) );
    addEntry( dictionary, "RANDOMWALKPERTURBERTABLENAME",          0, 0, 
                list_of( "RWPERTURBERTABLE" ) );
    addEntry( dictionary, "RANDOMWALKOUTPUTTABLENAME",             0, 0, 
                list_of( "RWOUTPUTTABLE" ) );    

    return dictionary;
}

//! Get dictionary for random walk database generator application mode.
tudat::input_output::dictionary::DictionaryPointer getRandomWalkDatabaseGeneratorDictionary( )
{
    DictionaryPointer dictionary = make_shared< Dictionary >( );

    // Add required parameters.
    addEntry( dictionary, "RANDOMWALKRUN",                         1, 0, list_of( "RWRUN" ) );
    addEntry( dictionary, "TESTPARTICLECASE",                      1, 0, list_of( "TPCASE" ) );
    addEntry( dictionary, "MONTECARLOPOPULATION",                  1, 0, list_of( "MCPOP" ) );
    addEntry( dictionary, "PERTURBERRINGNUMBERDENSITY",            1, 0, 
                list_of( "NRHOPERTURBERRING" ) );    
    addEntry( dictionary, "PERTURBERRINGMASS",                     1, 0, 
                list_of( "MPERTURBERRING" ) );
    addEntry( dictionary, "OBSERVATIONPERIOD",                     1, 0, list_of( "POBS" ) );
    addEntry( dictionary, "EPOCHWINDOWSIZE",                       1, 0, list_of( "EPOCHSIZE" ) );
    addEntry( dictionary, "NUMBEROFEPOCHWINDOWS",                  1, 0, list_of( "NUMEPOCHS" ) ); 

    // Add optional parameters.
    addEntry( dictionary, "TESTPARTICLECASETABLENAME",             
                0, 0, list_of( "TPCASETABLE" ) );
    addEntry( dictionary, "TESTPARTICLEINPUTTABLENAME",            
                0, 0, list_of( "TPINPUTTABLE" ) );
    addEntry( dictionary, "TESTPARTICLEKICKTABLENAME",             
                0, 0, list_of( "TPKICKTABLE" ) );
    addEntry( dictionary, "RANDOMWALKRUNTABLENAME",                0, 0, 
                list_of( "RWRUNTABLE" ) );
    addEntry( dictionary, "RANDOMWALKINPUTTABLENAME",              0, 0, 
                list_of( "RWINPUTTABLE" ) );
    addEntry( dictionary, "RANDOMWALKPERTURBERTABLENAME",          0, 0, 
                list_of( "RWPERTURBERTABLE" ) );
    addEntry( dictionary, "RANDOMWALKOUTPUTTABLENAME",             0, 0, 
                list_of( "RWOUTPUTTABLE" ) );              

    return dictionary;
}

//! Get dictionary for test particle simulator application mode.
DictionaryPointer getTestParticleSimulatorDictionary( )
{
    // Retrieve dictionary for test particle simulator.
    DictionaryPointer dictionary = make_shared< Dictionary >( );

    // Add required parameters.
    addEntry( dictionary, "TESTPARTICLECASE",                      1, 0, list_of( "TPCASE" ) );    

    // Add optional parameters.
    addEntry( dictionary, "NUMBEROFTHREADS",                       0, 0, list_of( "THREADS" ) );
    addEntry( dictionary, "OUTPUTMODE",                            0, 0, list_of( "OUTPUT" ) );
    addEntry( dictionary, "FILEOUTPUTDIRECTORY",                   0, 0, 
                list_of( "FILEOUTPUTDIR" ) );
    addEntry( dictionary, "OUTPUTINTERVAL",                        0, 0, list_of( "TOUTPUT" ) );  
    addEntry( dictionary, "TESTPARTICLESIMULATIONS",               0, 0, list_of( "TPSIMS" ) );
    addEntry( dictionary, "RANDOMWALKSIMULATIONPERIOD",            0, 0, list_of( "TRANDOM" ) );
    addEntry( dictionary, "CENTRALBODYGRAVITATIONALPARAMETER",     0, 0, list_of( "GRAVPARAM" ) );
    addEntry( dictionary, "PERTURBEDBODYRADIUS",                   0, 0, list_of( "RPERTURBED") );
    addEntry( dictionary, "PERTURBEDBODYBULKDENSITY",              0, 0, 
                list_of( "RHOPERTURBED" ) );
    addEntry( dictionary, "PERTURBEDBODYSEMIMAJORAXISATT0",        0, 0, list_of( "SMA0" ) );
    addEntry( dictionary, "PERTURBEDBODYECCENTRICITYATT0",         0, 0, list_of( "ECC0" ) );
    addEntry( dictionary, "PERTURBEDBODYINCLINATIONATT0",          0, 0, list_of( "INC0" ) );
    addEntry( dictionary, "PERTURBEDBODYARGUMENTOFPERIAPSISATT0",  0, 0, list_of( "AOP0" ) );
    addEntry( dictionary, "PERTURBEDBODYLONGITUDEOFASCENDINGNODEATT0",
                0, 0, list_of( "LAN0" ) );
    addEntry( dictionary, "PERTURBEDBODYTRUEANOMALYATT0",          0, 0, list_of( "TRAN0" ) ); 
    addEntry( dictionary, "SYNODICPERIODMAXIMUM",                  0, 0, 
                list_of( "TSYNODICMAX" ) );
    addEntry( dictionary, "STARTUPINTEGRATIONPERIOD",              0, 0, list_of( "TSTARTUP" ) );
    addEntry( dictionary, "CENTRALBODYJ2GRAVITYCOEFFICIENT",       0, 0, list_of( "J2" ) );
    addEntry( dictionary, "CENTRALBODYEQUATORIALRADIUS",           0, 0, list_of( "RCENTRAL" ) );
    addEntry( dictionary, "CONJUNCTIONEVENTDETECTIONDISTANCE",     0, 0, 
                list_of( "DCONJUNCTION") );
    addEntry( dictionary, "OPPOSITIONEVENTDETECTIONDISTANCE",      0, 0,
                list_of( "DOPPOSITION") );
    addEntry( dictionary, "NUMERICALINTEGRATORTYPE",               0, 0, list_of( "INTEGRATOR" ) );
    addEntry( dictionary, "INITIALSTEPSIZE",                       0, 0, list_of( "STEPSIZE" ) );
    addEntry( dictionary, "RUNGEKUTTARELATIVEERRORTOLERANCE",      0, 0, list_of( "RELTOL" ) );
    addEntry( dictionary, "RUNGEKUTTAABSOLUTEERRORTOLERANCE",      0, 0, list_of( "ABSTOL" ) );
    addEntry( dictionary, "TESTPARTICLECASETABLENAME",             
                0, 0, list_of( "TPCASETABLE" ) );
    addEntry( dictionary, "TESTPARTICLEINPUTTABLENAME",            
                0, 0, list_of( "TPINPUTTABLE" ) );
    addEntry( dictionary, "TESTPARTICLEKICKTABLENAME",             
                0, 0, list_of( "TPKICKTABLE" ) );
    addEntry( dictionary, "RANDOMWALKRUNTABLENAME",                0, 0, 
                list_of( "RWRUNTABLE" ) );
    addEntry( dictionary, "RANDOMWALKINPUTTABLENAME",              0, 0, 
                list_of( "RWINPUTTABLE" ) );
    addEntry( dictionary, "RANDOMWALKPERTURBERTABLENAME",          0, 0, 
                list_of( "RWPERTURBERTABLE" ) );
    addEntry( dictionary, "RANDOMWALKOUTPUTTABLENAME",             0, 0, 
                list_of( "RWOUTPUTTABLE" ) ); 

    return dictionary;
}

//! Get dictionary for random walk simulator application mode.
DictionaryPointer getRandomWalkSimulatorDictionary( )
{
    // Retrieve dictionary for random walk simulator.
    DictionaryPointer dictionary = make_shared< Dictionary >( );

    // Add required parameters.
    addEntry( dictionary, "RANDOMWALKRUN",                         1, 0, list_of( "RWRUN" ) );

    // Add optional parameters.
    addEntry( dictionary, "NUMBEROFTHREADS",                       0, 0, list_of( "THREADS" ) );
    addEntry( dictionary, "OUTPUTMODE",                            0, 0, list_of( "OUTPUT" ) );   
    addEntry( dictionary, "FILEOUTPUTDIRECTORY",                   0, 0, 
                list_of( "FILEOUTPUTDIR" ) ); 
    addEntry( dictionary, "RANDOMWALKSIMULATIONS",                 0, 0, list_of( "RWSIMS" ) );   
    addEntry( dictionary, "PERTURBERRINGNUMBERDENSITY",            1, 0, 
                list_of( "NRHOPERTURBERRING" ) );    
    addEntry( dictionary, "PERTURBERRINGMASS",                     1, 0, 
                list_of( "MPERTURBERRING" ) );
    addEntry( dictionary, "OBSERVATIONPERIOD",                     0, 0, list_of( "POBS" ) );
    addEntry( dictionary, "EPOCHWINDOWSIZE",                       0, 0, list_of( "EPOCHSIZE" ) );
    addEntry( dictionary, "NUMBEROFEPOCHWINDOWS",                  0, 0, list_of( "NUMEPOCHS" ) );
    addEntry( dictionary, "TESTPARTICLECASETABLENAME",             
                0, 0, list_of( "TPCASETABLE" ) );
    addEntry( dictionary, "TESTPARTICLEINPUTTABLENAME",            
                0, 0, list_of( "TPINPUTTABLE" ) );
    addEntry( dictionary, "TESTPARTICLEKICKTABLENAME",             
                0, 0, list_of( "TPKICKTABLE" ) );
    addEntry( dictionary, "RANDOMWALKRUNTABLENAME",                0, 0, 
                list_of( "RWRUNTABLE" ) );
    addEntry( dictionary, "RANDOMWALKINPUTTABLENAME",              0, 0, 
                list_of( "RWINPUTTABLE" ) );
    addEntry( dictionary, "RANDOMWALKPERTURBERTABLENAME",          0, 0, 
                list_of( "RWPERTURBERTABLE" ) );
    addEntry( dictionary, "RANDOMWALKOUTPUTTABLENAME",             0, 0, 
                list_of( "RWOUTPUTTABLE" ) );     

    return dictionary;
}

} // namespace input_output
} // namespace stomi
