/*    
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 *    All rights reserved.
 *    See http://bit.ly/12SHPLR for license details.
 */

#include <boost/assign/list_of.hpp>
#include <boost/make_shared.hpp>

#include "StoMi/InputOutput/dictionaries.h"

namespace stochastic_migration
{
namespace input_output
{

using namespace tudat::input_output::dictionary;
using boost::assign::list_of;
using boost::make_shared;

//! Get dictionary for test particle database generator application.
DictionaryPointer getTestParticleDatabaseGeneratorDictionary( )
{
    DictionaryPointer dictionary = make_shared< Dictionary >( );

    // Add required parameters.
    addEntry( dictionary, "DATABASEPATH",                          1, 0, 
                list_of( "DATABASE" )( "DB" ) );    
    addEntry( dictionary, "CASE",                                  1, 0 );
    addEntry( dictionary, "NUMBEROFSIMULATIONS",                   1, 0, list_of( "POPULATION" ) );
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
    addEntry( dictionary, "ECCENTRICITYDISTRIBUTIONANGLE",         0, 0, list_of( "ECCANGLE" ) );
    addEntry( dictionary, "ECCENTRICITYDISTRIBUTIONFULLWIDTHHALFMAXIMUM", 
                0, 0, list_of( "ECCFWHM" ) );
    addEntry( dictionary, "INCLINATIONDISTRIBUTIONMEAN",           0, 0, list_of( "INCMEAN" ) );
    addEntry( dictionary, "INCLINATIONDISTRIBUTIONANGLE",          0, 0, list_of( "INCANGLE" ) );
    addEntry( dictionary, "INCLINATIONDISTRIBUTIONFULLWIDTHHALFMAXIMUM",                       
                0, 0, list_of( "INCFWHM" ) );
    addEntry( dictionary, "NUMERICALINTEGRATORTYPE",               0, 0, list_of( "INTEGRATOR" ) );
    addEntry( dictionary, "INITIALSTEPSIZE",                       0, 0, list_of( "STEPSIZE0" ) );
    addEntry( dictionary, "RUNGEKUTTARELATIVEERRORTOLERANCE",      0, 0, list_of( "RELTOL" ) );
    addEntry( dictionary, "RUNGEKUTTAABSOLUTEERRORTOLERANCE",      0, 0, list_of( "ABSTOL" ) );
    addEntry( dictionary, "TESTPARTICLECASETABLENAME",             
                0, 0, list_of( "TPCASETABLE" ) );
    addEntry( dictionary, "TESTPARTICLEINPUTTABLENAME",            
                0, 0, list_of( "TPINPUTTABLE" ) );
    addEntry( dictionary, "TESTPARTICLEKICKTABLENAME",             
                0, 0, list_of( "TPKICKTABLE" ) );

    return dictionary;
}

//! Get dictionary for random walk database generator application.
tudat::input_output::dictionary::DictionaryPointer getRandomWalkDatabaseGeneratorDictionary( )
{
    DictionaryPointer dictionary = make_shared< Dictionary >( );

    // Add required parameters.
    addEntry( dictionary, "DATABASEPATH",                          1, 0, 
                list_of( "DATABASE" )( "DB" ) );    
    addEntry( dictionary, "CASE",                                  1, 0 );
    addEntry( dictionary, "TESTPARTICLECASE",                      1, 0, list_of( "TPCASE" ) );
    addEntry( dictionary, "MONTECARLOPOPULATION",                  1, 0, list_of( "MCPOP" ) );
    addEntry( dictionary, "PERTURBERDENSITY",                      1, 0, 
                list_of( "RHOPERTURBER" ) );    
    addEntry( dictionary, "PERTURBERRINGMASS",                     1, 0, 
                list_of( "MPERTURBERRING" ) );
    addEntry( dictionary, "OBSERVATIONPERIOD",                     1, 0, list_of( "POBS" ) );
    addEntry( dictionary, "EPOCHWINDOWSIZE",                       1, 0, list_of( "EPOCHSIZE" ) );
    addEntry( dictionary, "NUMBEROFEPOCHWINDOWS",                  1, 0, list_of( "NUMEPOCHS" ) ); 

    // Add optional parameters.
    addEntry( dictionary, "RANDOMWALKCASETABLENAME",               0, 0, 
                list_of( "RWCASETABLE" ) );
    addEntry( dictionary, "RANDOMWALKINPUTTABLENAME",              0, 0, 
                list_of( "RWINPUTTABLE" ) );
    addEntry( dictionary, "RANDOMWALKWINDOWSTABLENAME",            0, 0, 
                list_of( "RWWINDOWSTABLE" ) );
    addEntry( dictionary, "RANDOMWALKPERTURBERTABLENAME",          0, 0, 
                list_of( "RWPERTURBTABLE" ) );
    addEntry( dictionary, "RANDOMWALKOUTPUTTABLENAME",             0, 0, 
                list_of( "RWOUTPUTTABLE" ) );
    addEntry( dictionary, "TESTPARTICLECASETABLENAME",             0, 0, 
                list_of( "TPCASETABLE" ) );  
    addEntry( dictionary, "TESTPARTICLEINPUTTABLENAME",            0, 0, 
                list_of( "TPINPUTTABLE" ) );                  

    return dictionary;
}

//! Get dictionary for testParticleSimulator application.
DictionaryPointer getTestParticleSimulatorDictionary( )
{
    // Retrieve dictionary for test particle simulator.
    DictionaryPointer dictionary = make_shared< Dictionary >( );

    // Add required parameters.
    addEntry( dictionary, "CASE",                                  1, 0 );
    addEntry( dictionary, "DATABASE",                              1, 0, list_of( "DB" ) );

    // Add optional parameters.
    addEntry( dictionary, "NUMBEROFTHREADS",                       0, 0, list_of( "THREADS" ) );
    addEntry( dictionary, "OUTPUTMODE",                            0, 0, list_of( "OUTPUT" ) );
    addEntry( dictionary, "FILEOUTPUTDIRECTORY",                   0, 0, 
                list_of( "FILEOUTPUTDIR" ) );
    addEntry( dictionary, "OUTPUTINTERVAL",                        0, 0, list_of( "TOUTPUT" ) );  
    addEntry( dictionary, "SIMULATIONSTOEXECUTE",                  0, 0, 
                list_of( "SIMULATIONS" ) );
    addEntry( dictionary, "TESTPARTICLECASETABLENAME",             0, 0 );
    addEntry( dictionary, "TESTPARTICLEINPUTTABLENAME",            0, 0 );
    addEntry( dictionary, "TESTPARTICLEKICKTABLENAME",             0, 0 );
    addEntry( dictionary, "RANDOMWALKSIMULATIONPERIOD",            0, 0 ); 
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
    addEntry( dictionary, "INITIALSTEPSIZE",                       0, 0, list_of( "STEPSIZE0" ) );
    addEntry( dictionary, "RUNGEKUTTARELATIVEERRORTOLERANCE",      0, 0, list_of( "RELTOL" ) );
    addEntry( dictionary, "RUNGEKUTTAABSOLUTEERRORTOLERANCE",      0, 0, list_of( "ABSTOL" ) );

    return dictionary;
}

//! Get dictionary for randomWalkSimulator application.
DictionaryPointer getRandomWalkSimulatorDictionary( )
{
    // Retrieve dictionary for random walk simulator.
    DictionaryPointer dictionary = make_shared< Dictionary >( );

    // Add required parameters.
    addEntry( dictionary, "CASE",                                  1, 0 );
    addEntry( dictionary, "DATABASE",                              1, 0, list_of( "DB" ) );

    // Add optional parameters.
    addEntry( dictionary, "NUMBEROFTHREADS",                       0, 0, list_of( "THREADS" ) );
    addEntry( dictionary, "OUTPUTMODE",                            0, 0, list_of( "OUTPUT" ) );   
    addEntry( dictionary, "FILEOUTPUTDIRECTORY",                   0, 0, 
                list_of( "FILEOUTPUTDIR" ) ); 
    addEntry( dictionary, "MONTECARLORUNSTOEXECUTE",               0, 0, 
                list_of( "MONTECARLORUNS" ) );    
    addEntry( dictionary, "RANDOMWALKCASETABLENAME",               0, 0, 
                list_of( "RWCASETABLE" ) );
    addEntry( dictionary, "RANDOMWALKINPUTTABLENAME",              0, 0, 
                list_of( "RWINPUTTABLE" ) );
    addEntry( dictionary, "RANDOMWALKPERTURBERTABLENAME",          0, 0, 
                list_of( "RWPERTURBERTABLE" ) );    
    addEntry( dictionary, "RANDOMWALKOUTPUTTABLENAME",             0, 0, 
                list_of( "RWOUTPUTTABLE" ) );        
    addEntry( dictionary, "TESTPARTICLECASETABLENAME",             0, 0, list_of( "TPCASE" ) );
    addEntry( dictionary, "TESTPARTICLEKICKTABLENAME",             0, 0, list_of( "TPKICK" ) );
    addEntry( dictionary, "MONTECARLOPOPULATION",                  0, 0, list_of( "MCPOP" ) );
    addEntry( dictionary, "PERTURBERDENSITY",                      0, 0, 
                list_of( "RHOPERTURBER" ) );    
    addEntry( dictionary, "PERTURBERRINGMASS",                     0, 0, 
                list_of( "MPERTURBERRING" ) );        
    addEntry( dictionary, "OBSERVATIONPERIOD",                     0, 0, list_of( "POBS" ) );
    addEntry( dictionary, "EPOCHWINDOWSIZE",                       0, 0, list_of( "EPOCHSIZE" ) );
    addEntry( dictionary, "NUMBEROFEPOCHWINDOWS",                  0, 0, list_of( "NUMEPOCHS" ) );    


    return dictionary;
}

} // namespace input_output
} // namespace stochastic_migration