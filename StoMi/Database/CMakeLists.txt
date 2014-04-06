 #    Copyright (c) 2010-2014, Delft University of Technology
 #    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 #    All rights reserved.
 #    See http://bit.ly/12SHPLR for license details.

# Set source files.
set(DATABASE_SOURCES
  "${SRCROOT}${DATABASEDIR}/databaseReadFunctions.cpp"
  "${SRCROOT}${DATABASEDIR}/databaseWriteFunctions.cpp"
  "${SRCROOT}${DATABASEDIR}/randomWalkCase.cpp"
  "${SRCROOT}${DATABASEDIR}/randomWalkInput.cpp"
  "${SRCROOT}${DATABASEDIR}/testParticleCase.cpp"
  "${SRCROOT}${DATABASEDIR}/testParticleInput.cpp"
  "${SRCROOT}${DATABASEDIR}/testParticleKick.cpp"
)

# Set header files.
set(DATABASE_HEADERS
  "${SRCROOT}${DATABASEDIR}/databaseReadFunctions.h"
  "${SRCROOT}${DATABASEDIR}/databaseWriteFunctions.h"
  "${SRCROOT}${DATABASEDIR}/randomWalkCase.h"
  "${SRCROOT}${DATABASEDIR}/randomWalkInput.h"
  "${SRCROOT}${DATABASEDIR}/testParticleCase.h"
  "${SRCROOT}${DATABASEDIR}/testParticleInput.h"
  "${SRCROOT}${DATABASEDIR}/testParticleKick.h"
)

# Set unit test files.
set(DATABASE_UNIT_TESTS
  "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestDatabase.cpp"  
  "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestDatabaseReadFunctions.cpp"
  "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestDatabaseWriteFunctions.cpp"  
  "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestRandomWalkCase.cpp"    
  "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestRandomWalkInput.cpp"    
  "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestTestParticleCase.cpp"    
  "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestTestParticleInput.cpp"    
  "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestTestParticleKick.cpp"    
)

# Add static library.
add_library(stomi_database STATIC ${DATABASE_SOURCES} ${DATABASE_HEADERS})
setup_library_target(stomi_database)

# Add unit tests.
add_executable(test_Database ${DATABASE_UNIT_TESTS})
setup_unit_test_target(test_Database)
target_link_libraries(test_Database 
                      stomi_database stomi_input_output
                      tudat_input_output
                      ${TUDAT_CORE_LIBRARIES} 
                      ${SQLITECPP_LIBRARY} 
                      ${SQLITE3_LIBRARY}  
                      ${Boost_LIBRARIES})
