# Copyright (c) 2010-2014, Delft University of Technology
# Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
# All rights reserved.
# See http://bit.ly/12SHPLR for license details.

# Set source files.
set(DATABASE_SOURCES
  "${SOURCE_ROOT}${DATABASE_DIR}/databaseReadFunctions.cpp"
  "${SOURCE_ROOT}${DATABASE_DIR}/databaseWriteFunctions.cpp"
  "${SOURCE_ROOT}${DATABASE_DIR}/randomWalkInput.cpp"
  "${SOURCE_ROOT}${DATABASE_DIR}/randomWalkRun.cpp"  
  "${SOURCE_ROOT}${DATABASE_DIR}/testParticleCase.cpp"
  "${SOURCE_ROOT}${DATABASE_DIR}/testParticleInput.cpp"
  "${SOURCE_ROOT}${DATABASE_DIR}/testParticleKick.cpp"
)

# Set header files.
set(DATABASE_HEADERS
  "${SOURCE_ROOT}${DATABASE_DIR}/databaseReadFunctions.h"
  "${SOURCE_ROOT}${DATABASE_DIR}/databaseWriteFunctions.h"
  "${SOURCE_ROOT}${DATABASE_DIR}/randomWalkInput.h"
  "${SOURCE_ROOT}${DATABASE_DIR}/randomWalkRun.h"  
  "${SOURCE_ROOT}${DATABASE_DIR}/testParticleCase.h"
  "${SOURCE_ROOT}${DATABASE_DIR}/testParticleInput.h"
  "${SOURCE_ROOT}${DATABASE_DIR}/testParticleKick.h"
)

# Set unit test files.
set(DATABASE_UNIT_TESTS
  "${SOURCE_ROOT}${DATABASE_DIR}/UnitTests/unitTestDatabase.cpp"  
  "${SOURCE_ROOT}${DATABASE_DIR}/UnitTests/unitTestDatabaseReadFunctions.cpp"
  "${SOURCE_ROOT}${DATABASE_DIR}/UnitTests/unitTestDatabaseWriteFunctions.cpp"  
  "${SOURCE_ROOT}${DATABASE_DIR}/UnitTests/unitTestRandomWalkInput.cpp"  
  "${SOURCE_ROOT}${DATABASE_DIR}/UnitTests/unitTestRandomWalkRun.cpp"        
  "${SOURCE_ROOT}${DATABASE_DIR}/UnitTests/unitTestTestParticleCase.cpp"    
  "${SOURCE_ROOT}${DATABASE_DIR}/UnitTests/unitTestTestParticleInput.cpp"    
  "${SOURCE_ROOT}${DATABASE_DIR}/UnitTests/unitTestTestParticleKick.cpp"    
)

# Add static library.
add_library(stomi_database STATIC ${DATABASE_SOURCES} ${DATABASE_HEADERS})
setup_library_target(stomi_database)

# Add unit tests.
add_executable(stomi_test_Database ${DATABASE_UNIT_TESTS})
setup_unit_test_target(stomi_test_Database)
target_link_libraries(stomi_test_Database 
                      stomi_database stomi_input_output
                      tudat_input_output
                      ${TUDAT_CORE_LIBRARIES} 
                      ${SQLITECPP_LIBRARY} 
                      ${SQLITE3_LIBRARY}  
                      ${Boost_LIBRARIES})
