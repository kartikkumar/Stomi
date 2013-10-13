 #    Copyright (c) 2010-2013, Delft University of Technology
 #    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 #    All rights reserved.
 #    See http://bit.ly/12SHPLR for license details.

# Set source files.
set(DATABASE_SOURCES
  "${SRCROOT}${DATABASEDIR}/databaseReadFunctions.cpp"
  "${SRCROOT}${DATABASEDIR}/databaseWriteFunctions.cpp"
#  "${SRCROOT}${DATABASEDIR}/randomWalkMonteCarloRun.cpp"
#  "${SRCROOT}${DATABASEDIR}/randomWalkPerturber.cpp"
  "${SRCROOT}${DATABASEDIR}/randomWalkCase.cpp"
  "${SRCROOT}${DATABASEDIR}/testParticleCase.cpp"
  "${SRCROOT}${DATABASEDIR}/testParticleInput.cpp"
  "${SRCROOT}${DATABASEDIR}/testParticleKick.cpp"
)

# Set header files.
set(DATABASE_HEADERS
  "${SRCROOT}${DATABASEDIR}/databaseReadFunctions.h"
  "${SRCROOT}${DATABASEDIR}/databaseWriteFunctions.h"
#  "${SRCROOT}${DATABASEDIR}/randomWalkMonteCarloRun.h"
#  "${SRCROOT}${DATABASEDIR}/randomWalkPerturber.h"
  "${SRCROOT}${DATABASEDIR}/randomWalkCase.h"
  "${SRCROOT}${DATABASEDIR}/testParticleCase.h"
  "${SRCROOT}${DATABASEDIR}/testParticleInput.h"
  "${SRCROOT}${DATABASEDIR}/testParticleKick.h"
)

# Add static library.
add_library(stochastic_migration_database STATIC ${DATABASE_SOURCES} ${DATABASE_HEADERS})
setup_library_target(stochastic_migration_database "${SRCROOT}${DATABASEDIR}")

# Add unit tests.
add_executable(test_DatabaseReadFunctions 
               "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestDatabaseReadFunctions.cpp")
setup_unit_test_target(test_DatabaseReadFunctions "${SRCROOT}${DATABASEDIR}")
target_link_libraries(test_DatabaseReadFunctions 
                      stochastic_migration_database stochastic_migration_basics 
                      tudat_input_output
                      ${TUDAT_CORE_LIBRARIES} 
                      ${SQLITECPP_LIBRARY} 
                      ${SQLITE3_LIBRARY}  
                      ${Boost_LIBRARIES})

#add_executable(test_DatabaseWriteFunctions 
#               "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestDatabaseWriteFunctions.cpp")
#setup_unit_test_target(test_DatabaseWriteFunctions "${SRCROOT}${DATABASEDIR}")
#target_link_libraries(test_DatabaseWriteFunctions 
#                      stochastic_migration_database stochastic_migration_basics 
#                      ${TUDAT_CORE_LIBRARIES}
#                      ${SQLITECPP_LIBRARY} 
#                      ${SQLITE3_LIBRARY} 
#                      ${Boost_LIBRARIES})

# add_executable(test_RandomWalkMonteCarloRun
#                "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestRandomWalkMonteCarloRun.cpp")
# setup_unit_test_target(test_RandomWalkMonteCarloRun "${SRCROOT}${DATABASEDIR}")
# target_link_libraries(test_RandomWalkMonteCarloRun stochastic_migration_database assist_basics
#                       ${TUDAT_CORE_LIBRARIES} ${SQLITE3_LIBRARY} ${Boost_LIBRARIES})
                      
# add_executable(test_RandomWalkPerturber 
#                "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestRandomWalkPerturber.cpp")
# setup_unit_test_target(test_RandomWalkPerturber "${SRCROOT}${DATABASEDIR}")
# target_link_libraries(test_RandomWalkPerturber stochastic_migration_database assist_basics
#                       ${TUDAT_CORE_LIBRARIES} ${SQLITE3_LIBRARY} ${Boost_LIBRARIES})

add_executable(test_RandomWalkCase 
               "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestRandomWalkCase.cpp")
setup_unit_test_target(test_RandomWalkCase "${SRCROOT}${DATABASEDIR}")
target_link_libraries(test_RandomWalkCase
                      stochastic_migration_database 
                      assist_basics 
                      ${TUDAT_CORE_LIBRARIES} 
                      ${Boost_LIBRARIES})

add_executable(test_TestParticleCase 
               "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestTestParticleCase.cpp")
setup_unit_test_target(test_TestParticleCase "${SRCROOT}${DATABASEDIR}")
target_link_libraries(test_TestParticleCase 
                      stochastic_migration_database 
                      assist_basics 
                      tudat_basic_mathematics 
                      ${TUDAT_CORE_LIBRARIES} 
                      ${Boost_LIBRARIES})

add_executable(test_TestParticleInput 
               "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestTestParticleInput.cpp")
setup_unit_test_target(test_TestParticleInput "${SRCROOT}${DATABASEDIR}")
target_link_libraries(test_TestParticleInput 
                      stochastic_migration_database 
                      assist_basics
                      tudat_basic_mathematics
                      ${TUDAT_CORE_LIBRARIES} 
                      ${Boost_LIBRARIES})

add_executable(test_TestParticleKick 
               "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestTestParticleKick.cpp")
setup_unit_test_target(test_TestParticleKick "${SRCROOT}${DATABASEDIR}")
target_link_libraries(test_TestParticleKick 
                      stochastic_migration_database 
                      assist_basics 
                      ${TUDAT_CORE_LIBRARIES} 
                      ${Boost_LIBRARIES})
