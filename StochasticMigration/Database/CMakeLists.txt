 #    Copyright (c) 2010-2013, Delft University of Technology
 #    All rights reserved.
 #
 #    Redistribution and use in source and binary forms, with or without modification, are
 #    permitted provided that the following conditions are met:
 #      - Redistributions of source code must retain the above copyright notice, this list of
 #        conditions and the following disclaimer.
 #      - Redistributions in binary form must reproduce the above copyright notice, this list of
 #        conditions and the following disclaimer in the documentation and/or other materials
 #        provided with the distribution.
 #      - Neither the name of the Delft University of Technology nor the names of its contributors
 #        may be used to endorse or promote products derived from this software without specific
 #        prior written permission.
 #
 #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 #    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 #    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 #    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 #    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 #    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 #    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 #    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 #    OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      120214    K. Kumar          Adapted file to work with MabSimulations/includes/ sub-directory.
 #      120523    K. Kumar          Adapted for Database/ sub-directory.
 #      130217    K. Kumar          Updated static library prefix to stochastic_migration_.
 #
 #    References
 #
 #    Notes
 #

# Set source files.
set(DATABASE_SOURCES
  "${SRCROOT}${DATABASEDIR}/databaseHelpFunctions.cpp"
  "${SRCROOT}${DATABASEDIR}/databaseReadFunctions.cpp"
  "${SRCROOT}${DATABASEDIR}/databaseWriteFunctions.cpp"
  "${SRCROOT}${DATABASEDIR}/randomWalkMonteCarloRun.cpp"
  "${SRCROOT}${DATABASEDIR}/randomWalkPerturber.cpp"
  "${SRCROOT}${DATABASEDIR}/testParticleCase.cpp"
  "${SRCROOT}${DATABASEDIR}/testParticleInput.cpp"
  "${SRCROOT}${DATABASEDIR}/testParticleKick.cpp"
)

# Set header files.
set(DATABASE_HEADERS
  "${SRCROOT}${DATABASEDIR}/databaseHelpFunctions.h"
  "${SRCROOT}${DATABASEDIR}/databaseReadFunctions.h"
  "${SRCROOT}${DATABASEDIR}/databaseWriteFunctions.h"
  "${SRCROOT}${DATABASEDIR}/randomWalkMonteCarloRun.h"
  "${SRCROOT}${DATABASEDIR}/randomWalkPerturber.h"
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
target_link_libraries(test_DatabaseReadFunctions stochastic_migration_database 
                      stochastic_migration_basics ${SQLITECPP_LIBRARY} ${TUDAT_CORE_LIBRARIES} 
                      ${SQLITE3_LIBRARY}  ${Boost_LIBRARIES})

add_executable(test_DatabaseWriteFunctions 
               "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestDatabaseWriteFunctions.cpp")
setup_unit_test_target(test_DatabaseWriteFunctions "${SRCROOT}${DATABASEDIR}")
target_link_libraries(test_DatabaseWriteFunctions stochastic_migration_database 
                      stochastic_migration_basics ${SQLITECPP_LIBRARY} ${TUDAT_CORE_LIBRARIES}
                      ${SQLITE3_LIBRARY} ${Boost_LIBRARIES})

add_executable(test_RandomWalkMonteCarloRun
               "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestRandomWalkMonteCarloRun.cpp")
setup_unit_test_target(test_RandomWalkMonteCarloRun "${SRCROOT}${DATABASEDIR}")
target_link_libraries(test_RandomWalkMonteCarloRun stochastic_migration_database assist_basics
                      ${TUDAT_CORE_LIBRARIES} ${SQLITE3_LIBRARY} ${Boost_LIBRARIES})
                      
add_executable(test_RandomWalkPerturber 
               "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestRandomWalkPerturber.cpp")
setup_unit_test_target(test_RandomWalkPerturber "${SRCROOT}${DATABASEDIR}")
target_link_libraries(test_RandomWalkPerturber stochastic_migration_database assist_basics
                      ${TUDAT_CORE_LIBRARIES} ${SQLITE3_LIBRARY} ${Boost_LIBRARIES})

add_executable(test_TestParticleCase 
               "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestTestParticleCase.cpp")
setup_unit_test_target(test_TestParticleCase "${SRCROOT}${DATABASEDIR}")
target_link_libraries(test_TestParticleCase stochastic_migration_database assist_basics
                      ${TUDAT_CORE_LIBRARIES} ${SQLITE3_LIBRARY} ${Boost_LIBRARIES})

add_executable(test_TestParticleKick 
               "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestTestParticleKick.cpp")
setup_unit_test_target(test_TestParticleKick "${SRCROOT}${DATABASEDIR}")
target_link_libraries(test_TestParticleKick stochastic_migration_database assist_basics 
                      ${TUDAT_CORE_LIBRARIES} ${SQLITE3_LIBRARY} ${Boost_LIBRARIES})

add_executable(test_TestParticleInput 
               "${SRCROOT}${DATABASEDIR}/UnitTests/unitTestTestParticleInput.cpp")
setup_unit_test_target(test_TestParticleInput "${SRCROOT}${DATABASEDIR}")
target_link_libraries(test_TestParticleInput stochastic_migration_database assist_basics
                      ${TUDAT_CORE_LIBRARIES} ${SQLITE3_LIBRARY} ${Boost_LIBRARIES})
