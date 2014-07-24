# Copyright (c) 2010-2014, Delft University of Technology
# Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
# All rights reserved.
# See http://bit.ly/12SHPLR for license details.

# Set source files.
set(ASTRODYNAMICS_SOURCES
  "${SOURCE_ROOT}${ASTRODYNAMICS_DIR}/propagationDataPoint.cpp"
  "${SOURCE_ROOT}${ASTRODYNAMICS_DIR}/randomWalkFunctions.cpp"
)

# Set header files.
set(ASTRODYNAMICS_HEADERS
  "${SOURCE_ROOT}${ASTRODYNAMICS_DIR}/propagationDataPoint.h"
  "${SOURCE_ROOT}${ASTRODYNAMICS_DIR}/randomWalkFunctions.h"
)

# Set unit test files.
set(ASTRODYNAMICS_UNIT_TESTS
  "${SOURCE_ROOT}${ASTRODYNAMICS_DIR}/UnitTests/unitTestAstrodynamics.cpp"  
  "${SOURCE_ROOT}${ASTRODYNAMICS_DIR}/UnitTests/unitTestBody.cpp"
  "${SOURCE_ROOT}${ASTRODYNAMICS_DIR}/UnitTests/unitTestDataUpdater.cpp"
  "${SOURCE_ROOT}${ASTRODYNAMICS_DIR}/UnitTests/unitTestRandomWalkFunctions.cpp"  
)

# Add static library.
add_library(stomi_astrodynamics STATIC ${ASTRODYNAMICS_SOURCES} ${ASTRODYNAMICS_HEADERS})
setup_library_target(stomi_astrodynamics)

# Add unit tests.
add_executable(stomi_test_Astrodynamics ${ASTRODYNAMICS_UNIT_TESTS})
setup_unit_test_target(stomi_test_Astrodynamics)
target_link_libraries(stomi_test_Astrodynamics 
                      stomi_astrodynamics stomi_database
                      ${TUDAT_CORE_LIBRARIES} 
                      ${SQLITECPP_LIBRARY} 
                      ${SQLITE3_LIBRARY}  
                      ${Boost_LIBRARIES})
