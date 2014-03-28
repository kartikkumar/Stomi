 #    Copyright (c) 2010-2014, Delft University of Technology
 #    Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
 #    All rights reserved.
 #    See http://bit.ly/12SHPLR for license details.

# Set source files.
set(ASTRODYNAMICS_SOURCES
  "${SRCROOT}${ASTRODYNAMICSDIR}/propagationDataPoint.cpp"
  "${SRCROOT}${ASTRODYNAMICSDIR}/randomWalkFunctions.cpp"
)

# Set header files.
set(ASTRODYNAMICS_HEADERS
  "${SRCROOT}${ASTRODYNAMICSDIR}/hillSphere.h"
  "${SRCROOT}${ASTRODYNAMICSDIR}/propagationDataPoint.h"
  "${SRCROOT}${ASTRODYNAMICSDIR}/randomWalkFunctions.h"
)

# Add static library.
add_library(stomi_astrodynamics STATIC ${ASTRODYNAMICS_SOURCES} ${ASTRODYNAMICS_HEADERS})
setup_library_target(stomi_astrodynamics "${SRCROOT}${ASTRODYNAMICSDIR}")

# Add unit tests.
add_executable(test_RandomWalkFunctions 
               "${SRCROOT}${ASTRODYNAMICSDIR}/UnitTests/unitTestRandomWalkFunctions.cpp")
setup_unit_test_target(test_RandomWalkFunctions "${SRCROOT}${DATABASEDIR}")
target_link_libraries(test_RandomWalkFunctions 
                      stomi_astrodynamics stomi_database
                      assist_basics
                      ${TUDAT_CORE_LIBRARIES} 
                      ${SQLITECPP_LIBRARY} 
                      ${SQLITE3_LIBRARY}  
                      ${Boost_LIBRARIES})
