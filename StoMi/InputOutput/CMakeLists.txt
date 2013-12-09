 #    Copyright (c) 2010-2013, Delft University of Technology
 #    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 #    All rights reserved.
 #    See http://bit.ly/12SHPLR for license details.

# Set source files.
set(INPUTOUTPUT_SOURCES
#  "${SRCROOT}${INPUTOUTPUTDIR}/basicInputOutput.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/dictionaries.cpp"
)

# Set header files.
set(INPUTOUTPUT_HEADERS
#  "${SRCROOT}${INPUTOUTPUTDIR}/basicInputOutput.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/dictionaries.h"
)

# Add static library.
add_library(stomi_input_output STATIC ${INPUTOUTPUT_SOURCES} ${INPUTOUTPUT_HEADERS})
setup_library_target(stomi_input_output "${SRCROOT}${INPUTOUTPUTDIR}")

## Add unit tests.
#add_executable(test_BasicInputOutput 
#               "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestBasicInputOutput.cpp")
#setup_unit_test_target(test_BasicInputOutput "${SRCROOT}${INPUTOUTPUTDIR}")
#target_link_libraries(test_BasicInputOutput stomi_input_output 
#                      stomi_database ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
