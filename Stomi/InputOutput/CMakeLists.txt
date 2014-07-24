# Copyright (c) 2010-2014, Delft University of Technology
# Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
# All rights reserved.
# See http://bit.ly/12SHPLR for license details.

# Set source files.
set(INPUTOUTPUT_SOURCES
  "${SOURCE_ROOT}${INPUTOUTPUT_DIR}/dictionaries.cpp"
)

# Set header files.
set(INPUTOUTPUT_HEADERS
  "${SOURCE_ROOT}${INPUTOUTPUT_DIR}/dictionaries.h"
  "${SOURCE_ROOT}${INPUTOUTPUT_DIR}/rootPath.h"
)

# Set unit test files.
set(INPUTOUTPUT_UNIT_TESTS
    "${SOURCE_ROOT}${INPUTOUTPUT_DIR}/UnitTests/unitTestInputOutput.cpp"
    "${SOURCE_ROOT}${INPUTOUTPUT_DIR}/UnitTests/unitTestRootPath.cpp"
)

# Add static library.
add_library(stomi_input_output STATIC ${INPUTOUTPUT_SOURCES} ${INPUTOUTPUT_HEADERS})
setup_library_target(stomi_input_output)

# Add unit tests.
add_executable(stomi_test_InputOutput ${INPUTOUTPUT_UNIT_TESTS})
setup_unit_test_target(stomi_test_InputOutput)
target_link_libraries(stomi_test_InputOutput 
                      stomi_input_output 
                      ${TUDAT_CORE_LIBRARIES}
                      ${Boost_LIBRARIES})
