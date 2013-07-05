 #
 #    Copyright (c) 2010-2013, Delft University of Technology
 #    All rights reserved.
 #    See COPYING for license details.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      120214    K. Kumar          Adapted file to work with MabSimulations/includes/
 #                                  sub-directory.
 #      120523    K. Kumar          Adapted for InputOutput/ sub-directory.
 #      130217    K. Kumar          Updated static library prefix to stochastic_migration_.
 #
 #    References
 #
 #    Notes
 #

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
add_library(stochastic_migration_input_output STATIC ${INPUTOUTPUT_SOURCES} ${INPUTOUTPUT_HEADERS})
setup_library_target(stochastic_migration_input_output "${SRCROOT}${INPUTOUTPUTDIR}")

## Add unit tests.
#add_executable(test_BasicInputOutput 
#               "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestBasicInputOutput.cpp")
#setup_unit_test_target(test_BasicInputOutput "${SRCROOT}${INPUTOUTPUTDIR}")
#target_link_libraries(test_BasicInputOutput stochastic_migration_input_output 
#                      stochastic_migration_database ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
