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
 #      120214    K. Kumar          Adapted to work with Tudat and Tudat Core.
 #      130211    K. Kumar          Added GeneralTools library.
 #      130217    K. Kumar          Renamed project StochasticMigration.
 #
 #    References
 #
 #    Notes
 #

# Specify minimum CMake version required.
cmake_minimum_required(VERSION 2.6)

# Specify project name.
project( StochasticMigration )

# Set module path to local CMake scripts.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../CMakeModules)

# Set root-directory for code to current source directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Define the basic input (trunk).
set(SRCROOT "${CODEROOT}/StochasticMigration")

# Set lib and bin directories where static libraries, applications, and unit tests are built.
set(LIBROOT "${CODEROOT}/lib")
set(BINROOT "${CODEROOT}/bin")

# Set testing options based on platform.
enable_testing( )

# Set platform-specific compiler flags.
if (WIN32)
  if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /Ox /W3 /FC -D_SCL_SECURE_NO_WARNINGS")
    # Because we are using static boost libraries, with static runtime, we need to force MSVC to
    # also use static runtime: (from http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace)
    foreach(flag_var
          CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
          CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      # Find all dynamic runtime (MD) references and replace with static (MT)
      if(${flag_var} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
      endif(${flag_var} MATCHES "/MD")
    endforeach(flag_var)
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
  endif( )
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
endif( )

# ------------------------------------------------------------------------------------
# Don't define any new variables with SET that the user might change beyond this point.
#
# Load user specific settings. These can overwrite all previous settings.
if (NOT ASSIST_USER_SETTINGS AND NOT DEFINED ENV{ASSIST_USER_SETTINGS})
  # Load from default paths, least important first.
  INCLUDE( "../UserSettings.txt" OPTIONAL)
  INCLUDE( "UserSettings.txt" OPTIONAL)
else()
  MESSAGE(STATUS "Loading preset")
  # Load the user specified file.
  if ( ASSIST_USER_SETTINGS )
    INCLUDE( $ASSIST_USER_SETTINGS OPTIONAL)
  else()
    INCLUDE( $ENV{ASSIST_USER_SETTINGS} OPTIONAL)
  endif()
endif()

# Set CMake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release". NB: It is important this
# section is placed AFTER the "UserSettings.txt" file is included, to ensure that that the correct
# indiciation is given for which build type has been selected.
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    message(STATUS "WARNING: building release version!")
    set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "WARNING: building debug version!")
endif()

# Set the global macros for setting up targets.
macro(setup_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}/applications")
endmacro(setup_executable_target)

macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
endmacro(setup_library_target)

macro(setup_unit_test_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${BINROOT}/unit_tests/${target_name}")
endmacro(setup_unit_test_target)

# Set the main sub-directories.
#set(ASTRODYNAMICSDIR "/Astrodynamics")
set(BASICSDIR "/Basics")
set(DATABASEDIR "/Database")
set(INPUTOUTPUTDIR "/InputOutput")

# Include the top-level directories.
include_directories(AFTER
  "${CODEROOT}"
)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif( )

# Configure Boost libraries.
if(NOT Boost_USE_STATIC_LIBS)
  set(Boost_USE_STATIC_LIBS ON)
endif()
if(NOT Boost_USE_MULTITHREADED)
  set(Boost_USE_MULTITHREADED ON)
endif()
if(NOT Boost_USE_STATIC_RUNTIME)
  set(Boost_USE_STATIC_RUNTIME ON)
endif()

# Find Boost libraries on local system.
find_package(Boost 1.52.0 COMPONENTS thread date_time system unit_test_framework filesystem regex serialization REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif( )

# Find Tudat Core library on local system.
find_package(TudatCore 2 REQUIRED)

# Include TudatCore directories.
# Set CMake flag to suppress TudatCore warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${TUDAT_CORE_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_CORE_INCLUDE_DIR}\"")
endif( )

# Find Tudat library on local system.
find_package(Tudat 2 REQUIRED)

# Include Tudat directories.
# Set CMake flag to suppress Tudat warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif( )

# Find Assist library on local system.
find_package(Assist 1 REQUIRED)

# Include Assist directories.
# Set CMake flag to suppress Assist warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${ASSIST_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${ASSIST_INCLUDE_DIR}\"")
endif( )

# Find SQLite 3 library on local system.
find_package(sqlite3 REQUIRED)

# Include SQLite3 directories.
# Set CMake flag to suppress SQLite3 warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${SQLITE3_INCLUDE_DIRS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SQLITE3_INCLUDE_DIR}\"")
endif()

# Set flags to compile using openmp.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp")

# Set flags to profile using gprof.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

# Set sub-directories.
set(SUBDIRS ${BASICSDIR} ${DATABASEDIR} ${INPUTOUTPUTDIR})
#${ASTRODYNAMICSDIR}  ${DATABASEDIR} )

# Add sub-directories to CMake process.
foreach(CURRENT_SUBDIR ${SUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_SUBDIR}")
endforeach(CURRENT_SUBDIR)

# Add sandbox.
add_executable(sandbox "${SRCROOT}/Applications/sandbox.cpp")
setup_executable_target(sandbox "${SRCROOT}")
target_link_libraries(sandbox tudat_gravitation ${SQLITE3_LIBRARY} ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

## Add SQLite database generator for stochastic migration simulations.
#add_executable(stochasticMigrationDatabaseGenerator "${SRCROOT}/Applications/stochasticMigrationDatabaseGenerator.cpp")
#setup_executable_target(stochasticMigrationDatabaseGenerator "${SRCROOT}")
#target_link_libraries(stochasticMigrationDatabaseGenerator stochastic_migration_input_output stochastic_migration_database general_tools_input_output general_tools_database general_tools_basics tudat_gravitation tudat_input_output tudat_basic_astrodynamics ${SQLITE3_LIBRARY} ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

## Add test particle simulator (uses OpenMP).
#add_executable(testParticleSimulator "${SRCROOT}/Applications/testParticleSimulator.cpp")
#setup_executable_target(testParticleSimulator "${SRCROOT}")
#target_link_libraries(testParticleSimulator stochastic_migration_database stochastic_migration_input_output general_tools_astrodynamics general_tools_input_output tudat_state_derivative_models tudat_basic_mathematics tudat_basic_astrodynamics tudat_gravitation tudat_numerical_integrators tudat_input_output ${SQLITE3_LIBRARY} ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

## Add Mab random walk application (uses openmp).
#add_executable(randomWalkSimulator "${SRCROOT}/Applications/randomWalkSimulator.cpp")
#setup_executable_target(randomWalkSimulator "${SRCROOT}")
#target_link_libraries(randomWalkSimulator mab_database mab_input_output mab_simulation_includes tudat_basic_astrodynamics tudat_gravitation tudat_statistics tudat_numerical_integrators tudat_input_output ${SQLITE3_LIBRARY} ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

## Add verification random walk simulator (uses openmp).
#add_executable(randomWalkCheckSimulator "${SRCROOT}/Applications/randomWalkCheckSimulator.cpp")
#setup_executable_target(randomWalkCheckSimulator "${SRCROOT}")
#target_link_libraries(randomWalkCheckSimulator mab_database mab_input_output mab_simulation_includes tudat_basic_astrodynamics tudat_root_finding_methods tudat_gravitation tudat_statistics tudat_numerical_integrators ${SQLITE3_LIBRARY} ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
