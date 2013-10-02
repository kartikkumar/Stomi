 #    Copyright (c) 2010-2013, Delft University of Technology
 #    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 #    All rights reserved.
 #    See http://bit.ly/12SHPLR for license details.

# Specify minimum CMake version required.
cmake_minimum_required(VERSION 2.6)

# Specify project name.
project(StochasticMigration)

# Set module path to local CMake scripts.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../CMakeModules)

# Define the directory with the source code.
set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Define the code root directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Set testing options based on platform.
enable_testing()

# Set platform-specific compiler flags.
if (WIN32)
  if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /Ox /W3 /FC -D_SCL_SECURE_NO_WARNINGS")
    # Because we are using static boost libraries, with static runtime, we need to force MSVC to
    # also use static runtime: (from http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace).
    foreach(flag_var
          CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
          CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      # Find all dynamic runtime (MD) references and replace with static (MT)
      if(${flag_var} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
      endif(${flag_var} MATCHES "/MD")
    endforeach(flag_var)
  else()
    set(CMAKE_CXX_FLAGS 
        "${CMAKE_CXX_FLAGS} -Wall -g -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
  endif()
else()
  set(CMAKE_CXX_FLAGS 
      "${CMAKE_CXX_FLAGS} -Wall -g -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
endif()

# -------------------------------------------------------------------------------------------------
# DO NOT MANUALLY EDIT ANY SETTINGS AFTER THIS POINT IN THE FILE.
# 
# Pass any additional settings directly to CMake (e.g., via command-line) or by using a config file
# (e.g., UserSettings.txt).
# -------------------------------------------------------------------------------------------------

# Load user specific settings. These can overwrite all previous settings.
if (NOT STOCHASTIC_MIGRATION_USER_SETTINGS AND NOT DEFINED ENV{STOCHASTIC_MIGRATION_USER_SETTINGS})
  # Load from default paths, least important first.
  include( "../UserSettings.txt" OPTIONAL)
  include( "UserSettings.txt" OPTIONAL)
else()
  MESSAGE(STATUS "Loading preset...")
  # Load the user specified file.
  if (STOCHASTIC_MIGRATION_USER_SETTINGS)
    include($STOCHASTIC_MIGRATION_USER_SETTINGS OPTIONAL)
  else()
    include($ENV{STOCHASTIC_MIGRATION_USER_SETTINGS} OPTIONAL)
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

# Set lib and bin directories where static libraries and unit tests are built.
if(NOT LIB_ROOT)
  set(LIB_ROOT "${PROJECT_BINARY_DIR}/../lib")
endif()
if(NOT BIN_ROOT)
  set(BIN_ROOT "${PROJECT_BINARY_DIR}/../bin")
endif()

# Set the global macros for setting up targets.
macro(setup_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BIN_ROOT}/applications")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BIN_ROOT}/applications")
endmacro(setup_executable_target)

macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIB_ROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIB_ROOT}")
endmacro(setup_library_target)

macro(setup_unit_test_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BIN_ROOT}/unit_tests")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${BIN_ROOT}/unit_tests/${target_name}")
endmacro(setup_unit_test_target)

if(NOT TUDAT_BUNDLE_DISTRIBUTION_PATH)
    set(TUDAT_BUNDLE_DISTRIBUTION_PATH "${CODEROOT}")
endif(NOT TUDAT_BUNDLE_DISTRIBUTION_PATH)

if(NOT STOCHASTIC_MIGRATION_DISTRIBUTION_PATH)
    set(STOCHASTIC_MIGRATION_DISTRIBUTION_PATH 
        "${TUDAT_BUNDLE_DISTRIBUTION_PATH}/tudatApplications/stochasticMigration")
endif(NOT STOCHASTIC_MIGRATION_DISTRIBUTION_PATH)

# Install files.
install(DIRECTORY "${SRCROOT}/"
        DESTINATION "${STOCHASTIC_MIGRATION_DISTRIBUTION_PATH}/StochasticMigration"
        PATTERN ".DS_STORE" EXCLUDE
        PATTERN "CMakeLists.txt.user" EXCLUDE
        PATTERN ".svn" EXCLUDE
        PATTERN ".git" EXCLUDE
        PATTERN ".bzr" EXCLUDE
)

install(FILES ${CODEROOT}/LICENSE ${CODEROOT}/NOTICE ${CODEROOT}/README 
              ${CODEROOT}/UserSettings.txt.example
        DESTINATION "${STOCHASTIC_MIGRATION_DISTRIBUTION_PATH}")

# Set the main sub-directories.
set(ASTRODYNAMICSDIR "/Astrodynamics")
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
find_package(Boost 1.52.0 
             COMPONENTS thread date_time system unit_test_framework 
                        filesystem regex serialization REQUIRED)

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
find_package(Assist 1.1 REQUIRED)

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
  include_directories(SYSTEM AFTER "${SQLITE3_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SQLITE3_INCLUDE_DIR}\"")
endif()

# Find SQLiteCpp library on local system.
find_package(SQLiteCpp 0.050 REQUIRED)

# Include SQLiteCpp directories.
# Set CMake flag to suppress SQLiteCpp warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${SQLITECPP_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SQLITECPP_INCLUDE_DIR}\"")
endif( )

# If documentation should be built, find Doxygen package and setup config file.
if(BUILD_DOCUMENTATION)
  find_package(Doxygen)

  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR 
            "Doxygen is needed to build the documentation. Please install it correctly.")
  endif()

  # Configure the Template Doxyfile for our specific project.
  configure_file(${CODEROOT}/Doxyfile.in 
                 ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)

  # Add a custom target to run Doxygen when ever the project is built.
  # IF you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the above command.
  add_custom_target (docs ALL COMMAND ${DOXYGEN_EXECUTABLE} 
                     ${PROJECT_BINARY_DIR}/Doxyfile 
                     SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)

  # Install documentation.
  install(DIRECTORY "${CODEROOT}/Documentation" 
          DESTINATION "${STOCHASTIC_MIGRATION_DISTRIBUTION_PATH}")
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
set(SUBDIRS ${ASTRODYNAMICSDIR} ${BASICSDIR} ${DATABASEDIR} ${INPUTOUTPUTDIR})

# Add sub-directories to CMake process.
foreach(CURRENT_SUBDIR ${SUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_SUBDIR}")
endforeach(CURRENT_SUBDIR)

# Add SQLite database generator for test particle simulations.
add_executable(testParticleDatabaseGenerator 
               "${SRCROOT}/Applications/testParticleDatabaseGenerator.cpp")
setup_executable_target(testParticleDatabaseGenerator "${SRCROOT}")
target_link_libraries(testParticleDatabaseGenerator 
                      stochastic_migration_astrodynamics stochastic_migration_input_output 
                      assist_astrodynamics assist_input_output assist_mathematics  
                      tudat_basic_astrodynamics tudat_input_output
                      ${TUDAT_CORE_LIBRARIES} 
                      ${SQLITECPP_LIBRARY} 
                      ${SQLITE3_LIBRARY} 
                      ${Boost_LIBRARIES})

# Add test particle simulator (uses OpenMP).
add_executable(testParticleSimulator "${SRCROOT}/Applications/testParticleSimulator.cpp")
setup_executable_target(testParticleSimulator "${SRCROOT}")
target_link_libraries(testParticleSimulator 
                      stochastic_migration_astrodynamics stochastic_migration_basics 
                      stochastic_migration_database stochastic_migration_input_output 
                      assist_astrodynamics assist_basics assist_input_output
                      tudat_state_derivative_models tudat_basic_mathematics 
                      tudat_basic_astrodynamics tudat_gravitation tudat_numerical_integrators 
                      tudat_input_output 
                      ${TUDAT_CORE_LIBRARIES} 
                      ${SQLITECPP_LIBRARY} 
                      ${SQLITE3_LIBRARY} 
                      ${Boost_LIBRARIES})

# # Add Mab random walk application (uses OpenMP).
# add_executable(randomWalkSimulator "${SRCROOT}/Applications/randomWalkSimulator.cpp")
# setup_executable_target(randomWalkSimulator "${SRCROOT}")
# target_link_libraries(randomWalkSimulator
#                       stochastic_migration_astrodynamics stochastic_migration_database 
#                       stochastic_migration_input_output
#                       assist_astrodynamics assist_basics assist_input_output assist_mathematics
#                       tudat_input_output
#                       ${TUDAT_CORE_LIBRARIES} 
#                       ${SQLITECPP_LIBRARY} 
#                       ${SQLITE3_LIBRARY} 
#                       ${Boost_LIBRARIES})
