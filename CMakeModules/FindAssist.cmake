 #    Copyright (c) 2010-2013, Delft University of Technology
 #    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 #    All rights reserved.
 #    See http://bit.ly/12SHPLR for license details.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      12xxxx    B. Tong Minh      File created based on FindEigen3.cmake.
 #      130211    K. Kumar          Updated for General Tools project. 
 #      130323    K. Kumar          Renamed file and updated for Assist project.
 #
 #    References
 #      FindEigen3.cmake.
 #
 #    Notes
 #      This script tries to find the Assist library. This module supports requiring a minimum 
 #      version, e.g. you can do version, e.g. you can do find_package(Assist 3.1.2) to require i
 #      version 3.1.2 or newer of Assist.
 #
 #      Once done, this will define:
 #
 #          ASSIST_FOUND - system has Assist lib with correct version;
 #          ASSIST_INCLUDE_DIR - the Assist include directory.
 #
 #      Original copyright statements (from FindEigen3.cmake:
 #          Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
 #          Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
 #          Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
 #
 #      FindEigen3.cmake states that redistribution and use is allowed according to the terms of
 #      the 2-clause BSD license.

macro(_assist_check_version)
  file(READ "${ASSIST_INCLUDE_DIR}/Assist/assistVersion.h" _assist_header)

  string(REGEX MATCH "define[ \t]+ASSIST_VERSION_MAJOR[ \t]+([0-9]+)" _assist_major_version_match "${_assist_header}")
  set(ASSIST_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+ASSIST_VERSION_MINOR[ \t]+([0-9]+)" _assist_minor_version_match "${_assist_header}")
  set(ASSIST_MINOR_VERSION "${CMAKE_MATCH_1}")

  set(ASSIST_VERSION ${ASSIST_MAJOR_VERSION}.${ASSIST_MINOR_VERSION})
  if(${ASSIST_VERSION} VERSION_LESS ${Assist_FIND_VERSION})
    set(ASSIST_VERSION_OK FALSE)
  else(${ASSIST_VERSION} VERSION_LESS ${Assist_FIND_VERSION})
    set(ASSIST_VERSION_OK TRUE)
  endif(${ASSIST_VERSION} VERSION_LESS ${Assist_FIND_VERSION})

  if(NOT ASSIST_VERSION_OK)

    message(STATUS "Assist version ${ASSIST_VERSION} found in ${ASSIST_INCLUDE_DIR}, "
                   "but at least version ${Assist_FIND_VERSION} is required!")
  endif(NOT ASSIST_VERSION_OK)

  set(ASSIST_LIBRARIES "assist")
  link_directories(${ASSIST_LIBRARIES_DIR})
endmacro(_assist_check_version)

if (ASSIST_INCLUDE_DIR)

  # in cache already
  _assist_check_version()
  set(ASSIST_FOUND ${ASSIST_VERSION_OK})

else (ASSIST_INCLUDE_DIR)

  find_path(ASSIST_BASE_PATH NAMES assistVersion.h
      PATHS
      ${PROJECT_SOURCE_DIR}/../assist
      ${PROJECT_SOURCE_DIR}/../../assist
      ${PROJECT_SOURCE_DIR}/../../../assist
      ${PROJECT_SOURCE_DIR}/../../../../assist
      PATH_SUFFIXES Assist 
    )
  set(ASSIST_INCLUDE_DIR ${ASSIST_BASE_PATH}/../)
  set(ASSIST_LIBRARIES_DIR ${ASSIST_BASE_PATH}/../lib)

  if(ASSIST_INCLUDE_DIR)
    _assist_check_version()
  endif(ASSIST_INCLUDE_DIR)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Assist DEFAULT_MSG ASSIST_INCLUDE_DIR ASSIST_VERSION_OK)

  mark_as_advanced(ASSIST_INCLUDE_DIR)

endif(ASSIST_INCLUDE_DIR)
