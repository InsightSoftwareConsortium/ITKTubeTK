##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

set( ENV{TUBETK_RUN_MODEL} "Nightly" )

set( SCRIPT_NAME "BuildTest" )
set( SCRIPT_BINARY_SUBDIR "TubeTK-Build" )
set( SCRIPT_TubeTK_USE_SUPERBUILD ON )
include( ${TUBETK_SCRIPT_DIR}/cmakecache.cmake )
ctest_start( "$ENV{TUBETK_RUN_MODEL}" )

if( SITE_NIGHTLY_BUILD )
  ctest_empty_binary_directory( "${TUBETK_BINARY_DIR}" )
  ctest_update( SOURCE ${CTEST_SOURCE_DIRECTORY} )
  ctest_configure( BUILD "${CTEST_BINARY_DIRECTORY}/.." )
  ctest_read_custom_files( "${CTEST_BINARY_DIRECTORY}/.." )
  ctest_build( BUILD "${CTEST_BINARY_DIRECTORY}/.." )
else()
  ctest_read_custom_files( "${CTEST_BINARY_DIRECTORY}/.." )
endif()

if( SITE_NIGHTLY_TEST )
  ctest_test( BUILD "${CTEST_BINARY_DIRECTORY}" )
endif()

if( SITE_NIGHTLY_COVERAGE )
  ctest_coverage( BUILD "${CTEST_BINARY_DIRECTORY}" )
endif()

if( SITE_NIGHTLY_MEMORY )
  ctest_memcheck( BUILD "${CTEST_BINARY_DIRECTORY}" )
endif()

function( TubeTK_Package )
  execute_process(
    COMMAND ${CMAKE_COMMAND} --build ${CTEST_BINARY_DIRECTORY} --target package --config ${CTEST_BUILD_CONFIGURATION}
    WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_FILE CPackOutputFiles.txt
    )
endfunction( TubeTK_Package )

function( TubeTK_Upload )
  set(package_list)
  set(regexp ".*CPack: - package: (.*) generated\\.")
  set(raw_package_list)
  file(STRINGS ${CTEST_BINARY_DIRECTORY}/CPackOutputFiles.txt raw_package_list REGEX ${regexp})
  foreach(package ${raw_package_list})
    string(REGEX REPLACE ${regexp} "\\1" package_path "${package}" )
    list(APPEND package_list ${package_path})
  endforeach()
  ctest_upload( FILES ${package_list} )
endfunction( TubeTK_Upload )

if( SITE_NIGHTLY_PACKAGE )
  TubeTK_Package()
endif()

if( SITE_NIGHTLY_UPLOAD )
  TubeTK_Upload()
endif()

ctest_submit()

if( SITE_NIGHTLY_STYLE )
  include( "${TUBETK_SCRIPT_DIR}/style.cmake" )
endif()

set(CTEST_RUN_CURRENT_SCRIPT 0)
