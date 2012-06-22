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

set( ENV{TUBETK_RUN_MODEL} "Experimental" )

set( SCRIPT_NAME "BuildTest" )
include( ${TUBETK_SCRIPT_DIR}/cmakecache.cmake )
ctest_start( "$ENV{TUBETK_RUN_MODEL}" )

if( SITE_EXPERIMENTAL_BUILD )
  #ctest_update( SOURCE ${CTEST_SOURCE_DIRECTORY} )
  ctest_configure( BUILD "${TUBETK_BINARY_DIR}" )
  ctest_read_custom_files( "${TUBETK_BINARY_DIR}" )
  ctest_build( BUILD "${TUBETK_BINARY_DIR}" )
else()
  ctest_read_custom_files( "${TUBETK_BINARY_DIR}" )
endif()

if( SITE_EXPERIMENTAL_TEST )
  ctest_test( BUILD "${TUBETK_BINARY_DIR}/TubeTK-Build" )
endif()

if( SITE_EXPERIMENTAL_COVERAGE )
  ctest_coverage( BUILD "${TUBETK_BINARY_DIR}/TubeTK-Build" )
endif()

if( SITE_EXPERIMENTAL_MEMORY )
  ctest_memcheck( BUILD "${TUBETK_BINARY_DIR}/TubeTK-Build" )
endif()

function( TubeTK_Package )
  execute_process(
    COMMAND ${CMAKE_COMMAND} --build ${TUBETK_BINARY_DIR}/TubeTK-Build --target package --config ${TUBETK_BINARY_DIR}/TubeTK-Build
    WORKING_DIRECTORY ${TUBETK_BINARY_DIR}/TubeTK-Build
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_FILE CPackOutputFiles.txt
    )
endfunction( TubeTK_Package )

function( TubeTK_Upload )
  set(package_list)
  set(regexp ".*CPack: - package: (.*) generated\\.")
  set(raw_package_list)
  file(STRINGS ${TUBETK_BINARY_DIR}/TubeTK-Build/CPackOutputFiles.txt raw_package_list REGEX ${regexp})
  foreach(package ${raw_package_list})
    string(REGEX REPLACE ${regexp} "\\1" package_path "${package}" )
    list(APPEND package_list ${package_path})
  endforeach()
  ctest_upload( FILES ${package_list} )
endfunction( TubeTK_Upload )

if( SITE_EXPERIMENTAL_PACKAGE )
  TubeTK_Package()
endif()

if( SITE_EXPERIMENTAL_UPLOAD )
  TubeTK_Upload()
endif()

ctest_submit()

if( SITE_EXPERIMENTAL_STYLE )
  include( "${TUBETK_SCRIPT_DIR}/style.cmake" )
endif()

set(CTEST_RUN_CURRENT_SCRIPT 0)
