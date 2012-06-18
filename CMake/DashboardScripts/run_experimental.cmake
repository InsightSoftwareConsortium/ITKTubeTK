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
set( ENV{TUBETK_FORCE_BUILD} "0" )

set( SCRIPT_NAME "BuildTest" )
set( SCRIPT_BINARY_SUBDIR "TubeTK-Build" )
set( SCRIPT_TubeTK_USE_SUPERBUILD OFF )
include( ${TUBETK_SCRIPT_DIR}/cmakecache.cmake )
ctest_start( "$ENV{TUBETK_RUN_MODEL}" )

  ctest_read_custom_files( "${CTEST_BINARY_DIRECTORY}/.." )
  ctest_configure( BUILD "${CTEST_BINARY_DIRECTORY}/.." )
  ctest_submit( PARTS configure )

  if( SITE_EXPERIMENTAL_BUILD )
    ctest_build( BUILD "${CTEST_BINARY_DIRECTORY}/.." )
    ctest_submit( PARTS build )
  endif()

  if( SITE_EXPERIMENTAL_TEST )
    ctest_test( BUILD "${CTEST_BINARY_DIRECTORY}" )
    ctest_submit( PARTS test )
  endif()

  if( SITE_EXPERIMENTAL_COVERAGE )
    ctest_coverage( BUILD "${CTEST_BINARY_DIRECTORY}" )
    ctest_submit( PARTS Coverage )
  endif()

  if( SITE_EXPERIMENTAL_PACKAGE )
    EXECUTE_PROCESS( COMMAND make -C "${CTEST_BINARY_DIRECTORY}" package )
  endif()

  if( SITE_EXPERIMENTAL_UPLOAD )
    include( "${TUBETK_SCRIPT_DIR}/upload.cmake" )
  endif()

  if( SITE_EXPERIMENTAL_MEMORY )
    ctest_memcheck( BUILD "${CTEST_BINARY_DIRECTORY}" )
    ctest_submit( PARTS MemCheck )
  endif()

if( SITE_EXPERIMENTAL_STYLE )
  include( "${TUBETK_SCRIPT_DIR}/style.cmake" )
endif()

set(CTEST_RUN_CURRENT_SCRIPT 0)
