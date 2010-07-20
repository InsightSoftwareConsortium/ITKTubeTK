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

cmake_minimum_required( VERSION 2.6 )

include( ${CTEST_SCRIPT_DIRECTORY}/../../tubetk_config.cmake )
include( ${SITE_BINARY_DIR}/tubetk_runtime_config.cmake )

set( SCRIPT_NAME "BuildTest" )
set( SCRIPT_BINARY_SUBDIR "" )
set( SCRIPT_TubeTK_USE_SUPERBUILD ON )

include( ${CTEST_SCRIPT_DIRECTORY}/tubetk_cmakecache.cmake )

ctest_start( "${RUN_MODEL}" )
ctest_update( SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res )
if( res GREATER 0 OR res LESS 0 OR RUN_MODEL STREQUAL "Experimental" )

  message( "Changes detected. Rebuilding and testing..." )

  ctest_configure( BUILD "${CTEST_BINARY_DIRECTORY}" )
  ctest_read_custom_files( "${CTEST_BINARY_DIRECTORY}" )
  ctest_build( BUILD "${CTEST_BINARY_DIRECTORY}" )
  ctest_submit()
  ctest_test( BUILD "${CTEST_BINARY_DIRECTORY}/TubeTK-Build" )
  ctest_submit()
  set( res 0 )
  set( ENV{TUBETK_CONTINUOUS_UPDATE} "1" )

else()

  message( "No changes detected." )
  set( ENV{TUBETK_CONTINUOUS_UPDATE} "0" )

endif()

