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

cmake_minimum_required(VERSION 2.6)

include( ${CTEST_SCRIPT_DIRECTORY}/../../tubetk_config.cmake )

set( CTEST_SITE "${SITE_NAME}" )
set( CTEST_BUILD_NAME "${SITE_BUILD_NAME}" )
set( CTEST_BUILD_CONFIGURATION "${SITE_BUILD_TYPE}" )

set( CTEST_TEST_TIMEOUT 1500 )

set( CTEST_CMAKE_GENERATOR "Unix Makefiles" )
set( CTEST_NOTES_FILES "${SITE_SCRIPT_DIR}/${CTEST_SCRIPT_NAME}" )
set( CTEST_SOURCE_DIRECTORY "${SITE_SOURCE_DIR}" )
set( CTEST_BINARY_DIRECTORY "${SITE_BINARY_DIR}" )

set( ENV{DISPLAY} ":0" )

set( CTEST_BUILD_COMMAND "${SITE_MAKE_COMMAND}" )
set( CTEST_CMAKE_COMMAND "${SITE_CMAKE_COMMAND}" )
set( CTEST_QMAKE_COMMAND "${SITE_QMAKE_COMMAND}" )
set( CTEST_UPDATE_COMMAND "${SITE_UPDATE_COMMAND}" )
set( CTEST_COVERAGE_COMMAND "${SITE_COVERAGE_COMMAND}" )

set( CTEST_MEMORYCHECK_COMMAND "${SITE_MEMORYCHECK_COMMAND}" )
set( CTEST_MEMORYCHECK_SUPPRESSIONS_FILE
   "${CTEST_SCRIPT_DIRECTORY}/tubetk_valgrind_supression.txt" )
set( CTEST_MEMORYCHECK_COMMAND_OPTIONS "${MEMORYCHECK_OPTIONS}" )

message( "Running memory check..." )

ctest_start( Nightly )
ctest_read_custom_files( "${CTEST_BINARY_DIRECTORY}" )
ctest_memcheck( BUILD "${CTEST_BINARY_DIRECTORY}" )
ctest_submit()

