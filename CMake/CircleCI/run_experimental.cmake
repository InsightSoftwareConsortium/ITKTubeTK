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

set( CTEST_BUILD_NAME "${SITE_BUILD_NAME}-BuildTest-Experimental" )

if( SITE_EXPERIMENTAL_BUILD_DOCUMENTATION )
  set( BUILD_DOCUMENTATION ON )
  set( TubeTK_USE_DOXYGEN ON )
else( SITE_EXPERIMENTAL_BUILD_DOCUMENTATION )
  set( BUILD_DOCUMENTATION OFF )
  set( TubeTK_USE_DOXYGEN OFF )
endif( SITE_EXPERIMENTAL_BUILD_DOCUMENTATION )

if( SITE_EXPERIMENTAL_CPPCHECK )
  set( TubeTK_USE_CPPCHECK ON )
else( SITE_EXPERIMENTAL_CPPCHECK )
  set( TubeTK_USE_CPPCHECK OFF )
endif( SITE_EXPERIMENTAL_CPPCHECK )

if( SITE_EXPERIMENTAL_KWSTYLE )
  set( TubeTK_USE_KWSTYLE ON )
else( SITE_EXPERIMENTAL_KWSTYLE )
  set( TubeTK_USE_KWSTYLE OFF )
endif( SITE_EXPERIMENTAL_KWSTYLE )

if( SITE_EXPERIMENTAL_MEMORY )
  set( TubeTK_USE_VALGRIND ON )
else( SITE_EXPERIMENTAL_MEMORY )
  set( TubeTK_USE_VALGRIND OFF )
endif( SITE_EXPERIMENTAL_MEMORY )

set( BUILD_TESTING ON )

configure_file( ${TubeTK_SOURCE_DIR}/CMake/InitCMakeCache.cmake.in
  ${TubeTK_BINARY_DIR}/InitCMakeCache.cmake IMMEDIATE @ONLY )
set( CTEST_NOTES_FILES "${TubeTK_BINARY_DIR}/InitCMakeCache.cmake" )

ctest_start( "Experimental" )

if( SITE_EXPERIMENTAL_BUILD )
  ctest_configure( BUILD "${TubeTK_BINARY_DIR}"
    SOURCE "${TubeTK_SOURCE_DIR}"
    OPTIONS "-C${TubeTK_BINARY_DIR}/InitCMakeCache.cmake" )
  ctest_read_custom_files( "${TubeTK_BINARY_DIR}" )
  ctest_build( BUILD "${TubeTK_BINARY_DIR}" )
  ctest_submit( PARTS Notes Configure Build )
else( SITE_EXPERIMENTAL_BUILD )
  ctest_read_custom_files( "${TubeTK_BINARY_DIR}" )
  ctest_submit( PARTS Notes )
endif( SITE_EXPERIMENTAL_BUILD )

if( SITE_EXPERIMENTAL_TEST )
  ctest_test( BUILD "${TubeTK_BINARY_DIR}/TubeTK-build" )
  ctest_submit( PARTS Test )
endif( SITE_EXPERIMENTAL_TEST )

if( SITE_EXPERIMENTAL_COVERAGE )
  ctest_coverage( BUILD "${TubeTK_BINARY_DIR}/TubeTK-build" )
  ctest_submit( PARTS Coverage )
endif( SITE_EXPERIMENTAL_COVERAGE )

if( SITE_EXPERIMENTAL_MEMORY )
  ctest_memcheck( BUILD "${TubeTK_BINARY_DIR}/TubeTK-build" )
  ctest_submit( PARTS MemCheck )
endif( SITE_EXPERIMENTAL_MEMORY )

if( SITE_EXPERIMENTAL_PACKAGE )
  execute_process( COMMAND ${CMAKE_COMMAND}
      --build ${TubeTK_BINARY_DIR}/TubeTK-build
      --target package
    WORKING_DIRECTORY ${TubeTK_BINARY_DIR}/TubeTK-build
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_FILE CPackOutputFiles.txt )
endif( SITE_EXPERIMENTAL_PACKAGE )

if( SITE_EXPERIMENTAL_UPLOAD )
  set( package_list)
  set( regexp ".*CPack: - package: (.*) generated\\." )
  set( raw_package_list )
  file( STRINGS ${TubeTK_BINARY_DIR}/TubeTK-build/CPackOutputFiles.txt raw_package_list REGEX ${regexp} )
  foreach( package ${raw_package_list} )
    string( REGEX REPLACE ${regexp} "\\1" package_path "${package}" )
    list( APPEND package_list ${package_path} )
  endforeach( package ${raw_package_list} )
  ctest_upload( FILES ${package_list} )
  ctest_submit( PARTS Upload )
endif( SITE_EXPERIMENTAL_UPLOAD )

if( SITE_EXPERIMENTAL_KWSTYLE )
  set( CTEST_BUILD_NAME "${SITE_BUILD_NAME}-Style-Experimental" )
  configure_file( ${TubeTK_SOURCE_DIR}/CMake/InitCMakeCache.cmake.in
    ${TubeTK_BINARY_DIR}/InitCMakeCache.cmake IMMEDIATE @ONLY )
  set( CTEST_NOTES_FILES "${TubeTK_BINARY_DIR}/InitCMakeCache.cmake" )
  ctest_start( "Experimental" )
  ctest_configure( BUILD "${TubeTK_BINARY_DIR}"
    SOURCE "${TubeTK_SOURCE_DIR}"
    OPTIONS "-C${TubeTK_BINARY_DIR}/InitCMakeCache.cmake" )
  ctest_read_custom_files( "${TubeTK_BINARY_DIR}" )
  set( CTEST_BUILD_TARGET "StyleCheck" )
  execute_process( COMMAND ${CMAKE_COMMAND}
    --build ${TubeTK_BINARY_DIR}/TubeTK-build
    --target StyleCheck )
  ctest_submit()
endif( SITE_EXPERIMENTAL_KWSTYLE )

set( CTEST_RUN_CURRENT_SCRIPT 0 )
