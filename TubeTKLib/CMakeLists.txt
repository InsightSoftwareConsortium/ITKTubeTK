##############################################################################
#
# Library:   TubeTKLib
#
# Copyright Kitware Inc.
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
cmake_minimum_required(VERSION 3.10.2)

####
####
####
####
project( TubeTKLib )

list( APPEND CMAKE_MODULE_PATH ${TubeTKLib_SOURCE_DIR}/CMake )

####
# Default output directories
####
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" )

set( CMAKE_POSITION_INDEPENDENT_CODE ON )

####
# GPU ArrayFire setup
####
option( TubeTKLib_USE_ARRAYFIRE
  "Use the ArrayFire library to speedup filtering opertions using the GPU."
  OFF)

####
# Cppcheck setup.
####
option( TubeTKLib_USE_CPPCHECK "Use Cppcheck to enable additional tests."
  OFF )
mark_as_advanced( TubeTKLib_USE_CPPCHECK )

####
# KWStyle setup.
####
option( TubeTKLib_USE_KWSTYLE
  "Use KWStyle to check for coding standards violations." OFF )
mark_as_advanced( TubeTKLib_USE_KWSTYLE )

####
# LIBSVM setup.
####
option( TubeTKLib_USE_LIBSVM "Use LIBSVM for additional functionality."
  OFF )
mark_as_advanced( TubeTKLib_USE_LIBSVM )
mark_as_advanced( LIBSVM_DIR )

####
# RandomForest setup.
####
option( TubeTKLib_USE_RANDOMFOREST
  "Use RandomForest to enable additional functionality." OFF )
mark_as_advanced( TubeTKLib_USE_RANDOMFOREST )
mark_as_advanced( RandomForest_DIR )

####
# Provide option for Valgrind
####
option( TubeTKLib_USE_VALGRIND
  "Use Valgrind to enable additional functionality." OFF )
if( WIN32 )
  mark_as_advanced( TubeTKLib_USE_VALGRIND )
endif( WIN32 )

####
# ITK Setup.
####
if( NOT ITK_DIR )
  set( ITK_DIR "NOT_FOUND" CACHE PATH
	  "Path to ITK installation/build directory." FORCE )
endif( NOT ITK_DIR )

####
# VTK setup.
####
option( TubeTKLib_USE_VTK "Use VTK to enable additional functionality." ON )
if( TubeTKLib_USE_VTK )
  if( NOT VTK_DIR )
    set( VTK_DIR "NOT_FOUND" CACHE PATH
	    "Path to ITK installation/build directory." FORCE )
  endif( NOT VTK_DIR )
endif( TubeTKLib_USE_VTK )

####
#  Configure Project
####
set( TubeTKLib_CONFIG_BINARY_DIR ${TubeTKLib_BINARY_DIR} )
mark_as_advanced( TubeTKLib_CONFIG_BINARY_DIR )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )
mark_as_advanced( CLEAR ITK_DIR )
set( ITK_HAS_MINIMAL_PATH_EXTRACTION FALSE )
foreach( itk_module in LISTS ITK_MODULE_ENABLED )
  if( ${itk_module} MATCHES "MinimalPathExtraction" )
    set( ITK_HAS_MINIMAL_PATH_EXTRACTION TRUE )
  endif()
endforeach()

option( TubeTKLib_SUPERBUILD "Build TubeTKLib's subprojects" ON )
mark_as_advanced( TubeTKLib_SUPERBUILD )
if( NOT ITK_HAS_MINIMAL_PATH_EXTRACTION )
  if( TubeTKLib_SUPERBUILD )
    add_subdirectory( ${TubeTKLib_SOURCE_DIR}/SuperBuild )
    return()
  endif()
endif()

if( TubeTKLib_USE_VTK )
  find_package( VTK REQUIRED )
  include( ${VTK_USE_FILE} )
  mark_as_advanced( CLEAR VTK_DIR )
endif( TubeTKLib_USE_VTK )

if( TubeTKLib_USE_LIBSVM )
  find_package( LIBSVM REQUIRED )
endif( TubeTKLib_USE_LIBSVM )

if( TubeTKLib_USE_RANDOMFOREST )
  find_package( RandomForest REQUIRED )
endif( TubeTKLib_USE_RANDOMFOREST )

if( TubeTKLib_USE_ARRAYFIRE )
  find_package( ArrayFire REQUIRED )
  include_directories( ${ArrayFire_INCLUDE_DIRS} )
endif( TubeTKLib_USE_ARRAYFIRE )

if( TubeTKLib_USE_KWSTYLE )
  include( TubeTKLibKWStyleConfig )
endif( TubeTKLib_USE_KWSTYLE )

if( TubeTKLib_USE_CPPCHECK )
  list( APPEND CMAKE_MODULE_PATH "${TubeTKLib_SOURCE_DIR}/CMake/Cppcheck" )
  find_package( Cppcheck REQUIRED )
  set( CPPCHECK_PATHS_TO_CHECK
    "${TubeTKLib_SOURCE_DIR}" )
  set( CPPCHECK_SUPPRESSIONS_FILE
    "${TubeTKLib_SOURCE_DIR}/CMake/Cppcheck/TubeTKLibCppcheckSuppressions.txt"
    )
  set( CPPCHECK_COMMAND_PREFIX "${CPPCHECK_EXECUTABLE}" -q --xml
    --suppressions-list=${CPPCHECK_SUPPRESSIONS_FILE} )
  add_test( Cppcheck
    ${CPPCHECK_COMMAND_PREFIX} --enable=all ${CPPCHECK_PATHS_TO_CHECK} )
  set_tests_properties( Cppcheck
    PROPERTIES FAIL_REGULAR_EXPRESSION "error file=" )
endif( TubeTKLib_USE_CPPCHECK )

configure_file( "${TubeTKLib_SOURCE_DIR}/CMake/tubetklibConfigure.h.in"
  tubetklibConfigure.h @ONLY )

include_directories( ${TubeTKLib_BINARY_DIR} )

####
# Setup Launcher for Apps
####
if( WIN32 )
  set( TubeTKLib_EXECUTABLE_DIRS CACHE INTERNAL
    "Bin and Lib dirs for running apps." FORCE )
  foreach( _build_type "" Debug Release )
    list( APPEND TubeTKLib_EXECUTABLE_DIRS
      ${ITK_DIR}/bin/${_build_type}
      ${TubeTKLib_BINARY_DIR}/bin/${_build_type}
      ${TubeTKLib_BINARY_DIR}/lib/TubeTKLib/${_build_type} )
    if( TubeTKLib_USE_VTK )
      list( APPEND TubeTKLib_EXECUTABLE_DIRS ${VTK_DIR}/bin/${_build_type} )
    endif( TubeTKLib_USE_VTK )
  endforeach()
else( WIN32 )
  set( TubeTKLib_EXECUTABLE_DIRS
    ${ITK_DIR}/bin
    ${TubeTKLib_BINARY_DIR}/bin
    ${TubeTKLib_BINARY_DIR}/lib/TubeTKLib
    CACHE INTERNAL "Bin and Lib dirs for running apps." FORCE )
  if( TubeTKLib_USE_VTK )
    list( APPEND TubeTKLib_EXECUTABLE_DIRS ${VTK_DIR}/bin )
  endif( TubeTKLib_USE_VTK )
endif( WIN32 )

message( STATUS "Configuring Launcher script" )
set( TubeTKLib_LAUNCHER )
if( WIN32 )
  set( _launcher_platform "windows" )
  configure_file( ${TubeTKLib_SOURCE_DIR}/CMake/TubeTKLibLauncher.bat.in
    ${TubeTKLib_BINARY_DIR}/TubeTKLibLauncher.bat @ONLY )

  find_program( CMD_EXECUTABLE "cmd" )
  if( NOT CMD_EXECUTABLE )
    message( FATAL_ERROR
      "Could not find 'cmd' executable required to run tests and launcher" )
  endif()
  mark_as_advanced( CMD_EXECUTABLE )

  set( TubeTKLib_LAUNCHER ${TubeTKLib_BINARY_DIR}/TubeTKLibLauncher.bat )
elseif( UNIX )
  set( _launcher_platform "unix" )
  configure_file( ${TubeTKLib_SOURCE_DIR}/CMake/TubeTKLibLauncher.sh.in
    ${TubeTKLib_BINARY_DIR}/TubeTKLauncher.sh @ONLY )

  find_program( SH_EXECUTABLE "sh" )
  if( NOT SH_EXECUTABLE )
    message( FATAL_ERROR
      "Could not find 'sh' executable required to test using the launcher" )
  endif()
  set( TubeTKLib_LAUNCHER ${SH_EXECUTABLE}
    ${TubeTKLib_BINARY_DIR}/TubeTKLibLauncher.sh )
elseif( NOT UNIX )
  message( FATAL_ERROR
    "Configuring Launcher script - failed [unknown platform]" )
endif()
message( STATUS "Configuring Launcher script - ok [${_launcher_platform}]" )

add_subdirectory( Utilities )
add_subdirectory( Common )
add_subdirectory( CompareTools )
add_subdirectory( MetaIO )
add_subdirectory( Numerics )
add_subdirectory( ObjectDocuments )
add_subdirectory( Filtering )
add_subdirectory( Registration )
add_subdirectory( Segmentation )
add_subdirectory( IO )
add_subdirectory( USTK )


####
# Configure for building against a build tree
####
set( TubeTKLib_INCLUDE_DIRS
  ${TubeTKLib_BINARY_DIR}
  ${TubeTKLib_SOURCE_DIR}/Common
  ${TubeTKLib_SOURCE_DIR}/Filtering
  ${TubeTKLib_SOURCE_DIR}/IO
  ${TubeTKLib_SOURCE_DIR}/MetaIO
  ${TubeTKLib_SOURCE_DIR}/Numerics
  ${TubeTKLib_SOURCE_DIR}/ObjectDocuments
  ${TubeTKLib_SOURCE_DIR}/Registration
  ${TubeTKLib_SOURCE_DIR}/Segmentation
  ${TubeTKLib_SOURCE_DIR}/USTK )

set( TubeTKLib_LIBRARY_DIRS ${TubeTKLib_EXECUTABLE_DIRS} )

get_property( TubeTKLib_TARGETS GLOBAL PROPERTY TubeTKLib_TARGETS )
export( TARGETS ${TubeTKLib_TARGETS} FILE
  ${TubeTKLib_CONFIG_BINARY_DIR}/TubeTKLibTargets.cmake )

set( TubeTKLib_TARGETS_FILE
  ${TubeTKLib_CONFIG_BINARY_DIR}/TubeTKLibTargets.cmake )

set( TubeTKLib_USE_FILE
  ${TubeTKLib_CONFIG_BINARY_DIR}/UseTubeTKLib.cmake )

configure_file( ${TubeTKLib_SOURCE_DIR}/CMake/TubeTKLibConfig.cmake.in
  ${TubeTKLib_CONFIG_BINARY_DIR}/TubeTKLibConfig.cmake @ONLY )

configure_file( ${TubeTKLib_SOURCE_DIR}/CMake/UseTubeTKLib.cmake.in
  ${TubeTKLib_CONFIG_BINARY_DIR}/UseTubeTKLib.cmake @ONLY )

####
# Package
####
include( TubeTKLibCPackConfig )

####
# Enable access to remote testing data (must be placed last)
####
include( TubeTKLibExternalData )
if( DEFINED "ENV{ExternalData_OBJECT_STORES}" )
  set( ExternalData_OBJECT_STORES $ENV{ExternalData_OBJECT_STORES} )
endif()
include( ${CMAKE_ROOT}/Modules/ExternalData.cmake )
ExternalData_Add_Target( TubeTKLibData )
