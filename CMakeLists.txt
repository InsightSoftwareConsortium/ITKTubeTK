##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 ( the "License" );
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

project( TubeTK )

#-----------------------------------------------------------------------------
# Require CMake 2.8 for the superbuild functionality
cmake_minimum_required( VERSION 2.8 )

#-----------------------------------------------------------------------------
#  Set TubeTK
set( TubeTK_VERSION_MAJOR "0" )
set( TubeTK_VERSION_MINOR "9" )
include( CMake/kwsysDateStamp.cmake )
set( TubeTK_VERSION_BUILD
  ${KWSYS_DATE_STAMP_YEAR}${KWSYS_DATE_STAMP_MONTH}${KWSYS_DATE_STAMP_DAY} )
set( TubeTK_VERSION
  "${TubeTK_VERSION_MAJOR}.${TubeTK_VERSION_MINOR}.${TubeTK_VERSION_BUILD}" )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" )


#-----------------------------------------------------------------------------
# Dashboard setup
include( ${CMAKE_ROOT}/Modules/Dart.cmake )
mark_as_advanced( TCL_TCLSH DART_ROOT )
enable_testing()
if( BUILD_TESTING )
  set( BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard" )
  mark_as_advanced( BUILDNAME )
endif( BUILD_TESTING )
set( BUILD_WARNING_REPORT_LIMIT -1 )
set( BUILD_ERROR_REPORT_LIMIT -1 )
include( "CTestConfig.cmake" )
configure_file( ${TubeTK_SOURCE_DIR}/CMake/CTestCustom.cmake.in
  ${TubeTK_BINARY_DIR}/CTestCustom.cmake )


#-----------------------------------------------------------------------------
# Doxygen setup
add_subdirectory( CMake/Doxygen )

#-----------------------------------------------------------------------------
# KWStyle setup
option( TubeTK_USE_KWSTYLE
  "Run KWStyle in order to check for violations of the coding standard."
  OFF )
mark_as_advanced( TubeTK_USE_KWSTYLE )
if( TubeTK_USE_KWSTYLE )
  find_program( KWSTYLE_EXECUTABLE
    NAMES KWStyle
    PATHS /usr/local/bin
    )
  mark_as_advanced( KWSTYLE_EXECUTABLE )

  option( KWSTYLE_USE_VIM_FORMAT
    "Set KWStyle to generate errors with a VIM-compatible format."
    OFF )
  option( KWSTYLE_USE_MSVC_FORMAT
    "Set KWStyle to generate errors with a VisualStudio-compatible format."
    OFF )
  option( KWSTYLE_DASHBOARD_SUBMISSION
    "Set KWStyle arguments to generate report for dashboard submission."
    ON )

  if( KWSTYLE_DASHBOARD_SUBMISSION )
    set( KWSTYLE_ARGUMENTS -lesshtml
      -xml ${TubeTK_BINARY_DIR}/KWStyle.kws.xml
      -o ${TubeTK_BINARY_DIR}/KWStyle.Overwrite.txt
      -dart ${TubeTK_BINARY_DIR} -1 1
      -D ${TubeTK_BINARY_DIR}/KWStyle.Files.txt )
  else( KWSTYLE_DASHBOARD_SUBMISSION )
    set( KWSTYLE_ARGUMENTS -xml ${TubeTK_BINARY_DIR}/KWStyle.kws.xml
      -v -o ${TubeTK_BINARY_DIR}/KWStyle.Overwrite.txt
      -D ${TubeTK_BINARY_DIR}/KWStyle.Files.txt )
  endif( KWSTYLE_DASHBOARD_SUBMISSION )


  if( KWSTYLE_USE_VIM_FORMAT )
    set( KWSTYLE_ARGUMENTS -vim ${KWSTYLE_ARGUMENTS} )
  endif( KWSTYLE_USE_VIM_FORMAT )

  if( KWSTYLE_USE_MSVC_FORMAT )
    set( KWSTYLE_ARGUMENTS -msvc ${KWSTYLE_ARGUMENTS} )
  endif( KWSTYLE_USE_MSVC_FORMAT )

  configure_file( ${TubeTK_SOURCE_DIR}/CMake/KWStyle/KWStyle.kws.xml.in
    ${TubeTK_BINARY_DIR}/KWStyle.kws.xml )

  configure_file( ${TubeTK_SOURCE_DIR}/CMake/KWStyle/KWStyle.Files.txt.in
    ${TubeTK_BINARY_DIR}/KWStyle.Files.txt )

  configure_file( ${TubeTK_SOURCE_DIR}/CMake/KWStyle/KWStyle.Overwrite.txt.in
    ${TubeTK_BINARY_DIR}/KWStyle.Overwrite.txt )

  add_custom_target( StyleCheck
    COMMAND ${KWSTYLE_EXECUTABLE} ${KWSTYLE_ARGUMENTS}
    WORKING_DIRECTORY ${TubeTK_BINARY_DIR}
    COMMENT "Coding Style Checker" VERBATIM )
endif( TubeTK_USE_KWSTYLE )

#-----------------------------------------------------------------------------
# MSVC setup
if( WIN32 )
  add_definitions( -D_SCL_SECURE_NO_WARNINGS )
  if( NOT CMAKE_FIRST_PASS_COMPLETED )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /MP"
         CACHE STRING "Flags used by all build types" FORCE )
  endif( NOT CMAKE_FIRST_PASS_COMPLETED )
  option( DISABLE_MSVC8_DEPRECATED_WARNINGS
         "Disable Visual Studio 8 deprecated warnings" ON )
  mark_as_advanced( FORCE DISABLE_MSVC8_DEPRECATED_WARNINGS )
  if( DISABLE_MSVC8_DEPRECATED_WARNINGS )
    add_definitions( -D_CRT_SECURE_NO_DEPRECATE )
  endif( DISABLE_MSVC8_DEPRECATED_WARNINGS )
endif( WIN32 )

#-----------------------------------------------------------------------------
# Mac setup
if( APPLE )
  if( NOT CMAKE_FIRST_PASS_COMPLETED )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-common"
         CACHE STRING "Flags used by all build types" FORCE )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-common"
         CACHE STRING "Flags used by all build types" FORCE )
  endif( NOT CMAKE_FIRST_PASS_COMPLETED )
endif( APPLE )

#-----------------------------------------------------------------------------
# Linux setup
if( UNIX )
  if( NOT CMAKE_FIRST_PASS_COMPLETED )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing"
         CACHE STRING "Flags used by all build types" FORCE )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-strict-aliasing"
         CACHE STRING "Flags used by all build types" FORCE )
    if( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC"
           CACHE STRING "Flags used by all build types" FORCE )
      set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC"
           CACHE STRING "Flags used by all build types" FORCE )
    endif( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
  endif( NOT CMAKE_FIRST_PASS_COMPLETED )
endif( UNIX )

#-----------------------------------------------------------------------------
# ITK Setup
option( USE_SYSTEM_ITK
  "Exclude ITK from SuperBuild and use an existing build instead." OFF )
if( USE_SYSTEM_ITK )
  find_package( ITK REQUIRED )
  include( ${ITK_USE_FILE} )
endif( USE_SYSTEM_ITK )

#-----------------------------------------------------------------------------
# ITKv4 Switch
option( TubeTK_USE_ITKV4
  "Compile with ITKv4. Only supported applications will be built."
  OFF )
mark_as_advanced( TubeTK_USE_ITKV4 )

#-----------------------------------------------------------------------------
# QT Setup
option( TubeTK_USE_QT
  "Use QT to compile TubeTK, enables additional applications and libraries."
  ON )

if( TubeTK_USE_QT )
  find_package( Qt4 )
  set( minimum_required_qt_version "4.6" )
  if( QT4_FOUND )
    if( "${QT_VERSION_MAJOR}.${QT_VERSION_MINOR}"
        VERSION_LESS "${minimum_required_qt_version}" )
      message( FATAL_ERROR "Error: TubeTK requires Qt >= 4.6" )
    endif()
    include( ${QT_USE_FILE} )
  else()
    message( SEND_ERROR "Error: Qt not found." )
  endif()
endif()

#-----------------------------------------------------------------------------
# VTK Setup
option( TubeTK_USE_VTK
  "Use VTK to compile TubeTK, enables additional applications and libraries."
  ON )
if( TubeTK_USE_VTK )
  #-----------------------------------------------------------------------------
  # VTK Superbuild or System
  option( USE_SYSTEM_VTK
    "Exclude VTK from SuperBuild and use an existing build instead." OFF )
  if( USE_SYSTEM_VTK )
    find_package( VTK REQUIRED )
    include( ${VTK_USE_FILE} )
  endif( USE_SYSTEM_VTK )
endif( TubeTK_USE_VTK )

#-----------------------------------------------------------------------------
# CTK Setup
if( TubeTK_USE_QT )
  option( TubeTK_USE_CTK
    "Use CTK to provide Qt-based applications and build Slicer loadable modules."
    ON )
  #-----------------------------------------------------------------------------
  # CTK Superbuild or System
  option( USE_SYSTEM_CTK
    "Exclude CTK from SuperBuild and use an existing build instead." OFF )
  if( USE_SYSTEM_CTK )
    find_package( CTK REQUIRED )
    include( ${CTK_USE_FILE} )
  endif( USE_SYSTEM_CTK )
endif( TubeTK_USE_QT )

#-----------------------------------------------------------------------------
# Indicate that this cmakelists.txt file has already been parsed once.
# Keeps append commands from having an effect every time cmake is run.
set( CMAKE_FIRST_PASS_COMPLETED ON
  CACHE BOOL "Modify the cache vars only on first pass thru cmake file" FORCE )
mark_as_advanced( CMAKE_FIRST_PASS_COMPLETED )

#-----------------------------------------------------------------------------
# Superbuild setup
option( TubeTK_USE_SUPERBUILD
  "Build TubeTK and the projects it depends on via SuperBuild.cmake." ON )

if( TubeTK_USE_SUPERBUILD )

  include( "${CMAKE_CURRENT_SOURCE_DIR}/CMake/SuperBuild.cmake" )
  return()

else( TubeTK_USE_SUPERBUILD )

  find_package( ITK REQUIRED )
  include( ${ITK_USE_FILE} )

  if( TubeTK_USE_VTK )
    find_package( VTK REQUIRED )
    include( ${VTK_USE_FILE} )
  endif( TubeTK_USE_VTK )

  find_package( GenerateCLP REQUIRED )
  include( ${GenerateCLP_USE_FILE} )

  if( TubeTK_USE_QT )
    if( TubeTK_USE_CTK )
      find_package( CTK REQUIRED )
      include( ${CTK_USE_FILE} )
    endif( TubeTK_USE_CTK )
  endif( TubeTK_USE_QT )

endif( TubeTK_USE_SUPERBUILD )


#-----------------------------------------------------------------------------
# Misc.
include( ${CMAKE_ROOT}/Modules/TestBigEndian.cmake )
test_big_endian( CMAKE_WORDS_BIGENDIAN )

#-----------------------------------------------------------------------------
# Setup to conform to Slicer3 paths and define plugin copying macros
if( NOT TubeTK_USE_ITKV4 )
include( CMake/Slicer3Config.cmake )

#-----------------------------------------------------------------------------
# Build Slicer3 cli applications essential to TubeTK
svnGetSlicerBase( Resources ${CLI_SOURCE_DIR} )
svnGetSlicerCLI( Templates ${CLI_SOURCE_DIR} )
svnGetSlicerCLI( RegisterImagesModule ${CLI_SOURCE_DIR} )

set( CXX_TEST_PATH  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR} )
set( LAUNCH_EXE  "" )

set( SUPPORT_2D_IMAGES ON )
add_subdirectory( ${CLI_SOURCE_DIR}/RegisterImagesModule
  ${CLI_BINARY_DIR}/RegisterImagesModule )
endif( NOT TubeTK_USE_ITKV4 )

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Build TubeTK
if( NOT TubeTK_USE_ITKV4 )
  add_subdirectory( Utilities )
endif( NOT TubeTK_USE_ITKV4 )
add_subdirectory( Base )
add_subdirectory( Applications )

#-----------------------------------------------------------------------------
# Create the directory where tests store temp results
make_directory( ${TubeTK_BINARY_DIR}/Temporary )

install(
  FILES README.txt
  DESTINATION doc
  )

#-----------------------------------------------------------------------------
# Configure vars to export

set( TubeTK_INCLUDE_DIRS
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/IO
  ${TubeTK_SOURCE_DIR}/Base/Segmentation
  ${TubeTK_SOURCE_DIR}/Base/Registration )

if( WIN32 )
 set( TubeTK_EXECUTABLE_DIRS
   ${GenerateCLP_DIR}/${CMAKE_BUILD_TYPE}
   ${ITK_DIR}/bin/${CMAKE_BUILD_TYPE}
   ${VTK_DIR}/bin/${CMAKE_BUILD_TYPE}
   ${TubeTK_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE}
   ${TubeTK_BINARY_DIR}/lib/TubeTK/Plugins/${CMAKE_BUILD_TYPE} )
else()
 set( TubeTK_EXECUTABLE_DIRS
   ${GenerateCLP_DIR}
   ${ITK_DIR}/bin
   ${VTK_DIR}/bin
   ${TubeTK_BINARY_DIR}/bin
   ${TubeTK_BINARY_DIR}/lib/TubeTK/Plugins )
endif()

set( TubeTK_LIBRARY_DIRS
  ${TubeTK_EXECUTABLE_DIRS} )

set( TubeTK_LIBRARIES
  TubeIO
  TubeNumerics )

set( TubeTK_CMAKE_EXTENSIONS_DIR
  ${TubeTK_SOURCE_DIR}/CMake )

set( TubeTK_USE_FILE
  ${TubeTK_BINARY_DIR}/UseTubeTK.cmake )

#-----------------------------------------------------------------------------
# Configure TubeTKConfig.cmake
configure_file( ${TubeTK_SOURCE_DIR}/CMake/TubeTKConfig.cmake.in
  ${TubeTK_BINARY_DIR}/TubeTKConfig.cmake @ONLY )

#-----------------------------------------------------------------------------
# Configure UseTubeTK.cmake
configure_file( ${TubeTK_SOURCE_DIR}/CMake/UseTubeTK.cmake.in
  ${TubeTK_BINARY_DIR}/UseTubeTK.cmake @ONLY )

#-----------------------------------------------------------------------------
# Configure a launcher for running tubetk methods from the cmdline in windows
configure_file( ${TubeTK_SOURCE_DIR}/CMake/tubeLauncher.bat.in
  ${TubeTK_BINARY_DIR}/tubeLauncher.bat @ONLY )

#-----------------------------------------------------------------------------
# Configure a launcher for running tubetk methods from the cmdline in linux
configure_file( ${TubeTK_SOURCE_DIR}/CMake/tubeLauncher.sh.in
  ${TubeTK_BINARY_DIR}/tubeLauncher.sh @ONLY )
