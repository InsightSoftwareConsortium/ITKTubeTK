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
set( TubeTK_VERSION "${TubeTK_VERSION_MAJOR}.${TubeTK_VERSION_MINOR}" )

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
configure_file( ${TubeTK_SOURCE_DIR}/CTestCustom.cmake.in
                ${TubeTK_BINARY_DIR}/CTestCustom.cmake )


#-----------------------------------------------------------------------------
# Doxygen setup
add_subdirectory( Doxygen )

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
                        -xml ${TubeTK_BINARY_DIR}/KWStyle-TubeTK.kws.xml
                        -o ${TubeTK_BINARY_DIR}/KWStyle-TubeTKOverwrite.txt
                        -dart ${TubeTK_BINARY_DIR} -1 1
                         -D ${TubeTK_BINARY_DIR}/KWStyle-TubeTKFiles.txt )
  else( KWSTYLE_DASHBOARD_SUBMISSION )
    set( KWSTYLE_ARGUMENTS -xml ${TubeTK_BINARY_DIR}/KWStyle-TubeTK.kws.xml
                        -v -o ${TubeTK_BINARY_DIR}/KWStyle-TubeTKOverwrite.txt
                         -D ${TubeTK_BINARY_DIR}/KWStyle-TubeTKFiles.txt )
  endif( KWSTYLE_DASHBOARD_SUBMISSION )


  if( KWSTYLE_USE_VIM_FORMAT )
    set( KWSTYLE_ARGUMENTS -vim ${KWSTYLE_ARGUMENTS} )
  endif( KWSTYLE_USE_VIM_FORMAT )

  if( KWSTYLE_USE_MSVC_FORMAT )
    set( KWSTYLE_ARGUMENTS -msvc ${KWSTYLE_ARGUMENTS} )
  endif( KWSTYLE_USE_MSVC_FORMAT )

  configure_file( ${TubeTK_SOURCE_DIR}/KWStyle/TubeTK.kws.xml.in
                  ${TubeTK_BINARY_DIR}/KWStyle-TubeTK.kws.xml )

  configure_file( ${TubeTK_SOURCE_DIR}/KWStyle/TubeTKFiles.txt.in
                  ${TubeTK_BINARY_DIR}/KWStyle-TubeTKFiles.txt )

  configure_file( ${TubeTK_SOURCE_DIR}/KWStyle/TubeTKOverwrite.txt.in
                  ${TubeTK_BINARY_DIR}/KWStyle-TubeTKOverwrite.txt )

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
# QT Setup
option( TubeTK_USE_QT
  "Use QT to compile TubeTK, enables additional applications and libraries."
  ON )

if( TubeTK_USE_QT )
  set( QT_MIN_VERSION "4.6.0" )
  set( QT_OFFICIAL_VERSION "4.6" )
  set( QT_REQUIRED TRUE )
  find_package( Qt4 )
  if( NOT QT4_FOUND )
    message( SEND_ERROR
      "QT_QMAKE_EXECUTABLE must be qmake version 4.6.0 or greater." )
  endif( NOT QT4_FOUND )
endif( TubeTK_USE_QT )

#-----------------------------------------------------------------------------
# VTK Setup
option( TubeTK_USE_VTK
  "Use VTK to compile TubeTK, enables additional applications and libraries."
  ON )
if( TubeTK_USE_VTK )
  #-----------------------------------------------------------------------------
  # VTK Superbuild Setup
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

  include( "${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake" )
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
# Define locations for applications and libraries
#  We use Slicer3 variable names to simplify porting to Slicer3
set( Slicer3_INSTALL_BIN_DIR "bin" )
set( Slicer3_INSTALL_LIB_DIR "lib/TubeTK" )
set( Slicer3_INSTALL_INCLUDE_DIR "include/TubeTK" )
set( Slicer3_INSTALL_SHARE_DIR "share/TubeTK" )
set( Slicer3_INSTALL_ITKFACTORIES_DIR "${Slicer3_INSTALL_LIB_DIR}/ITKFactories" )

set( Slicer3_INSTALL_PLUGINS_BIN_DIR "${Slicer3_INSTALL_LIB_DIR}/Plugins" )
set( Slicer3_INSTALL_PLUGINS_LIB_DIR "${Slicer3_INSTALL_LIB_DIR}/Plugins" )
set( Slicer3_INSTALL_PLUGINS_INCLUDE_DIR
  "${Slicer3_INSTALL_INCLUDE_DIR}/Plugins" )
set( Slicer3_INSTALL_PLUGINS_SHARE_DIR "${Slicer3_INSTALL_SHARE_DIR}/Plugins" )
set( Slicer3_INSTALL_PLUGINS_CACHE_DIR "${Slicer3_INSTALL_LIB_DIR}/PluginsCache" )

set( Slicer3_HOME "${Slicer3_BINARY_DIR}" )

include( ${TubeTK_SOURCE_DIR}/CMake/Slicer3PluginsMacros.cmake )

# CLI_SOURCE_DIR is a Slicer3 variable.  Needs to point to the top of
#   the Slicer3 CLI modules.
set( CLI_SOURCE_DIR "${TubeTK_BINARY_DIR}/Slicer3" )
set( CLI_BINARY_DIR "${TubeTK_BINARY_DIR}/Slicer3-Build" )
make_directory( ${CLI_SOURCE_DIR} )
make_directory( ${CLI_BINARY_DIR} )
set( Slicer3_SOURCE_DIR ${CLI_SOURCE_DIR} )

# Define MACROs for downloading Slicer3 applications and libs used by TubeTK
if( NOT SVNCOMMAND )
  find_package( Subversion REQUIRED )
  set( SVNCOMMAND ${Subversion_SVN_EXECUTABLE} )
endif( NOT SVNCOMMAND )

macro( svnGetSlicerCLI FILENAME )
  set( svnCmd "co" )
  set( svnSite "http://svn.slicer.org/Slicer3/trunk/Applications/CLI" )
  execute_process(
    COMMAND "${SVNCOMMAND}" ${svnCmd} ${svnSite}/${FILENAME} ${CLI_SOURCE_DIR}/${FILENAME}
    OUTPUT_VARIABLE svnOutput
    ERROR_VARIABLE svnError
    RESULT_VARIABLE svnResult
    )
  if( svnError )
    message( "Error while fetching CLI module from Slicer3 SVN..." )
    message( "Svn = ${SVNCOMMAND}" )
    message( "Command = ${svnCmd}" )
    message( "Site = ${svnSite}" )
    message( "File = ${FILENAME}" )
    message( "Dir = ${CLI_SOURCE_DIR}" )
    message( "Output = ${svnOutput}" )
    message( FATAL "Result = ${svnResult}" )
  endif( svnError )
endmacro( svnGetSlicerCLI FILENAME )

macro( svnGetSlicerBase FILENAME )
  set( svnCmd "co" )
  set( svnSite "http://svn.slicer.org/Slicer3/trunk/" )
  execute_process(
    COMMAND "${SVNCOMMAND}" ${svnCmd} ${svnSite}/${FILENAME} ${CLI_SOURCE_DIR}/${FILENAME}
    OUTPUT_VARIABLE svnOutput
    ERROR_VARIABLE svnError
    RESULT_VARIABLE svnResult
    )
  if( svnError )
    message( "Error while fetching module from Slicer3 SVN..." )
    message( "Svn = ${SVNCOMMAND}" )
    message( "Command = ${svnCmd}" )
    message( "Site = ${svnSite}" )
    message( "File = ${FILENAME}" )
    message( "Dir = ${CLI_SOURCE_DIR}" )
    message( "Output = ${svnOutput}" )
    message( FATAL "Result = ${svnResult}" )
  endif( svnError )
endmacro( svnGetSlicerBase FILENAME )

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

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Build TubeTK
add_subdirectory( Utilities )
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
# Configure a launcher for running tubetk methods from the cmdline in windows
configure_file( ${TubeTK_SOURCE_DIR}/CMake/tubeLauncher.bat.in
                ${TubeTK_BINARY_DIR}/tubeLauncher.bat @ONLY )

#-----------------------------------------------------------------------------
# Configure a launcher for running tubetk methods from the cmdline in linux
configure_file( ${TubeTK_SOURCE_DIR}/CMake/tubeLauncher.sh.in
                ${TubeTK_BINARY_DIR}/tubeLauncher.sh @ONLY )

