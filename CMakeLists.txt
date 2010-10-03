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

project( TubeTK )

#-----------------------------------------------------------------------------
# Require CMake 2.8 for the superbuild functionality
cmake_minimum_required( VERSION 2.8)

#-----------------------------------------------------------------------------
#  Set TubeTK 
set(TubeTK_VERSION_MAJOR "0")
set(TubeTK_VERSION_MINOR "9")
set(TubeTK_VERSION "${TubeTK_VERSION_MAJOR}.${TubeTK_VERSION_MINOR}")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")


#-----------------------------------------------------------------------------
# Dashboard setup
INCLUDE( ${CMAKE_ROOT}/Modules/Dart.cmake )
MARK_AS_ADVANCED( TCL_TCLSH DART_ROOT )
ENABLE_TESTING( )
IF( BUILD_TESTING )
  SET( BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard" )
  MARK_AS_ADVANCED( BUILDNAME )
ENDIF( BUILD_TESTING )
SET( BUILD_WARNING_REPORT_LIMIT -1)
SET( BUILD_ERROR_REPORT_LIMIT -1)
INCLUDE( "CTestConfig.cmake" )
CONFIGURE_FILE( ${TubeTK_SOURCE_DIR}/CTestCustom.cmake.in
                ${TubeTK_BINARY_DIR}/CTestCustom.cmake )


#-----------------------------------------------------------------------------
# Doxygen setup
add_subdirectory( Doxygen )

#-----------------------------------------------------------------------------
# KWStyle setup
OPTION( TubeTK_USE_KWSTYLE 
  "Run KWStyle in order to check for violations of the coding standard." 
  OFF )
MARK_AS_ADVANCED( TubeTK_USE_KWSTYLE )
IF( TubeTK_USE_KWSTYLE )
  FIND_PROGRAM( KWSTYLE_EXECUTABLE
    NAMES KWStyle 
    PATHS
    /usr/local/bin
    )
  MARK_AS_ADVANCED( KWSTYLE_EXECUTABLE )

  OPTION( KWSTYLE_USE_VIM_FORMAT 
    "Set KWStyle to generate errors with a VIM-compatible format." 
    OFF )
  OPTION( KWSTYLE_USE_MSVC_FORMAT 
    "Set KWStyle to generate errors with a VisualStudio-compatible format." 
    OFF )
  OPTION( KWSTYLE_DASHBOARD_SUBMISSION 
    "Set KWStyle arguments to generate report for dashboard submission." 
    ON )

  IF( KWSTYLE_DASHBOARD_SUBMISSION )
    SET(KWSTYLE_ARGUMENTS -lesshtml 
                        -xml ${TubeTK_BINARY_DIR}/KWStyle-TubeTK.kws.xml
                        -o ${TubeTK_BINARY_DIR}/KWStyle-TubeTKOverwrite.txt 
                        -dart ${TubeTK_BINARY_DIR} -1 1
                         -D ${TubeTK_BINARY_DIR}/KWStyle-TubeTKFiles.txt)
  ELSE( KWSTYLE_DASHBOARD_SUBMISSION )
    SET(KWSTYLE_ARGUMENTS -xml ${TubeTK_BINARY_DIR}/KWStyle-TubeTK.kws.xml
                        -v -o ${TubeTK_BINARY_DIR}/KWStyle-TubeTKOverwrite.txt 
                         -D ${TubeTK_BINARY_DIR}/KWStyle-TubeTKFiles.txt)
  ENDIF( KWSTYLE_DASHBOARD_SUBMISSION )


  IF( KWSTYLE_USE_VIM_FORMAT )
    SET( KWSTYLE_ARGUMENTS -vim ${KWSTYLE_ARGUMENTS} )
  ENDIF( KWSTYLE_USE_VIM_FORMAT )

  IF( KWSTYLE_USE_MSVC_FORMAT )
    SET( KWSTYLE_ARGUMENTS -msvc ${KWSTYLE_ARGUMENTS} )
  ENDIF( KWSTYLE_USE_MSVC_FORMAT )

  CONFIGURE_FILE( ${TubeTK_SOURCE_DIR}/KWStyle/TubeTK.kws.xml.in
                  ${TubeTK_BINARY_DIR}/KWStyle-TubeTK.kws.xml )

  CONFIGURE_FILE( ${TubeTK_SOURCE_DIR}/KWStyle/TubeTKFiles.txt.in
                  ${TubeTK_BINARY_DIR}/KWStyle-TubeTKFiles.txt )

  CONFIGURE_FILE( ${TubeTK_SOURCE_DIR}/KWStyle/TubeTKOverwrite.txt.in
                  ${TubeTK_BINARY_DIR}/KWStyle-TubeTKOverwrite.txt )

  ADD_CUSTOM_TARGET( StyleCheck 
    COMMAND ${KWSTYLE_EXECUTABLE} ${KWSTYLE_ARGUMENTS}
    WORKING_DIRECTORY ${TubeTK_BINARY_DIR}
    COMMENT "Coding Style Checker" VERBATIM )
ENDIF( TubeTK_USE_KWSTYLE )

#-----------------------------------------------------------------------------
# MSVC setup
IF( WIN32 )
  ADD_DEFINITIONS( -D_SCL_SECURE_NO_WARNINGS )
  IF( NOT CMAKE_FIRST_PASS_COMPLETED )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /MP" 
         CACHE STRING "Flags used by all build types" FORCE )
  ENDIF( NOT CMAKE_FIRST_PASS_COMPLETED )
  OPTION( DISABLE_MSVC8_DEPRECATED_WARNINGS 
         "Disable Visual Studio 8 deprecated warnings" ON )
  MARK_AS_ADVANCED( FORCE DISABLE_MSVC8_DEPRECATED_WARNINGS )  
  IF( DISABLE_MSVC8_DEPRECATED_WARNINGS )
    ADD_DEFINITIONS( -D_CRT_SECURE_NO_DEPRECATE )
  ENDIF( DISABLE_MSVC8_DEPRECATED_WARNINGS )
ENDIF( WIN32 )

#-----------------------------------------------------------------------------
# Mac setup
IF( APPLE )
  IF( NOT CMAKE_FIRST_PASS_COMPLETED )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-common" 
         CACHE STRING "Flags used by all build types" FORCE )
    SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-common" 
         CACHE STRING "Flags used by all build types" FORCE )
  ENDIF( NOT CMAKE_FIRST_PASS_COMPLETED )
ENDIF( APPLE )

#-----------------------------------------------------------------------------
# Linux setup
IF( UNIX )
  IF( NOT CMAKE_FIRST_PASS_COMPLETED )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing" 
         CACHE STRING "Flags used by all build types" FORCE ) 
    SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-strict-aliasing" 
         CACHE STRING "Flags used by all build types" FORCE )
    IF( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
      SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC" 
           CACHE STRING "Flags used by all build types" FORCE )
      SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC" 
           CACHE STRING "Flags used by all build types" FORCE )
    ENDIF( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
  ENDIF( NOT CMAKE_FIRST_PASS_COMPLETED )
ENDIF( UNIX ) 

#-----------------------------------------------------------------------------
# ITK Setup
option( USE_SYSTEM_ITK 
  "Exclude ITK from SuperBuild and use an existing build instead." OFF )
if( USE_SYSTEM_ITK )
  find_package( ITK REQUIRED )
  include(${ITK_USE_FILE})
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
   MESSAGE(SEND_ERROR 
     "QT_QMAKE_EXECUTABLE must be qmake version 4.6.0 or greater." )
  endif( NOT QT4_FOUND )

  option( TubeTK_USE_QTOPENGL
    "Use QtOpenGL."
    OFF )

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
    include(${VTK_USE_FILE})
  endif( USE_SYSTEM_VTK )
endif( TubeTK_USE_VTK )
  
#-----------------------------------------------------------------------------
# OpenIGTLink Setup
option( TubeTK_USE_OpenIGTLink 
  "Use OpenIGTLink to provide communications with surgical systems." ON )

#-----------------------------------------------------------------------------
# FANN Setup
option( TubeTK_USE_FANN
  "Use fast artificial neural network library (FANN)." OFF )

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
SET( CMAKE_FIRST_PASS_COMPLETED ON 
  CACHE BOOL "Modify the cache vars only on first pass thru cmake file" FORCE )
MARK_AS_ADVANCED( CMAKE_FIRST_PASS_COMPLETED )

#-----------------------------------------------------------------------------
# Superbuild setup
option( TubeTK_USE_SUPERBUILD
  "Build TubeTK and the projects it depends on via SuperBuild.cmake." ON )

if( TubeTK_USE_SUPERBUILD )
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
else( TubeTK_USE_SUPERBUILD )

  find_package( ITK REQUIRED )
  include( ${ITK_USE_FILE} )

  find_package( GenerateCLP REQUIRED )
  include( ${GenerateCLP_USE_FILE} )

  if( TubeTK_USE_OpenIGTLink )
    find_package( OpenIGTLink REQUIRED )
    include( ${OpenIGTLink_USE_FILE} )
  endif( TubeTK_USE_OpenIGTLink )

  INCLUDE( ${CMAKE_ROOT}/Modules/FindPkgConfig.cmake )
  
  if (TubeTK_USE_FANN )
    pkg_check_modules(FANN REQUIRED fann)
  endif()

  if( TubeTK_USE_CTK )
    find_package( CTK REQUIRED )
    include( ${CTK_USE_FILE} )
  endif( TubeTK_USE_CTK )

endif( TubeTK_USE_SUPERBUILD )


#-----------------------------------------------------------------------------
# Misc.
include(${CMAKE_ROOT}/Modules/TestBigEndian.cmake)
test_big_endian(CMAKE_WORDS_BIGENDIAN)

#-----------------------------------------------------------------------------
# Define locations for applications and libraries
#  We use Slicer3 variable names to simplify porting to Slicer3
set(Slicer3_INSTALL_BIN_DIR "bin")
set(Slicer3_INSTALL_LIB_DIR "lib/TubeTK")
set(Slicer3_INSTALL_INCLUDE_DIR "include/TubeTK")
set(Slicer3_INSTALL_SHARE_DIR "share/TubeTK")
set(Slicer3_INSTALL_ITKFACTORIES_DIR "${Slicer3_INSTALL_LIB_DIR}/ITKFactories")

set(Slicer3_INSTALL_PLUGINS_BIN_DIR "${Slicer3_INSTALL_LIB_DIR}/Plugins")
set(Slicer3_INSTALL_PLUGINS_LIB_DIR "${Slicer3_INSTALL_LIB_DIR}/Plugins")
set(Slicer3_INSTALL_PLUGINS_INCLUDE_DIR 
  "${Slicer3_INSTALL_INCLUDE_DIR}/Plugins")
set(Slicer3_INSTALL_PLUGINS_SHARE_DIR "${Slicer3_INSTALL_SHARE_DIR}/Plugins")
set(Slicer3_INSTALL_PLUGINS_CACHE_DIR "${Slicer3_INSTALL_LIB_DIR}/PluginsCache")

set(Slicer3_HOME "${Slicer3_BINARY_DIR}")

include( ${TubeTK_SOURCE_DIR}/CMake/Slicer3PluginsMacros.cmake )

# CLI_SOURCE_DIR is a Slicer3 variable.  Needs to point to the top of
#   the Slicer3 CLI modules.
set( CLI_SOURCE_DIR "${TubeTK_BINARY_DIR}/Slicer3" )
set( CLI_BINARY_DIR "${TubeTK_BINARY_DIR}/Slicer3-Build" )
make_directory( ${CLI_SOURCE_DIR} )
make_directory( ${CLI_BINARY_DIR} )
set( Slicer3_SOURCE_DIR ${CLI_SOURCE_DIR} )

# Define MACROs for downloading Slicer3 applications and libs used by TubeTK
IF( NOT SVNCOMMAND )
  FIND_PACKAGE( Subversion REQUIRED )
  SET( SVNCOMMAND ${Subversion_SVN_EXECUTABLE} )
ENDIF( NOT SVNCOMMAND )

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

set( SUPPORT_2D_IMAGES 1 )
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

