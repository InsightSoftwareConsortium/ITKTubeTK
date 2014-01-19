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

# Require CMake 2.8 for its superbuild functionality.
cmake_minimum_required( VERSION 2.8.7 )
if( NOT ${CMAKE_VERSION} VERSION_LESS 2.8.13 )
  cmake_policy( SET CMP0026 OLD )
endif( NOT ${CMAKE_VERSION} VERSION_LESS 2.8.13 )

# Update CMake module path
set( TubeTK_CMAKE_EXTENSIONS_DIR ${TubeTK_SOURCE_DIR}/CMake )
set( CMAKE_MODULE_PATH ${TubeTK_CMAKE_EXTENSIONS_DIR} ${CMAKE_MODULE_PATH} )

# Set TubeTK version and output directories.
set( TubeTK_VERSION_MAJOR "0" )
set( TubeTK_VERSION_MINOR "9" )
set( TubeTK_VERSION_PATCH "0" )

include( TubeTKFunctionToday )
TubeTKFunctionToday( TubeTK_VERSION_BUILD )

set( TubeTK_VERSION
  "${TubeTK_VERSION_MAJOR}.${TubeTK_VERSION_MINOR}.${TubeTK_VERSION_PATCH}.${TubeTK_VERSION_BUILD}" )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/TubeTK" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/TubeTK" )

include( TubeTKBlockSetCMakeOSXVariables )
include( TubeTKMacroEmptyExternalProject )

# Setup for building 3D Slicer extensions.

option( TubeTK_BUILD_SLICER_EXTENSION
  "Build TubeTK as a 3D Slicer extension." OFF )
if( TubeTK_BUILD_SLICER_EXTENSION )
  set( EXTENSION_NAME TubeTK )
  set( EXTENSION_HOMEPAGE "http://www.tubetk.org" )
  set( EXTENSION_CATEGORY "TubeTK" )
  set( EXTENSION_CONTRIBUTORS "TubeTK Team at Kitware, Inc." )
  set( EXTENSION_DESCRIPTION
    "Using local image geometry for segmentation, registration, and filtering" )
  set( EXTENSION_ICONURL
    "http://public.kitware.com/Wiki/images/a/ad/TubeTK_Icon.jpg" )
  set( EXTENSION_SCREENSHOTURLS
    "http://public.kitware.com/Wiki/images/f/fd/VesselExtraction.jpg" )
  set( EXTENSION_README_FILE ${TubeTK_SOURCE_DIR}/README.md )
  set( EXTENSION_LICENSE_FILE ${TubeTK_SOURCE_DIR}/LICENSE.md )
  set( EXTENSION_BUILD_SUBDIRECTORY "TubeTK-build")
endif( TubeTK_BUILD_SLICER_EXTENSION )

# CDash dashboard setup.
include( CTest )
if( BUILD_TESTING )
  set( BUILDNAME "${BUILDNAME}" CACHE STRING
    "Name of the build on the CDash dashboard." )
  mark_as_advanced( BUILDNAME )
endif( BUILD_TESTING )
set( BUILD_WARNING_REPORT_LIMIT -1 )
set( BUILD_ERROR_REPORT_LIMIT -1 )
include( CTestConfig.cmake )
configure_file( ${TubeTK_SOURCE_DIR}/CMake/CTestCustom.cmake.in
  ${TubeTK_BINARY_DIR}/CTestCustom.cmake )

# Configure for the download of CTest unit testing data from a Midas Platform
# server.
set( MIDAS_USER "" CACHE STRING "Midas Platform user email address." )
mark_as_advanced( MIDAS_USER )
set( MIDAS_DEFAULT_API_KEY "" CACHE STRING "Midas Platform default API key." )
mark_as_advanced( MIDAS_DEFAULT_API_KEY )
set( MIDAS_DATA_DIR ${TubeTK_BINARY_DIR}/MIDAS_Data )

# Doxygen documentation setup.
option( BUILD_DOCUMENTATION
  "Build the software guide and source code documentation." OFF )
if( BUILD_DOCUMENTATION )
  option( TubeTK_USE_DOXYGEN
    "Use Doxygen to generate the source code documentation." ON )
  mark_as_advanced( TubeTK_USE_DOXYGEN )
  if( TubeTK_USE_DOXYGEN )
    add_subdirectory( CMake/Doxygen )
  endif( TubeTK_USE_DOXYGEN )
endif( BUILD_DOCUMENTATION )

# Image viewer setup.
option( TubeTK_BUILD_IMAGE_VIEWER "Use ImageViewer." OFF )

# Microsoft Visual Studio setup.
if( MSVC )
  add_definitions( -D_SCL_SECURE_NO_WARNINGS )
  if( NOT CMAKE_FIRST_PASS_COMPLETED AND "${CMAKE_CXX_COMPILER}" MATCHES "cl" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /MP"
         CACHE STRING "Flags used by all build types." FORCE )
  endif( NOT CMAKE_FIRST_PASS_COMPLETED AND "${CMAKE_CXX_COMPILER}" MATCHES "cl" )
  option( DISABLE_MSVC8_DEPRECATED_WARNINGS
    "Disable Microsoft Visual Studio 8 deprecated warnings." ON )
  mark_as_advanced( FORCE DISABLE_MSVC8_DEPRECATED_WARNINGS )
  if( DISABLE_MSVC8_DEPRECATED_WARNINGS )
    add_definitions( -D_CRT_SECURE_NO_DEPRECATE )
  endif( DISABLE_MSVC8_DEPRECATED_WARNINGS )
endif( MSVC )

# Apple OS X setup.
if( APPLE )
  if( NOT CMAKE_FIRST_PASS_COMPLETED )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-common"
      CACHE STRING "Flags used by all build types." FORCE )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-common"
      CACHE STRING "Flags used by all build types." FORCE )
  endif( NOT CMAKE_FIRST_PASS_COMPLETED )
endif( APPLE )

# Linux or UNIX setup.
if( UNIX )
  if( NOT CMAKE_FIRST_PASS_COMPLETED )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing"
      CACHE STRING "Flags used by all build types." FORCE )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-strict-aliasing"
      CACHE STRING "Flags used by all build types." FORCE )
    if( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC"
           CACHE STRING "Flags used by all build types." FORCE )
      set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC"
           CACHE STRING "Flags used by all build types." FORCE )
    endif( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
  endif( NOT CMAKE_FIRST_PASS_COMPLETED )
endif( UNIX )

# Implementation agnostic fPIC
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

# Slicer library setup.
# If available, Slicer can be used to retrieve external libraries and modules
# used by TubeTK.
option( TubeTK_USE_SLICER "Use Slicer to extract useful commands and CLIs."
  OFF )

# Slicer must be used if TubeTK is built as an extension.
if( TubeTK_BUILD_SLICER_EXTENSION )
  set( TubeTK_USE_SLICER ON CACHE BOOL
    "Use 3D Slicer to extract useful commands and CLIs." FORCE )
  mark_as_advanced( TubeTK_USE_SLICER )
else( TubeTK_BUILD_SLICER_EXTENSION )
  mark_as_advanced( CLEAR TubeTK_USE_SLICER )
endif( TubeTK_BUILD_SLICER_EXTENSION )

if( TubeTK_USE_SLICER )
  mark_as_advanced( TubeTK_USE_ITK TubeTK_USE_VTK TubeTK_USE_QT
    TubeTK_USE_CTK )

  unset( QT_QMAKE_EXECUTABLE )
  # Find Slicer to set the [ITK|VTK|CTK|...]_DIR variables
  find_package( Slicer REQUIRED )
  if( NOT Slicer_FOUND )
    set( Slicer_DIR "NOT_FOUND" CACHE PATH
      "Path to the Slicer-<buildtype>/Slicer-build directory." FORCE )
    mark_as_advanced( CLEAR Slicer_DIR )
  endif( NOT Slicer_FOUND )
  include( ${Slicer_USE_FILE} )
  # TubeTK_SLICER_SOURCE_DIR is needed to get ExpertAutomatedRegistration CLI
  if( Slicer_SOURCE_DIR )
    set( TubeTK_SLICER_SOURCE_DIR ${Slicer_SOURCE_DIR} )
  else( Slicer_SOURCE_DIR )
    set( TubeTK_SLICER_SOURCE_DIR ${Slicer_CMAKE_DIR}/.. )
  endif( Slicer_SOURCE_DIR )
  # All those dependencies are already provided by Slicer
  set( USE_SYSTEM_ITK ON CACHE BOOL
    "Exclude ITK from superbuild and use an existing build." FORCE )
  set( TubeTK_USE_VTK ON CACHE BOOL
    "Use VTK to enable additional functionality." FORCE )
  set( USE_SYSTEM_VTK ON CACHE BOOL
    "Exclude VTK from superbuild and use an existing build." FORCE )
  set( TubeTK_USE_QT ON CACHE BOOL
    "Use the Qt framework to enable additional functionality." FORCE )
  set( TubeTK_USE_CTK ON CACHE BOOL
    "Use CTK to enable additional functionality." FORCE )
  set( USE_SYSTEM_CTK ON CACHE BOOL
    "Exclude CTK from superbuild and use an existing build." FORCE )
  set( USE_SYSTEM_SLICER_EXECUTION_MODEL ON CACHE BOOL
    "Exclude SlicerExecutionModel from superbuild and use an existing build."
    FORCE )
  mark_as_advanced( CTK_DIR DCMTK_DIR ITK_DIR LibArchive_INCLUDE_DIR
    LibArchive_LIBRARY OpenIGTLink_DIR QT_QMAKE_EXECUTABLE QtTesting_DIR
    SlicerExecutionModel_DIR Slicer_EXTENSIONS_TRACK_QUALIFIER Teem_DIR
    VTK_DIR ZLIB_INCLUDE_DIR ZLIB_LIBRARY ZLIB_ROOT jqPlot_DIR qRestAPI_DIR
    MIDAS_PACKAGE_API_KEY MIDAS_PACKAGE_EMAIL MIDAS_PACKAGE_URL )
else( TubeTK_USE_SLICER )
  mark_as_advanced( CLEAR TubeTK_USE_ITK TubeTK_USE_VTK TubeTK_USE_QT )
  if( TubeTK_USE_QT )
    mark_as_advanced( CLEAR TubeTK_USE_CTK )
  endif( TubeTK_USE_QT )
  mark_as_advanced( Slicer_DIR )
endif( TubeTK_USE_SLICER )

# Boost libraries setup.
option( TubeTK_USE_BOOST
  "Use the Boost libraries to enable additional functionality." OFF )
if( TubeTK_USE_BOOST )
  find_package( Boost REQUIRED )
  if( NOT Boost_FOUND )
    set( BOOST_ROOT "NOT_FOUND" CACHE PATH
      "Path to the top-level directory of Boost." FORCE )
    mark_as_advanced( CLEAR BOOST_ROOT )
  endif( NOT Boost_FOUND )
endif( TubeTK_USE_BOOST )

#  Cppcheck setup.
set( TubeTK_USE_CPPCHECK_DEFAULT ON )
set( USE_SYSTEM_CPPCHECK_DEFAULT OFF )
if( WIN32 )
  set( TubeTK_USE_CPPCHECK_DEFAULT OFF )
  set( USE_SYSTEM_CPPCHECK_DEFAULT OFF )
endif( WIN32 )
option( TubeTK_USE_CPPCHECK "Use Cppcheck to enable additional tests."
  ${TubeTK_USE_CPPCHECK_DEFAULT} )
option( USE_SYSTEM_CPPCHECK
  "Exclude Cppcheck from superbuild and use an existing build."
  ${USE_SYSTEM_CPPCHECK_DEFAULT} )
mark_as_advanced( USE_SYSTEM_CPPCHECK )

set( TubeTK_USE_PYTHON_DEFAULT ON )
if( WIN32 )
  set( TubeTK_USE_PYTHON_DEFAULT OFF )
endif( WIN32 )
option( TubeTK_USE_PYTHON "Use Python to enable additional functionality."
  ${TubeTK_USE_PYTHON_DEFAULT} )
if( TubeTK_USE_PYTHON )
  # NumPy conversion support.
  option( TubeTK_USE_NUMPY "Use NumPy to enable additional functionality." ON )
  mark_as_advanced( CLEAR TubeTK_USE_NUMPY )

  # Test TubeTK examples via Python virtualenv.
  option( TubeTK_USE_NOTEBOOKS
    "Test IPython notebooks in a Python virtualenv." ON )
  mark_as_advanced( CLEAR TubeTK_USE_NOTEBOOKS )

  # PyQtGraph bridge support.
  option( TubeTK_USE_PYQTGRAPH
    "Use PyQtGraph to enable additional functionality." ON )
  mark_as_advanced( CLEAR TubeTK_USE_PYQTGRAPH )
else( TubeTK_USE_PYTHON )
  # NumPy conversion support.
  set( TubeTK_USE_NUMPY OFF CACHE BOOL
    "Use NumPy to enable additional functionality." FORCE )
  mark_as_advanced( TubeTK_USE_NUMPY )

  # Test TubeTK examples via Python virtualenv.
  set( TubeTK_USE_NOTEBOOKS OFF CACHE BOOL
    "Test IPython notebooks in a Python virtualenv." FORCE )
  mark_as_advanced( TubeTK_USE_NOTEBOOKS )

  # PyQtGraph bridge support.
  set( TubeTK_USE_PYQTGRAPH OFF CACHE BOOL
    "Use PyQtGraph to enable additional functionality." FORCE )
  mark_as_advanced( TubeTK_USE_PYQTGRAPH )
endif( TubeTK_USE_PYTHON )

# ITK setup.
option( USE_SYSTEM_ITK "Exclude ITK from superbuild and use an existing build."
  OFF )
mark_as_advanced( USE_SYSTEM_ITK )

# JsonCpp library setup.
option( USE_SYSTEM_JSONCPP
  "Exclude JsonCpp from superbuild and use an existing build." OFF )
mark_as_advanced( USE_SYSTEM_JSONCPP )

# KWStyle setup.
option( TubeTK_USE_KWSTYLE
  "Use KWStyle to check for coding standards violations." ON )
option( USE_SYSTEM_KWSTYLE
  "Exclude KWStyle from superbuild and use an existing build." OFF )
mark_as_advanced( USE_SYSTEM_KWSTYLE )

# LIBSVM setup.
option( TubeTK_USE_LIBSVM "Use LIBSVM to enable additional functionality." ON )
option( USE_SYSTEM_LIBSVM
  "Exclude LIBSVM from superbuild and use an existing build." OFF )
mark_as_advanced( USE_SYSTEM_LIBSVM )

# Parameter serializer setup.
option( USE_SYSTEM_SLICER_EXECUTION_MODEL
  "Exclude slicer execution model from superbuild and use an existing build."
  OFF )
mark_as_advanced( USE_SYSTEM_SLICER_EXECUTION_MODEL )
option( USE_SYSTEM_PARAMETER_SERIALIZER
  "Exclude parameter serializer from superbuild and use an existing build."
  OFF )
mark_as_advanced( USE_SYSTEM_PARAMETER_SERIALIZER )

# Qt framework setup.
option( TubeTK_USE_QT
  "Use the Qt framework to enable additional applications and libraries." ON )
if( TubeTK_USE_QT )
  # CTK setup.
  option( TubeTK_USE_CTK "Use CTK to enable additional functionality." ON )
  if( NOT TubeTK_USE_SLICER )
    mark_as_advanced( CLEAR TubeTK_USE_CTK )
  endif( NOT TubeTK_USE_SLICER )
  option( USE_SYSTEM_CTK
    "Exclude CTK from superbuild and use an existing build." OFF )
  if( NOT TubeTK_USE_SLICER )
    mark_as_advanced( USE_SYSTEM_CTK )
  endif( NOT TubeTK_USE_SLICER )

  find_package( Qt4 )
  # 3D Slicer only tests against a specific Qt 4 version. To
  # override the version check, set the TubeTK_REQUIRED_QT_VERSION CMake
  # variable.
  if( NOT DEFINED TubeTK_REQUIRED_QT_VERSION )
    set( TubeTK_REQUIRED_QT_VERSION "4.7.4" CACHE STRING
      "Exact Qt4 framework version required." FORCE )
    mark_as_advanced( TubeTK_REQUIRED_QT_VERSION )
  endif( NOT DEFINED TubeTK_REQUIRED_QT_VERSION )
  if( QT4_FOUND )
    if( "${QT_VERSION_MAJOR}.${QT_VERSION_MINOR}.${QT_VERSION_PATCH}"
        VERSION_LESS "${TubeTK_REQUIRED_QT_VERSION}" )
      message(FATAL_ERROR "error: TubeTK requires at least Qt ${TubeTK_REQUIRED_QT_VERSION} "
                          "-- you cannot use Qt ${QT_VERSION_MAJOR}.${QT_VERSION_MINOR}.${QT_VERSION_PATCH}. ${extra_error_message}")
    endif( "${QT_VERSION_MAJOR}.${QT_VERSION_MINOR}.${QT_VERSION_PATCH}"
      VERSION_LESS "${TubeTK_REQUIRED_QT_VERSION}" )
    include( ${QT_USE_FILE} )
  else( QT4_FOUND )
    message( SEND_ERROR "Error: Qt4 framework not found." )
  endif( QT4_FOUND )
else( TubeTK_USE_QT )
  option( TubeTK_USE_CTK "Use CTK to enable additional functionality." OFF )
  mark_as_advanced( TubeTK_USE_CTK )
  option( USE_SYSTEM_CTK
    "Exclude CTK from superbuild and use an existing build." OFF )
  mark_as_advanced( USE_SYSTEM_CTK )
endif( TubeTK_USE_QT )

# Slicer Execution model Setup
option( USE_SYSTEM_SLICER_EXECUTION_MODEL
  "Exclude SlicerExecutionModel from superbuild and use an existing build."
  OFF )
mark_as_advanced( USE_SYSTEM_SLICER_EXECUTION_MODEL )

option( TubeTK_USE_VALGRIND "Use Valgrind to enable additional functionality."
  OFF )

# VTK setup.
option( TubeTK_USE_VTK "Use VTK to enable additional functionality." ON )
option( USE_SYSTEM_VTK "Exclude VTK from superbuild and use an existing build."
  OFF )
mark_as_advanced( USE_SYSTEM_VTK )

# Modules
option( TubeTK_BUILD_MODULES
  "Build modules. If not, only TubeTK library is built" ON )

# Indicate that this CMakeLists.txt file has already been parsed once. Keeps
# append commands from having an effect every time CMake is run.
set( CMAKE_FIRST_PASS_COMPLETED ON CACHE BOOL
  "Modify the cache variables only on the first pass through the CMake file."
  FORCE )
mark_as_advanced( CMAKE_FIRST_PASS_COMPLETED )


option( TubeTK_USE_SUPERBUILD
  "Build TubeTK and the projects it depends on using a superbuild." ON )
mark_as_advanced( TubeTK_USE_SUPERBUILD )

# If TubeTK_CONFIG_BINARY_DIR isn't defined, it means TubeTK is *NOT* build
# using Superbuild. In that specific case, TubeTK_CONFIG_BINARY_DIR should
# default to TubeTK_BINARY_DIR
if( NOT DEFINED TubeTK_CONFIG_BINARY_DIR )
  set( TubeTK_CONFIG_BINARY_DIR ${CMAKE_BINARY_DIR} )
endif( NOT DEFINED TubeTK_CONFIG_BINARY_DIR )

if( TubeTK_USE_SUPERBUILD )
  if( TubeTK_BUILD_MODULES )
    # Expose the TubeTK_BUILD_* variables to Superbuild
    file( GLOB_RECURSE TubeTK_modules_files "TubeTKModules.cmake" )
    foreach( f ${TubeTK_modules_files} )
      include( ${f} )
    endforeach( f ${TubeTK_modules_files} )
  endif( TubeTK_BUILD_MODULES )

  add_subdirectory( "${CMAKE_CURRENT_SOURCE_DIR}/CMake/Superbuild" )
  return()
endif( TubeTK_USE_SUPERBUILD )

set( TubeTK_INSTALL_BIN_DIR "bin" )
set( TubeTK_INSTALL_LIB_DIR "lib/TubeTK" )
set( TubeTK_INSTALL_INCLUDE_DIR "include/TubeTK" )
set( TubeTK_INSTALL_SHARE_DIR "share/TubeTK" )

set( TubeTK_INSTALL_PLUGINS_BIN_DIR "${TubeTK_INSTALL_BIN_DIR}" )
set( TubeTK_INSTALL_PLUGINS_CACHE_DIR "${TubeTK_INSTALL_LIB_DIR}/CACHE" )
set( TubeTK_INSTALL_PLUGINS_INCLUDE_DIR "${TubeTK_INSTALL_INCLUDE_DIR}" )
set( TubeTK_INSTALL_PLUGINS_LIB_DIR "${TubeTK_INSTALL_LIB_DIR}" )
set( TubeTK_INSTALL_PLUGINS_SHARE_DIR "${TubeTK_INSTALL_SHARE_DIR}" )

if( TubeTK_USE_CPPCHECK )
  list( APPEND CMAKE_MODULE_PATH "${TubeTK_SOURCE_DIR}/CMake/Cppcheck" )
  find_package( Cppcheck REQUIRED )
  set( CPPCHECK_PATHS_TO_CHECK "${TubeTK_SOURCE_DIR}/Applications"
    "${TubeTK_SOURCE_DIR}/Base" )
  set( CPPCHECK_SUPPRESSIONS_FILE
    "${TubeTK_SOURCE_DIR}/CMake/Cppcheck/TubeTKCppcheckSuppressions.txt" )
  set( CPPCHECK_COMMAND_PREFIX "${CPPCHECK_EXECUTABLE}" -q --xml --suppressions
    "${CPPCHECK_SUPPRESSIONS_FILE}" )
  add_test( Cppcheck
    ${CPPCHECK_COMMAND_PREFIX} --enable=all ${CPPCHECK_PATHS_TO_CHECK} )
  set_tests_properties( Cppcheck
    PROPERTIES FAIL_REGULAR_EXPRESSION "error file=" )
endif( TubeTK_USE_CPPCHECK )

if( TubeTK_USE_QT )
  if( TubeTK_USE_CTK )
    find_package( CTK REQUIRED )
    include( ${CTK_USE_FILE} )
  endif( TubeTK_USE_CTK )
endif( TubeTK_USE_QT )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

find_package( JsonCpp REQUIRED )
include_directories( ${JsonCpp_INCLUDE_DIRS} )

if( TubeTK_USE_KWSTYLE )
  find_package( KWStyle REQUIRED )
  find_package( Git )
  if( GIT_FOUND )
    execute_process( COMMAND ${GIT_EXECUTABLE} config hooks.KWStyle.path
      "${KWSTYLE_EXECUTABLE}"
      WORKING_DIRECTORY ${TubeTK_SOURCE_DIR} )
  endif( GIT_FOUND )

  option( KWSTYLE_USE_VIM_FORMAT
  "Set KWStyle to generate errors in a Vim compatible format." OFF )
  mark_as_advanced( KWSTYLE_USE_VIM_FORMAT )

  option( KWSTYLE_USE_MSVC_FORMAT
    "Set KWStyle to generate errors in a Microsoft Visual Studio compatible format."
    OFF )
  mark_as_advanced( KWSTYLE_USE_MSVC_FORMAT )

  option( KWSTYLE_DASHBOARD_SUBMISSION
    "Set the KWStyle arguments to generate a report for CDash dashboard submission."
    OFF )
  mark_as_advanced( KWSTYLE_DASHBOARD_SUBMISSION )

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
    COMMENT "Coding style checker" VERBATIM )
endif( TubeTK_USE_KWSTYLE )

if( TubeTK_USE_LIBSVM )
  find_package( LIBSVM REQUIRED )
endif( TubeTK_USE_LIBSVM )

find_package( ParameterSerializer REQUIRED )
include_directories( ${ParameterSerializer_INCLUDE_DIRS} )

set( SlicerExecutionModel_CLI_RUNTIME_INSTALL_DIRECTORY
  ${TubeTK_INSTALL_BIN_DIR} )
set( SlicerExecutionModel_CLI_LIBRARY_INSTALL_DIRECTORY
  ${TubeTK_INSTALL_LIB_DIR} )
set( SlicerExecutionModel_CLI_INCLUDE_INSTALL_DIRECTORY
  ${TubeTK_INSTALL_INCLUDE_DIR} )
set( SlicerExecutionModel_CLI_SHARE_INSTALL_DIRECTORY
  ${TubeTK_INSTALL_SHARE_DIR} )
set( SlicerExecutionModel_CLI_ARCHIVE_INSTALL_DIRECTORY
  ${TubeTK_INSTALL_LIB_DIR} )

set( SlicerExecutionModel_CLI_RUNTIME_OUTPUT_DIRECTORY
  "${CMAKE_BINARY_DIR}/bin" )
set( SlicerExecutionModel_CLI_LIBRARY_OUTPUT_DIRECTORY
  "${CMAKE_BINARY_DIR}/lib/TubeTK" )
set( SlicerExecutionModel_CLI_INCLUDE_OUTPUT_DIRECTORY
  "${CMAKE_BINARY_DIR}/include/TubeTK" )
set( SlicerExecutionModel_CLI_SHARE_OUTPUT_DIRECTORY
  "${CMAKE_BINARY_DIR}/share/TubeTK" )
set( SlicerExecutionModel_CLI_ARCHIVE_OUTPUT_DIRECTORY
  "${CMAKE_BINARY_DIR}/lib/TubeTK" )

find_package( SlicerExecutionModel REQUIRED GenerateCLP )
include( ${GenerateCLP_USE_FILE} )

if( TubeTK_USE_VTK )
  find_package( VTK REQUIRED )
  include( ${VTK_USE_FILE} )
endif( TubeTK_USE_VTK )

# Check endianness.
include( TestBigEndian )
test_big_endian( CMAKE_WORDS_BIGENDIAN )

# Set TubeTK_EXECUTABLE_DIRS useful to configure launcher scripts
# and set TubeTK_LIBRARY_DIRS
if( MSVC )
  set(TubeTK_EXECUTABLE_DIRS )
  foreach(_build_type "" Debug Release)
    list( APPEND TubeTK_EXECUTABLE_DIRS
      ${JsonCpp_DIR}/bin/${_build_type}
      ${ParameterSerializer_DIR}/bin/${_build_type}
      ${GenerateCLP_DIR}/${_build_type}
      ${ITK_DIR}/bin/${_build_type}
      ${TubeTK_BINARY_DIR}/bin/${_build_type}
      ${TubeTK_BINARY_DIR}/lib/TubeTK/Plugins/${_build_type} )
    if( TubeTK_USE_VTK )
      list( APPEND TubeTK_EXECUTABLE_DIRS
        ${VTK_DIR}/bin/${_build_type} )
    endif( TubeTK_USE_VTK )
    if( TubeTK_USE_SLICER )
      list( APPEND TubeTK_EXECUTABLE_DIRS
        ${Slicer_DIR}/bin/${_build_type} )
    endif( TubeTK_USE_SLICER )
  endforeach()
else( MSVC )
  set( TubeTK_EXECUTABLE_DIRS
    ${JsonCpp_DIR}/bin
    ${ParameterSerializer_DIR}/bin
    ${GenerateCLP_DIR}
    ${ITK_DIR}/bin
    ${TubeTK_BINARY_DIR}/bin
    ${TubeTK_BINARY_DIR}/lib
    ${TubeTK_BINARY_DIR}/lib/TubeTK/Plugins )
  if( TubeTK_USE_VTK )
    list( APPEND ${VTK_DIR}/bin )
  endif( TubeTK_USE_VTK )
endif( MSVC )


# Configure a launcher for running TubeTK methods from the command line
message(STATUS "Configuring Launcher script")
set(TubeTK_LAUNCHER )
if( WIN32 )
  # Microsoft Windows.
  set(_launcher_platform "windows")
  configure_file( ${TubeTK_SOURCE_DIR}/CMake/TubeTKLauncher.bat.in
    ${TubeTK_BINARY_DIR}/TubeTKLauncher.bat @ONLY )

  find_program(CMD_EXECUTABLE "cmd")
  if( NOT CMD_EXECUTABLE )
    message(FATAL_ERROR "Could not find 'cmd' executable required to run test using the launcher")
  endif()

  set(TubeTK_LAUNCHER ${TubeTK_BINARY_DIR}/TubeTKLauncher.bat)
elseif( UNIX AND NOT TubeTK_USE_SLICER )
  # Apple OS X, Linux, or UNIX.
  set(_launcher_platform "unix")
  configure_file( ${TubeTK_SOURCE_DIR}/CMake/TubeTKLauncher.sh.in
    ${TubeTK_BINARY_DIR}/TubeTKLauncher.sh @ONLY )

  find_program(SH_EXECUTABLE "sh")
  if( NOT SH_EXECUTABLE )
    message(FATAL_ERROR "Could not find 'sh' executable required to run test using the launcher")
  endif()
  set(TubeTK_LAUNCHER ${SH_EXECUTABLE} ${TubeTK_BINARY_DIR}/TubeTKLauncher.sh)
elseif( NOT UNIX )
  message(FATAL_ERROR "Configuring Launcher script - failed [unknown platform]")
endif()
message(STATUS "Configuring Launcher script - ok [${_launcher_platform}]")

# Setup to conform to 3D Slicer paths and define plugin copying macros.
if( TubeTK_BUILD_MODULES AND NOT TubeTK_USE_SLICER )
  include( Slicer4MacroGetSource )

  Slicer4MacroGetSource( Resources )
  Slicer4MacroGetSource( Base/CLI )
  Slicer4MacroGetSource( Libs/ITKFactoryRegistration )

  Slicer4MacroGetCLISource( ExpertAutomatedRegistration )

  # Build ITKFactorRegistration
  include( Slicer4Config )
  set( CXX_TEST_PATH ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR} )
  set( LAUNCH_EXE "" )
  set( SUPPORT_2D_IMAGES on )

  set( ITKFactoryRegistration_INSTALL_BIN_DIR ${TubeTK_INSTALL_BIN_DIR} )
  set( ITKFactoryRegistration_INSTALL_LIB_DIR ${TubeTK_INSTALL_LIB_DIR} )
  add_subdirectory( ${TubeTK_SLICER_SOURCE_DIR}/Libs/ITKFactoryRegistration
    ${TubeTK_SLICER_BINARY_DIR}/Libs/ITKFactoryRegistration )
  include_directories( ${Slicer_BINARY_DIR}/Libs/itkFactoryRegistration )

  # Build SlicerBaseCLI
  add_subdirectory( ${TubeTK_SLICER_SOURCE_DIR}/Base/CLI
    ${TubeTK_SLICER_BINARY_DIR}/Base/CLI )

  # Config SlicerExecutionModel
  set( extras ITKFactoryRegistration )
  if( TubeTK_USE_VTK )
    list( APPEND extras SlicerBaseCLI )
  endif( TubeTK_USE_VTK )
  set( SlicerExecutionModel_EXTRA_EXECUTABLE_TARGET_LIBRARIES
    ${extras}
    CACHE INTERNAL "3D Slicer execution model extra libraries." FORCE )

endif( TubeTK_BUILD_MODULES AND NOT TubeTK_USE_SLICER )

include( TubeTKMacroInstallPlugins )

# Setup Python VirtualEnv for testing.
if( BUILD_TESTING AND TubeTK_USE_PYTHON AND NOT TubeTK_USE_SLICER )
  include( CMake/TubeTKVirtualEnvSetup.cmake )
endif( BUILD_TESTING AND TubeTK_USE_PYTHON AND NOT TubeTK_USE_SLICER )

# Test TubeTK IPython examples.
if( TubeTK_USE_NOTEBOOKS )
  add_subdirectory( Examples )
endif( TubeTK_USE_NOTEBOOKS )

add_subdirectory( Utilities )
add_subdirectory( Base )
if( TubeTK_BUILD_MODULES )
  add_subdirectory( Applications )
endif( TubeTK_BUILD_MODULES )

# Optionally build 3D Slicer modules.
if( TubeTK_BUILD_MODULES AND TubeTK_USE_SLICER )
  add_subdirectory( SlicerModules )
endif( TubeTK_BUILD_MODULES AND TubeTK_USE_SLICER )

# Create the directory where CTest unit tests store temporary results.
make_directory( ${TubeTK_BINARY_DIR}/Temporary )

# Configure variables to export.
set( TubeTK_INCLUDE_DIRS
  ${TubeTK_SOURCE_DIR}/Base/CLI
  ${TubeTK_SOURCE_DIR}/Base/Common
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/IO
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/ObjectDocuments
  ${TubeTK_SOURCE_DIR}/Base/Registration
  ${TubeTK_SOURCE_DIR}/Base/Segmentation
  ${TubeTK_SOURCE_DIR}/Base/USTK )

set( TubeTK_LIBRARY_DIRS
  ${TubeTK_EXECUTABLE_DIRS} )

# Configure TubeTKTargets.cmake.
get_property( TubeTK_TARGETS GLOBAL PROPERTY TubeTK_TARGETS )
export( TARGETS ${TubeTK_TARGETS} FILE ${TubeTK_CONFIG_BINARY_DIR}/TubeTKTargets.cmake )

set( TubeTK_TARGETS_FILE
  ${TubeTK_CONFIG_BINARY_DIR}/TubeTKTargets.cmake )

set( TubeTK_USE_FILE
  ${TubeTK_CONFIG_BINARY_DIR}/UseTubeTK.cmake )

# Configure TubeTKConfig.cmake.
# The variables set above are for building against a TubeTK build tree.
# A different TubeTKConfig.cmake will need to be generated in order to build
# against a TubeTK installation.
configure_file( ${TubeTK_SOURCE_DIR}/CMake/TubeTKConfig.cmake.in
  ${TubeTK_CONFIG_BINARY_DIR}/TubeTKConfig.cmake @ONLY )

# Configure UseTubeTK.cmake.
configure_file( ${TubeTK_SOURCE_DIR}/CMake/UseTubeTK.cmake.in
  ${TubeTK_CONFIG_BINARY_DIR}/UseTubeTK.cmake @ONLY )

# Configure CPack packaging.
if( NOT TubeTK_USE_SLICER )
  # Record ITK, VTK, CTK, and GenerateCLP libraries used by TubeTK for
  # packaging.
  file( GLOB_RECURSE ITK_BINARIES "${ITK_DIR}/lib/lib*" )
  if( TubeTK_USE_VTK )
    file( GLOB_RECURSE VTK_BINARIES "${VTK_DIR}/bin/lib*" )
  endif( TubeTK_USE_VTK )
  if( TubeTK_USE_CTK )
    file( GLOB_RECURSE CTK_BINARIES "${CTK_DIR}/CTK-build/bin/lib*" )
  endif( TubeTK_USE_CTK )
  file( GLOB_RECURSE SEM_BINARIES "${GenerateCLP_DIR}/bin/lib*" )
  file( GLOB_RECURSE JsonCpp_BINARIES "${JsonCpp_DIR}/lib/lib*" )
  list( APPEND TubeTK_SUPPORT_BINARIES
    ${SEM_BINARIES} ${CTK_BINARIES} ${VTK_BINARIES} ${ITK_BINARIES}
    ${JsonCpp_BINARIES} )
  install( FILES ${TubeTK_SUPPORT_BINARIES}
    DESTINATION lib )

  # Record the environment setup files for packaging.
  set( TubeTK_LINUX_SETUP
    ${TubeTK_SOURCE_DIR}/CMake/TubeTKSourceMeToSetupLinux.sh )
  set( TubeTK_WINDOWS_SETUP
    ${TubeTK_SOURCE_DIR}/CMake/TubeTKRunMeToSetupWindows.bat )

  include( InstallRequiredSystemLibraries )

  set( CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "TubeTK for the analysis of geometric structures in images." )
  set( CPACK_PACKAGE_VENDOR "Kitware, Inc." )
  set( CPACK_PACKAGE_DESCRIPTION_FILE "${TubeTK_SOURCE_DIR}/CMake/CPack/DESCRIPTION.txt" )
  set( CPACK_RESOURCE_FILE_LICENSE "${TubeTK_SOURCE_DIR}/CMake/CPack/LICENSE.txt" )

  set( CPACK_PACKAGE_VERSION_MAJOR ${TubeTK_VERSION_MAJOR} )
  set( CPACK_PACKAGE_VERSION_MINOR ${TubeTK_VERSION_MINOR} )
  set( CPACK_PACKAGE_VERSION_PATCH ${TubeTK_VERSION_PATCH} )
  set( CPACK_PACKAGE_VERSION ${TubeTK_VERSION_BUILD} )

  set( CPACK_PACKAGE_CONTACT "stephen.aylward@kitware.com" )

  set( CPACK_PACKAGE_INSTALL_DIRECTORY
    "TubeTK_${TubeTK_VERSION_MAJOR}.${TubeTK_VERSION_MINOR}" )

  if( WIN32 AND NOT UNIX )
    # There is a bug in NSI that does not handle full UNIX paths properly. Make
    # sure there is at least one set of four backslashes.
    set( CPACK_PACKAGE_ICON
      "${TubeTK_SOURCE_DIR}\\\\Icon.png" )
    set( CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\TubeTK.exe" )
    set( CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} TubeTK" )
    set( CPACK_NSIS_HELP_LINK "http:\\\\\\\\tubetk.org" )
    set( CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\tubetk.org" )
    set( CPACK_NSIS_CONTACT "stephen.aylward@kitware.com" )
    set( CPACK_NSIS_MODIFY_PATH ON )
  else( WIN32 AND NOT UNIX )
    set( CPACK_STRIP_FILES true )
    set( CPACK_SOURCE_STRIP_FILES true )
  endif( WIN32 AND NOT UNIX )

  set( CMAKE_INSTALL_RPATH "$ORIGIN/../..;." )

  # Standard CPack setup.
  option( CMAKE_INSTALL_DEBUG_LIBRARIES
    "Install the Microsoft runtime debug libraries with CMake." FALSE )
  mark_as_advanced( CMAKE_INSTALL_DEBUG_LIBRARIES )

  # By default, do not warn when built on machines using only Microsoft Visual
  # Studio Express.
  if( NOT DEFINED CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS )
    set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS ON )
  endif( NOT DEFINED CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS )

  # Make this explicit here, rather than accepting the CPack default value, so
  # we can refer to it.
  set( CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}" )

  # Installers for 32- vs. 64-bit CMake:
  #   - Root install directory (displayed to end user at installer runtime).
  #   - NSIS package/display name (text used in the installer GUI).
  #   - Registry key used to store information about the installation.
  if( CMAKE_CL_64 )
    set( CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64" )
    set( CPACK_NSIS_PACKAGE_NAME
      "${CPACK_PACKAGE_INSTALL_DIRECTORY} ( Win64 )" )
    set( CPACK_PACKAGE_INSTALL_REGISTRY_KEY
      "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION} ( Win64 )" )
  else( CMAKE_CL_64 )
    set( CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES" )
    set( CPACK_NSIS_PACKAGE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}" )
    set( CPACK_PACKAGE_INSTALL_REGISTRY_KEY
      "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION}" )
  endif( CMAKE_CL_64 )

  if( NOT DEFINED CPACK_SYSTEM_NAME )
    set( CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR} )
  endif( NOT DEFINED CPACK_SYSTEM_NAME )

  if( ${CPACK_SYSTEM_NAME} MATCHES Windows )
    if( CMAKE_CL_64 )
      set( CPACK_SYSTEM_NAME win64-x64 )
    else( CMAKE_CL_64 )
      set( CPACK_SYSTEM_NAME win32-x86 )
    endif( CMAKE_CL_64 )
  endif( ${CPACK_SYSTEM_NAME} MATCHES Windows )

  set( CPACK_SOURCE_PACKAGE_FILE_NAME "TubeTK-${TubeTK_VERSION_BUILD}" )
  set( CPACK_PACKAGE_FILE_NAME
    "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}" )

  # Include CPack model once all variables are set.
  include( CPack )
else( NOT TubeTK_USE_SLICER )
  include( ${Slicer_EXTENSION_CPACK} )
endif( NOT TubeTK_USE_SLICER )
