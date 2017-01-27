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

cmake_minimum_required( VERSION 3.0.0 )

# Update CMake module path
set( TubeTK_CMAKE_EXTENSIONS_DIR ${TubeTK_SOURCE_DIR}/CMake )
set( CMAKE_MODULE_PATH ${TubeTK_CMAKE_EXTENSIONS_DIR} ${CMAKE_MODULE_PATH} )

include( ExternalProject )
include( TubeTKExternalProjectDependency )
include( CMakeDependentOption )
include( TubeTKCheckPythonLibrary )

# Set TubeTK version and output directories.
#
set( TubeTK_VERSION_MAJOR "0" )
set( TubeTK_VERSION_MINOR "9" )
set( TubeTK_VERSION_PATCH "0" )

option( TubeTK_INSTALL_DEVELOPMENT
  "Add install rules for TubeTK SDK dependencies (headers, libraries, ...)" ON )
mark_as_advanced( TubeTK_INSTALL_DEVELOPMENT )

mark_as_advanced( CMAKE_BACKWARDS_COMPATIBILITY )
mark_as_advanced( EXECUTABLE_OUTPUT_PATH )
mark_as_advanced( LIBRARY_OUTPUT_PATH )
mark_as_advanced( CMAKE_INSTALL_PREFIX )
mark_as_advanced( BUILD_TESTING )
mark_as_advanced( BUILD_DOCUMENTATION )

include( TubeTKFunctionToday )
TubeTKFunctionToday( TubeTK_VERSION_BUILD )

set( TubeTK_VERSION
  "${TubeTK_VERSION_MAJOR}.${TubeTK_VERSION_MINOR}.${TubeTK_VERSION_PATCH}.${TubeTK_VERSION_BUILD}" )

# Set a variable if not already defined
macro( tubetk_set_if_not_defined var defaultvalue )
  if( NOT DEFINED ${var} )
    set( ${var} "${defaultvalue}" )
  endif()
endmacro()

#
# TubeTK install directories
#
tubetk_set_if_not_defined( CMAKE_BUILD_TYPE "Release" )
tubetk_set_if_not_defined( TubeTK_INSTALL_BIN_DIR "bin" )
tubetk_set_if_not_defined( TubeTK_INSTALL_LIB_DIR "lib" )
tubetk_set_if_not_defined( TubeTK_INSTALL_INCLUDE_DIR "include" )
tubetk_set_if_not_defined( TubeTK_INSTALL_SHARE_DIR "share/TubeTK" )

# Pass value to TubeTK inner build
mark_as_superbuild( VARS
  TubeTK_INSTALL_BIN_DIR
  TubeTK_INSTALL_LIB_DIR
  TubeTK_INSTALL_INCLUDE_DIR
  TubeTK_INSTALL_SHARE_DIR )

#
# TubeTK applications install directories
#
tubetk_set_if_not_defined( SlicerExecutionModel_CLI_INSTALL_RUNTIME_DESTINATION
  ${TubeTK_INSTALL_BIN_DIR} )
tubetk_set_if_not_defined( SlicerExecutionModel_CLI_INSTALL_LIBRARY_DESTINATION
  ${TubeTK_INSTALL_LIB_DIR} )
tubetk_set_if_not_defined( SlicerExecutionModel_CLI_INSTALL_INCLUDE_DESTINATION
  ${TubeTK_INSTALL_INCLUDE_DIR} )
tubetk_set_if_not_defined( SlicerExecutionModel_CLI_INSTALL_SHARE_DESTINATION
  ${TubeTK_INSTALL_SHARE_DIR} )
tubetk_set_if_not_defined( SlicerExecutionModel_CLI_INSTALL_ARCHIVE_DESTINATION
  ${TubeTK_INSTALL_LIB_DIR} )

# Pass value to TubeTK inner build
mark_as_superbuild( VARS
  SlicerExecutionModel_CLI_INSTALL_RUNTIME_DESTINATION
  SlicerExecutionModel_CLI_INSTALL_LIBRARY_DESTINATION
  SlicerExecutionModel_CLI_INSTALL_INCLUDE_DESTINATION
  SlicerExecutionModel_CLI_INSTALL_SHARE_DESTINATION
  SlicerExecutionModel_CLI_INSTALL_ARCHIVE_DESTINATION )

#
# Default output directories
#
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY
  "${CMAKE_BINARY_DIR}/${TubeTK_INSTALL_BIN_DIR}" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY
  "${CMAKE_BINARY_DIR}/${TubeTK_INSTALL_LIB_DIR}" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  "${CMAKE_BINARY_DIR}/${TubeTK_INSTALL_LIB_DIR}" )

#
# CDash dashboard setup.
#
include( CTest )
if( BUILD_TESTING )
  if( NOT TubeTK_SUPPORT_2D_IMAGES )
    message( STATUS
     "Forcing TubeTK_SUPPORT_2D_IMAGES to ON because BUILD_TESTING is ON" )
    set( TubeTK_SUPPORT_2D_IMAGES ON CACHE BOOL
      "Compile CLI with support for 2D images" FORCE )
  endif(  NOT TubeTK_SUPPORT_2D_IMAGES  )
  set( BUILDNAME "${BUILDNAME}" CACHE STRING
    "Name of the build on the CDash dashboard." )
  mark_as_advanced( BUILDNAME )
endif( BUILD_TESTING )

set( BUILD_WARNING_REPORT_LIMIT -1 )
set( BUILD_ERROR_REPORT_LIMIT -1 )

include( CTestConfig.cmake )

configure_file( ${TubeTK_SOURCE_DIR}/CMake/CTestCustom.cmake.in
  ${TubeTK_BINARY_DIR}/CTestCustom.cmake )

set( TubeTK_BUILD_TESTING ${BUILD_TESTING} )

#
# Configure for the download of CTest unit testing data from Midas
#
set( MIDAS_USER "" CACHE STRING "Midas Platform user email address." )
mark_as_advanced( MIDAS_USER )
set( MIDAS_DEFAULT_API_KEY "" CACHE STRING "Midas Platform default API key." )
mark_as_advanced( MIDAS_DEFAULT_API_KEY )
set( MIDAS_DATA_DIR ${TubeTK_BINARY_DIR}/MIDAS_Data )

#
# Doxygen documentation setup.
#
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

# Append an item to a string if it does not already contain it.  Any additional
# arguments are treated as extra "set" arguments
macro( tubetk_append_if_absent VARNAME VALUE )
  if( NOT ${VARNAME} MATCHES "(^| )${VALUE}($| )" )
    if( ${VARNAME} )
      set( ${VARNAME} "${${VARNAME}} ${VALUE}" ${ARGN} )
    else()
      set( ${VARNAME} "${VALUE}" ${ARGN} )
    endif()
  endif()
endmacro()

#
# Microsoft Visual Studio setup.
#
if( MSVC )
  add_definitions( -D_SCL_SECURE_NO_WARNINGS )
  if( "${CMAKE_CXX_COMPILER}" MATCHES "cl" )
    tubetk_append_if_absent( CMAKE_CXX_FLAGS "/bigobj" )
    tubetk_append_if_absent( CMAKE_CXX_FLAGS "/MP" )
  endif( "${CMAKE_CXX_COMPILER}" MATCHES "cl" )
  option( DISABLE_MSVC8_DEPRECATED_WARNINGS
    "Disable Microsoft Visual Studio 8 deprecated warnings." ON )
  mark_as_advanced( FORCE DISABLE_MSVC8_DEPRECATED_WARNINGS )
  if( DISABLE_MSVC8_DEPRECATED_WARNINGS )
    add_definitions( -D_CRT_SECURE_NO_DEPRECATE )
  endif( DISABLE_MSVC8_DEPRECATED_WARNINGS )
endif( MSVC )

#
# Apple OS X setup.
#
if( APPLE )
  include( TubeTKBlockSetCMakeOSXVariables )
  if( NOT ${CMAKE_VERSION} VERSION_LESS 3.1 )
    cmake_policy( SET CMP0042 OLD )
  endif( NOT ${CMAKE_VERSION} VERSION_LESS 3.1 )
  tubetk_append_if_absent( CMAKE_CXX_FLAGS "-fno-common" )
  tubetk_append_if_absent( CMAKE_C_FLAGS "-fno-common" )
endif( APPLE )

#
# Linux or UNIX setup.
#
if( UNIX )
  tubetk_append_if_absent( CMAKE_CXX_FLAGS "-fno-strict-aliasing" )
  tubetk_append_if_absent( CMAKE_C_FLAGS "-fno-strict-aliasing" )
  if( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
    tubetk_append_if_absent( CMAKE_CXX_FLAGS "-fPIC" )
    tubetk_append_if_absent( CMAKE_C_FLAGS "-fPIC" )
  endif( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
endif( UNIX )

#
# Gcc setup.
#
if( CMAKE_COMPILER_IS_GNUCXX )
  include( CheckCCompilerFlag )
  check_c_compiler_flag(-Wno-unused-local-typedefs
    HAS_NO_UNUSED_LOCAL_TYPEDEFS)
  if( HAS_NO_UNUSED_LOCAL_TYPEDEFS )
    tubetk_append_if_absent( CMAKE_CXX_FLAGS "-Wno-unused-local-typedefs" )
    tubetk_append_if_absent( CMAKE_C_FLAGS "-Wno-unused-local-typedefs" )
  endif( HAS_NO_UNUSED_LOCAL_TYPEDEFS )
endif( CMAKE_COMPILER_IS_GNUCXX )

#
# Implementation agnostic fPIC
#
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

#
# Build TubeTK to operate as a Slicer module
#
option( TubeTK_BUILD_USING_SLICER
  "Use Slicer to extract useful commands and CLIs." OFF )
option( TubeTK_BUILD_WITHIN_SLICER
  "Whether TubeTK's CLIs modules are configured to be run from Slicer or alone."
  OFF )
mark_as_advanced( TubeTK_BUILD_WITHIN_SLICER )
# Slicer's extension build tool defines TubeTK_BUILD_SLICER_EXTENSION
# to indicate that this project has to be compiled using Slicer, within Slicer.
if( TubeTK_BUILD_SLICER_EXTENSION )
  set(TubeTK_BUILD_USING_SLICER TRUE)
  set(TubeTK_BUILD_WITHIN_SLICER TRUE)
endif()
if( TubeTK_BUILD_USING_SLICER )
  include( TubeTKSlicer4ExtensionConfig )

  mark_as_advanced( TubeTK_USE_ITK TubeTK_USE_VTK )

  #set( TubeTK_BUILD_WITHIN_SLICER ${TubeTK_BUILD_USING_SLICER} )
  # Do not unset these variables if project is built as
  # a Slicer extension within Slicer.
  if( NOT TubeTK_BUILD_WITHIN_SLICER)
    unset( PYTHON_EXECUTABLE CACHE)
    unset( PYTHON_INCLUDE_DIR CACHE)
    unset( PYTHON_LIBRARY CACHE)
    unset( PYTHON_DEBUG_LIBRARY CACHE)
  endif()

  # Find Slicer to set the [ITK|VTK|...]_DIR variables
  find_package( Slicer REQUIRED )
  if( NOT Slicer_USE_FILE )
    message( FATAL_ERROR
      "Slicer_USE_FILE not defined. Find Slicer silently failed." )
  endif( NOT Slicer_USE_FILE )
  include( ${Slicer_USE_FILE} )

  mark_as_advanced( CLEAR Subversion_SVN_EXECUTABLE )
  find_package( Subversion REQUIRED )
  mark_as_advanced( CLEAR Subversion_SVN_EXECUTABLE )

  if( Slicer_SOURCE_DIR )
    set( TubeTK_SLICER_SOURCE_DIR "${Slicer_SOURCE_DIR}" )
  else( Slicer_SOURCE_DIR )
    set( TubeTK_SLICER_SOURCE_DIR "${Slicer_CMAKE_DIR}/.." )
  endif( Slicer_SOURCE_DIR )

  if( Slicer_BINARY_DIR )
    set( TubeTK_SLICER_BINARY_DIR "${Slicer_BINARY_DIR}" )
  else( Slicer_BINARY_DIR )
    set( TubeTK_SLICER_BINARY_DIR "${Slicer_HOME}" )
  endif( Slicer_BINARY_DIR )

  # All those dependencies are already provided by Slicer
  set( USE_SYSTEM_ITK ON CACHE BOOL
    "Exclude ITK from superbuild and use an existing build." FORCE )
  set( TubeTK_USE_VTK ON CACHE BOOL
    "Use VTK to enable additional functionality." FORCE )
  set( USE_SYSTEM_VTK ON CACHE BOOL
    "Exclude VTK from superbuild and use an existing build." FORCE )
  set( USE_SYSTEM_SLICER_EXECUTION_MODEL ON CACHE BOOL
    "Exclude SlicerExecutionModel from superbuild and use an existing build."
    FORCE )
  mark_as_advanced( CTK_DIR DCMTK_DIR ITK_DIR LibArchive_INCLUDE_DIR
    LibArchive_LIBRARY OpenIGTLink_DIR QT_QMAKE_EXECUTABLE QtTesting_DIR
    SlicerExecutionModel_DIR Slicer_EXTENSIONS_TRACK_QUALIFIER Teem_DIR
    VTK_DIR ZLIB_INCLUDE_DIR ZLIB_LIBRARY ZLIB_ROOT jqPlot_DIR qRestAPI_DIR
    CURL_INCLUDE_DIR CURL_LIBRARY ZLIB_LIBRARY ZLIB_ROOT ZLIB_INCLUDE_DIR
    jpPlot_DIR MIDAS_PACKAGE_API_KEY MIDAS_PACKAGE_EMAIL MIDAS_PACKAGE_URL
    SimpleITK_DIR GENERATECLP_EXE OPENSSL_CRYPTO_LIBRARY OPENSSL_INCLUDE_DIR
    OPENSSL_SSL_LIBRARY HIDE_CONSOLE NUMPY_INCLUDE_DIR SH_EXECUTABLE
    LIB_EAY_DEBUG LIB_EAY_RELEASE SSL_EAY_DEBUG
    SSL_EAY_RELEASE CMAKE_CXX_COMPILER CMAKE_C_COMPILER CMAKE_CXX_MP_FLAG
    CMAKE_CXX_MP_NUM_PROCESSORS ParameterSerializer_DIR )
else( TubeTK_BUILD_USING_SLICER )
  mark_as_advanced( CLEAR TubeTK_USE_ITK TubeTK_USE_VTK )
  mark_as_advanced( Slicer_DIR )
endif( TubeTK_BUILD_USING_SLICER )

include( TubeTKSlicer4Config )

#
# 2D image support.
#
option( TubeTK_SUPPORT_2D_IMAGES "Compile CLI with support for 2D images" ON )
mark_as_advanced( TubeTK_SUPPORT_2D_IMAGES )
if( TubeTK_SUPPORT_2D_IMAGES )
  add_definitions( -DSUPPORT_2D_IMAGES )
endif()

#
# GPU ArrayFire setup
#
option( TubeTK_USE_GPU_ARRAYFIRE
  "Use the ArrayFire library to speedup filtering opertions using the GPU." OFF)

#
# Boost libraries setup.
#
option( TubeTK_USE_BOOST
  "Use the Boost libraries to enable additional functionality." OFF )
if( TubeTK_USE_BOOST )
  set( Boost_USE_STATIC_LIBS ON )
  set( Boost_USE_MULTITHREADED ON )
  set( Boost_USE_STATIC_RUNTIME OFF )
  find_package( Boost REQUIRED )
  if( NOT Boost_FOUND )
    set( BOOST_ROOT "NOT_FOUND" CACHE PATH
      "Path to the top-level directory of Boost." FORCE )
    mark_as_advanced( CLEAR BOOST_ROOT )
  endif( NOT Boost_FOUND )
endif( TubeTK_USE_BOOST )

#
#  Cppcheck setup.
#
option( TubeTK_USE_CPPCHECK "Use Cppcheck to enable additional tests." ON )
option( USE_SYSTEM_CPPCHECK
  "Exclude Cppcheck from superbuild and use an existing build." OFF )
mark_as_advanced( USE_SYSTEM_CPPCHECK )

#
# Setup python
#
CMAKE_DEPENDENT_OPTION( TubeTK_USE_PYTHON
  "Use Python to enable additional functionality." ON
  "NOT TubeTK_BUILD_USING_SLICER" OFF )

CMAKE_DEPENDENT_OPTION( TubeTK_USE_NUMPY_STACK
  "Use NumPy to enable additional functionality." ON
  "TubeTK_USE_PYTHON" OFF )

CMAKE_DEPENDENT_OPTION( TubeTK_USE_PYQTGRAPH
  "Use PyQtGraph to enable additional functionality." ON
  "TubeTK_USE_PYTHON" OFF )

if( TubeTK_USE_PYTHON )
  find_package( PythonInterp REQUIRED )
  find_package( PythonLibs REQUIRED )

  # NumPy conversion support.
  if( TubeTK_USE_NUMPY_STACK )
    TubeTKCheckPythonLibraries( REQUIRED LIBRARIES numpy scipy
      ERROR_MESSAGE "Set TubeTK_USE_NUMPY_STACK to OFF" )
  endif()

  # PyQtGraph bridge support.
  if( TubeTK_USE_PYQTGRAPH )
    TubeTKCheckPythonLibraries( REQUIRED LIBRARIES pyqtgraph
      ERROR_MESSAGE "Set TubeTK_USE_PYQTGRAPH to OFF" )
  endif()

  set( NOTEBOOK_TEST_DRIVER
    ${TubeTK_SOURCE_DIR}/Utilities/Python/EvaluateIPythonNotebook.py
    CACHE INTERNAL "Test driver command for IPython Notebooks." FORCE )
endif( TubeTK_USE_PYTHON )

#
# ITK setup.
#
option( USE_SYSTEM_ITK "Exclude ITK from superbuild and use an existing build."
  OFF )
mark_as_advanced( USE_SYSTEM_ITK )

#
# JsonCpp library setup.
#
option( TubeTK_USE_JsonCpp "Use JsonCpp to enable additional functionality."
  ON )
option( USE_SYSTEM_JSONCPP
  "Exclude JsonCpp from superbuild and use an existing build." OFF )
mark_as_advanced( USE_SYSTEM_JSONCPP )
mark_as_advanced( JsonCpp_DIR )

#
# KWStyle setup.
#
option( TubeTK_USE_KWSTYLE
  "Use KWStyle to check for coding standards violations." ON )
option( USE_SYSTEM_KWSTYLE
  "Exclude KWStyle from superbuild and use an existing build." OFF )
mark_as_advanced( USE_SYSTEM_KWSTYLE )

#
# LIBSVM setup.
#
option( TubeTK_USE_LIBSVM "Use LIBSVM to enable additional functionality." OFF )
option( USE_SYSTEM_LIBSVM
  "Exclude LIBSVM from superbuild and use an existing build." OFF )
mark_as_advanced( USE_SYSTEM_LIBSVM )
mark_as_advanced( LIBSVM_DIR )

#
# RandomForest setup.
#
option( TubeTK_USE_RandomForest
  "Use RandomForest to enable additional functionality." OFF )
option( USE_SYSTEM_RandomForest
  "Exclude RandomForest from superbuild and use an existing build." OFF )
mark_as_advanced( USE_SYSTEM_RandomForest )
mark_as_advanced( RandomForest_DIR )


#
# Slicer Execution model Setup
#
option( USE_SYSTEM_SLICER_EXECUTION_MODEL
  "Exclude SlicerExecutionModel from superbuild and use an existing build."
  OFF )
mark_as_advanced( USE_SYSTEM_SLICER_EXECUTION_MODEL )

#
# Provide option for Valgrind
#
option( TubeTK_USE_VALGRIND "Use Valgrind to enable additional functionality."
  OFF )

#
# VTK setup.
#
option( TubeTK_USE_VTK "Use VTK to enable additional functionality." ON )
option( USE_SYSTEM_VTK "Exclude VTK from superbuild and use an existing build."
  OFF )
mark_as_advanced( USE_SYSTEM_VTK )

#
# Modules
#
option( TubeTK_BUILD_APPLICATIONS
  "Build applications. If not, only TubeTK library is built" ON )
mark_as_advanced( TubeTK_BUILD_APPLICATIONS )

# If TubeTK_CONFIG_BINARY_DIR isn't defined, it means TubeTK is *NOT* being
# built using Superbuild. In that specific case, TubeTK_CONFIG_BINARY_DIR
# should default to TubeTK_BINARY_DIR
if( NOT DEFINED TubeTK_CONFIG_BINARY_DIR )
  set( TubeTK_CONFIG_BINARY_DIR ${CMAKE_BINARY_DIR} )
endif( NOT DEFINED TubeTK_CONFIG_BINARY_DIR )
mark_as_advanced( TubeTK_CONFIG_BINARY_DIR )

# Setup Python VirtualEnv for testing.
if( BUILD_TESTING AND TubeTK_USE_PYTHON AND NOT TubeTK_BUILD_USING_SLICER )
  option( BUILD_TESTING_VIRTUAL_ENV
    "Build a Python virtual environment to run tests." OFF )
  if( BUILD_TESTING_VIRTUAL_ENV )
    include( TubeTKVirtualEnvSetup )
  else()
    find_package( PythonInterp REQUIRED )
    find_package( PythonLibs REQUIRED )
    set( PYTHON_TESTING_EXECUTABLE ${PYTHON_EXECUTABLE} )
  endif()
endif( BUILD_TESTING AND TubeTK_USE_PYTHON AND NOT TubeTK_BUILD_USING_SLICER )

#
# Testing using Examples
#
option( TubeTK_USE_EXAMPLES_AS_TESTS
 "Use TubeTK examples as (long running) tests." OFF )

#
# Setup Superbuild
#
option( TubeTK_USE_SUPERBUILD
  "Build TubeTK and the projects it depends on using a superbuild." ON )
mark_as_advanced( TubeTK_USE_SUPERBUILD )

if( TubeTK_USE_SUPERBUILD )
  add_subdirectory( "${CMAKE_CURRENT_SOURCE_DIR}/CMake/Superbuild" )
  return()
endif( TubeTK_USE_SUPERBUILD )

#
# Once superbuild is done, build TubeTK
#

set( TubeTK_INSTALL_PLUGINS_BIN_DIR "${TubeTK_INSTALL_BIN_DIR}" )
set( TubeTK_INSTALL_PLUGINS_CACHE_DIR "${TubeTK_INSTALL_LIB_DIR}/CACHE" )
set( TubeTK_INSTALL_PLUGINS_INCLUDE_DIR "${TubeTK_INSTALL_INCLUDE_DIR}" )
set( TubeTK_INSTALL_PLUGINS_LIB_DIR "${TubeTK_INSTALL_LIB_DIR}" )
set( TubeTK_INSTALL_PLUGINS_SHARE_DIR "${TubeTK_INSTALL_SHARE_DIR}" )

find_package( JsonCpp REQUIRED )
include_directories( ${JsonCpp_INCLUDE_DIRS} )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )
if( TubeTK_USE_VTK )
  find_package( VTK REQUIRED )
  include( ${VTK_USE_FILE} )
endif( TubeTK_USE_VTK )

if( TubeTK_USE_LIBSVM )
  find_package( LIBSVM REQUIRED )
endif( TubeTK_USE_LIBSVM )

if( TubeTK_USE_RandomForest )
  find_package( RandomForest REQUIRED )
endif( TubeTK_USE_RandomForest )

if( TubeTK_USE_GPU_ARRAYFIRE )
  find_package( ArrayFire REQUIRED )
  include_directories( ${ArrayFire_INCLUDE_DIRS} )
endif( TubeTK_USE_GPU_ARRAYFIRE )

if( TubeTK_USE_BOOST )
  set( Boost_USE_STATIC_LIBS ON )
  set( Boost_USE_MULTITHREADED ON )
  set( Boost_USE_STATIC_RUNTIME OFF )
  find_package( Boost REQUIRED )
endif( TubeTK_USE_BOOST )

if( TubeTK_USE_KWSTYLE )
  include( TubeTKKWStyleConfig )
endif( TubeTK_USE_KWSTYLE )

if( TubeTK_USE_CPPCHECK )
  list( APPEND CMAKE_MODULE_PATH "${TubeTK_SOURCE_DIR}/CMake/Cppcheck" )
  find_package( Cppcheck REQUIRED )
  set( CPPCHECK_PATHS_TO_CHECK "${TubeTK_SOURCE_DIR}/Applications"
    "${TubeTK_SOURCE_DIR}/Base"
    "${TubeTK_SOURCE_DIR}/ITKModules" )
  set( CPPCHECK_SUPPRESSIONS_FILE
    "${TubeTK_SOURCE_DIR}/CMake/Cppcheck/TubeTKCppcheckSuppressions.txt" )
  set( CPPCHECK_COMMAND_PREFIX "${CPPCHECK_EXECUTABLE}" -q --xml
    --suppressions-list=${CPPCHECK_SUPPRESSIONS_FILE} )
  add_test( Cppcheck
    ${CPPCHECK_COMMAND_PREFIX} --enable=all ${CPPCHECK_PATHS_TO_CHECK} )
  set_tests_properties( Cppcheck
    PROPERTIES FAIL_REGULAR_EXPRESSION "error file=" )
endif( TubeTK_USE_CPPCHECK )

# Set SlicerExecutionModel output directories absolute path using install
#   directory relative path.
foreach( type RUNTIME LIBRARY INCLUDE SHARE ARCHIVE )
  set( SlicerExecutionModel_CLI_${type}_OUTPUT_DIRECTORY
    ${CMAKE_BINARY_DIR}/${SlicerExecutionModel_CLI_INSTALL_${type}_DESTINATION}
    )
endforeach()

configure_file ( "${TubeTK_SOURCE_DIR}/CMake/tubetkConfigure.h.in"
  tubetkConfigure.h @ONLY )

include_directories( ${TubeTK_BINARY_DIR} )

# Check endianness.
include( TestBigEndian )
test_big_endian( CMAKE_WORDS_BIGENDIAN )

#
# Setup Launcher
#

# Set TubeTK_EXECUTABLE_DIRS useful to configure launcher scripts
# and set TubeTK_LIBRARY_DIRS
if( WIN32 )
  set( TubeTK_EXECUTABLE_DIRS CACHE INTERNAL
    "Bin and Lib dirs for running apps." FORCE )
  foreach( _build_type "" Debug Release )
    list( APPEND TubeTK_EXECUTABLE_DIRS
      ${GenerateCLP_DIR}/${_build_type}
      ${ITK_DIR}/bin/${_build_type}
      ${TubeTK_BINARY_DIR}/bin/${_build_type}
      ${TubeTK_BINARY_DIR}/lib/TubeTK/${_build_type} )
    if( TubeTK_USE_VTK )
      list( APPEND TubeTK_EXECUTABLE_DIRS ${VTK_DIR}/bin/${_build_type} )
    endif( TubeTK_USE_VTK )
    if( TubeTK_BUILD_USING_SLICER )
      list( APPEND TubeTK_EXECUTABLE_DIRS ${Slicer_DIR}/bin/${_build_type} )
    endif( TubeTK_BUILD_USING_SLICER )
  endforeach()
else( WIN32 )
  set( TubeTK_EXECUTABLE_DIRS
    ${GenerateCLP_DIR}
    ${ITK_DIR}/bin
    ${TubeTK_BINARY_DIR}/bin
    ${TubeTK_BINARY_DIR}/lib/TubeTK
    CACHE INTERNAL "Bin and Lib dirs for running apps." FORCE )
  if( TubeTK_USE_VTK )
    list( APPEND TubeTK_EXECUTABLE_DIRS ${VTK_DIR}/bin )
  endif( TubeTK_USE_VTK )
  if( TubeTK_BUILD_USING_SLICER )
    list( APPEND TubeTK_EXECUTABLE_DIRS ${Slicer_DIR}/bin/${_build_type} )
  endif( TubeTK_BUILD_USING_SLICER )
endif( WIN32 )

# Boost may be dynamically linked, so...
if( TubeTK_USE_BOOST )
  list( APPEND TubeTK_EXECUTABLE_DIRS ${Boost_LIBRARY_DIR} )
endif( TubeTK_USE_BOOST )

# TubeTK uses a static lib for JsonCpp, but others (e.g., Slicer) may use
# a dynamic linked version, so we add it to our launcher path.
if( USE_SYSTEM_SLICER_EXECUTION_MODEL )
  get_filename_component( JsonCpp_LIBRARY_DIR ${JsonCpp_LIBRARY} DIRECTORY )
  list( APPEND TubeTK_EXECUTABLE_DIRS ${JsonCpp_LIBRARY_DIR} )
endif( USE_SYSTEM_SLICER_EXECUTION_MODEL )

#
# Configure a launcher for running TubeTK methods from the command line
#
message( STATUS "Configuring Launcher script" )
set( TubeTK_LAUNCHER )
if( WIN32 )
  # Microsoft Windows.
  set( _launcher_platform "windows" )
  configure_file( ${TubeTK_SOURCE_DIR}/CMake/TubeTKLauncher.bat.in
    ${TubeTK_BINARY_DIR}/TubeTKLauncher.bat @ONLY )

  find_program( CMD_EXECUTABLE "cmd" )
  if( NOT CMD_EXECUTABLE )
    message( FATAL_ERROR
      "Could not find 'cmd' executable required to run test with the launcher" )
  endif()

  set( TubeTK_LAUNCHER ${TubeTK_BINARY_DIR}/TubeTKLauncher.bat )
elseif( UNIX AND TubeTK_BUILD_USING_SLICER )
  # Apple OS X, Linux, or UNIX.
  set( _launcher_platform "unix" )
  configure_file( ${TubeTK_SOURCE_DIR}/CMake/TubeTKLauncher.sh.in
    ${TubeTK_BINARY_DIR}/TubeTKLauncher.sh @ONLY )

  find_program( SH_EXECUTABLE "sh" )
  if( NOT SH_EXECUTABLE )
    message( FATAL_ERROR
      "Could not find 'sh' executable required to run test using the launcher" )
  endif()
  set( TubeTK_LAUNCHER ${SH_EXECUTABLE} ${TubeTK_BINARY_DIR}/TubeTKLauncher.sh )
elseif( NOT UNIX )
  message( FATAL_ERROR
    "Configuring Launcher script - failed [unknown platform]" )
endif()
message( STATUS "Configuring Launcher script - ok [${_launcher_platform}]" )

include( TubeTKMacroInstallPlugins )

# Test TubeTK IPython examples.
if( TubeTK_USE_EXAMPLES_AS_TESTS )
  add_subdirectory( Examples )
endif( TubeTK_USE_EXAMPLES_AS_TESTS )

add_subdirectory( Utilities )
add_subdirectory( Base )
add_subdirectory( ITKModules )

if( TubeTK_BUILD_APPLICATIONS )
  add_subdirectory( Applications )
endif( TubeTK_BUILD_APPLICATIONS )

if( TubeTK_BUILD_USING_SLICER )
  add_subdirectory( SlicerModules )
endif( TubeTK_BUILD_USING_SLICER )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

# Create the directory where CTest unit tests store temporary results.
make_directory( ${TubeTK_BINARY_DIR}/Temporary )

# Configure variables to export.
set( TubeTK_INCLUDE_DIRS
  ${TubeTK_BINARY_DIR}
  ${TubeTK_SOURCE_DIR}/Base/CLI
  ${TubeTK_SOURCE_DIR}/Base/Common
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/IO
  ${TubeTK_SOURCE_DIR}/Base/MetaIO
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/ObjectDocuments
  ${TubeTK_SOURCE_DIR}/Base/Registration
  ${TubeTK_SOURCE_DIR}/Base/Segmentation
  ${TubeTK_SOURCE_DIR}/Base/USTK )

set( TubeTK_LIBRARY_DIRS
  ${TubeTK_EXECUTABLE_DIRS} )

# Configure TubeTKTargets.cmake.
get_property( TubeTK_TARGETS GLOBAL PROPERTY TubeTK_TARGETS )
export( TARGETS ${TubeTK_TARGETS} FILE
  ${TubeTK_CONFIG_BINARY_DIR}/TubeTKTargets.cmake )

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
if( NOT TubeTK_BUILD_USING_SLICER )
  include( TubeTKCPackConfig )
 else( NOT TubeTK_BUILD_USING_SLICER )
  include( ${Slicer_EXTENSION_CPACK} )
endif( NOT TubeTK_BUILD_USING_SLICER )
