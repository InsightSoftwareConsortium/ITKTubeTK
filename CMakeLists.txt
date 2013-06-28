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

include( ${CMAKE_CURRENT_SOURCE_DIR}/CMake/SlicerBlockSetCMakeOSXVariables.cmake )

project( TubeTK )

set( Slicer_SKIP_PROJECT_COMMAND ON )
#-------------------------------------------------------------------------------
# Require CMake 2.8 for its superbuild functionality.
cmake_minimum_required( VERSION 2.8.8 )
#-------------------------------------------------------------------------------
# Set TubeTK version and output directories.
set( TubeTK_VERSION_MAJOR "0" )
set( TubeTK_VERSION_MINOR "9" )
set( TubeTK_VERSION_PATCH "0" )
include( CMake/CMakeFunctionToday.cmake )
TODAY( TubeTK_VERSION_BUILD )
set( TubeTK_VERSION
  "${TubeTK_VERSION_MAJOR}.${TubeTK_VERSION_MINOR}.${TubeTK_VERSION_BUILD}" )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/TubeTK" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/TubeTK" )
#-------------------------------------------------------------------------------
# Setup for building 3D Slicer extensions.
set( SlicerExecutionModel_CLI_RUNTIME_INSTALL_DIRECTORY "bin" )
set( SlicerExecutionModel_CLI_LIBRARY_INSTALL_DIRECTORY "lib/TubeTK" )
set( SlicerExecutionModel_CLI_INCLUDE_INSTALL_DIRECTORY "include/TubeTK" )
set( SlicerExecutionModel_CLI_SHARE_INSTALL_DIRECTORY "share/TubeTK" )
set( SlicerExecutionModel_CLI_ARHCHIVE_INSTALL_DIRECTORY "lib/TubeTK" )

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
endif( TubeTK_BUILD_SLICER_EXTENSION )
#-------------------------------------------------------------------------------
# CDash dashboard setup.
include( ${CMAKE_ROOT}/Modules/Dart.cmake )
mark_as_advanced( TCL_TCLSH DART_ROOT )
enable_testing()
if( BUILD_TESTING )
  set( BUILDNAME "${BUILDNAME}" CACHE STRING
    "Name of the build on the CDash dashboard." )
  mark_as_advanced( BUILDNAME )
endif( BUILD_TESTING )
set( BUILD_WARNING_REPORT_LIMIT -1 )
set( BUILD_ERROR_REPORT_LIMIT -1 )
include( "CTestConfig.cmake" )
configure_file( ${TubeTK_SOURCE_DIR}/CMake/CTestCustom.cmake.in
  ${TubeTK_BINARY_DIR}/CTestCustom.cmake )
#-------------------------------------------------------------------------------
# Configure for the download of CTest unit testing data from a Midas Platform
# server.
set( MIDAS_USER "" CACHE STRING "Midas Platform user email address." )
mark_as_advanced( MIDAS_USER )
set( MIDAS_DEFAULT_API_KEY "" CACHE STRING "Midas Platform default API key." )
mark_as_advanced( MIDAS_DEFAULT_API_KEY )
set( MIDAS_DATA_DIR ${TubeTK_BINARY_DIR}/MIDAS_Data )
#-------------------------------------------------------------------------------
# Doxygen documentation setup.
option( BUILD_DOCUMENTATION
  "Build the software guide and source code documentation."
  OFF )
if( BUILD_DOCUMENTATION )
  option( TubeTK_USE_DOXYGEN
    "Use Doxygen to generate the source code documentation."
    ON )
  mark_as_advanced( TubeTK_USE_DOXYGEN )
  if( TubeTK_USE_DOXYGEN )
    add_subdirectory( CMake/Doxygen )
  endif( TubeTK_USE_DOXYGEN )
endif( BUILD_DOCUMENTATION )
#-------------------------------------------------------------------------------
#  Cppcheck setup.
option( TubeTK_USE_CPPCHECK
  "Use Cppcheck to perform a static analysis of the source code."
  OFF )
mark_as_advanced( TubeTK_USE_CPPCHECK )
if( TubeTK_USE_CPPCHECK )
  list( APPEND CMAKE_MODULE_PATH "${TubeTK_SOURCE_DIR}/CMake/Cppcheck" )
  include( FindCppcheck )
  find_package( Cppcheck REQUIRED )
  set( cppcheck_paths_to_check
    "${TubeTK_SOURCE_DIR}/Applications"
    "${TubeTK_SOURCE_DIR}/Base" )
  set( cppcheck_suppressions_file
    "${TubeTK_SOURCE_DIR}/CMake/Cppcheck/TubeTKCppcheck.supp" )
  set( cppcheck_cmd_prefix
    "${CPPCHECK_EXECUTABLE}" -q --xml --suppressions "${cppcheck_suppressions_file}" )
  add_test( CppcheckStyle
    ${cppcheck_cmd_prefix} --enable=style ${cppcheck_paths_to_check} )
  add_test( CppcheckPerformance
    ${cppcheck_cmd_prefix} --enable=performance ${cppcheck_paths_to_check} )
  add_test( CppcheckPortability
    ${cppcheck_cmd_prefix} --enable=portability ${cppcheck_paths_to_check} )
  add_test( CppcheckInformation
    ${cppcheck_cmd_prefix} --enable=information ${cppcheck_paths_to_check} )
  add_test( CppcheckUnusedFunction
    ${cppcheck_cmd_prefix} --enable=unusedFunction ${cppcheck_paths_to_check} )
  set_tests_properties( CppcheckStyle CppcheckPerformance CppcheckPortability
    CppcheckInformation CppcheckUnusedFunction
    PROPERTIES FAIL_REGULAR_EXPRESSION "error file=" )
endif( TubeTK_USE_CPPCHECK )
#-------------------------------------------------------------------------------
# KWStyle setup.
option( TubeTK_USE_KWSTYLE
  "Use KWStyle to check for coding standards violations."
  OFF )
mark_as_advanced( TubeTK_USE_KWSTYLE )
if( TubeTK_USE_KWSTYLE )
  find_program( KWSTYLE_EXECUTABLE
    NAMES KWStyle
    PATHS /usr/local/bin )
  mark_as_advanced( KWSTYLE_EXECUTABLE )

  option( KWSTYLE_USE_VIM_FORMAT
    "Set KWStyle to generate errors in a Vim compatible format."
    OFF )
  mark_as_advanced( KWSTYLE_USE_VIM_FORMAT )
  option( KWSTYLE_USE_MSVC_FORMAT
    "Set KWStyle to generate errors in a Microsoft Visual Studio compatible format."
    OFF )
  mark_as_advanced( KWSTYLE_USE_MSVC_FORMAT )
  option( KWSTYLE_DASHBOARD_SUBMISSION
    "Set the KWStyle arguments to generate a report for CDash dashboard submission."
    ON )
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
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
# Apple OS X setup.
if( APPLE )
  if( NOT CMAKE_FIRST_PASS_COMPLETED )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-common"
      CACHE STRING "Flags used by all build types." FORCE )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-common"
      CACHE STRING "Flags used by all build types." FORCE )
  endif( NOT CMAKE_FIRST_PASS_COMPLETED )
endif( APPLE )
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
# Implementation agnostic fPIC.
set( CMAKE_POSITION_INDEPENDENT_CODE ON )
#-------------------------------------------------------------------------------
# JsonCpp library setup.
option( USE_SYSTEM_JsonCpp
  "Exclude JsonCpp from the superbuild and use an existing build instead." OFF )
mark_as_advanced( USE_SYSTEM_JsonCpp )
if( USE_SYSTEM_JsonCpp )
  find_package( JsonCpp REQUIRED )
  include_directories( ${JsonCpp_INCLUDE_DIRS} )
endif( USE_SYSTEM_JsonCpp )
#-------------------------------------------------------------------------------
# LIBSVM library setup.
option( TubeTK_USE_LIBSVM
  "Use LIBSVM for support vector machine support to enable additional applications."
  OFF )
if( TubeTK_USE_LIBSVM )
  list( APPEND CMAKE_MODULE_PATH "${TubeTK_SOURCE_DIR}/CMake/LibSVM" )
  include( FindLibSVM )
  find_package( LibSVM REQUIRED )
endif( TubeTK_USE_LIBSVM )
#-------------------------------------------------------------------------------
# Boost libraries setup.
option( TubeTK_USE_Boost
  "Use the Boost libraries to enable additional applications." OFF )
#-------------------------------------------------------------------------------
# Parameter serializer setup.
option( USE_SYSTEM_ParameterSerializer
  "Exclude the parameter serializer from the superbuild and use an existing build instead."
  OFF )
mark_as_advanced( USE_SYSTEM_ParameterSerializer )
if( USE_SYSTEM_ParameterSerializer )
  find_package( ParameterSerializer REQUIRED )
  include_directories( ${ParameterSerializer_INCLUDE_DIRS} )
endif( USE_SYSTEM_ParameterSerializer )

if( NOT TubeTK_BUILD_SLICER_EXTENSION ) #-------------------------------------------------------------------------------
  # ITK setup.
  option( USE_SYSTEM_ITK
    "Exclude ITK from the superbuild and use an existing build instead." OFF )
  mark_as_advanced( USE_SYSTEM_ITK )
  if( USE_SYSTEM_ITK )
    find_package( ITK REQUIRED )
    include( ${ITK_USE_FILE} )
  endif( USE_SYSTEM_ITK )
  option( TubeTK_USE_SimpleITK
    "Use SimpleITK to enable additional applications and libraries."
    OFF ) #-------------------------------------------------------------------------------
  # Qt framework setup.
  option( TubeTK_USE_QT
    "Use the Qt framework to enable additional applications and libraries."
    ON )

  if( TubeTK_USE_QT )
    find_package( Qt4 )
    # 3D Slicer only tests against a specific Qt 4 version. To
    # override the version check, set the TubeTK_REQUIRED_QT_VERSION CMake
    # variable.
    if( NOT DEFINED TubeTK_REQUIRED_QT_VERSION )
      set( TubeTK_REQUIRED_QT_VERSION
        "4.7.4" CACHE STRING "Exact Qt version required." FORCE )
      mark_as_advanced( TubeTK_REQUIRED_QT_VERSION )
    endif( NOT DEFINED TubeTK_REQUIRED_QT_VERSION )
    if( QT4_FOUND )
      if( NOT "${QT_VERSION_MAJOR}.${QT_VERSION_MINOR}.${QT_VERSION_PATCH}"
        VERSION_EQUAL "${TubeTK_REQUIRED_QT_VERSION}" )
        message( FATAL_ERROR
          "Error: TubeTK requires Qt == ${TubeTK_REQUIRED_QT_VERSION}" )
      endif( NOT "${QT_VERSION_MAJOR}.${QT_VERSION_MINOR}.${QT_VERSION_PATCH}"
        VERSION_EQUAL "${TubeTK_REQUIRED_QT_VERSION}" )
      include( ${QT_USE_FILE} )
    else( QT4_FOUND )
      message( SEND_ERROR "Error: Qt not found." )
    endif( QT4_FOUND )
  endif( TubeTK_USE_QT ) #-------------------------------------------------------------------------------
  # VTK setup.
  option( TubeTK_USE_VTK
    "Use VTK to enable additional applications and libraries."
    ON )
  if( TubeTK_USE_VTK )
    option( USE_SYSTEM_VTK
      "Exclude VTK from the superbuild and use an existing build instead." OFF )
    mark_as_advanced( USE_SYSTEM_VTK )
    if( USE_SYSTEM_VTK )
      find_package( VTK REQUIRED )
      include( ${VTK_USE_FILE} )
    endif( USE_SYSTEM_VTK )
  endif( TubeTK_USE_VTK ) #-------------------------------------------------------------------------------
  # CTK setup.
  if( TubeTK_USE_QT )
    option( TubeTK_USE_CTK
      "Use CTK to enable additional applications and build the 3D Slicer modules."
      ON )
    option( USE_SYSTEM_CTK
      "Exclude CTK from the superbuild and use an existing build instead." OFF )
    mark_as_advanced( USE_SYSTEM_CTK )
    if( USE_SYSTEM_CTK )
      find_package( CTK REQUIRED )
      include( ${CTK_USE_FILE} )
    endif( USE_SYSTEM_CTK )
  endif( TubeTK_USE_QT )
else( NOT TubeTK_BUILD_SLICER_EXTENSION )
  unset( ITK_DIR )
  unset( VTK_DIR )
  unset( CTK_DIR )
  unset( QT_QMAKE_EXECUTABLE )
  unset( SlicerExecutionModel_DIR )
  unset( TubeTK_USE_SUPERBUILD )
endif( NOT TubeTK_BUILD_SLICER_EXTENSION )
#-------------------------------------------------------------------------------
# Image viewer setup.
option( TubeTK_BUILD_ImageViewer "Build the lightweight image viewer." OFF )
#-------------------------------------------------------------------------------
# Indicate that this CMakeLists.txt file has already been parsed once. Keeps
# append commands from having an effect every time CMake is run.
set( CMAKE_FIRST_PASS_COMPLETED ON
  CACHE BOOL
  "Modify the cache variables only on the first pass through the CMake file."
  FORCE )
mark_as_advanced( CMAKE_FIRST_PASS_COMPLETED )
#-------------------------------------------------------------------------------
# Superbuild setup.
set( superbuild_enabled_by_default ON )
if( TubeTK_BUILD_SLICER_EXTENSION )
  set( superbuild_enabled_by_default OFF )
endif( TubeTK_BUILD_SLICER_EXTENSION )

option( TubeTK_USE_SUPERBUILD
  "Build TubeTK and the projects it depends on using a superbuild."
  ${superbuild_enabled_by_default} )
mark_as_advanced( TubeTK_USE_SUPERBUILD )
if( TubeTK_USE_SUPERBUILD )
  if( TubeTK_BUILD_SLICER_EXTENSION )
    find_package( Slicer REQUIRED )
    include( ${Slicer_USE_FILE} )
  endif( TubeTK_BUILD_SLICER_EXTENSION )

  include( "${CMAKE_CURRENT_SOURCE_DIR}/CMake/SuperBuild.cmake" )
  return()
else( TubeTK_USE_SUPERBUILD )
  find_package( JsonCpp REQUIRED )
  include_directories( ${JsonCpp_INCLUDE_DIRS} )

  find_package( ParameterSerializer )
  include_directories( ${ParameterSerializer_INCLUDE_DIRS} )

  if( TubeTK_BUILD_SLICER_EXTENSION )
    find_package( Slicer REQUIRED )
    include( ${Slicer_USE_FILE} )

    set( TubeTK_USE_VTK ON )
    set( TubeTK_USE_QT ON )
    set( TubeTK_USE_CTK ON )
  else( TubeTK_BUILD_SLICER_EXTENSION )
    find_package( ITK REQUIRED )
    include( ${ITK_USE_FILE} )

    if( TubeTK_USE_VTK )
      find_package( VTK REQUIRED )
      include( ${VTK_USE_FILE} )
    endif( TubeTK_USE_VTK )

    find_package( SlicerExecutionModel REQUIRED GenerateCLP )
    include( ${GenerateCLP_USE_FILE} )

    if( TubeTK_USE_QT )
      if( TubeTK_USE_CTK )
        find_package( CTK REQUIRED )
        include( ${CTK_USE_FILE} )
      endif( TubeTK_USE_CTK )
    endif( TubeTK_USE_QT )
  endif( TubeTK_BUILD_SLICER_EXTENSION )
endif( TubeTK_USE_SUPERBUILD )
#-------------------------------------------------------------------------------
# Check endianness.
include( ${CMAKE_ROOT}/Modules/TestBigEndian.cmake )
test_big_endian( CMAKE_WORDS_BIGENDIAN )
#-------------------------------------------------------------------------------
# Setup to conform to 3D Slicer paths and define plugin copying macros.
include( CMake/Slicer4Config.cmake )
#-------------------------------------------------------------------------------
# Build 3D Slicer command line applications essential to TubeTK.
svnGetSlicerBase( Resources )
svnGetSlicerBase( Base/CLI )
svnGetSlicerBase( Libs/ITKFactoryRegistration )
svnGetSlicerCLI( ExpertAutomatedRegistration )

set( CXX_TEST_PATH  ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR} )
set( LAUNCH_EXE  "" )
set( SUPPORT_2D_IMAGES on )

if( NOT TubeTK_BUILD_SLICER_EXTENSION )
  add_subdirectory( ${Slicer_SOURCE_DIR}/Libs/ITKFactoryRegistration
    ${Slicer_BINARY_DIR}/Libs/ITKFactoryRegistration )
  include_directories( ${Slicer_BINARY_DIR}/Libs/itkFactoryRegistration )
  add_subdirectory( ${Slicer_SOURCE_DIR}/Base/CLI
    ${Slicer_BINARY_DIR}/Base/CLI )

  set( SlicerExecutionModel_EXTRA_EXECUTABLE_TARGET_LIBRARIES
    SlicerBaseCLI
    ITKFactoryRegistration
    CACHE INTERNAL "3D Slicer execution model extra libraries." FORCE )
endif( NOT TubeTK_BUILD_SLICER_EXTENSION )

#-------------------------------------------------------------------------------
# Build TubeTK.
add_subdirectory( Utilities )
add_subdirectory( Base )
add_subdirectory( Applications )
#-------------------------------------------------------------------------------
# Optionally build 3D Slicer modules.
if( TubeTK_BUILD_SLICER_EXTENSION )
  add_subdirectory( SlicerModules )
endif( TubeTK_BUILD_SLICER_EXTENSION )
#-------------------------------------------------------------------------------
# Create the directory where CTest unit tests store temporary results.
make_directory( ${TubeTK_BINARY_DIR}/Temporary )
#-------------------------------------------------------------------------------
# Configure variables to export.
set( TubeTK_INCLUDE_DIRS
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/IO
  ${TubeTK_SOURCE_DIR}/Base/Segmentation
  ${TubeTK_SOURCE_DIR}/Base/Registration
  ${TubeTK_SOURCE_DIR}/Base/USTK )

if( MSVC )
  set( TubeTK_EXECUTABLE_DIRS
    ${GenerateCLP_DIR}/${CMAKE_BUILD_TYPE}
    ${ITK_DIR}/bin/${CMAKE_BUILD_TYPE}
    ${VTK_DIR}/bin/${CMAKE_BUILD_TYPE}
    ${TubeTK_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE}
    ${TubeTK_BINARY_DIR}/lib/TubeTK/Plugins/${CMAKE_BUILD_TYPE} )
else( MSVC )
  set( TubeTK_EXECUTABLE_DIRS
    ${GenerateCLP_DIR}
    ${ITK_DIR}/bin
    ${VTK_DIR}/bin
    ${TubeTK_BINARY_DIR}/bin
    ${TubeTK_BINARY_DIR}/lib
    ${TubeTK_BINARY_DIR}/lib/TubeTK/Plugins )
endif( MSVC )

set( TubeTK_LIBRARY_DIRS
  ${TubeTK_EXECUTABLE_DIRS} )

set( TubeTK_LIBS
  TubeIO
  TubeNumerics )
set( TubeTK_LIBRARIES )
foreach( library ${TubeTK_LIBS} )
  get_target_property( library_location ${library} LOCATION )
  list( APPEND TubeTK_LIBRARIES ${library_location} )
endforeach( library ${TubeTK_LIBS} )

set( TubeTK_CMAKE_EXTENSIONS_DIR
  ${TubeTK_SOURCE_DIR}/CMake )

set( TubeTK_USE_FILE
  ${TubeTK_BINARY_DIR}/UseTubeTK.cmake )
#-------------------------------------------------------------------------------
# Configure TubeTKConfig.cmake.
# The variables set above are for building against a TubeTK build tree.
# A different TubeTKConfig.cmake will need to be generated in order to build
# against a TubeTK installation.
configure_file( ${TubeTK_SOURCE_DIR}/CMake/TubeTKConfig.cmake.in
  ${TubeTK_BINARY_DIR}/TubeTKConfig.cmake @ONLY )
#-------------------------------------------------------------------------------
# Configure UseTubeTK.cmake.
configure_file( ${TubeTK_SOURCE_DIR}/CMake/UseTubeTK.cmake.in
  ${TubeTK_BINARY_DIR}/UseTubeTK.cmake @ONLY )
#-------------------------------------------------------------------------------
# Configure a launcher for running TubeTK methods from the command line in
# Microsoft Windows.
configure_file( ${TubeTK_SOURCE_DIR}/CMake/tubeLauncher.bat.in
  ${TubeTK_BINARY_DIR}/tubeLauncher.bat @ONLY )
#-------------------------------------------------------------------------------
# Configure a launcher for running TubeTK methods from the command line in
# Apple OS X, Linux, or UNIX.
configure_file( ${TubeTK_SOURCE_DIR}/CMake/tubeLauncher.sh.in
  ${TubeTK_BINARY_DIR}/tubeLauncher.sh @ONLY )
#-------------------------------------------------------------------------------
# Configure CPack packaging.
if( NOT TubeTK_BUILD_SLICER_EXTENSION )
  # Record ITK, VTK, CTK, and GenerateCLP libraries used by TubeTK for
  # packaging.
  file( GLOB_RECURSE ITK_BINARIES "${ITK_DIR}/bin/lib*" )
  file( GLOB_RECURSE VTK_BINARIES "${VTK_DIR}/bin/lib*" )
  file( GLOB_RECURSE CTK_BINARIES "${CTK_DIR}/CTK-build/bin/lib*" )
  file( GLOB_RECURSE SEM_BINARIES "${GenerateCLP_DIR}/bin/lib*" )
  list( APPEND TubeTK_SUPPORT_BINARIES
    ${SEM_BINARIES} ${CTK_BINARIES} ${VTK_BINARIES} ${ITK_BINARIES} )
  install( FILES ${TubeTK_SUPPORT_BINARIES}
    DESTINATION lib )

  # Record the environment setup files for packaging.
  set( TubeTK_LINUX_SETUP
    ${TubeTK_SOURCE_DIR}/CMake/tubeSourceMeToSetupLinux.sh )
  set( TubeTK_WINDOWS_SETUP
    ${TubeTK_SOURCE_DIR}/CMake/tubeRunMeToSetupWindows.bat )

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
      "${TubeTK_SOURCE_DIR}/Documentation\\\\TubeTK_Icon.bmp" )
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
else( NOT TubeTK_BUILD_SLICER_EXTENSION )
  include( ${Slicer_EXTENSION_CPACK} )
endif( NOT TubeTK_BUILD_SLICER_EXTENSION )
