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

include( ${CMAKE_CURRENT_SOURCE_DIR}/CMake/SlicerBlockSetCMakeOSXVariables.cmake )

project( TubeTK )

set( Slicer_SKIP_PROJECT_COMMAND ON )

#-----------------------------------------------------------------------------
# Require CMake 2.8 for the superbuild functionality
cmake_minimum_required( VERSION 2.8.8 )

#-----------------------------------------------------------------------------
#  Set TubeTK
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

#-----------------------------------------------------------------------------
#  Setup for building extensions
set( SlicerExecutionModel_CLI_RUNTIME_INSTALL_DIRECTORY "bin" )
set( SlicerExecutionModel_CLI_LIBRARY_INSTALL_DIRECTORY "lib/TubeTK" )
set( SlicerExecutionModel_CLI_INCLUDE_INSTALL_DIRECTORY "include/TubeTK" )
set( SlicerExecutionModel_CLI_SHARE_INSTALL_DIRECTORY "share/TubeTK" )
set( SlicerExecutionModel_CLI_ARHCHIVE_INSTALL_DIRECTORY "lib/TubeTK" )

set( SlicerExecutionModel_CLI_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
set( SlicerExecutionModel_CLI_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/TubeTK" )
set( SlicerExecutionModel_CLI_INCLUDE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/include/TubeTK" )
set( SlicerExecutionModel_CLI_SHARE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/share/TubeTK" )
set( SlicerExecutionModel_CLI_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/TubeTK" )

option( TubeTK_BUILD_SLICER_EXTENSION
        "Build TubeTK as a Slicer extension" OFF )
if( TubeTK_BUILD_SLICER_EXTENSION )
  set( EXTENSION_NAME TubeTK )
  set( EXTENSION_HOMEPAGE "http://www.tubetk.org" )
  set( EXTENSION_CATEGORY "TubeTK" )
  set( EXTENSION_CONTRIBUTORS "TubeTK Team @ Kitware" )
  set( EXTENSION_DESCRIPTION
    "Using local image geometry for segmentation, registration, and filtering" )
  set( EXTENSION_ICONURL
    "http://public.kitware.com/Wiki/images/a/ad/TubeTK_Icon.jpg" )
  set( EXTENSION_SCREENSHOTURLS
    "http://public.kitware.com/Wiki/images/f/fd/VesselExtraction.jpg" )
endif( TubeTK_BUILD_SLICER_EXTENSION )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/Modules")

#-----------------------------------------------------------------------------
# Dashboard setup
include( ${CMAKE_ROOT}/Modules/Dart.cmake )
mark_as_advanced( TCL_TCLSH DART_ROOT )
enable_testing()
if( BUILD_TESTING )
  set( BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard" )
  mark_as_advanced( BUILDNAME )
endif()
set( BUILD_WARNING_REPORT_LIMIT -1 )
set( BUILD_ERROR_REPORT_LIMIT -1 )
include( "CTestConfig.cmake" )
configure_file( ${TubeTK_SOURCE_DIR}/CMake/CTestCustom.cmake.in
  ${TubeTK_BINARY_DIR}/CTestCustom.cmake )

#-----------------------------------------------------------------------------
# Configure MIDAS_USER and API_KEY to download testing data from MIDAS
set( MIDAS_USER "" CACHE STRING "MIDAS User email address" )
mark_as_advanced( MIDAS_USER )
set( MIDAS_DEFAULT_API_KEY "" CACHE STRING "MIDAS DEFAULT API Key" )
mark_as_advanced( MIDAS_DEFAULT_API_KEY )
set( MIDAS_DATA_DIR ${TubeTK_BINARY_DIR}/MIDAS_Data )

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
    PATHS /usr/local/bin )
  mark_as_advanced( KWSTYLE_EXECUTABLE )

  option( KWSTYLE_USE_VIM_FORMAT
    "Set KWStyle to generate errors with a VIM-compatible format."
    OFF )
  mark_as_advanced( KWSTYLE_USE_VIM_FORMAT )
  option( KWSTYLE_USE_MSVC_FORMAT
    "Set KWStyle to generate errors with a VisualStudio-compatible format."
    OFF )
  mark_as_advanced( KWSTYLE_USE_MSVC_FORMAT )
  option( KWSTYLE_DASHBOARD_SUBMISSION
    "Set KWStyle arguments to generate report for dashboard submission."
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
  endif()


  if( KWSTYLE_USE_VIM_FORMAT )
    set( KWSTYLE_ARGUMENTS -vim ${KWSTYLE_ARGUMENTS} )
  endif()

  if( KWSTYLE_USE_MSVC_FORMAT )
    set( KWSTYLE_ARGUMENTS -msvc ${KWSTYLE_ARGUMENTS} )
  endif()

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
endif()

#-----------------------------------------------------------------------------
# MSVC setup
if( MSVC )
  add_definitions( -D_SCL_SECURE_NO_WARNINGS )
  if( NOT CMAKE_FIRST_PASS_COMPLETED AND "${CMAKE_CXX_COMPILER}" MATCHES "cl" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /MP"
         CACHE STRING "Flags used by all build types" FORCE )
  endif()
  option( DISABLE_MSVC8_DEPRECATED_WARNINGS
         "Disable Visual Studio 8 deprecated warnings" ON )
  mark_as_advanced( FORCE DISABLE_MSVC8_DEPRECATED_WARNINGS )
  if( DISABLE_MSVC8_DEPRECATED_WARNINGS )
    add_definitions( -D_CRT_SECURE_NO_DEPRECATE )
  endif()
endif( MSVC )

#-----------------------------------------------------------------------------
# Mac setup
if( APPLE )
  if( NOT CMAKE_FIRST_PASS_COMPLETED )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-common"
         CACHE STRING "Flags used by all build types" FORCE )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-common"
         CACHE STRING "Flags used by all build types" FORCE )
  endif()
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
    endif()
  endif()
endif()

#-----------------------------------------------------------------------------
# JsonCpp Setup
option( USE_SYSTEM_JsonCpp
  "Exclude JsonCpp from SuperBuild and use and existing build instead" OFF )
mark_as_advanced( USE_SYSTEM_JsonCpp )
if( USE_SYSTEM_JsonCpp )
  find_package( JsonCpp REQUIRED )
  include_directories( ${JsonCpp_INCLUDE_DIRS} )
endif()

#-----------------------------------------------------------------------------
# libSVM Setup
option( TubeTK_USE_LIBSVM
  "Build libSVM for Support Vector Machine support."
  OFF )

#-----------------------------------------------------------------------------
# Boost setup
option( TubeTK_USE_Boost
  "Use Boost libraries."
  OFF )

if( NOT TubeTK_BUILD_SLICER_EXTENSION )
  #-----------------------------------------------------------------------------
  # ITK Setup
  option( USE_SYSTEM_ITK
    "Exclude ITK from SuperBuild and use an existing build instead." OFF )
  mark_as_advanced( USE_SYSTEM_ITK )
  if( USE_SYSTEM_ITK )
    find_package( ITK REQUIRED )
    include( ${ITK_USE_FILE} )
  endif()

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
    # Slicer only tests against a specific Qt version.  To override the
    # version check, set the TubeTK_REQUIRED_QT_VERSION CMake variable.
    if( NOT DEFINED TubeTK_REQUIRED_QT_VERSION )
      set( TubeTK_REQUIRED_QT_VERSION
        "4.7.4" CACHE STRING "Exact required Qt version" FORCE)
    mark_as_advanced( TubeTK_REQUIRED_QT_VERSION )
    endif()
    if( QT4_FOUND )
      if( NOT "${QT_VERSION_MAJOR}.${QT_VERSION_MINOR}.${QT_VERSION_PATCH}"
          VERSION_EQUAL "${TubeTK_REQUIRED_QT_VERSION}" )
        message( FATAL_ERROR
          "Error: TubeTK requires Qt == ${TubeTK_REQUIRED_QT_VERSION}" )
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
    mark_as_advanced( USE_SYSTEM_VTK )
    if( USE_SYSTEM_VTK )
      find_package( VTK REQUIRED )
      include( ${VTK_USE_FILE} )
    endif()
  endif()

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
    mark_as_advanced( USE_SYSTEM_CTK )
    if( USE_SYSTEM_CTK )
      find_package( CTK REQUIRED )
      include( ${CTK_USE_FILE} )
    endif()
  endif()

else( NOT TubeTK_BUILD_SLICER_EXTENSION )

  unset( ITK_DIR )
  unset( VTK_DIR )
  unset( CTK_DIR )
  unset( QT_QMAKE_EXECUTABLE )
  unset( SlicerExecutionModel_DIR )
  unset( TubeTK_USE_SUPERBUILD )

endif( NOT TubeTK_BUILD_SLICER_EXTENSION )

#-----------------------------------------------------------------------------
# Indicate that this cmakelists.txt file has already been parsed once.
# Keeps append commands from having an effect every time cmake is run.
set( CMAKE_FIRST_PASS_COMPLETED ON
  CACHE BOOL "Modify the cache vars only on first pass thru cmake file" FORCE )
mark_as_advanced( CMAKE_FIRST_PASS_COMPLETED )

#-----------------------------------------------------------------------------
# Superbuild setup

set( superbuild_enabled_by_default ON )
if( TubeTK_BUILD_SLICER_EXTENSION )
  set( superbuild_enabled_by_default OFF )
endif()

option( TubeTK_USE_SUPERBUILD
  "Build TubeTK and the projects it depends on via SuperBuild.cmake."
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

  if( TubeTK_BUILD_SLICER_EXTENSION )

    find_package( Slicer REQUIRED )
    include( ${Slicer_USE_FILE} )

    set( TubeTK_USE_ITKV4 OFF )
    if( ${ITK_VERSION_MAJOR} STREQUAL "4" )
      set( TubeTK_USE_ITKV4 ON )
    endif()

    set( TubeTK_USE_VTK ON )
    set( TubeTK_USE_QT ON )
    set( TubeTK_USE_CTK ON )

  else( TubeTK_BUILD_SLICER_EXTENSION )

    find_package( ITK REQUIRED )
    include( ${ITK_USE_FILE} )

    if( TubeTK_USE_VTK )
      find_package( VTK REQUIRED )
      include( ${VTK_USE_FILE} )
    endif()

    find_package( SlicerExecutionModel REQUIRED GenerateCLP )
    include( ${GenerateCLP_USE_FILE} )

    if( TubeTK_USE_QT )
      if( TubeTK_USE_CTK )
        find_package( CTK REQUIRED )
        include( ${CTK_USE_FILE} )
      endif()
    endif()
  endif( TubeTK_BUILD_SLICER_EXTENSION )

endif()


#-----------------------------------------------------------------------------
# Misc.
include( ${CMAKE_ROOT}/Modules/TestBigEndian.cmake )
test_big_endian( CMAKE_WORDS_BIGENDIAN )

#-----------------------------------------------------------------------------
# Setup to conform to Slicer3 paths and define plugin copying macros
if( NOT TubeTK_USE_ITKV4 )
  include( CMake/Slicer3Config.cmake )

  #---------------------------------------------------------------------------
  # Build Slicer3 cli applications essential to TubeTK
  svnGetSlicerBase( Resources ${CLI_SOURCE_DIR} )
  svnGetSlicerCLI( Templates ${CLI_SOURCE_DIR} )
  svnGetSlicerCLI( RegisterImagesModule ${CLI_SOURCE_DIR} )

  set( CXX_TEST_PATH  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR} )
  set( LAUNCH_EXE  "" )

  set( SUPPORT_2D_IMAGES ON )
  add_subdirectory( ${CLI_SOURCE_DIR}/RegisterImagesModule
    ${CLI_BINARY_DIR}/RegisterImagesModule )
endif()

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Build TubeTK
add_subdirectory( Utilities )
add_subdirectory( Base )
add_subdirectory( Applications )

#-----------------------------------------------------------------------------
# Optionally build Slicer-specific modules
if( TubeTK_BUILD_SLICER_EXTENSION )
  add_subdirectory( SlicerModules )
endif( TubeTK_BUILD_SLICER_EXTENSION )

#-----------------------------------------------------------------------------
# Create the directory where tests store temp results
make_directory( ${TubeTK_BINARY_DIR}/Temporary )

#-----------------------------------------------------------------------------
# Configure vars to export

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
else()
 set( TubeTK_EXECUTABLE_DIRS
   ${GenerateCLP_DIR}
   ${ITK_DIR}/bin
   ${VTK_DIR}/bin
   ${TubeTK_BINARY_DIR}/bin
   ${TubeTK_BINARY_DIR}/lib
   ${TubeTK_BINARY_DIR}/lib/TubeTK/Plugins )
endif()

set( TubeTK_LIBRARY_DIRS
  ${TubeTK_EXECUTABLE_DIRS} )

set( TubeTK_LIBS
  TubeIO
  TubeNumerics )
set( TubeTK_LIBRARIES )
foreach( library ${TubeTK_LIBS} )
  get_target_property( library_location ${library} LOCATION )
  list( APPEND TubeTK_LIBRARIES ${library_location} )
endforeach()

set( TubeTK_CMAKE_EXTENSIONS_DIR
  ${TubeTK_SOURCE_DIR}/CMake )

set( TubeTK_USE_FILE
  ${TubeTK_BINARY_DIR}/UseTubeTK.cmake )

#-----------------------------------------------------------------------------
# Configure TubeTKConfig.cmake
# The variables set above are for building against a TubeTK build tree.
# A different TubeTKConfig.cmake will need to be generated in order to build
# against a TubeTK installation.

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

#-----------------------------------------------------------------------------
# Configure packaging

if( NOT TubeTK_BUILD_SLICER_EXTENSION )

  # Record itk, vtk, ctk, and GenerateCLP libs used by TubeTK for packaging
  file( GLOB_RECURSE ITK_BINARIES "${ITK_DIR}/bin/lib*" )
  file( GLOB_RECURSE VTK_BINARIES "${VTK_DIR}/bin/lib*" )
  file( GLOB_RECURSE CTK_BINARIES "${CTK_DIR}/CTK-build/bin/lib*" )
  file( GLOB_RECURSE SEM_BINARIES "${GenerateCLP_DIR}/bin/lib*" )
  list( APPEND TubeTK_SUPPORT_BINARIES
    ${SEM_BINARIES} ${CTK_BINARIES} ${VTK_BINARIES} ${ITK_BINARIES} )
  install( FILES ${TubeTK_SUPPORT_BINARIES}
    DESTINATION lib )

  # Record the environment setup files for packaging
  set( TubeTK_LINUX_SETUP
    ${TubeTK_SOURCE_DIR}/CMake/tubeSourceMeToSetupLinux.sh )
  set( TubeTK_WINDOWS_SETUP
    ${TubeTK_SOURCE_DIR}/CMake/tubeRunMeToSetupWindows.bat )

  include( InstallRequiredSystemLibraries )

  set( CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "TubeTK for the analysis of geometric structures in images." )
  set( CPACK_PACKAGE_VENDOR "Kitware, Inc." )
  set( CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt" )
  set( CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt" )

  set( CPACK_PACKAGE_VERSION_MAJOR ${TubeTK_VERSION_MAJOR} )
  set( CPACK_PACKAGE_VERSION_MINOR ${TubeTK_VERSION_MINOR} )
  set( CPACK_PACKAGE_VERSION_PATCH ${TubeTK_VERSION_PATCH} )
  set( CPACK_PACKAGE_VERSION ${TubeTK_VERSION_BUILD} )

  set( CPACK_PACKAGE_CONTACT "stephen.aylward@kitware.com" )

  set( CPACK_PACKAGE_INSTALL_DIRECTORY
    "TubeTK_${TubeTK_VERSION_MAJOR}.${TubeTK_VERSION_MINOR}" )

  if( WIN32 AND NOT UNIX )
    # There is a bug in NSI that does not handle full unix paths properly. Make
    # sure there is at least one set of four ( 4 ) backlasshes.
    set( CPACK_PACKAGE_ICON
      "${TubeTK_SOURCE_DIR}/Documentation\\\\TubeTK_Icon.bmp" )
    set( CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\TubeTK.exe" )
    set( CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} TubeTK" )
    set( CPACK_NSIS_HELP_LINK "http:\\\\\\\\tubetk.org" )
    set( CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\tubetk.org" )
    set( CPACK_NSIS_CONTACT "stephen.aylward@kitware.com" )
    set( CPACK_NSIS_MODIFY_PATH ON )
  else()
    set( CPACK_STRIP_FILES true )
    set( CPACK_SOURCE_STRIP_FILES true )
  endif()

  set( CMAKE_INSTALL_RPATH "$ORIGIN/../..;." )

  #
  # Standard CPack setup
  OPTION( CMAKE_INSTALL_DEBUG_LIBRARIES
    "Install Microsoft runtime debug libraries with CMake." FALSE )
  MARK_AS_ADVANCED( CMAKE_INSTALL_DEBUG_LIBRARIES )

  # By default, do not warn when built on machines using only VS Express:
  if( NOT DEFINED CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS )
    set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS ON )
  endif()

  # Make this explicit here, rather than accepting the CPack default value,
  # so we can refer to it:
  set( CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}" )

  # Installers for 32- vs. 64-bit CMake:
  #  - Root install directory ( displayed to end user at installer-run time )
  #  - "NSIS package/display name" ( text used in the installer GUI )
  #  - Registry key used to store info about the installation
  if( CMAKE_CL_64 )
    set( CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64" )
    set( CPACK_NSIS_PACKAGE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} ( Win64 )" )
    set( CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION} ( Win64 )" )
  else()
    set( CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES" )
    set( CPACK_NSIS_PACKAGE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}" )
    set( CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION}" )
  endif()

  if( NOT DEFINED CPACK_SYSTEM_NAME )
    set( CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR} )
  endif()

  if( ${CPACK_SYSTEM_NAME} MATCHES Windows )
    if( CMAKE_CL_64 )
      set( CPACK_SYSTEM_NAME win64-x64 )
    else()
      set( CPACK_SYSTEM_NAME win32-x86 )
    endif()
  endif()

  set( CPACK_SOURCE_PACKAGE_FILE_NAME "TubeTK-${TubeTK_VERSION_BUILD}" )
  set( CPACK_PACKAGE_FILE_NAME
    "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}" )

  # include CPack model once all variables are set
  INCLUDE( CPack )
else( NOT TubeTK_BUILD_SLICER_EXTENSION )
  include( ${Slicer_EXTENSION_CPACK} )
endif( NOT TubeTK_BUILD_SLICER_EXTENSION )
