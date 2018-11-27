##############################################################################
#
# Library:   TubeTK
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

cmake_policy( SET CMP0070 NEW )
cmake_policy( SET CMP0071 NEW )

# Include TubeTKLib
if( NOT ITK_SOURCE_DIR )
  include( itk-module-init.cmake )
endif()

####
####
####
####
project( TubeTK )

list( APPEND CMAKE_MODULE_PATH ${TubeTK_SOURCE_DIR}/CMake )

include( TubeTKExternalData )
if( DEFINED "ENV{ExternalData_OBJECT_STORES}" )
  set( ExternalData_OBJECT_STORES $ENV{ExternalData_OBJECT_STORES} )
endif()


####
# Default output directories
####
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" )

mark_as_advanced( BUILD_SHARED_LIBS )
mark_as_advanced( CMAKE_BACKWARDS_COMPATIBILITY )
mark_as_advanced( EXECUTABLE_OUTPUT_PATH )
mark_as_advanced( LIBRARY_OUTPUT_PATH )
mark_as_advanced( CMAKE_INSTALL_PREFIX )
mark_as_advanced( BUILD_TESTING )
mark_as_advanced( BUILD_DOCUMENTATION )

####
# CDash dashboard setup.
####
include( CTest )

make_directory( ${TubeTK_BINARY_DIR}/Temporary )

set( BUILDNAME "${BUILDNAME}" CACHE STRING
  "Name of the build on the CDash dashboard." )
mark_as_advanced( BUILDNAME )

set( BUILD_WARNING_REPORT_LIMIT -1 )
set( BUILD_ERROR_REPORT_LIMIT -1 )

include( CTestConfig.cmake )

configure_file( ${TubeTK_SOURCE_DIR}/CMake/CTestCustom.cmake.in
  ${TubeTK_BINARY_DIR}/CTestCustom.cmake )

set( TubeTK_BUILD_TESTING ${BUILD_TESTING} )

include( TubeTKExternalData )
if( DEFINED "ENV{ExternalData_OBJECT_STORES}" )
  set( ExternalData_OBJECT_STORES $ENV{ExternalData_OBJECT_STORES} )
endif()

####
# Setup TubeTKLib definition since UseTubeTKLib.cmake won't be available
# until after superbuild of TubeTKLib has completed.
####
set( TubeTKLib_DIR
  ${TubeTK_BINARY_DIR}/TubeTKLib
  )
list( APPEND CMAKE_PREFIX_PATH ${TubeTKLib_DIR} )

set( TubeTKLib_LIBRARIES
  TubeTKLibMetaIO
  TubeTKLibCommon
  TubeTKLibNumerics
  )

set( TubeTKLib_LIB_DIRS
  ${TubeTK_BINARY_DIR}/TubeTKLib/Lib
  )

set( TubeTKLib_CompareTextFiles_EXE
  ${TubeTK_BINARY_DIR}/TubeTKLib/bin/CompareTextFiles
  )

set( TubeTKLib_CompareImages_EXE
  ${TubeTK_BINARY_DIR}/TubeTKLib/bin/CompareImages
  )

set( TubeTKLib_INCLUDE_DIRS
  ${TubeTK_BINARY_DIR}/TubeTKLib
  ${TubeTK_SOURCE_DIR}/TubeTKLib/Common
  ${TubeTK_SOURCE_DIR}/TubeTKLib/Filtering
  ${TubeTK_SOURCE_DIR}/TubeTKLib/IO
  ${TubeTK_SOURCE_DIR}/TubeTKLib/MetaIO
  ${TubeTK_SOURCE_DIR}/TubeTKLib/Numerics
  ${TubeTK_SOURCE_DIR}/TubeTKLib/ObjectDocuments
  ${TubeTK_SOURCE_DIR}/TubeTKLib/Registration
  ${TubeTK_SOURCE_DIR}/TubeTKLib/Segmentation
  ${TubeTK_BINARY_DIR}/TubeTKLib/MinimalPathExtraction/include
  )

####
# TubeTK includes TubeTKLib's libraries and include directories
####
set( TubeTK_LIBRARIES
  TubeTK
  ${TubeTKLib_LIBRARIES}
  )

####
# When this module is loaded by an app, load TubeTKLib too.
####
set( TubeTK_EXPORT_CODE_INSTALL "
  set( TubeTKLib_DIR \"${TubeTKLib_DIR}\" )
  find_package( TubeTKLib REQUIRED HINTS \"${TubeTK_DIR}/TubeTKLib\" )
  ")

set( TubeTK_EXPORT_CODE_BUILD "
  if( NOT ITK_BINARY_DIR )
    set( TubeTKLib_DIR \"${TubeTKLib_DIR}\" )
    find_package( TubeTKLib REQUIRED HINTS \"${TubeTK_DIR}/TubeTKLib\" )
  endif()
  ")

####
# 
####
set( TubeTK_SYSTEM_INCLUDE_DIRS ${TubeTKLib_INCLUDE_DIRS} )
set( TubeTK_SYSTEM_LIBRARY_DIRS ${TubeTKLib_LIB_DIRS} )

if( NOT ITK_SOURCE_DIR )

  find_package( ITK REQUIRED )
  include( ${ITK_USE_FILE} )
  list( APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR} )

  include( ITKModuleExternal )

else()

  set( ITK_DIR ${CMAKE_BINARY_DIR} )
  itk_module_impl()

endif()

add_subdirectory( python )
add_subdirectory( apps )

#
# Enable access to remote testing data (must be placed last)
#
include( ${CMAKE_ROOT}/Modules/ExternalData.cmake )
ExternalData_Add_Target( TubeTKData )
