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
cmake_minimum_required(VERSION 3.16.0)

project( TubeTK )

####
# Setup CMake
####
list( APPEND CMAKE_MODULE_PATH ${TubeTK_SOURCE_DIR}/CMake )

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # if using Linux gcc in debug mode, enable gcov computation
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
endif()

####
# Setup External Data Cache
####
include( TubeTKExternalData )
if( DEFINED "ENV{ExternalData_OBJECT_STORES}" )
  set( ExternalData_OBJECT_STORES $ENV{ExternalData_OBJECT_STORES} )
endif()

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

include( TubeTKExternalData )
if( DEFINED "ENV{ExternalData_OBJECT_STORES}" )
  set( ExternalData_OBJECT_STORES $ENV{ExternalData_OBJECT_STORES} )
endif()


####
# Find Python
####
if( ITK_WRAP_PYTHON )
  find_package( PythonInterp )
  find_package( PythonLibs REQUIRED )
  find_package( NumPy REQUIRED )
  mark_as_advanced( NUMPY_INCLUDE_DIR )
endif()

####
# VTK setup.
####
option( TubeTK_USE_VTK "Use VTK to enable additional functionality." OFF )
if( TubeTK_USE_VTK )
  find_package( VTK REQUIRED )
  mark_as_advanced( CLEAR VTK_DIR )
  if( VTK_USE_FILE )
    include( ${VTK_USE_FILE} )
  endif()
endif()

####
# ArrayFire setup
####
option( TubeTK_USE_ARRAYFIRE
  "Use the ArrayFire library to speedup filtering opertions using the GPU."
  OFF)
if( TubeTK_USE_ARRAYFIRE )
  find_package( ArrayFire REQUIRED )
  mark_as_advanced( CLEAR ArrayFire_DIR )
endif()

####
# Find ITK and pre-built modules
####
set( TubeTK_ITK_COMPONENTS
  MinimalPathExtraction
  ITKAnisotropicSmoothing
  ITKBinaryMathematicalMorphology
  ITKCommon
  ITKDistanceMap
  ITKFFT
  ITKHDF5
  ITKImageFunction
  ITKImageIntensity
  ITKImageIO
  ITKIO
  ITKIOCSV
  ITKIOImageBase
  ITKIOMeta
  ITKIOSpatialObjects
  ITKLabelVoting
  ITKOptimizers
  ITKPDEDeformableRegistration
  ITKRegionGrowing
  ITKRegistrationCommon
  ITKSmoothing
  ITKSpatialObjects
  ITKTestKernel
  ITKTransform
  ITKTransformIO
  ITKStatistics
  )

if( NOT ITK_SOURCE_DIR )
  find_package( ITK COMPONENTS ${TubeTK_ITK_COMPONENTS} REQUIRED )
  mark_as_advanced( CLEAR ITK_DIR )
  if( NOT ITK_FOUND )
    find_package( ITK )
    if( ITK_FOUND )
      message( FATAL_ERROR
        "TubeTK requires ITK with MinimalPathExtraction module enabled.")
    endif()
    set( ITK_DIR "NOT_FOUND" CACHE PATH
	    "Path to ITK installation/build directory." FORCE )
  endif()
  include( ${ITK_USE_FILE} )
endif()

####
# Create TubeTK Configuration file (to pass ArrayFire and other flags)
####
configure_file( "${TubeTK_SOURCE_DIR}/CMake/tubetkConfigure.h.in"
  ${TubeTK_BINARY_DIR}/include/tubetkConfigure.h @ONLY )

####
# Setup include directories
####
set( TubeTK_INCLUDE_DIRS
  $<INSTALL_INTERFACE:include>
  $<INSTALL_INTERFACE:include/Common>
  $<INSTALL_INTERFACE:include/Filtering>
  $<INSTALL_INTERFACE:include/IO>
  $<INSTALL_INTERFACE:include/MetaIO>
  $<INSTALL_INTERFACE:include/Numerics>
  $<INSTALL_INTERFACE:include/ObjectDocuments>
  $<INSTALL_INTERFACE:include/Registration>
  $<INSTALL_INTERFACE:include/Segmentation> )

set( TubeTK_LIBRARIES TubeTK )

if( TubeTK_USE_ARRAYFIRE )
  list( APPEND TubeTKLib_INCLUDE_DIRS
    ${ArrayFire_INCLUDE_DIRS} )
endif()

####
# Build base library and optional components
####
set( TubeTK_SYSTEM_INCLUDE_DIRS
  ${TubeTK_BINARY_DIR}/include
  ${TubeTK_SOURCE_DIR}/src
  ${TubeTK_SOURCE_DIR}/src/Common
  ${TubeTK_SOURCE_DIR}/src/Filtering
  ${TubeTK_SOURCE_DIR}/src/IO
  ${TubeTK_SOURCE_DIR}/src/MetaIO
  ${TubeTK_SOURCE_DIR}/src/Numerics
  ${TubeTK_SOURCE_DIR}/src/ObjectDocuments
  ${TubeTK_SOURCE_DIR}/src/Registration
  ${TubeTK_SOURCE_DIR}/src/Segmentation )

include_directories( ${TubeTK_SYSTEM_INCLUDE_DIRS} )

####
# Build TubeTK module
####
if( NOT ITK_SOURCE_DIR )
  list( APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR} )
  include( ITKModuleExternal )
else()
  set(ITK_DIR ${CMAKE_BINARY_DIR})
  itk_module_impl()
endif()

####
# Valgrind / MemCheck Options
####
set( MEMORYCHECK_SUPPRESSIONS_FILE
 "${TubeTK_SOURCE_DIR}/CMake/TubeTKValgrind.supp" )

####
# Build Examples?
####
itk_module_examples()
