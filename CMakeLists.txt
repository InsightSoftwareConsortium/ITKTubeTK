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

####
####
####
####
project( TubeTK )

####
# Setup CMake
####
list( APPEND CMAKE_MODULE_PATH ${TubeTK_SOURCE_DIR}/CMake )

####
# Setup External Data Cache
####
include( TubeTKExternalData )
if( DEFINED "ENV{ExternalData_OBJECT_STORES}" )
  set( ExternalData_OBJECT_STORES $ENV{ExternalData_OBJECT_STORES} )
endif()


####
# Default output directories
####
make_directory( "${TubeTK_BINARY_DIR}/bin" )
make_directory( "${TubeTK_BINARY_DIR}/lib" )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${TubeTK_BINARY_DIR}/bin" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${TubeTK_BINARY_DIR}/lib" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${TubeTK_BINARY_DIR}/lib" )

mark_as_advanced( BUILD_SHARED_LIBS )
mark_as_advanced( CMAKE_INSTALL_PREFIX )

####
# Setup testing
####
if( NOT BUILD_TESTING )
  option( BUILD_TESTING "Build TubeTK's tests" on )
  mark_as_advanced( BUILD_TESTING )
endif()
set( TubeTK_DATA_ROOT ${TubeTK_SOURCE_DIR}/data_keys)

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
find_package( PythonInterp )
find_package( PythonLibs REQUIRED )
mark_as_advanced( NUMPY_INCLUDE_DIR )

####
# VTK setup.
####
option( TubeTK_USE_VTK "Use VTK to enable additional functionality." ON )
if( TubeTK_USE_VTK )
  find_package( VTK REQUIRED )
  mark_as_advanced( CLEAR VTK_DIR )
  include( ${VTK_USE_FILE} )
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
  include_directories( ${ArrayFire_INCLUDE_DIRS} )
endif()


####
# LIBSVM setup.
####
option( TubeTK_USE_LIBSVM "Use LIBSVM for additional functionality."
  OFF )
mark_as_advanced( TubeTK_USE_LIBSVM )
if( TubeTK_USE_LIBSVM )
  find_package( LIBSVM REQUIRED )
  mark_as_advanced( CLEAR LIBSVM_DIR )
  include_directories( ${LIBSVM_INCLUDE_DIRS} )
endif()

####
# RandomForest setup.
####
option( TubeTK_USE_RANDOMFOREST
  "Use RandomForest to enable additional functionality." OFF )
mark_as_advanced( TubeTK_USE_RANDOMFOREST )
if( TubeTK_USE_RANDOMFOREST )
  find_package( RandomForest REQUIRED )
  mark_as_advanced( CLEAR RandomForest_DIR )
  include_directories( ${RandomForest_INCLUDE_DIRS} )
endif()

####
# Create TubeTK Configuration file (to pass ArrayFire and other flags)
####
configure_file( "${TubeTK_SOURCE_DIR}/CMake/tubetkConfigure.h.in"
  tubetkConfigure.h @ONLY )
include_directories( ${TubeTK_BINARY_DIR} )

####
# CompareImages and CompareTextFiles for app tests
####
set( TubeTKLib_CompareTextFiles_EXE
  ${TubeTK_BINARY_DIR}/bin/CompareTextFiles
  )
set( TubeTKLib_CompareImages_EXE
  ${TubeTK_BINARY_DIR}/bin/CompareImages
  )

####
# Configure for building against a build tree
####
set( TubeTKLib_LIBRARIES
  TubeTKLibCommon
  TubeTKLibMetaIO
  TubeTKLibNumerics
  TubeTKLibObjectDocuments
  TubeTKLibUSTK )

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
  ${TubeTK_SOURCE_DIR}/TubeTKLib/USTK )

set( TubeTKLib_LIBRARY_DIRS
  ${TubeTK_BINARY_DIR}/bin
  ${TubeTK_BINARY_DIR}/lib )

####
# TubeTK includes TubeTKLib's libraries and include directories
####
set( TubeTK_LIBRARIES
  TubeTK
  ${TubeTKLib_LIBRARIES} )

set( TubeTK_INCLUDE_DIRS
  ${ITK_INCLUDE_DIRS}
  ${TubeTK_SOURCE_DIR}/include
  ${TubeTKLib_INCLUDE_DIRS} )

set( TubeTK_LIBRARY_DIRS
  ${TubeTK_BINARY_DIR}/bin 
  ${TubeTK_BINARY_DIR}/lib
  ${TubeTK_EXECUTABLE_DIRS} )

set( TubeTK_SYSTEM_INCLUDE_DIRS ${TubeTK_INCLUDE_DIRS} )
set( TubeTK_SYSTEM_LIBRARY_DIRS ${TubeTK_LIBRARY_DIRS} )

####
# Find ITK and pre-built modules
####
set( TubeTK_ITK_COMPONENTS
  MinimalPathExtraction
  ITKImageIO
  ITKTransformIO
  ITKIOSpatialObjects
  ITKRegistrationCommon
  ITKBinaryMathematicalMorphology
  ITKFFT
  ITKCommon
  ITKHDF5
  ITKIOCSV
  ITKRegionGrowing
  ITKAnisotropicSmoothing
  ITKLabelVoting
  ITKPDEDeformableRegistration
  ITKIOImageBase
  ITKImageFunction
  ITKImageIntensity
  ITKDistanceMap
  ITKTestKernel
  ITKOptimizers
  ITKSmoothing
  ITKSpatialObjects
  ITKStatistics
  ITKTransform
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
  list( APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR} )
  include( ITKModuleExternal )
else()
  set( ITK_DIR ${CMAKE_BINARY_DIR} )
  itk_module_impl()
endif()

add_subdirectory( TubeTKLib )
add_subdirectory( python )
add_subdirectory( apps )

####
# Enable access to remote testing data (must be placed last)
####
include( ${CMAKE_ROOT}/Modules/ExternalData.cmake )
ExternalData_Add_Target( TubeTKData )
