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

project( TubeTK )

if(NOT ITK_SOURCE_DIR)
  include(itk-module-init.cmake)
endif()

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
option(TubeTK_USE_VTK "Build TubeTK filters that require VTK." ON)
if(TubeTK_USE_VTK)
  find_package(VTK REQUIRED
               HINTS ./VTK-build ../VTK-build ../../VTK-build)

  # Add compiler flags needed to use VTK.
  set(CMAKE_C_FLAGS
      "${CMAKE_C_FLAGS} ${VTK_REQUIRED_C_FLAGS}")
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} ${VTK_REQUIRED_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} ${VTK_REQUIRED_EXE_LINKER_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS
      "${CMAKE_SHARED_LINKER_FLAGS} ${VTK_REQUIRED_SHARED_LINKER_FLAGS}")
  set(CMAKE_MODULE_LINKER_FLAGS
      "${CMAKE_MODULE_LINKER_FLAGS} ${VTK_REQUIRED_MODULE_LINKER_FLAGS}")

  # The VTK DICOMParser and vtkmetaio includes conflict with the ITK
  # versions. Here we remove them from the include directories.
  #
  string(REGEX REPLACE "[^;]*MetaIO;"
           "" TubeTK_VTK_INCLUDE_DIRS "${VTK_INCLUDE_DIRS}")
  string(REGEX REPLACE "[^;]*vtkmetaio;"
           "" TubeTK_VTK_INCLUDE_DIRS "${VTK_INCLUDE_DIRS}")
  string(REGEX REPLACE "[^;]*DICOMParser;"
           "" TubeTK_VTK_INCLUDE_DIRS "${VTK_INCLUDE_DIRS}")

  set(TubeTK_EXPORT_CODE_INSTALL "
  set(VTK_DIR \"${VTK_DIR}\")
  find_package(VTK REQUIRED)
  ")

  set(TubeTK_EXPORT_CODE_BUILD "
  if(NOT ITK_BINARY_DIR)
    set(VTK_DIR \"${VTK_DIR}\")
    find_package(VTK REQUIRED)
  endif()
  ")
endif() # ITKTubeTK_USE_VTK

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
  ITKVTK
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
endif()

####
# Create TubeTK Configuration file (to pass flags)
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
if( TubeTK_USE_VTK )
  list( APPEND TubeTK_LIBRARIES ${VTK_LIBRARIES} )
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

if( TubeTK_USE_VTK )
  list( APPEND TubeTK_SYSTEM_INCLUDE_DIRS ${TubeTK_VTK_INCLUDE_DIRS} )
endif()

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
# Build Applications?
####
option( TubeTK_BUILD_APPLICATIONS "Build applications based on TubeTK." OFF )

itk_module_examples()
