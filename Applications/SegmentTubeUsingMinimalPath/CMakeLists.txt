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
set( MODULE_NAME SegmentTubeUsingMinimalPath )
project( ${MODULE_NAME} )

#Check if MinimalPathExtraction is built in ITK as a remote module
list( FIND ITK_MODULES_ENABLED MinimalPathExtraction
  MinimalPathExtraction_ITK )

if( NOT( MinimalPathExtraction_DIR ) AND
  ( ${MinimalPathExtraction_ITK} EQUAL "-1" ) )

  if( NOT MinimalPathExtraction_SOURCE_DIR )

    message( FATAL_ERROR "MinimalPathExtraction_SOURCE_DIR not found" )

  endif()

  # XXX Only clean build will ensure expected version is fetched.
  if( NOT EXISTS ${MinimalPathExtraction_SOURCE_DIR} )

    set( GIT_CLONE_ARGS "clone" "${MinimalPathExtraction_LOCAL_BUILD_URL}"
      "${MinimalPathExtraction_SOURCE_DIR}" )
    message(STATUS "Cloning MinimalPathExtraction into '${MinimalPathExtraction_SOURCE_DIR}'")
    execute_process( COMMAND ${GIT_EXECUTABLE} ${GIT_CLONE_ARGS}
      WORKING_DIRECTORY ${TubeTK_BINARY_DIR}
      OUTPUT_QUIET
      ERROR_VARIABLE _error
      RESULT_VARIABLE _result )
    if( NOT _result EQUAL 0 )
      message(FATAL_ERROR "Failed to clone '${MinimalPathExtraction_LOCAL_BUILD_URL}': ${_error}")
    endif()

    set( GIT_HASH_ARGS "checkout"
      "${MinimalPathExtraction_LOCAL_BUILD_HASH_OR_TAG}" )
    message(STATUS "Checking out MinimalPathExtraction '${MinimalPathExtraction_LOCAL_BUILD_HASH_OR_TAG}'")
    execute_process( COMMAND ${GIT_EXECUTABLE} ${GIT_HASH_ARGS}
      WORKING_DIRECTORY ${MinimalPathExtraction_SOURCE_DIR}
      OUTPUT_QUIET
      ERROR_VARIABLE _error
      RESULT_VARIABLE _result )
    if( NOT _result EQUAL 0 )
      message(FATAL_ERROR "Failed to checkout MinimalPathExtraction '${MinimalPathExtraction_LOCAL_BUILD_HASH_OR_TAG}': ${_error}")
    endif()

  endif( NOT EXISTS ${MinimalPathExtraction_SOURCE_DIR} )

  # Find ITK
  find_package( ITK REQUIRED )
  if( TubeTK_BUILD_WITHIN_SLICER )
    set( ITK_NO_IO_FACTORY_REGISTER_MANAGER 1 )
  endif( TubeTK_BUILD_WITHIN_SLICER )
  include( ${ITK_USE_FILE} )

  #Find MinimalPathExtraction additional sources
  set( ITKMinimalPathExtraction_SRCS
    "${MinimalPathExtraction_SOURCE_DIR}/src/itkIterateNeighborhoodOptimizer.cxx" )

  SEMMacroBuildCLI(
    NAME ${MODULE_NAME}
    LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
    TARGET_LIBRARIES
      ${ITK_LIBRARIES} TubeTKIO ITKIOMeta ITKIOSpatialObjects ITKIOImageBase
        TubeCLI TubeTKNumerics TubeTKCommon
    INCLUDE_DIRECTORIES
      ${MinimalPathExtraction_SOURCE_DIR}/include
    ADDITIONAL_SRCS
      ${ITKMinimalPathExtraction_SRCS} )

else()

  if( NOT TubeTK_SOURCE_DIR )
    find_package( TubeTK REQUIRED )
    include( ${TubeTK_USE_FILE} )
  endif( NOT TubeTK_SOURCE_DIR )

  find_package( SlicerExecutionModel REQUIRED )
  include( ${SlicerExecutionModel_USE_FILE} )

  find_package( ITK REQUIRED COMPONENTS MinimalPathExtraction )
  include( ${ITK_USE_FILE} )
  set( MinimalPathExtraction_REQUIRED_LIB ${ITK_LIBRARIES})

  find_package( ITK REQUIRED )
  if( TubeTK_BUILD_WITHIN_SLICER )
    set( ITK_NO_IO_FACTORY_REGISTER_MANAGER 1 )
    # See Libs/ITKFactoryRegistration/CMakeLists.txt
  endif( TubeTK_BUILD_WITHIN_SLICER )
  include( ${ITK_USE_FILE} )
  set( ITK_LIBRARIES ${MinimalPathExtraction_REQUIRED_LIB} ${ITK_LIBRARIES} )
  list( REMOVE_DUPLICATES ITK_LIBRARIES )

  SEMMacroBuildCLI(
    NAME ${MODULE_NAME}
    LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
    TARGET_LIBRARIES
      ${ITK_LIBRARIES} ITKIOMeta ITKIOSpatialObjects ITKIOImageBase
      TubeCLI TubeTKCommon TubeTKSegmentation )

endif()

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
