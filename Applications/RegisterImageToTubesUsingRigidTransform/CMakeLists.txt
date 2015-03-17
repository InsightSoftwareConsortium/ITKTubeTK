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

set( MODULE_NAME RegisterImageToTubesUsingRigidTransform )
project( ${MODULE_NAME} )

# To silence the warning about getting LOCATION (Tuner only targets non-VS for
# now)
if( NOT ${CMAKE_VERSION} VERSION_LESS 3.0 )
  cmake_policy(SET CMP0026 OLD)
endif( NOT ${CMAKE_VERSION} VERSION_LESS 3.0 )

if( NOT TubeTK_SOURCE_DIR )
  find_package( TubeTK REQUIRED )
  include( ${TubeTK_USE_FILE} )
endif( NOT TubeTK_SOURCE_DIR )

find_package( ITK REQUIRED )
if( TubeTK_RUN_MODULES_IN_SLICER )
  set( ITK_NO_IO_FACTORY_REGISTER_MANAGER 1 ) # See Libs/ITKFactoryRegistration/CMakeLists.txt
endif( TubeTK_RUN_MODULES_IN_SLICER )
include( ${ITK_USE_FILE} )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

if( SlicerExecutionModel_USE_SERIALIZER )
  add_definitions( -DSlicerExecutionModel_USE_SERIALIZER )
endif()

set( ${MODULE_NAME}_INCLUDE_DIRS
  ${TubeTK_SOURCE_DIR}/Base/CLI
  ${TubeTK_SOURCE_DIR}/Base/Common
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Registration
  ${TubeTK_SOURCE_DIR}/Base/USTK )

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES} ITKIOSpatialObjects ITKIOMeta
  INCLUDE_DIRECTORIES ${${MODULE_NAME}_INCLUDE_DIRS} )

SEMMacroBuildCLI(
  NAME ComputeImageToTubeRigidMetricImage
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES} ITKIOSpatialObjects ITKIOMeta
  INCLUDE_DIRECTORIES ${${MODULE_NAME}_INCLUDE_DIRS} )

SEMMacroBuildCLI(
  NAME ExtractMetricImageSlice
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES} ITKIOSpatialObjects ITKIOMeta
  INCLUDE_DIRECTORIES ${${MODULE_NAME}_INCLUDE_DIRS} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
