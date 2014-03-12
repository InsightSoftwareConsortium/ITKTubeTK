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

find_package( ParameterSerializer REQUIRED )
include_directories( ${ParameterSerializer_INCLUDE_DIRS} )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )
if( SlicerExecutionModel_USE_SERIALIZER )
  add_definitions( -DSlicerExecutionModel_USE_SERIALIZER )
endif()

find_package( ITK REQUIRED )
include( ${USE_ITK_FILE} )

set(tubetk_base_includes
  ${TubeTK_SOURCE_DIR}/Base/CLI
  ${TubeTK_SOURCE_DIR}/Base/Common
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer/Filtering
  ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer/Registration
  ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer/USTK
  ${TubeTK_SOURCE_DIR}/Base/Registration
  ${TubeTK_SOURCE_DIR}/Base/USTK )

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES} ${ParameterSerializer_LIBRARIES}
  INCLUDE_DIRECTORIES ${tubetk_base_includes} )

SEMMacroBuildCLI(
  NAME ComputeImageToTubeRigidMetricImage
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES} ${ParameterSerializer_LIBRARIES}
  INCLUDE_DIRECTORIES ${tubetk_base_includes} )

SEMMacroBuildCLI(
  NAME ExtractMetricImageSlice
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES} ${ParameterSerializer_LIBRARIES}
  INCLUDE_DIRECTORIES ${tubetk_base_includes} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
