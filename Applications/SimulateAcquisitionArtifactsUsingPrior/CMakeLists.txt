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

set( MODULE_NAME SimulateAcquisitionArtifactsUsingPrior )
project( ${MODULE_NAME} )

set( SimulateAcquisitionArtifactsUsingPrior_H_Files
  tubeCompareImageWithPrior.h )

set( SimulateAcquisitionArtifactsUsingPrior_HXX_Files
  tubeCompareImageWithPrior.hxx )

add_custom_target( SimulateAcquisitionArtifactsUsingPriorInclude SOURCES
  ${SimulateAcquisitionArtifactsUsingPrior_H_Files}
  ${SimulateAcquisitionArtifactsUsingPrior_HXX_Files} )

if( NOT TubeTK_SOURCE_DIR )
  find_package( TubeTK REQUIRED )
  include( ${TubeTK_USE_FILE} )
endif( NOT TubeTK_SOURCE_DIR )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

find_package( ITK REQUIRED )
include( ${USE_ITK_FILE} )

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES} TubeNumerics
  INCLUDE_DIRECTORIES
    ${TubeTK_SOURCE_DIR}/Base/CLI
    ${TubeTK_SOURCE_DIR}/Base/Common
    ${TubeTK_SOURCE_DIR}/Base/Numerics
    ${TubeTK_SOURCE_DIR}/Base/Filtering
    ${TubeTK_SOURCE_DIR}/Base/Segmentation
    ${TubeTK_SLICER_SOURCE_DIR}/Modules/CLI/ExpertAutomatedRegistration
    ${TubeTK_SLICER_SOURCE_DIR}/Modules/CLI/ExpertAutomatedRegistration/ITKRegistrationHelper )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
