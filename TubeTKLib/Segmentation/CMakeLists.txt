##############################################################################
#
# Library:   TubeTKLib
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

project( TubeTKLibSegmentation )

if( TubeTKLib_USE_LIBSVM )
  find_package( LIBSVM )
  include_directories( ${LIBSVM_INCLUDE_DIRS} )
endif( TubeTKLib_USE_LIBSVM )

if( TubeTKLib_USE_RANDOMFOREST )
  find_package( RandomForest )
  include_directories( ${RandomForest_INCLUDE_DIR} )
endif( TubeTKLib_USE_RANDOMFOREST )

set( TubeTKLib_Base_Segmentation_H_Files
  itktubeComputeSegmentTubesParameters.h
  itktubePDFSegmenterBase.h
  itktubePDFSegmenterParzen.h
  itktubeRadiusExtractor2.h
  itktubeRidgeExtractor.h
  itktubeSegmentTubes.h
  itktubeSegmentBinaryImageSkeleton.h
  itktubeSegmentBinaryImageSkeleton3D.h
  itktubeSegmentTubesUsingMinimalPathFilter.h
  itktubeTubeExtractor.h
  itktubeRidgeSeedFilter.h
  itktubeComputeTrainingMask.h)
if( TubeTKLib_USE_LIBSVM )
  list( APPEND TubeTKLib_Base_Segmentation_H_Files
    itktubePDFSegmenterSVM.h )
endif( TubeTKLib_USE_LIBSVM )
if( TubeTKLib_USE_RANDOMFOREST )
  list( APPEND TubeTKLib_Base_Segmentation_H_Files
    itktubePDFSegmenterRandomForest.h )
endif( TubeTKLib_USE_RANDOMFOREST )

set( TubeTKLib_Base_Segmentation_HXX_Files
  itktubeComputeSegmentTubesParameters.hxx
  itktubePDFSegmenterBase.hxx
  itktubePDFSegmenterParzen.hxx
  itktubePDFSegmenterRandomForest.hxx
  itktubeRadiusExtractor2.hxx
  itktubeRidgeExtractor.hxx
  itktubeSegmentTubes.hxx
  itktubeSegmentBinaryImageSkeleton.hxx
  itktubeSegmentBinaryImageSkeleton3D.hxx
  itktubeSegmentTubesUsingMinimalPathFilter.hxx
  itktubeTubeExtractor.hxx
  itktubeRidgeSeedFilter.hxx
  itktubeComputeTrainingMask.hxx)
if( TubeTKLib_USE_LIBSVM )
  list( APPEND TubeTKLib_Base_Segmentation_HXX_Files
    itktubePDFSegmenterSVM.hxx )
endif( TubeTKLib_USE_LIBSVM )
if( TubeTKLib_USE_RANDOMFOREST )
  list( APPEND TubeTKLib_Base_Segmentation_HXX_Files
    itktubePDFSegmenterRandomForest.hxx )
endif( TubeTKLib_USE_RANDOMFOREST )

add_library( ${PROJECT_NAME} INTERFACE )

target_include_directories( ${PROJECT_NAME} INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR} )

target_link_libraries( ${PROJECT_NAME} INTERFACE
   TubeTKLibMetaIO
   TubeTKLibFiltering
   TubeTKLibNumerics
   MinimalPathExtraction )

if( TubeTKLib_BUILD_TESTING )
  add_subdirectory( Testing )
endif( TubeTKLib_BUILD_TESTING )

if( TubeTKLib_INSTALL_DEVELOPMENT )
  install( FILES
    ${TubeTKLib_Base_Segmentation_H_Files}
    ${TubeTKLib_Base_Segmentation_HXX_Files}
    DESTINATION ${TubeTKLib_INSTALL_INCLUDE_DIR}
    COMPONENT Development )

  TubeTKLibMacroInstallPlugins( ${PROJECT_NAME} )
endif()

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS ${PROJECT_NAME} )
