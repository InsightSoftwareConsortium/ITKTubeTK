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

set( TubeTKLib_Base_Segmentation_H_Files
  itktubeComputeSegmentTubesParameters.h
  itktubePDFSegmenterBase.h
  itktubePDFSegmenterParzen.h
  itktubeRadiusExtractor2.h
  itktubeRidgeExtractor.h
  itktubeSegmentBinaryImageSkeleton.h
  itktubeSegmentBinaryImageSkeleton3D.h
  itktubeSegmentTubeUsingMinimalPathFilter.h
  itktubeTubeExtractor.h
  itktubeRidgeSeedFilter.h
  itktubeComputeTrainingMask.h)
if( TubeTK_USE_LIBSVM )
  list( APPEND TubeTKLib_Base_Segmentation_H_Files
    itktubePDFSegmenterSVM.h )
endif()
if( TubeTK_USE_RANDOMFOREST )
  list( APPEND TubeTKLib_Base_Segmentation_H_Files
    itktubePDFSegmenterRandomForest.h )
endif()

set( TubeTKLib_Base_Segmentation_HXX_Files
  itktubeComputeSegmentTubesParameters.hxx
  itktubePDFSegmenterBase.hxx
  itktubePDFSegmenterParzen.hxx
  itktubePDFSegmenterRandomForest.hxx
  itktubeRadiusExtractor2.hxx
  itktubeRidgeExtractor.hxx
  itktubeSegmentBinaryImageSkeleton.hxx
  itktubeSegmentBinaryImageSkeleton3D.hxx
  itktubeSegmentTubeUsingMinimalPathFilter.hxx
  itktubeTubeExtractor.hxx
  itktubeRidgeSeedFilter.hxx
  itktubeComputeTrainingMask.hxx)

if( TubeTK_USE_LIBSVM )
  list( APPEND TubeTKLib_Base_Segmentation_HXX_Files
    itktubePDFSegmenterSVM.hxx )
endif()

if( TubeTK_USE_RANDOMFOREST )
  list( APPEND TubeTKLib_Base_Segmentation_HXX_Files
    itktubePDFSegmenterRandomForest.hxx )
endif()

add_library( TubeTKLibSegmentation INTERFACE )

target_include_directories( TubeTKLibSegmentation INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories(
  ${ITK_INCLUDE_DIRS} )

target_link_libraries( TubeTKLibSegmentation INTERFACE
   TubeTKLibIO
   TubeTKLibMetaIO
   TubeTKLibFiltering
   TubeTKLibNumerics
   MinimalPathExtraction
   ${ITK_LIBRARIES} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTKLib_Base_Segmentation_H_Files}
  ${TubeTKLib_Base_Segmentation_HXX_Files}
  DESTINATION include
  COMPONENT Development )

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS
  TubeTKLibSegmentation )
