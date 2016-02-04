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

project( TubeTKSegmentation )

set( TubeTK_Base_Segmentation_H_Files
  itktubeCVTImageFilter.h
  itktubeLabelOverlapMeasuresImageFilter.h
  itktubePDFSegmenter.h
  itktubeRadiusExtractor2.h
  itktubeRidgeExtractor.h
  itktubeRidgeSeedFilter.h
  itktubeTubeExtractor.h )

set( TubeTK_Base_Segmentation_HXX_Files
  itktubeCVTImageFilter.hxx
  itktubeLabelOverlapMeasuresImageFilter.hxx
  itktubePDFSegmenter.hxx
  itktubeRadiusExtractor2.hxx
  itktubeRidgeExtractor.hxx
  itktubeRidgeSeedFilter.hxx
  itktubeTubeExtractor.hxx )

add_library( TubeTKSegmentation INTERFACE )

target_include_directories( TubeTKSegmentation INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR} )

target_link_libraries( TubeTKSegmentation INTERFACE
   TubeTKFiltering # Required because of itktubePDFSegmenter.hxx including itktubeSmoothingRecursiveGaussianImageFilter.h
   TubeTKNumerics )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_Segmentation_H_Files}
  ${TubeTK_Base_Segmentation_HXX_Files}
  DESTINATION include )

TubeTKMacroInstallPlugins( TubeTKSegmentation )

set_property( GLOBAL APPEND PROPERTY TubeTK_TARGETS TubeTKSegmentation )
