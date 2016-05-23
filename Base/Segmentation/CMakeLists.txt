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

if( TubeTK_USE_LIBSVM )
  find_package( LIBSVM )
  include_directories( ${LIBSVM_INCLUDE_DIRS} )
endif( TubeTK_USE_LIBSVM )

set( TubeTK_Base_Segmentation_H_Files
  itktubePDFSegmenterBase.h
  itktubePDFSegmenterParzen.h
  itktubeRadiusExtractor2.h
  itktubeRidgeExtractor.h
  itktubeSegmentTubesUsingMinimalPathFilter.h
  itktubeTubeExtractor.h )
if( TubeTK_USE_LIBSVM )
  list( APPEND TubeTK_Base_Segmentation_H_Files
    itktubePDFSegmenterSVM.h
    itktubeRidgeSeedFilter.h )
endif( TubeTK_USE_LIBSVM )

set( TubeTK_Base_Segmentation_HXX_Files
  itktubePDFSegmenterBase.hxx
  itktubePDFSegmenterParzen.hxx
  itktubeRadiusExtractor2.hxx
  itktubeRidgeExtractor.hxx
  itktubeSegmentTubesUsingMinimalPathFilter.hxx
  itktubeTubeExtractor.hxx )
if( TubeTK_USE_LIBSVM )
  list( APPEND TubeTK_Base_Segmentation_HXX_Files
    itktubePDFSegmenterSVM.hxx
    itktubeRidgeSeedFilter.hxx )
endif( TubeTK_USE_LIBSVM )

add_library( ${PROJECT_NAME} INTERFACE )

target_include_directories( ${PROJECT_NAME} INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR} )

target_link_libraries( ${PROJECT_NAME} INTERFACE
   TubeTKFiltering
   TubeTKNumerics
   MinimalPathExtraction )

if( TubeTK_BUILD_TESTING )
  add_subdirectory( Testing )
endif( TubeTK_BUILD_TESTING )

if( TubeTK_INSTALL_DEVELOPMENT )
  install( FILES
    ${TubeTK_Base_Segmentation_H_Files}
    ${TubeTK_Base_Segmentation_HXX_Files}
    DESTINATION ${TubeTK_INSTALL_INCLUDE_DIR}
    COMPONENT Development )

  TubeTKMacroInstallPlugins( ${PROJECT_NAME} )
endif()

set_property( GLOBAL APPEND PROPERTY TubeTK_TARGETS ${PROJECT_NAME} )
