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

project( TubeTKLibIO )

set( TubeTKLib_Base_IO_H_Files
  itktubePDFSegmenterParzenIO.h
  itktubeRidgeSeedFilterIO.h
  itktubeTubeExtractorIO.h
  itktubeTubeXIO.h )
if( TubeTKLib_USE_LIBSVM )
  find_package( LIBSVM REQUIRED )
  include_directories( ${LIBSVM_INCLUDE_DIRS} )
  list( APPEND TubeTKLib_Base_IO_H_Files
    itktubePDFSegmenterSVMIO.h )
endif( TubeTKLib_USE_LIBSVM )
if( TubeTKLib_USE_RANDOMFOREST )
  find_package( RandomForest REQUIRED )
  include_directories( ${RandomForest_INCLUDE_DIR} )
  list( APPEND TubeTKLib_Base_IO_H_Files
    itktubePDFSegmenterRandomForestIO.h )
endif( TubeTKLib_USE_RANDOMFOREST )

set( TubeTKLib_Base_IO_HXX_Files
  itktubePDFSegmenterParzenIO.hxx
  itktubeRidgeSeedFilterIO.hxx
  itktubeTubeExtractorIO.hxx
  itktubeTubeXIO.hxx )
if( TubeTKLib_USE_LIBSVM )
  list( APPEND TubeTKLib_Base_IO_HXX_Files
    itktubePDFSegmenterSVMIO.hxx )
endif( TubeTKLib_USE_LIBSVM )
if( TubeTKLib_USE_RANDOMFOREST )
  list( APPEND TubeTKLib_Base_IO_HXX_Files
    itktubePDFSegmenterRandomForestIO.hxx )
endif( TubeTKLib_USE_RANDOMFOREST )

add_library( ${PROJECT_NAME} INTERFACE )

target_link_libraries( ${PROJECT_NAME} INTERFACE
  TubeTKLibMetaIO TubeTKLibCommon TubeTKLibNumerics TubeTKLibSegmentation )

target_include_directories( ${PROJECT_NAME} INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" )

if( TubeTKLib_BUILD_TESTING )
  add_subdirectory( Testing )
endif( TubeTKLib_BUILD_TESTING )

if( TubeTKLib_INSTALL_DEVELOPMENT )
  install( FILES
    ${TubeTKLib_Base_IO_H_Files}
    ${TubeTKLib_Base_IO_HXX_Files}
    DESTINATION ${TubeTKLib_INSTALL_INCLUDE_DIR}
    COMPONENT Development )
endif()

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS ${PROJECT_NAME} )
