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

set( TubeTKLib_Base_IO_H_Files
  itktubePDFSegmenterParzenIO.h
  itktubeRidgeSeedFilterIO.h
  itktubeTubeExtractorIO.h
  itktubeTubeXIO.h )

if( TubeTK_USE_LIBSVM )
  list( APPEND TubeTKLib_Base_IO_H_Files
    itktubePDFSegmenterSVMIO.h )
endif()

if( TubeTK_USE_RANDOMFOREST )
  list( APPEND TubeTKLib_Base_IO_H_Files
    itktubePDFSegmenterRandomForestIO.h )
endif()

set( TubeTKLib_Base_IO_HXX_Files
  itktubePDFSegmenterParzenIO.hxx
  itktubeRidgeSeedFilterIO.hxx
  itktubeTubeExtractorIO.hxx
  itktubeTubeXIO.hxx )

if( TubeTK_USE_LIBSVM )
  list( APPEND TubeTKLib_Base_IO_HXX_Files
    itktubePDFSegmenterSVMIO.hxx )
endif()

if( TubeTK_USE_RANDOMFOREST )
  list( APPEND TubeTKLib_Base_IO_HXX_Files
    itktubePDFSegmenterRandomForestIO.hxx )
endif()

add_library( TubeTKLibIO INTERFACE )

target_link_libraries( TubeTKLibIO INTERFACE
  TubeTKLibMetaIO
  TubeTKLibCommon
  TubeTKLibNumerics
  TubeTKLibSegmentation )

target_include_directories( TubeTKLibIO INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTKLib_Base_IO_H_Files}
  ${TubeTKLib_Base_IO_HXX_Files}
  DESTINATION include
  COMPONENT Development )

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS TubeTKLibIO )
