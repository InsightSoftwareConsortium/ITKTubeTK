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

project( TubeTKIO )

set( TubeTK_Base_IO_H_Files
  itktubePDFSegmenterParzenIO.h
  itktubeRidgeSeedFilterIO.h
  itktubeTubeExtractorIO.h
  itktubeTubeXIO.h )
if( TubeTK_USE_LIBSVM )
  find_package( LIBSVM REQUIRED )
  include_directories( ${LIBSVM_INCLUDE_DIRS} )
  list( APPEND TubeTK_Base_IO_H_Files
    itktubePDFSegmenterSVMIO.h )
endif( TubeTK_USE_LIBSVM )
if( TubeTK_USE_RandomForest )
  find_package( RandomForest REQUIRED )
  include_directories( ${RandomForest_INCLUDE_DIR} )
  list( APPEND TubeTK_Base_IO_H_Files
    itktubePDFSegmenterRandomForestIO.h )
endif( TubeTK_USE_RandomForest )

set( TubeTK_Base_IO_HXX_Files
  itktubePDFSegmenterParzenIO.hxx
  itktubeRidgeSeedFilterIO.hxx
  itktubeTubeExtractorIO.hxx
  itktubeTubeXIO.hxx )
if( TubeTK_USE_LIBSVM )
  list( APPEND TubeTK_Base_IO_HXX_Files
    itktubePDFSegmenterSVMIO.hxx )
endif( TubeTK_USE_LIBSVM )
if( TubeTK_USE_RandomForest )
  list( APPEND TubeTK_Base_IO_HXX_Files
    itktubePDFSegmenterRandomForestIO.hxx )
endif( TubeTK_USE_RandomForest )

add_library( ${PROJECT_NAME} INTERFACE )

target_link_libraries( ${PROJECT_NAME} INTERFACE
  TubeTKMetaIO TubeTKCommon TubeTKNumerics TubeTKSegmentation )

target_include_directories( ${PROJECT_NAME} INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" )

if( TubeTK_BUILD_TESTING )
  add_subdirectory( Testing )
endif( TubeTK_BUILD_TESTING )

if( TubeTK_INSTALL_DEVELOPMENT )
  install( FILES
    ${TubeTK_Base_IO_H_Files}
    ${TubeTK_Base_IO_HXX_Files}
    DESTINATION ${TubeTK_INSTALL_INCLUDE_DIR}
    COMPONENT Development )

  TubeTKMacroInstallPlugins( ${PROJECT_NAME} )
endif()

set_property( GLOBAL APPEND PROPERTY TubeTK_TARGETS ${PROJECT_NAME} )
