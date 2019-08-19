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

set( TubeTKLib_Base_USTK_H_Files
  itkAngleOfIncidenceImageFilter.h
  itkLabelMapToAcousticImpedanceFunctor.h
  itkLabelMapToAcousticImpedanceImageFilter.h
  itkUltrasoundProbeGeometryCalculator.h
  SyncRecord.h
  SyncRecordManager.h
  itktubeInnerOpticToPlusImageReader.h
  itktubeMarkDuplicateFramesInvalidImageFilter.h )

set( TubeTKLib_Base_USTK_HXX_Files
  itkAngleOfIncidenceImageFilter.hxx
  itkLabelMapToAcousticImpedanceFunctor.hxx
  itkLabelMapToAcousticImpedanceImageFilter.hxx
  itkUltrasoundProbeGeometryCalculator.hxx
  itktubeMarkDuplicateFramesInvalidImageFilter.hxx )

set( TubeTKLib_Base_USTK_SRCS
  SyncRecord.cpp
  SyncRecordManager.cpp
  itktubeInnerOpticToPlusImageReader.cxx )

add_library( TubeTKLibUSTK STATIC
  ${TubeTKLib_Base_USTK_H_Files}
  ${TubeTKLib_Base_USTK_HXX_Files}
  ${TubeTKLib_Base_USTK_SRCS} )

target_link_libraries( TubeTKLibUSTK PUBLIC TubeTKLibSegmentation )

target_include_directories( TubeTKLibUSTK PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTKLib_Base_USTK_H_Files}
  ${TubeTKLib_Base_USTK_HXX_Files}
  DESTINATION include
  COMPONENT Development )

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS TubeTKLibUSTK )
