##############################################################################
#
# Library:   TubeTK
#
# Copyright 2012 Kitware Inc. 28 Corporate Drive,
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

project( TubeUSTK )

set( TubeTK_Base_USTK_H_Files
  itkAcousticImpulseResponseImageFilter.h
  itkAngleOfIncidenceImageFilter.h
  itkGradientBasedAngleOfIncidenceImageFilter.h
  itkLabelMapToAcousticImpedanceFunctor.h
  itkLabelMapToAcousticImpedanceImageFilter.h
  itkUltrasoundProbeGeometryCalculator.h
  SyncRecord.h
  SyncRecordManager.h
  itktubeInnerOpticToPlusImageReader.h
  itktubeMarkDuplicateFramesInvalidImageFilter.h )

set( TubeTK_Base_USTK_HXX_Files
  itkAcousticImpulseResponseImageFilter.hxx
  itkAngleOfIncidenceImageFilter.hxx
  itkGradientBasedAngleOfIncidenceImageFilter.hxx
  itkLabelMapToAcousticImpedanceFunctor.hxx
  itkLabelMapToAcousticImpedanceImageFilter.hxx
  itkUltrasoundProbeGeometryCalculator.hxx
  itktubeMarkDuplicateFramesInvalidImageFilter.hxx )

set( TubeTK_Base_USTK_SRCS
  SyncRecord.cpp
  SyncRecordManager.cpp
  itktubeInnerOpticToPlusImageReader.cxx )

add_library( TubeTKUSTK STATIC
  ${TubeTK_Base_USTK_H_Files}
  ${TubeTK_Base_USTK_HXX_Files}
  ${TubeTK_Base_USTK_SRCS} )

target_link_libraries( TubeTKUSTK PUBLIC TubeTKSegmentation )

target_include_directories( TubeTKUSTK PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR} )

set( TARGETS TubeTKUSTK )
TubeTKMacroInstallPlugins( ${TARGETS} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_USTK_H_Files}
  ${TubeTK_Base_USTK_HXX_Files}
  DESTINATION include )

set_property( GLOBAL APPEND PROPERTY TubeTK_TARGETS TubeTKUSTK )
