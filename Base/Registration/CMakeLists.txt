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

project( TubeTKRegistration )

set( TubeTK_Base_Registration_H_Files
  itktubeAnisotropicDiffusiveRegistrationFunction.h
  itktubeDiffusiveRegistrationFilter.h
  itktubeDiffusiveRegistrationFilterUtils.h
  itktubeImageToTubeRigidMetric.h
  itktubeImageToTubeRigidRegistration.h
  itktubeMeanSquareRegistrationFunction.h
  itktubeTubeExponentialResolutionWeightFunction.h
  itktubeTubeParametricExponentialResolutionWeightFunction.h
  itktubeTubeParametricExponentialWithBoundsResolutionWeightFunction.h
  itktubeTubeToTubeTransformFilter.h )

if( TubeTK_USE_VTK )
  list( APPEND TubeTK_Base_Registration_H_Files
    itktubeAnisotropicDiffusiveSparseRegistrationFilter.h
    itktubeAnisotropicDiffusiveRegistrationFilter.h )
endif( TubeTK_USE_VTK )

set( TubeTK_Base_Registration_HXX_Files
  itktubeAnisotropicDiffusiveRegistrationFunction.hxx
  itktubeDiffusiveRegistrationFilter.hxx
  itktubeDiffusiveRegistrationFilterUtils.hxx
  itktubeImageToTubeRigidMetric.hxx
  itktubeImageToTubeRigidRegistration.hxx
  itktubeMeanSquareRegistrationFunction.hxx
  itktubeTubeToTubeTransformFilter.hxx )

if( TubeTK_USE_VTK )
  list( APPEND TubeTK_Base_Registration_HXX_Files
    itktubeAnisotropicDiffusiveSparseRegistrationFilter.hxx
    itktubeAnisotropicDiffusiveRegistrationFilter.hxx )
endif( TubeTK_USE_VTK )

add_library( TubeTKRegistration INTERFACE )

target_include_directories( TubeTKRegistration INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR} )

target_link_libraries( TubeTKRegistration INTERFACE
  TubeTKFiltering
  TubeTKNumerics # Required by AnisotropicDiffusiveSparseRegistrationFilter including itkTubeMath.h
  )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_Segmentation_H_Files}
  ${TubeTK_Base_Segmentation_HXX_Files}
  DESTINATION include )

TubeTKMacroInstallPlugins( TubeTKRegistration )

set_property( GLOBAL APPEND PROPERTY TubeTK_TARGETS TubeTKRegistration )
