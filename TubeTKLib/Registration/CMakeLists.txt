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

set( TubeTKLib_Base_Registration_H_Files
  itktubeAnisotropicDiffusiveRegistrationFunction.h
  itktubeDiffusiveRegistrationFilter.h
  itktubeDiffusiveRegistrationFilterUtils.h
  itktubeImageToTubeRigidMetric.h
  itktubeImageToTubeRigidRegistration.h
  itktubeMeanSquareRegistrationFunction.h
  itktubeMergeAdjacentImagesFilter.h
  itktubeTubeExponentialResolutionWeightFunction.h
  itktubeTubeParametricExponentialResolutionWeightFunction.h
  itktubeTubeParametricExponentialWithBoundsResolutionWeightFunction.h
  itktubeTubeToTubeTransformFilter.h )

if( TubeTK_USE_VTK )
  list( APPEND TubeTKLib_Base_Registration_H_Files
    itktubeAnisotropicDiffusiveSparseRegistrationFilter.h
    itktubeAnisotropicDiffusiveRegistrationFilter.h )
endif()

set( TubeTKLib_Base_Registration_HXX_Files
  itktubeAnisotropicDiffusiveRegistrationFunction.hxx
  itktubeDiffusiveRegistrationFilter.hxx
  itktubeDiffusiveRegistrationFilterUtils.hxx
  itktubeImageToTubeRigidMetric.hxx
  itktubeImageToTubeRigidRegistration.hxx
  itktubeMeanSquareRegistrationFunction.hxx
  itktubeMergeAdjacentImagesFilter.hxx
  itktubeTubeToTubeTransformFilter.hxx )

if( TubeTK_USE_VTK )
  list( APPEND TubeTKLib_Base_Registration_HXX_Files
    itktubeAnisotropicDiffusiveSparseRegistrationFilter.hxx
    itktubeAnisotropicDiffusiveRegistrationFilter.hxx )
endif()

add_library( TubeTKLibRegistration INTERFACE )

target_include_directories( TubeTKLibRegistration INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" )

target_link_libraries( TubeTKLibRegistration INTERFACE
  TubeTKLibFiltering
  TubeTKLibNumerics
  )

include_directories(
  ${ITK_INCLUDE_DIRS} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTKLib_Base_Segmentation_H_Files}
  ${TubeTKLib_Base_Segmentation_HXX_Files}
  DESTINATION ${TubeTK_INSTALL_INCLUDE_DIR}
  COMPONENT Development )

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS
  TubeTKLibRegistration )
