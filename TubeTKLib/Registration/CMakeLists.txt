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

project( TubeTKLibRegistration )

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

if( TubeTKLib_USE_VTK )
  list( APPEND TubeTKLib_Base_Registration_H_Files
    itktubeAnisotropicDiffusiveSparseRegistrationFilter.h
    itktubeAnisotropicDiffusiveRegistrationFilter.h )
endif( TubeTKLib_USE_VTK )

set( TubeTKLib_Base_Registration_HXX_Files
  itktubeAnisotropicDiffusiveRegistrationFunction.hxx
  itktubeDiffusiveRegistrationFilter.hxx
  itktubeDiffusiveRegistrationFilterUtils.hxx
  itktubeImageToTubeRigidMetric.hxx
  itktubeImageToTubeRigidRegistration.hxx
  itktubeMeanSquareRegistrationFunction.hxx
  itktubeMergeAdjacentImagesFilter.hxx
  itktubeTubeToTubeTransformFilter.hxx )

if( TubeTKLib_USE_VTK )
  list( APPEND TubeTKLib_Base_Registration_HXX_Files
    itktubeAnisotropicDiffusiveSparseRegistrationFilter.hxx
    itktubeAnisotropicDiffusiveRegistrationFilter.hxx )
endif( TubeTKLib_USE_VTK )

add_library( ${PROJECT_NAME} INTERFACE )

target_include_directories( ${PROJECT_NAME} INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR} )

target_link_libraries( ${PROJECT_NAME} INTERFACE
  TubeTKLibFiltering
  TubeTKLibNumerics
  )

if( TubeTKLib_BUILD_TESTING )
  add_subdirectory( Testing )
endif( TubeTKLib_BUILD_TESTING )

if( TubeTKLib_INSTALL_DEVELOPMENT )
  install( FILES
    ${TubeTKLib_Base_Segmentation_H_Files}
    ${TubeTKLib_Base_Segmentation_HXX_Files}
    DESTINATION ${TubeTKLib_INSTALL_INCLUDE_DIR}
    COMPONENT Development )
endif()

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS ${PROJECT_NAME} )
