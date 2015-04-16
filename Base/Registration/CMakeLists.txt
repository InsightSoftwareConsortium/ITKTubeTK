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
  itktubeAnisotropicDiffusiveRegistrationFilter.h
  itktubeAnisotropicDiffusiveRegistrationFunction.h
  itktubeAnisotropicDiffusiveSparseRegistrationFilter.h
  itktubeDiffusiveRegistrationFilter.h
  itktubeDiffusiveRegistrationFilterUtils.h
  itktubeImageToTubeRigidMetric.h
  itktubeImageToTubeRigidRegistration.h
  itktubeMeanSquareRegistrationFunction.h
  itktubeTubeExponentialResolutionWeightFunction.h
  itktubeTubeParametricExponentialResolutionWeightFunction.h
  itktubeTubeParametricExponentialWithBoundsResolutionWeightFunction.h
  itktubeTubeToTubeTransformFilter.h )

set( TubeTK_Base_Registration_HXX_Files
  itktubeAnisotropicDiffusiveRegistrationFilter.hxx
  itktubeAnisotropicDiffusiveRegistrationFunction.hxx
  itktubeAnisotropicDiffusiveSparseRegistrationFilter.hxx
  itktubeDiffusiveRegistrationFilter.hxx
  itktubeDiffusiveRegistrationFilterUtils.hxx
  itktubeImageToTubeRigidMetric.hxx
  itktubeImageToTubeRigidRegistration.hxx
  itktubeMeanSquareRegistrationFunction.hxx
  itktubeTubeToTubeTransformFilter.hxx )

add_library( TubeTKRegistration INTERFACE )

target_link_libraries( TubeTKRegistration INTERFACE
  TubeTKCommon TubeTKNumerics TubeTKFiltering
  ${ITK_LIBRARIES}
)

TubeTKMacroInstallPlugins( TubeTKRegistration )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_Segmentation_H_Files}
  ${TubeTK_Base_Segmentation_HXX_Files}
  DESTINATION include )

# Export target
set_property( GLOBAL APPEND PROPERTY TubeTK_TARGETS TubeTKRegistration )
