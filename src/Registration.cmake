##############################################################################
#
# Library:   TubeTK
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

set( TubeTK_Registration_H_Files
  Registration/itktubeAnisotropicDiffusiveRegistrationFunction.h
  Registration/itktubeDiffusiveRegistrationFilter.h
  Registration/itktubeDiffusiveRegistrationFilterUtils.h
  Registration/itktubeImageToTubeRigidMetric.h
  Registration/itktubeImageToTubeRigidRegistration.h
  Registration/itktubeMeanSquareRegistrationFunction.h
  Registration/itktubeMergeAdjacentImagesFilter.h
  Registration/itktubeTubeExponentialResolutionWeightFunction.h
  Registration/itktubeTubeParametricExponentialResolutionWeightFunction.h
  Registration/itktubeTubeParametricExponentialWithBoundsResolutionWeightFunction.h
  Registration/itktubeTubeToTubeTransformFilter.h )

if( TubeTK_USE_VTK )
  list( APPEND TubeTK_Registration_H_Files
    Registration/itktubeAnisotropicDiffusiveSparseRegistrationFilter.h
    Registration/itktubeAnisotropicDiffusiveRegistrationFilter.h )
endif()

set( TubeTK_Registration_HXX_Files
  Registration/itktubeAnisotropicDiffusiveRegistrationFunction.hxx
  Registration/itktubeDiffusiveRegistrationFilter.hxx
  Registration/itktubeDiffusiveRegistrationFilterUtils.hxx
  Registration/itktubeImageToTubeRigidMetric.hxx
  Registration/itktubeImageToTubeRigidRegistration.hxx
  Registration/itktubeMeanSquareRegistrationFunction.hxx
  Registration/itktubeMergeAdjacentImagesFilter.hxx
  Registration/itktubeTubeToTubeTransformFilter.hxx )

if( TubeTK_USE_VTK )
  list( APPEND TubeTK_Registration_HXX_Files
    Registration/itktubeAnisotropicDiffusiveSparseRegistrationFilter.hxx
    Registration/itktubeAnisotropicDiffusiveRegistrationFilter.hxx )
endif()

set( TubeTK_Registration_CXX_Files )

list( APPEND TubeTK_SRCS
  ${TubeTK_Registration_H_Files}
  ${TubeTK_Registration_HXX_Files}
  ${TubeTK_Registration_CXX_Files} )
