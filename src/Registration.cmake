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
  Registration/itkAffineImageToImageRegistrationMethod.h
  Registration/itkAnisotropicSimilarity3DTransform.h
  Registration/itkAnisotropicSimilarityLandmarkBasedTransformInitializer.h
  Registration/itkBSplineImageToImageRegistrationMethod.h
  Registration/itkImageRegionMomentsCalculator.h
  Registration/itkImageRegionSplitter.h
  Registration/itkImageToImageRegistrationHelper.h
  Registration/itkImageToImageRegistrationMethod.h
  Registration/itkInitialImageToImageRegistrationMethod.h
  Registration/itkOptimizedImageToImageRegistrationMethod.h
  Registration/itkRigidImageToImageRegistrationMethod.h
  Registration/itkScaleSkewAngle2DImageToImageRegistrationMethod.h
  Registration/itkScaleSkewAngle2DTransform.h
  Registration/itkScaleSkewVersor3DImageToImageRegistrationMethod.h
  Registration/itkSimilarity2DTransform.h
  Registration/itktubeAffineSpatialObjectToImageRegistrationMethod.h
  Registration/itktubeAnisotropicDiffusiveRegistrationFunction.h
  Registration/itktubeDiffusiveRegistrationFilter.h
  Registration/itktubeDiffusiveRegistrationFilterUtils.h
  Registration/itktubeInitialSpatialObjectToImageRegistrationMethod.h
  Registration/itktubeMeanSquareRegistrationFunction.h
  Registration/itktubeMergeAdjacentImagesFilter.h
  Registration/itktubeOptimizedSpatialObjectToImageRegistrationMethod.h
  Registration/itktubePointBasedSpatialObjectToImageMetric.h
  Registration/itktubePointBasedSpatialObjectTransformFilter.h
  Registration/itktubeRigidSpatialObjectToImageRegistrationMethod.h
  Registration/itktubeScaleSkewAngle2DSpatialObjectToImageRegistrationMethod.h
  Registration/itktubeScaleSkewVersor3DSpatialObjectToImageRegistrationMethod.h
  Registration/itktubeSpatialObjectToImageMetric.h
  Registration/itktubeSpatialObjectToImageRegistrationHelper.h
  Registration/itktubeSpatialObjectToImageRegistrationMethod.h
  )

if( TubeTK_USE_VTK )
  list( APPEND TubeTK_Registration_H_Files
    Registration/itktubeAnisotropicDiffusiveSparseRegistrationFilter.h
    Registration/itktubeAnisotropicDiffusiveRegistrationFilter.h )
endif()

set( TubeTK_Registration_HXX_Files
  Registration/itkAffineImageToImageRegistrationMethod.hxx
  Registration/itkAnisotropicSimilarity3DTransform.hxx
  Registration/itkAnisotropicSimilarityLandmarkBasedTransformInitializer.hxx
  Registration/itkBSplineImageToImageRegistrationMethod.hxx
  Registration/itkImageRegionMomentsCalculator.hxx
  Registration/itkImageRegionSplitter.hxx
  Registration/itkImageToImageRegistrationHelper.hxx
  Registration/itkImageToImageRegistrationMethod.hxx
  Registration/itkInitialImageToImageRegistrationMethod.hxx
  Registration/itkOptimizedImageToImageRegistrationMethod.hxx
  Registration/itkRigidImageToImageRegistrationMethod.hxx
  Registration/itkScaleSkewAngle2DImageToImageRegistrationMethod.hxx
  Registration/itkScaleSkewAngle2DTransform.hxx
  Registration/itkScaleSkewVersor3DImageToImageRegistrationMethod.hxx
  Registration/itkSimilarity2DTransform.hxx
  Registration/itktubeAffineSpatialObjectToImageRegistrationMethod.hxx
  Registration/itktubeAnisotropicDiffusiveRegistrationFunction.hxx
  Registration/itktubeDiffusiveRegistrationFilter.hxx
  Registration/itktubeDiffusiveRegistrationFilterUtils.hxx
  Registration/itktubeInitialSpatialObjectToImageRegistrationMethod.hxx
  Registration/itktubeMeanSquareRegistrationFunction.hxx
  Registration/itktubeMergeAdjacentImagesFilter.hxx
  Registration/itktubeOptimizedSpatialObjectToImageRegistrationMethod.hxx
  Registration/itktubePointBasedSpatialObjectToImageMetric.hxx
  Registration/itktubePointBasedSpatialObjectTransformFilter.hxx
  Registration/itktubeRigidSpatialObjectToImageRegistrationMethod.hxx
  Registration/itktubeScaleSkewAngle2DSpatialObjectToImageRegistrationMethod.hxx
  Registration/itktubeScaleSkewVersor3DSpatialObjectToImageRegistrationMethod.hxx
  Registration/itktubeSpatialObjectToImageRegistrationHelper.hxx
  Registration/itktubeSpatialObjectToImageRegistrationMethod.hxx
  )

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
