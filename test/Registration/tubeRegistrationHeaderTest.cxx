/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include <itkImage.h>
#include <itkTubeSpatialObject.h>

#include "itkAffineImageToImageRegistrationMethod.h"
#include "itkAnisotropicSimilarity3DTransform.h"
#include "itkAnisotropicSimilarityLandmarkBasedTransformInitializer.h"
#include "itkBSplineImageToImageRegistrationMethod.h"
#include "itkImageRegionMomentsCalculator.h"
#include "itkImageRegionSplitter.h"
#include "itkImageToImageRegistrationHelper.h"
#include "itkImageToImageRegistrationMethod.h"
#include "itkInitialImageToImageRegistrationMethod.h"
#include "itkOptimizedImageToImageRegistrationMethod.h"
#include "itkRigidImageToImageRegistrationMethod.h"
#include "itkScaleSkewAngle2DImageToImageRegistrationMethod.h"
#include "itkScaleSkewAngle2DTransform.h"
#include "itkScaleSkewVersor3DImageToImageRegistrationMethod.h"
#include "itkSimilarity2DTransform.h"
//#include "itktubeAnisotropicDiffusiveRegistrationFilter.h"
#include "itktubeAnisotropicDiffusiveRegistrationFunction.h"
//#include "itktubeAnisotropicDiffusiveSparseRegistrationFilter.h"
#include "itktubeDiffusiveRegistrationFilter.h"
#include "itktubeDiffusiveRegistrationFilterUtils.h"
#include "itktubeInitialSpatialObjectToImageRegistrationMethod.h"
#include "itktubeMeanSquareRegistrationFunction.h"
#include "itktubeMergeAdjacentImagesFilter.h"
#include "itktubeOptimizedSpatialObjectToImageRegistrationMethod.h"
#include "itktubePointBasedSpatialObjectToImageMetric.h"
#include "itktubePointBasedSpatialObjectTransformFilter.h"
#include "itktubeRigidSpatialObjectToImageRegistrationMethod.h"
#include "itktubeScaleSkewAngle2DSpatialObjectToImageRegistrationMethod.h"
#include "itktubeScaleSkewVersor3DSpatialObjectToImageRegistrationMethod.h"
#include "itktubeSpatialObjectToImageMetric.h"
#include "itktubeSpatialObjectToImageRegistrationHelper.h"
#include "itktubeSpatialObjectToImageRegistrationMethod.h"

int tubeRegistrationHeaderTest( int itkNotUsed( argc ), char * itkNotUsed(
  argv )[] )
{
  return EXIT_SUCCESS;
}
