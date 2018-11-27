/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include "tubeTestMain.h"

#include "tubetklibConfigure.h"

#include <iostream>

void RegisterTests( void )
{
  REGISTER_TEST( tubeBaseFilteringPrintTest );
  REGISTER_TEST( itktubeCVTImageFilterTest );
  REGISTER_TEST( itktubeExtractTubePointsSpatialObjectFilterTest );
  REGISTER_TEST( itktubeFFTGaussianDerivativeIFFTFilterTest );
  REGISTER_TEST( itktubeRidgeFFTFilterTest );
  REGISTER_TEST( itktubeSubSampleTubeSpatialObjectFilterTest );
  REGISTER_TEST( itktubeSubSampleTubeTreeSpatialObjectFilterTest );
  REGISTER_TEST( itktubeStructureTensorRecursiveGaussianImageFilterTest );
  REGISTER_TEST( itktubeStructureTensorRecursiveGaussianImageFilterTestNew );
  REGISTER_TEST( itktubeTubeEnhancingDiffusion2DImageFilterTest );
  REGISTER_TEST( itktubeSheetnessMeasureImageFilterTest );
  REGISTER_TEST( itktubeSheetnessMeasureImageFilterTest2 );
  REGISTER_TEST( itktubeShrinkWithBlendingImageFilterTest );
  REGISTER_TEST( itktubeAnisotropicHybridDiffusionImageFilterTest );
  REGISTER_TEST( itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest );
  REGISTER_TEST( itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest );
  REGISTER_TEST( itktubeTortuositySpatialObjectFilterTest );

  #if defined( TubeTKLib_USE_ARRAYFIRE )
    REGISTER_TEST( itktubeGPUArrayFireGaussianDerivativeFilterTest );
  #endif
}
