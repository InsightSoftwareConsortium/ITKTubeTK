/*=========================================================================

Library:   TubeTKLib

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

#include "itktubeBasisFeatureVectorGenerator.h"
#include "itktubeBlurImageFunction.h"
#include "itktubeComputeImageSimilarityMetrics.h"
#include "itktubeFeatureVectorGenerator.h"
#include "itktubeImageRegionMomentsCalculator.h"
#include "itktubeJointHistogramImageFunction.h"
#include "itktubeNJetFeatureVectorGenerator.h"
#include "itktubeNJetImageFunction.h"
#include "itktubeRidgeFFTFeatureVectorGenerator.h"
#include "itktubeSingleValuedCostFunctionImageSource.h"
#include "itktubeVectorImageToListGenerator.h"
#include "itktubeVotingResampleImageFunction.h"
#include "tubeBrentOptimizer1D.h"
#include "tubeGoldenMeanOptimizer1D.h"
#include "tubeMatrixMath.h"
#include "tubeOptimizer1D.h"
#include "tubeOptimizerND.h"
#include "tubeParabolicFitOptimizer1D.h"
#include "tubeSpline1D.h"
#include "tubeSplineApproximation1D.h"
#include "tubeSplineND.h"
#include "tubeUserFunction.h"

#include <iostream>

int tubeNumericsHeaderTest ( int itkNotUsed( argc ), char * itkNotUsed( argv )[] )
{
  return EXIT_SUCCESS;
}
