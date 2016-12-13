/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#include "tubeBaseNumericsTestsCLP.h"
#include "tubeTestMain.h"

// HACK Redefine 'main' to fix missing symbol error
#undef main
#define main ModuleEntryPoint

#include <iostream>

void RegisterTests( void )
{
  REGISTER_TEST( tubeBaseNumericsPrintTest );
  REGISTER_TEST( itktubeBlurImageFunctionTest );
  REGISTER_TEST( itktubeImageRegionMomentsCalculatorTest );
  REGISTER_TEST( itktubeJointHistogramImageFunctionTest );
  REGISTER_TEST( itktubeNJetBasisFeatureVectorGeneratorTest );
  REGISTER_TEST( itktubeNJetFeatureVectorGeneratorTest );
  REGISTER_TEST( itktubeNJetImageFunctionTest );
  REGISTER_TEST( itktubeSingleValuedCostFunctionImageSourceTest );
  REGISTER_TEST( itktubeRecordOptimizationParameterProgressionCommandTest );
  REGISTER_TEST( itktubeRidgeBasisFeatureVectorGeneratorTest );
  REGISTER_TEST( itktubeRidgeFFTFeatureVectorGeneratorTest );
  REGISTER_TEST( itktubeVotingResampleImageFunctionTest );
  REGISTER_TEST( tubeBrentOptimizer1DTest );
  REGISTER_TEST( tubeBrentOptimizerNDTest );
  REGISTER_TEST( tubeGoldenMeanOptimizer1DTest );
  REGISTER_TEST( tubeMatrixMathTest );
  REGISTER_TEST( tubeParabolicFitOptimizer1DTest );
  REGISTER_TEST( tubeSplineApproximation1DTest );
  REGISTER_TEST( tubeSplineNDTest );
  REGISTER_TEST( tubeTubeMathTest );
  REGISTER_TEST( tubeUserFunctionTest );
}
