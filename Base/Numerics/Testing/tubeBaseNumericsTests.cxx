/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include <iostream>
#include "tubeTestMain.h"

void RegisterTests( void )
{
  REGISTER_TEST( tubeBaseNumericsPrintTest );
  REGISTER_TEST( itkImageRegionMomentsCalculatorTest );
  REGISTER_TEST( itkJointHistogramImageFunctionTest );
  REGISTER_TEST( itkTubeNJetImageFunctionTest );
  REGISTER_TEST( itkTubeBlurImageFunctionTest );
  REGISTER_TEST( itkTubeNJetFeatureVectorGeneratorTest );
  REGISTER_TEST( itkTubeNJetBasisFeatureVectorGeneratorTest );
  REGISTER_TEST( itkTubeRidgeFeatureVectorGeneratorTest );
  REGISTER_TEST( itkTubeRidgeBasisFeatureVectorGeneratorTest );
  REGISTER_TEST( itkVotingResampleImageFunctionTest );
  REGISTER_TEST( tubeMatrixMathTest );
  REGISTER_TEST( tubeOptBrentNDTest );
  REGISTER_TEST( tubeOptBrentTest );
  REGISTER_TEST( tubeOptGoldenMeanTest );
  REGISTER_TEST( tubeOptParabolicFitTest );
  REGISTER_TEST( tubeSplineApprox1DTest );
  REGISTER_TEST( tubeSplineNDTest );
  REGISTER_TEST( tubeUserFuncTest );
}
