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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <iostream>

#include "itkTestMain.h" 

void RegisterTests()
{
  REGISTER_TEST( tubeBaseSegmentationVTreePrintTest );
  REGISTER_TEST( itkBlur3DImageFunctionTest );
  REGISTER_TEST( itkMatrixMathTest );
  REGISTER_TEST( itkOptBrentTest );
  REGISTER_TEST( itkOptGoldenMeanTest );
  REGISTER_TEST( itkOptParabolicFitTest );
  REGISTER_TEST( itkOptBrentNDTest );
  REGISTER_TEST( itkRidgeExtractorTest );
  REGISTER_TEST( itkRidgeExtractorTest2 );
  REGISTER_TEST( itkSplineApprox1DTest );
  REGISTER_TEST( itkSplineNDTest );
  REGISTER_TEST( itkUserFuncTest );
}

