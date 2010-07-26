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

#include "itkBlur3DImageFunction.h"
#include "itkMatrixMath.h"
#include "itkOptBrent1D.h"
#include "itkOptGoldenMean1D.h"
#include "itkOptimizer1D.h"
#include "itkOptimizerND.h"
#include "itkOptParabolicFit1D.h"
#include "itkRadiusExtractor.h"
#include "itkRidgeExtractor.h"
#include "itkSpline1D.h"
#include "itkSplineND.h"
#include "itkSplineApproximation1D.h"
#include "itkTubeExtractor.h"
#include "itkTubeNetExtractor.h"
#include "UserFunc.h"

int main ( int , char ** )
{
  return EXIT_SUCCESS;
}

