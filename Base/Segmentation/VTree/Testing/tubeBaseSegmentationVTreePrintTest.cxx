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
#define ITK_LEAN_AND_MEAN

#include "itkOrientedImage.h"

#include "itkTubeBlurImageFunction.h"
#include "itkTubeRadiusExtractor.h"
#include "itkTubeRidgeExtractor.h"

int tubeBaseSegmentationVTreePrintTest( int, char* [] )
{
  typedef itk::OrientedImage< float, 3 > ImageType;

  itk::tube::BlurImageFunction< ImageType >::Pointer
    blurObject =
    itk::tube::BlurImageFunction< ImageType >::New();
  std::cout << "-------------itkTubeBlurImageFunction"
            << blurObject
            << std::endl;

  itk::tube::RadiusExtractor< ImageType >::Pointer
    radiusObject =
    itk::tube::RadiusExtractor< ImageType >::New();
  std::cout << "-------------itkTubeRadiusExtractor"
            << radiusObject
            << std::endl;

  itk::tube::RidgeExtractor< ImageType >::Pointer
    ridgeObject =
    itk::tube::RidgeExtractor< ImageType >::New();
  std::cout << "-------------itkTubeRidgeExtractor"
            << ridgeObject
            << std::endl;

  return EXIT_SUCCESS;
}
