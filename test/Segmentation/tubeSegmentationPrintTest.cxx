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

#include "tubetkConfigure.h"

#include "itktubePDFSegmenterParzen.h"
#include "itktubeRadiusExtractor2.h"
#include "itktubeRidgeExtractor.h"
#include "itktubeRidgeSeedFilter.h"
#include "itktubeTubeExtractor.h"

#include <itkImage.h>

int tubeSegmentationPrintTest( int itkNotUsed( argc ),
  char * itkNotUsed( argv )[] )
{
  typedef itk::Image< float, 2 >         ImageType;
  typedef itk::Image< unsigned char, 2 > CharImageType;

  itk::tube::PDFSegmenterParzen< ImageType, ImageType >::Pointer
    pdfSegmenterParzen = itk::tube::PDFSegmenterParzen< ImageType,
      ImageType >::New();
  std::cout << "-------------itktubePDFSegmenterParzen" << pdfSegmenterParzen
    << std::endl;

  itk::tube::RadiusExtractor2< ImageType >::Pointer
    radius2Object = itk::tube::RadiusExtractor2< ImageType >::New();
  std::cout << "-------------itktubeRadiusExtractor2" << radius2Object
    << std::endl;

  itk::tube::RidgeExtractor< ImageType >::Pointer
    ridgeObject = itk::tube::RidgeExtractor< ImageType >::New();
  std::cout << "-------------itktubeRidgeExtractor" << ridgeObject
    << std::endl;

  itk::tube::RidgeSeedFilter< ImageType, CharImageType >::Pointer
    seedObject = itk::tube::RidgeSeedFilter< ImageType,
    CharImageType >::New();
  std::cout << "-------------itktubeRidgeSeedFilter" << seedObject
    << std::endl;

  itk::tube::TubeExtractor< ImageType >::Pointer tubeObject =
    itk::tube::TubeExtractor< ImageType >::New();
  std::cout << "-------------itktubeTubeExtractor" << tubeObject <<
  std::endl;

  return EXIT_SUCCESS;
}
