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

#include "tubetkConfigure.h"

#include "itktubeMetaLDA.h"
#include "itktubeMetaNJetLDA.h"
#include "itktubeMetaClassPDF.h"
#include "itktubeMetaRidgeSeed.h"
#include "itktubeMetaTubeExtractor.h"
#include "itktubePDFSegmenterParzenIO.h"
#ifdef TubeTK_USE_LIBSVM
#  include "itktubeRidgeSeedFilterIO.h"
#endif
#include "itktubeTubeExtractorIO.h"
#include "itktubeTubeXIO.h"

int tubeBaseIOPrintTest( int tubeNotUsed( argc ), char * tubeNotUsed( argv )[] )
{
  itk::tube::MetaLDA metaLDA;
  std::cout << "-------------metaLDA" << std::endl;
  metaLDA.PrintInfo();

  itk::tube::MetaNJetLDA metaNJetLDA;
  std::cout << "-------------metaNJetLDA" << std::endl;
  metaNJetLDA.PrintInfo();

  itk::tube::MetaClassPDF metaClassPDF;
  std::cout << "-------------metaClassPDF" << std::endl;
  metaClassPDF.PrintInfo();

  itk::tube::MetaRidgeSeed metaRidgeSeed;
  std::cout << "-------------metaRidgeSeed" << std::endl;
  metaRidgeSeed.PrintInfo();

  itk::tube::MetaTubeExtractor metaTubeExtractor;
  std::cout << "-------------metaTubeExtractor" << std::endl;
  metaTubeExtractor.PrintInfo();

  typedef itk::Image< float, 3 > ImageType;
  itk::tube::PDFSegmenterParzenIO< ImageType,
    ImageType > pdfSegmenterParzenIO;
  std::cout << "-------------pdfSegmenterParzenIO" << std::endl;
  pdfSegmenterParzenIO.PrintInfo();

#ifdef TubeTK_USE_LIBSVM
  itk::tube::RidgeSeedFilterIO< ImageType, ImageType > ridgeSeedFilterIO;
  std::cout << "-------------ridgeSeedFilterIO" << std::endl;
  ridgeSeedFilterIO.PrintInfo();
#endif

  itk::tube::TubeExtractorIO< ImageType > tubeExtractorIO;
  std::cout << "-------------tubeExtractorIO" << std::endl;
  tubeExtractorIO.PrintInfo();

  itk::tube::TubeXIO< 3 >::Pointer tubeTubeXIO;
  std::cout << "-------------tubeTubeXIO" << tubeTubeXIO << std::endl;

  return EXIT_SUCCESS;
}
