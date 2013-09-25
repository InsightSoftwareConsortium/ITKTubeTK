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

#include "itktubeMetaPDF.h"

int itktubeMetaPDFTest( int argc, char * argv[] )
{
  if( argc != 2 )
    {
    std::cout << "Usage: testname <tempfilename>" << std::endl;
    return EXIT_FAILURE;
    }

  itk::tube::MetaPDF pdf1;

  std::vector< int > dimSize(2);
  dimSize[0] = 10;
  dimSize[1] = 10;
  std::vector< float > binMin(2);
  binMin[0] = -5;
  binMin[1] = 20;
  std::vector< float > binSize(2);
  binSize[0] = 10;
  binSize[1] = 5;
  float data[100];
  for( unsigned int i = 0; i < 100; i++ )
    {
    data[i] = i;
    }
  pdf1.InitializeEssential( 2, dimSize, binMin, binSize, data );
  pdf1.Write( argv[1] );

  itk::tube::MetaPDF pdf2( argv[1] );
  if( pdf1.GetNumberOfFeatures() != pdf2.GetNumberOfFeatures() ||
    pdf1.GetNumberOfBinsPerFeature()[0] !=
      pdf2.GetNumberOfBinsPerFeature()[0] ||
    pdf1.GetNumberOfBinsPerFeature()[1] !=
      pdf2.GetNumberOfBinsPerFeature()[1] ||
    pdf1.GetBinMin()[0] != pdf2.GetBinMin()[0] ||
    pdf1.GetBinMin()[1] != pdf2.GetBinMin()[1] ||
    pdf1.GetBinSize()[0] != pdf2.GetBinSize()[0] ||
    pdf1.GetBinSize()[1] != pdf2.GetBinSize()[1] )
    {
    std::cout << pdf2.GetBinMin()[0] << ", " << pdf2.GetBinMin()[1]
      << std::endl;
    std::cout << pdf2.GetBinSize()[0] << ", " << pdf2.GetBinSize()[1]
      << std::endl;
    std::cout << "Written file and read file do not match" << std::endl;
    return EXIT_FAILURE;
    }
  for( unsigned int i = 0; i < 100; i++ )
    {
    if( pdf2.GetPDF()[i] != i )
      {
      std::cout << "Written and read data does not match" << std::endl;
      return EXIT_FAILURE;
      }
    }

  itk::tube::MetaPDF pdf3( 2, dimSize, binMin, binSize, data );
  pdf3.Write( argv[1] );
  dimSize[1] = 20;
  pdf3.SetNumberOfBinsPerFeature( dimSize );
  pdf3.Read( argv[1] );
  if( pdf1.GetNumberOfFeatures() != pdf3.GetNumberOfFeatures() ||
    pdf1.GetNumberOfBinsPerFeature()[0] !=
      pdf3.GetNumberOfBinsPerFeature()[0] ||
    pdf1.GetNumberOfBinsPerFeature()[1] !=
      pdf3.GetNumberOfBinsPerFeature()[1] ||
    pdf1.GetBinMin()[0] != pdf3.GetBinMin()[0] ||
    pdf1.GetBinMin()[1] != pdf3.GetBinMin()[1] ||
    pdf1.GetBinSize()[0] != pdf3.GetBinSize()[0] ||
    pdf1.GetBinSize()[1] != pdf3.GetBinSize()[1] )
    {
    std::cout << "Re-written file and read file do not match" << std::endl;
    return EXIT_FAILURE;
    }
  for( unsigned int i = 0; i < 100; i++ )
    {
    if( pdf3.GetPDF()[i] != i )
      {
      std::cout << "Re-written and read data does not match" << std::endl;
      return EXIT_FAILURE;
      }
    }


  return EXIT_SUCCESS;
}
