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

// ITK includes
#include <itkImageFileReader.h>

// TubeTKITK includes
#include "tubeComputeImageSimilarityMetrics.h"

// TubeTK includes
#include "tubeMessage.h"

#include "ComputeImageSimilarityMetricsCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // typedefs
  typedef TPixel                                              PixelType;
  typedef itk::Image< PixelType, VDimension >                 ImageType;
  typedef itk::ImageFileReader< ImageType >                   ReaderType;
  typedef tube::ComputeImageSimilarityMetrics< ImageType >    FilterType;

  // read input image 1
  typename ReaderType::Pointer reader1 = ReaderType::New();

  try
    {
    reader1->SetFileName( inputVolume1.c_str() );
    reader1->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error reading input image 1: "
                        + std::string( err.GetDescription() ) );
    return EXIT_FAILURE;
    }

  // read input image 2
  typename ReaderType::Pointer reader2 = ReaderType::New();

  try
    {
    reader2->SetFileName( inputVolume2.c_str() );
    reader2->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error reading input image 2: "
                        + std::string( err.GetDescription() ) );
    return EXIT_FAILURE;
    }

  // compute image similarity
  typename FilterType::Pointer similarityCalculator = FilterType::New();

  similarityCalculator->SetInput1( reader1->GetOutput() );
  similarityCalculator->SetInput2( reader2->GetOutput() );
  similarityCalculator->SetSamplingRate( samplingRate );
  similarityCalculator->SetUseCorrelation( correlation );

  std::cout << similarityCalculator->GetOutput() << std::endl;

  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( inputVolume1, argc, argv );
}
