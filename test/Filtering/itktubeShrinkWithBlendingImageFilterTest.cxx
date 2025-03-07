/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#include "itktubeShrinkWithBlendingImageFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>

int itktubeShrinkWithBlendingImageFilterTest( int argc, char * argv[] )
{
  if( argc != 4 )
    {
    std::cout << "Missing arguments." << std::endl;
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0] << " inputImage outputImage pointImage"
      << std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  enum { Dimension = 2 };

  // Define the pixel type
  using PixelType = float;

  // Declare the types of the images
  using ImageType = itk::Image< PixelType, Dimension >;

  // Declare the reader and writer
  using ReaderType = itk::ImageFileReader< ImageType >;
  using WriterType = itk::ImageFileWriter< ImageType >;

  using PointPixelType = itk::Vector< float, Dimension >;
  using PointImageType = itk::Image< PointPixelType, Dimension >;
  using PointImageWriterType = itk::ImageFileWriter< PointImageType >;

  // Declare the type for the Filter
  using FilterType = itk::tube::ShrinkWithBlendingImageFilter< ImageType, ImageType >;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception caught during input read:" << std::endl << e
      << std::endl;
    return EXIT_FAILURE;
    }
  ImageType::Pointer inputImage = reader->GetOutput();

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inputImage );
  FilterType::ShrinkFactorsType shrinkFactors;
  shrinkFactors.Fill( 4 );
  filter->SetShrinkFactors( shrinkFactors );
  FilterType::InputIndexType overlap;
  overlap.Fill( 1 );
  filter->SetOverlap( overlap );
  filter->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetUseCompression( true );
  writer->SetInput( filter->GetOutput() );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception caught during write:" << std::endl << e
      << std::endl;
    return EXIT_FAILURE;
    }

  PointImageWriterType::Pointer pointImageWriter =
    PointImageWriterType::New();
  pointImageWriter->SetFileName( argv[3] );
  pointImageWriter->SetUseCompression( true );
  pointImageWriter->SetInput( filter->GetOutputMipPointImage() );
  try
    {
    pointImageWriter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception caught during image index write:" << std::endl
      << e << std::endl;
    return EXIT_FAILURE;
    }

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
