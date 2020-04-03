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

#include "itktubeRidgeFFTFeatureVectorGenerator.h"

int itktubeRidgeFFTFeatureVectorGeneratorTest( int argc, char * argv[] )
{
  if( argc != 4 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
      << " inputImage outputF0Image outputF1Image"
      << std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  enum { Dimension = 2 };

  // Define the pixel type
  typedef float PixelType;

  // Declare the types of the images
  typedef itk::Image<PixelType, Dimension>  ImageType;

  // Declare the reader and writer
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;


  // Declare the type for the Filter
  typedef itk::tube::RidgeFFTFeatureVectorGenerator< ImageType > FilterType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception caught during input read:" << std::endl << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer inputImage = reader->GetOutput();

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inputImage );

  FilterType::RidgeScalesType scales( 2 );
  scales[0] = 0.4;
  scales[1] = 1.0;
  filter->SetScales( scales );
  std::cout << filter << std::endl;

  filter->SetUpdateWhitenStatisticsOnUpdate( true );
  filter->Update();

  WriterType::Pointer imageFeature0Writer = WriterType::New();
  imageFeature0Writer->SetFileName( argv[2] );
  imageFeature0Writer->SetUseCompression( true );
  imageFeature0Writer->SetInput( filter->GetFeatureImage( 0 ) );
  try
    {
    imageFeature0Writer->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception caught during write:" << std::endl << e;
    return EXIT_FAILURE;
    }

  WriterType::Pointer imageFeature1Writer = WriterType::New();
  imageFeature1Writer->SetFileName( argv[3] );
  imageFeature1Writer->SetUseCompression( true );
  imageFeature1Writer->SetInput( filter->GetFeatureImage( 6 ) );
  try
    {
    imageFeature1Writer->Update();
    }
  catch ( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during write:" << std::endl << e;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
