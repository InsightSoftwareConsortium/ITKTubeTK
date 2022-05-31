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

#include "itktubeNJetFeatureVectorGenerator.h"

int itktubeNJetFeatureVectorGeneratorTest( int argc, char * argv[] )
{
  if( argc != 5 )
    {
    std::cout << "Missing arguments." << std::endl;
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0]
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
  typedef itk::tube::NJetFeatureVectorGenerator< ImageType > FilterType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception caught during input read:" << std::endl << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer inputImage = reader->GetOutput();

  // Create the reader and writer
  ReaderType::Pointer maskReader = ReaderType::New();
  maskReader->SetFileName( argv[2] );
  try
    {
    maskReader->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception caught during input read:" << std::endl << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer maskImage = maskReader->GetOutput();

  FilterType::NJetScalesType scales( 2 );
  scales[0] = 2;
  scales[1] = 4;
  FilterType::NJetScalesType scales2( 1 );
  scales2[0] = 2;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inputImage );
  filter->SetZeroScales( scales );
  filter->SetFirstScales( scales );
  filter->SetSecondScales( scales2 );
  filter->SetRidgeScales( scales2 );
  std::cout << filter << std::endl;

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetUseCompression( true );
  writer->SetInput( filter->GetFeatureImage( 0 ) );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception caught during write:" << std::endl << e;
    return EXIT_FAILURE;
    }

  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetFileName( argv[4] );
  writer2->SetUseCompression( true );
  writer2->SetInput( filter->GetFeatureImage( 1 ) );
  try
    {
    writer2->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception caught during write:" << std::endl << e;
    return EXIT_FAILURE;
    }

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
