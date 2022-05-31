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

#include "itktubeRidgeFFTFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int itktubeRidgeFFTFilterTest( int argc, char * argv[] )
{
  if( argc != 7 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " scale inputImage ridgenessImage roundnessImage curvatureImage levelnessImage" << std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  enum { Dimension = 3 };

  // Define the pixel type
  typedef float PixelType;

  // Declare the types of the images
  typedef itk::Image<PixelType, Dimension>  ImageType;

  // Declare the reader and writer
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  // Declare the type for the Filter
  typedef itk::tube::RidgeFFTFilter< ImageType > FunctionType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during read:\n"  << e;
    return EXIT_FAILURE;
    }

  ImageType::Pointer inputImage = reader->GetOutput();

  FunctionType::Pointer func = FunctionType::New();
  func->SetInput( inputImage );

  func->SetScale( atof( argv[1] ) );
  func->Update();

  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName( argv[3] );
  writer->SetInput( func->GetRidgeness() );
  writer->SetUseCompression( true );
  writer->Update();

  writer->SetFileName( argv[4] );
  writer->SetInput( func->GetRoundness() );
  writer->SetUseCompression( true );
  writer->Update();

  writer->SetFileName( argv[5] );
  writer->SetInput( func->GetCurvature() );
  writer->SetUseCompression( true );
  writer->Update();

  writer->SetFileName( argv[6] );
  writer->SetInput( func->GetLevelness() );
  writer->SetUseCompression( true );
  writer->Update();

  func->SetScale( atof( argv[1] ) + 1 );
  func->Update();
  std::string tmpStr = argv[3];
  tmpStr = tmpStr + "2.mha";
  writer->SetFileName( tmpStr.c_str() );
  writer->SetInput( func->GetRidgeness() );
  writer->SetUseCompression( true );
  writer->Update();


  return EXIT_SUCCESS;
}
