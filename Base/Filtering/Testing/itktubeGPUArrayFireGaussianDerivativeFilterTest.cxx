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

#include "itktubeGPUArrayFireGaussianDerivativeFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int itktubeGPUArrayFireGaussianDerivativeFilterTest( int argc, char * argv[] )
{
  if( argc < 7 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " orderX orderY scaleX scaleY inputImage outputImage"
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
  typedef itk::tube::GPUArrayFireGaussianDerivativeFilter< ImageType, ImageType >
    FunctionType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[5] );
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

  FunctionType::OrdersType orders;
  orders[0] = atoi( argv[1] );
  orders[1] = atoi( argv[2] );
  func->SetOrders( orders );
  FunctionType::SigmasType sigmas;
  sigmas[0] = atof( argv[3] );
  sigmas[1] = atof( argv[4] );
  func->SetSigmas( sigmas );
  func->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[6] );
  writer->SetInput( func->GetOutput() );
  writer->SetUseCompression( true );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during write:\n"  << e;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
