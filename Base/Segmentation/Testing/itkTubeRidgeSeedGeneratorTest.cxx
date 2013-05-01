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

#include <itkImage.h>
#include "itkFilterWatcher.h"
#include <itkExceptionObject.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkTubeRidgeSeedGenerator.h>

int itkTubeRidgeSeedGeneratorTest(int argc, char* argv[] )
{
  if( argc != 4 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImage maskImage outputImage"
      << std::endl;
    return EXIT_FAILURE;
    }

  bool error = false;

  // Define the dimension of the images
  const unsigned int Dimension = 3;

  // Define the pixel type
  typedef float                                   PixelType;
  // Declare the type of the images
  typedef itk::Image< PixelType, Dimension >      ImageType;

  // Define the mask pixel type
  typedef unsigned char                           MaskPixelType;
  // Declare the type of the mask image
  typedef itk::Image< MaskPixelType, Dimension >  MaskImageType;

  // Declare the reader and writer
  typedef itk::ImageFileReader< ImageType >       ReaderType;

  typedef itk::ImageFileWriter< ImageType >       WriterType;

  // Declare the mask reader and writer
  typedef itk::ImageFileReader< MaskImageType >   MaskReaderType;

  // Declare the type for the Filter
  typedef itk::tube::RidgeSeedGenerator< ImageType, MaskImageType >
                                                  FunctionType;

  // Create the reader
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch(itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during read:\n"  << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer inputImage = reader->GetOutput();

  // Create the mask reader
  MaskReaderType::Pointer maskReader = MaskReaderType::New();
  maskReader->SetFileName( argv[2] );
  try
    {
    maskReader->Update();
    }
  catch(itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during mask read:\n"  << e;
    return EXIT_FAILURE;
    }
  MaskImageType::Pointer maskImage = maskReader->GetOutput();

  FunctionType::Pointer func = FunctionType::New();
  func->SetRidgeImage( inputImage );

  func->SetLabelmap( maskImage );
  func->SetObjectId( 255 );
  func->AddObjectId( 127 );

  FunctionType::RidgeScalesType scales(5);
  scales[0] = 0.3;
  scales[1] = 0.5;
  scales[2] = 0.75;
  scales[3] = 1.0;
  scales[4] = 1.5;
  func->SetScales( scales );

  func->Update();

  func->SetLabelmap( NULL );
  func->UpdateLDAImages();

  char filename[4096];
  WriterType::Pointer writer = WriterType::New();
  writer->SetUseCompression( true );
  std::cout << "Number of LDA = " << func->GetNumberOfLDA() << std::endl;
  for( unsigned int i=0; i<3; i++ )
    {
    sprintf( filename, "%s.lda%02u.mha", argv[3], i );
    writer->SetFileName( filename );
    writer->SetInput( func->GetLDAImage(i) );
    std::cout << "LDA" << i << " = " << func->GetLDAValue(i) << " : "
      << func->GetLDAVector(i) << std::endl;
    try
      {
      writer->Update();
      }
    catch(itk::ExceptionObject& e)
      {
      std::cerr << "Exception caught during write:\n"  << e;
      error = true;
      }
    }

  for( unsigned int i=0; i<func->GetNumberOfFeatures(); i++ )
    {
    sprintf( filename, "%s.f%02u.mha", argv[3], i );
    writer->SetFileName( filename );
    writer->SetInput( func->GetFeatureImage(i) );
    try
      {
      writer->Update();
      }
    catch(itk::ExceptionObject& e)
      {
      std::cerr << "Exception caught during write:\n"  << e;
      error = true;
      }
    }


  // All objects should be automatically destroyed at this point
  if( !error )
    {
    return EXIT_SUCCESS;
    }
  else
    {
    return EXIT_FAILURE;
    }
}
