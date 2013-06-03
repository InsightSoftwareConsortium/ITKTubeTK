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

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <itkImage.h>
#include "itkFilterWatcher.h"
#include <itkMacro.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <itkTubePDFSegmenter.h>

int itkTubePDFSegmenterTest(int argc, char* argv[] )
{
  if( argc != 7 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
      << " inputImage inputMask force outputProbImage0 outputProbImage1 outputMask"
      << std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  const unsigned int Dimension = 2;

  // Define the pixel type
  typedef float PixelType;

  // Declare the types of the images
  typedef itk::Image<PixelType, Dimension>  ImageType;

  // Declare the reader and writer
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;


  // Declare the type for the Filter
  typedef itk::tube::PDFSegmenter< ImageType, 1, ImageType > FilterType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch(itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during input read:\n"  << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer inputImage = reader->GetOutput();

  // Create the reader and writer
  ReaderType::Pointer maskReader = ReaderType::New();
  maskReader->SetFileName( argv[3] );
  try
    {
    maskReader->Update();
    }
  catch(itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during input read:\n"  << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer maskImage = maskReader->GetOutput();

  FilterType::Pointer filter = FilterType::New();
  filter->SetInputVolume( 0, inputImage );
  filter->SetLabelmap( maskImage );
  filter->AddObjectId( 255 );
  filter->AddObjectId( 127 );
  filter->SetVoidId( 0 );
  filter->SetErodeRadius( 1 );
  filter->SetHoleFillIterations( 4 );
  filter->SetProbabilitySmoothingStandardDeviation( 1 );
  filter->SetFprWeight( 1.0 );
  filter->SetDraft( false );
  filter->SetReclassifyObjectMask( true );
  filter->SetReclassifyNotObjectMask( true );
  if( argv[2][0] == 't' || argv[2][0] == 'T' || argv[2][0] == '1' )
    {
    filter->SetForceClassification( true );
    }
  else
    {
    filter->SetForceClassification( false );
    }
  filter->Update();
  filter->ClassifyImages();

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetUseCompression( true );
  writer->SetInput( filter->GetClassProbabilityVolume(0) );
  try
    {
    writer->Update();
    }
  catch(itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during write:\n"  << e;
    return EXIT_FAILURE;
    }

  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetFileName( argv[5] );
  writer2->SetUseCompression( true );
  writer2->SetInput( filter->GetClassProbabilityVolume(1) );
  try
    {
    writer2->Update();
    }
  catch(itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during write:\n"  << e;
    return EXIT_FAILURE;
    }

  WriterType::Pointer maskWriter = WriterType::New();
  maskWriter->SetFileName( argv[6] );
  maskWriter->SetUseCompression( true );
  maskWriter->SetInput( filter->GetLabelmap() );
  try
    {
    maskWriter->Update();
    }
  catch(itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during write:\n"  << e;
    return EXIT_FAILURE;
    }

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
