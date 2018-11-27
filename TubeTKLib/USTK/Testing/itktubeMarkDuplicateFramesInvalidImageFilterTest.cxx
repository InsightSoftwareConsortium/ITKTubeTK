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

#include <itkImageFileWriter.h>
#include <itkArchetypeSeriesFileNames.h>
#include <itkTestingMacros.h>
#include <itkRGBToLuminanceImageFilter.h>
#include <itkMetaImageIO.h>

#include "itktubeInnerOpticToPlusImageReader.h"
#include "itktubeMarkDuplicateFramesInvalidImageFilter.h"

int itktubeMarkDuplicateFramesInvalidImageFilterTest( int argc, char * argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: "
              << argv[0]
              << "innerOpticMetadata "
              << "outputImage.mha "
              << std::endl;
    return EXIT_FAILURE;
    }
  const char * innerOpticMetadata = argv[1];
  const char * outputImageFile = argv[2];

  typedef itk::tube::InnerOpticToPlusImageReader ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( innerOpticMetadata );
  // Make sure the MetaDataDictionary is populated.
  TRY_EXPECT_NO_EXCEPTION( reader->Update() );
  typedef ReaderType::OutputImageType RGBImageType;
  RGBImageType::Pointer inputImage = reader->GetOutput();
  inputImage->DisconnectPipeline();

  typedef itk::Image< ReaderType::PixelComponentType, ReaderType::ImageDimension >
    OutputImageType;

  typedef itk::RGBToLuminanceImageFilter< RGBImageType, OutputImageType >
    LuminanceFilterType;
  LuminanceFilterType::Pointer luminanceFilter =
    LuminanceFilterType::New();
  luminanceFilter->SetInput( inputImage );

  typedef itk::tube::MarkDuplicateFramesInvalidImageFilter< OutputImageType >
    DuplicateFilterType;
  DuplicateFilterType::Pointer duplicateFilter = DuplicateFilterType::New();
  duplicateFilter->SetInput( luminanceFilter->GetOutput() );
  duplicateFilter->SetTolerance( 3 );
  TEST_EXPECT_EQUAL( duplicateFilter->GetTolerance(), 3 );
  duplicateFilter->SetFractionalThreshold( 0.7 );
  TEST_EXPECT_EQUAL( duplicateFilter->GetFractionalThreshold(), 0.7 );
  duplicateFilter->SetInputMetaDataDictionary(
    &( inputImage->GetMetaDataDictionary() ) );
  duplicateFilter->DebugOn();
  TRY_EXPECT_NO_EXCEPTION( duplicateFilter->Update() );

  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImageFile );
  writer->SetInput( luminanceFilter->GetOutput() );
  writer->SetUseInputMetaDataDictionary( false );
  typedef itk::MetaImageIO ImageIOType;
  ImageIOType::Pointer imageIO = ImageIOType::New();
  imageIO->SetMetaDataDictionary(
    duplicateFilter->GetOutputMetaDataDictionary() );
  writer->SetImageIO( imageIO );
  writer->SetUseCompression( true );
  TRY_EXPECT_NO_EXCEPTION( writer->Update() );

  TEST_EXPECT_EQUAL( inputImage
    ->GetMetaDataDictionary().GetKeys().size(), 9 );
  TEST_EXPECT_EQUAL( duplicateFilter->GetOutputMetaDataDictionary()
    .GetKeys().size(), 9 );

  return EXIT_SUCCESS;
}
