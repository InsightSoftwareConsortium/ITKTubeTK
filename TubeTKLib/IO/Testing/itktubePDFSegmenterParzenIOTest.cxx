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

#include "itktubePDFSegmenterParzenIO.h"

int itktubePDFSegmenterParzenIOTest( int argc, char * argv[] )
{
  if( argc != 8 )
    {
    std::cout << "Missing arguments." << std::endl;
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0]
      << " inputImage1 inputImage2 inputLabelMap outputLabelMap"
      << " pdfFile outputLabelMap2 pdfFile2"
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
  typedef itk::tube::PDFSegmenterParzen< ImageType, ImageType >
    FilterType;

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

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  try
    {
    reader2->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception caught during input image2 read:" << std::endl
      << e << std::endl;
    return EXIT_FAILURE;
    }
  ImageType::Pointer inputImage2 = reader2->GetOutput();

  // Create the reader and writer
  ReaderType::Pointer labelmapReader = ReaderType::New();
  labelmapReader->SetFileName( argv[3] );
  try
    {
    labelmapReader->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception caught during input read:" << std::endl << e
      << std::endl;
    return EXIT_FAILURE;
    }
  ImageType::Pointer labelmapImage = labelmapReader->GetOutput();

  FilterType::FeatureVectorGeneratorType::Pointer fvGen =
    FilterType::FeatureVectorGeneratorType::New();
  fvGen->SetInput( inputImage );
  fvGen->AddInput( inputImage2 );

  FilterType::Pointer filter = FilterType::New();
  filter->SetFeatureVectorGenerator( fvGen );
  filter->SetLabelMap( labelmapImage );
  filter->SetObjectId( 255 );
  filter->AddObjectId( 127 );
  filter->SetVoidId( 0 );
  filter->SetErodeDilateRadius( 0 );
  filter->SetHoleFillIterations( 5 );
  filter->SetProbabilityImageSmoothingStandardDeviation( 1 );
  filter->SetHistogramSmoothingStandardDeviation( 2 );
  filter->SetOutlierRejectPortion( 0.1 );
  filter->SetObjectPDFWeight( 0, 1.5 );
  filter->SetReclassifyObjectLabels( true );
  filter->SetReclassifyNotObjectLabels( true );
  filter->SetForceClassification( true );
  filter->Update();
  std::cout << "*** Filter 1 ***" << std::endl << filter << std::endl;
  filter->ClassifyImages();

  WriterType::Pointer labelmapWriter = WriterType::New();
  labelmapWriter->SetFileName( argv[4] );
  labelmapWriter->SetUseCompression( true );
  labelmapWriter->SetInput( filter->GetLabelMap() );
  try
    {
    labelmapWriter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception caught during write:" << std::endl << e
      << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cout << "Exception caught during label write." << std::endl;
    return EXIT_FAILURE;
    }

  itk::tube::PDFSegmenterParzenIO< ImageType, ImageType > PDFIO( filter );
  std::cout << "*** Writing Filter 1 ***" << std::endl;
  std::cout << "filename = " << argv[5] << std::endl;
  std::cout << "*** PDFIO ***" << std::endl;
  PDFIO.PrintInfo();
  try
    {
    PDFIO.Write( argv[5] );
    }
  catch( ... )
    {
    std::cout << "Exception caught during pdf write." << std::endl;
    return EXIT_FAILURE;
    }

  FilterType::Pointer filter2 = FilterType::New();
  filter2->SetFeatureVectorGenerator( fvGen );

  itk::tube::PDFSegmenterParzenIO< ImageType, ImageType > PDFIO2( filter2 );
  try
    {
    if( !PDFIO2.Read( argv[5] ) )
      {
      std::cout << "Error in reading PDF file." << std::endl;
      return EXIT_FAILURE;
      }
    }
  catch( ... )
    {
    std::cout << "Exception caught during pdf write." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "*** PDFIO2 ***" << std::endl;
  PDFIO2.PrintInfo();

  std::cout << "*** Filter 2 ***" << std::endl << filter2 << std::endl;
  filter2->ClassifyImages();

  WriterType::Pointer labelmapWriter2 = WriterType::New();
  labelmapWriter2->SetFileName( argv[6] );
  labelmapWriter2->SetUseCompression( true );
  labelmapWriter2->SetInput( filter2->GetLabelMap() );
  try
    {
    labelmapWriter2->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception caught during label2 write:" << std::endl << e
      << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cout << "Exception caught during label2 write." << std::endl;
    return EXIT_FAILURE;
    }

  itk::tube::PDFSegmenterParzenIO< ImageType, ImageType > PDFIO3( filter2 );
  PDFIO.Write( argv[7] );

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
