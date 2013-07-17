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

#include "itktubePDFSegmenter.h"

#include <itkFilterWatcher.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>

int itktubePDFSegmenterTest( int argc, char * argv[] )
{
  if( argc != 12 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
      << " inputImage1 inputImage2 inputLabelMap force blur outputProbImage0"
      << " outputPDF0 outputProbImage1 outputPDF1 outputLabelMap"
      << " labeledFeatureSpace"
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
  typedef itk::tube::PDFSegmenter< ImageType, 2, ImageType > FilterType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during input read:" << std::endl << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer inputImage = reader->GetOutput();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  try
    {
    reader2->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during input image2 read:" << std::endl
      << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer inputImage2 = reader2->GetOutput();

  // Create the reader and writer
  ReaderType::Pointer labelmapReader = ReaderType::New();
  labelmapReader->SetFileName( argv[5] );
  try
    {
    labelmapReader->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during input read:\n"  << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer labelmapImage = labelmapReader->GetOutput();

  FilterType::Pointer filter = FilterType::New();
  filter->SetInputVolume( 0, inputImage );
  filter->SetInputVolume( 1, inputImage2 );
  filter->SetLabelMap( labelmapImage );
  filter->SetObjectId( 255 );
  filter->AddObjectId( 127 );
  filter->SetVoidId( 0 );
  filter->SetErodeRadius( 0 );
  filter->SetHoleFillIterations( 5 );
  float blur = atof( argv[4] );
  filter->SetProbabilityImageSmoothingStandardDeviation( blur );
  filter->SetHistogramSmoothingStandardDeviation( 2 );
  filter->SetOutlierRejectPortion( 0.1 );
  filter->SetObjectPDFWeight( 0, 1.5 );
  filter->SetDraft( false );
  if( argv[3][0] == 't' || argv[3][0] == 'T' || argv[3][0] == '1' )
    {
    filter->SetReclassifyObjectLabels( true );
    filter->SetReclassifyNotObjectLabels( true );
    filter->SetForceClassification( true );
    }
  else
    {
    filter->SetReclassifyObjectLabels( false );
    filter->SetReclassifyNotObjectLabels( false );
    filter->SetForceClassification( false );
    }
  filter->Update();
  filter->ClassifyImages();

  WriterType::Pointer probWriter0 = WriterType::New();
  probWriter0->SetFileName( argv[6] );
  probWriter0->SetUseCompression( true );
  probWriter0->SetInput( filter->GetClassProbabilityForInputVolume( 0 ) );
  try
    {
    probWriter0->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during write:\n"  << e;
    return EXIT_FAILURE;
    }

  WriterType::Pointer pdfWriter0 = WriterType::New();
  pdfWriter0->SetFileName( argv[7] );
  pdfWriter0->SetUseCompression( true );
  pdfWriter0->SetInput( filter->GetClassPDFImage( 0 ) );
  try
    {
    pdfWriter0->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during write:\n"  << e;
    return EXIT_FAILURE;
    }

  WriterType::Pointer probWriter1 = WriterType::New();
  probWriter1->SetFileName( argv[8] );
  probWriter1->SetUseCompression( true );
  probWriter1->SetInput( filter->GetClassProbabilityForInputVolume( 1 ) );
  try
    {
    probWriter1->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during write:\n"  << e;
    return EXIT_FAILURE;
    }

  WriterType::Pointer pdfWriter1 = WriterType::New();
  pdfWriter1->SetFileName( argv[9] );
  pdfWriter1->SetUseCompression( true );
  pdfWriter1->SetInput( filter->GetClassPDFImage( 1 ) );
  try
    {
    pdfWriter1->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during write:\n"  << e;
    return EXIT_FAILURE;
    }

  WriterType::Pointer labelmapWriter = WriterType::New();
  labelmapWriter->SetFileName( argv[10] );
  labelmapWriter->SetUseCompression( true );
  labelmapWriter->SetInput( filter->GetLabelMap() );
  try
    {
    labelmapWriter->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during write:\n"  << e;
    return EXIT_FAILURE;
    }

  filter->GenerateLabeledFeatureSpace();

  WriterType::Pointer labeledFeatureSpaceWriter = WriterType::New();
  labeledFeatureSpaceWriter->SetFileName( argv[11] );
  labeledFeatureSpaceWriter->SetUseCompression( true );
  labeledFeatureSpaceWriter->SetInput( filter->GetLabeledFeatureSpace() );
  try
    {
    labeledFeatureSpaceWriter->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during write:\n"  << e;
    return EXIT_FAILURE;
    }

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
