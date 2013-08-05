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

#include "itktubeRidgeSeedFilter.h"

int itktubeRidgeSeedFilterTest( int argc, char * argv[] )
{
  if( argc != 8 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
      << " inputImage labelmapImage objId bkgId outputClass0PDF outputFeature0Image outputImage"
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

  typedef itk::Image< unsigned char, Dimension >    LabelMapType;
  typedef itk::ImageFileReader< LabelMapType >      LabelMapReaderType;
  typedef itk::ImageFileWriter< LabelMapType >      LabelMapWriterType;

  typedef itk::Image< float, 3 >                    PDFImageType;
  typedef itk::ImageFileWriter< PDFImageType >      PDFImageWriterType;

  typedef itk::Image< float, Dimension >            FeatureImageType;
  typedef itk::ImageFileWriter< FeatureImageType >  FeatureImageWriterType;

  // Declare the type for the Filter
  typedef itk::tube::RidgeSeedFilter< ImageType, LabelMapType >
    FilterType;

  // Create the reader
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

  // Create the mask reader
  LabelMapReaderType::Pointer mReader = LabelMapReaderType::New();
  mReader->SetFileName( argv[2] );
  try
    {
    mReader->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during input mask read:" << std::endl
      << e;
    return EXIT_FAILURE;
    }
  LabelMapType::Pointer labelmapImage = mReader->GetOutput();

  FilterType::RidgeScalesType scales(2);
  scales[0] = 0.15;
  scales[1] = 0.5;

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inputImage );
  filter->SetLabelMap( labelmapImage );
  filter->SetScales( scales );
  int objId = atoi( argv[3] );
  int bkgId = atoi( argv[4] );
  filter->SetObjectId( objId );
  filter->AddObjectId( bkgId );
  filter->AddObjectId( 0 );
  std::cout << filter << std::endl;
  filter->Update();
  std::cout << "Update done." << std::endl;

  filter->ClassifyImages();
  std::cout << "Classification done." << std::endl;

  PDFImageWriterType::Pointer pdfImageWriter = PDFImageWriterType::New();
  pdfImageWriter->SetFileName( argv[5] );
  pdfImageWriter->SetUseCompression( true );
  pdfImageWriter->SetInput( filter->GetPDFSegmenter()->
    GetClassPDFImage( 0 ) );
  try
    {
    pdfImageWriter->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during write:" << std::endl << e;
    return EXIT_FAILURE;
    }

  FeatureImageWriterType::Pointer feature2ImageWriter =
    FeatureImageWriterType::New();
  feature2ImageWriter->SetFileName( argv[6] );
  feature2ImageWriter->SetUseCompression( true );
  feature2ImageWriter->SetInput( filter->GetRidgeFeatureGenerator()->
    GetFeatureImage( 0 ) );
  try
    {
    feature2ImageWriter->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during write:" << std::endl << e;
    return EXIT_FAILURE;
    }

  LabelMapWriterType::Pointer labelmapWriter = LabelMapWriterType::New();
  labelmapWriter->SetFileName( argv[7] );
  labelmapWriter->SetUseCompression( true );
  labelmapWriter->SetInput( filter->GetOutput() );
  try
    {
    labelmapWriter->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during write:" << std::endl << e;
    return EXIT_FAILURE;
    }

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
