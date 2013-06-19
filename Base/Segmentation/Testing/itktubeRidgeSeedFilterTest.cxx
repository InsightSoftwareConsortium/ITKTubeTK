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

#include <itkFilterWatcher.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int itktubeRidgeSeedFilterTest( int argc, char * argv[] )
{
  if( argc != 7 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
      << " inputImage labelmapImage objId bkgId outputClass0PDF outputImage"
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

  typedef itk::Image< unsigned char, Dimension >    LabelmapType;
  typedef itk::ImageFileReader< LabelmapType >      LabelmapReaderType;
  typedef itk::ImageFileWriter< LabelmapType >      LabelmapWriterType;

  typedef itk::Image< float, 3 >                    PDFImageType;
  typedef itk::ImageFileWriter< PDFImageType >      PDFImageWriterType;

  // Declare the type for the Filter
  typedef itk::tube::RidgeSeedFilter< ImageType, LabelmapType >
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
  LabelmapReaderType::Pointer mReader = LabelmapReaderType::New();
  mReader->SetFileName( argv[2] );
  try
    {
    mReader->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during input mask read:" << std::endl << e;
    return EXIT_FAILURE;
    }
  LabelmapType::Pointer labelmapImage = mReader->GetOutput();

  FilterType::RidgeScalesType scales(3);
  scales[0] = 0.4;
  scales[1] = 0.8;
  scales[2] = 1.6;

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inputImage );
  filter->SetLabelmap( labelmapImage );
  filter->SetScales( scales );
  int objId = atoi( argv[3] );
  int bkgId = atoi( argv[4] );
  filter->SetObjectId( objId );
  filter->AddObjectId( bkgId );
  std::cout << filter << std::endl;
  filter->Update();
  std::cout << "Update done." << std::endl;

  filter->ClassifyImages();
  std::cout << "Classification done." << std::endl;

  PDFImageWriterType::Pointer pdfImageWriter = PDFImageWriterType::New();
  pdfImageWriter->SetFileName( argv[5] );
  pdfImageWriter->SetUseCompression( true );
  pdfImageWriter->SetInput( filter->GetPDFSegmenter()->GetClassPDFImage( 0 ) );
  try
    {
    pdfImageWriter->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during write:" << std::endl << e;
    return EXIT_FAILURE;
    }

  LabelmapWriterType::Pointer labelmapWriter = LabelmapWriterType::New();
  labelmapWriter->SetFileName( argv[6] );
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
