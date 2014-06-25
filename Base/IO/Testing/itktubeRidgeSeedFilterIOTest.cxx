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

#include "itktubeRidgeSeedFilterIO.h"

int itktubeRidgeSeedFilterIOTest( int argc, char * argv[] )
{
  if( argc != 6 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
      << " inputImage labelmapImage outputImage fileIO outputImage2"
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

  typedef itk::Image< unsigned char, Dimension >    LabelMapType;
  typedef itk::ImageFileReader< LabelMapType >      LabelMapReaderType;
  typedef itk::ImageFileWriter< LabelMapType >      LabelMapWriterType;

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
  filter->SetRidgeId( 255 );
  filter->SetBackgroundId( 127 );
  filter->SetUnknownId( 0 );
  filter->SetTrainClassifier( true );
  filter->Update();
  std::cout << "Update done." << std::endl;

  filter->ClassifyImages();
  std::cout << "Classification done." << std::endl;

  LabelMapWriterType::Pointer labelmapWriter = LabelMapWriterType::New();
  labelmapWriter->SetFileName( argv[3] );
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

  itk::tube::RidgeSeedFilterIO< ImageType, LabelMapType > filterIO(
    filter );
  filterIO.Write( argv[4] );

  FilterType::Pointer filter2 = FilterType::New();
  filter2->SetInput( inputImage );

  itk::tube::RidgeSeedFilterIO< ImageType, LabelMapType > filterIO2(
    filter2 );
  filterIO2.Read( argv[4] );

  filter2->SetTrainClassifier( false );
  filter2->Update();
  filter2->ClassifyImages();

  LabelMapWriterType::Pointer labelmapWriter2 = LabelMapWriterType::New();
  labelmapWriter2->SetFileName( argv[5] );
  labelmapWriter2->SetUseCompression( true );
  labelmapWriter2->SetInput( filter2->GetOutput() );
  try
    {
    labelmapWriter2->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during write2:" << std::endl << e;
    return EXIT_FAILURE;
    }

  char out3[255];
  sprintf( out3, "%s.mrs", argv[5] );
  itk::tube::RidgeSeedFilterIO< ImageType, LabelMapType > filterIO3(
    filter2 );
  filterIO3.Write( out3 );

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
