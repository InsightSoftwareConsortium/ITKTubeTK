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

#include "itktubeRidgeSeedFilterIO.h"

int itktubeRidgeSeedFilterIOTest( int argc, char * argv[] )
{
  itk::TimeProbesCollectorBase timeCollector;

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

  timeCollector.Start( "Reader" );
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
  timeCollector.Stop( "Reader" );

  FilterType::RidgeScalesType scales( 3 );
  scales[0] = 0.35;
  scales[1] = 1.05;
  scales[2] = 1.4;

  timeCollector.Start( "Filter Setup" );
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inputImage );
  filter->SetLabelMap( labelmapImage );
  filter->SetScales( scales );
  filter->SetRidgeId( 255 );
  filter->SetBackgroundId( 127 );
  filter->SetUnknownId( 0 );
  filter->SetTrainClassifier( true );
  timeCollector.Stop( "Filter Setup" );
  timeCollector.Start( "Filter Update" );
  filter->Update();
  std::cout << "Update done." << std::endl;
  timeCollector.Stop( "Filter Update" );

  timeCollector.Start( "Filter ClassifyImages" );
  filter->ClassifyImages();
  std::cout << "Classification done." << std::endl;
  timeCollector.Stop( "Filter ClassifyImages" );

  timeCollector.Start( "Labelmap Writer" );
  LabelMapWriterType::Pointer labelmapWriter = LabelMapWriterType::New();
  labelmapWriter->SetFileName( argv[3] );
  labelmapWriter->SetUseCompression( true );
  labelmapWriter->SetInput( filter->GetOutput() );
  try
    {
    labelmapWriter->Update();
    }
  catch ( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during write:" << std::endl << e;
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Labelmap Writer" );
  std::cout << "Labelmap Written" << std::endl;

  timeCollector.Start( "Filter copy and write" );
  itk::tube::RidgeSeedFilterIO< ImageType, LabelMapType > filterIO(
    filter );
  filterIO.Write( argv[4] );
  timeCollector.Stop( "Filter copy and write" );
  std::cout << "Filter copy and write done." << std::endl;

  timeCollector.Start( "Filter read and process 2" );
  FilterType::Pointer filter2 = FilterType::New();
  filter2->SetInput( inputImage );
  std::cout << "Filter2 setinput." << std::endl;

  itk::tube::RidgeSeedFilterIO< ImageType, LabelMapType > filterIO2(
    filter2 );
  std::cout << "FilterIO2 Init." << std::endl;
  filterIO2.Read( argv[4] );
  std::cout << "FilterIO2 Read." << std::endl;
  filter2->ClassifyImages();
  timeCollector.Stop( "Filter read and process 2" );
  std::cout << "Filter2 ClassifyImages." << std::endl;

  timeCollector.Start( "Labelmap Writer2" );
  LabelMapWriterType::Pointer labelmapWriter2 = LabelMapWriterType::New();
  labelmapWriter2->SetFileName( argv[5] );
  labelmapWriter2->SetUseCompression( true );
  labelmapWriter2->SetInput( filter2->GetOutput() );
  try
    {
    labelmapWriter2->Update();
    }
  catch ( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during write2:" << std::endl << e;
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Labelmap Writer2" );
  std::cout << "Labelmap2 Written" << std::endl;

  timeCollector.Start( "Filter Writer3" );
  char out3[255];
  snprintf( out3, 254, "%s.mrs", argv[5] );
  itk::tube::RidgeSeedFilterIO< ImageType, LabelMapType > filterIO3(
    filter2 );
  filterIO3.Write( out3 );
  timeCollector.Stop( "Filter Writer3" );
  std::cout << "Filter Writer3" << std::endl;

  timeCollector.Report();

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
