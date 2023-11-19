/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include "itktubeRidgeSeedFilter.h"

int itktubeRidgeSeedFilterTest( int argc, char * argv[] )
{
  if( argc != 9 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
      << " inputImage labelmapImage objId bkgId pdfMethod"
      << " outputFeature0Image outputImage maxScaleImage"
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

  FilterType::RidgeScalesType scales( 3 );
  scales[0] = 0.35;
  scales[1] = 0.7;
  scales[2] = 1.05;

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inputImage );
  filter->SetLabelMap( labelmapImage );
  filter->SetScales( scales );
  int objId = atoi( argv[3] );
  int bkgId = atoi( argv[4] );
  if( argv[5][0] == '1' )
    {
    std::cerr << "LIBSVM not enabled." << std::endl;
    return EXIT_FAILURE;
    }
  else if( argv[5][0] == '2' )
    {
    std::cerr << "LIBRandomForest not enabled." << std::endl;
    return EXIT_FAILURE;
    }
  filter->SetRidgeId( objId );
  filter->SetBackgroundId( bkgId );
  filter->SetUnknownId( 0 );
  filter->SetTrainClassifier( true );
  std::cout << "Update started." << std::endl;
  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Error in RidgeSeedFilter update." << std::endl;
    std::cout << e << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cout << "Error in RidgeSeedFilter update." << std::endl;
    return EXIT_FAILURE;
    }
  try
    {
    filter->ClassifyImages();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Error in RidgeSeedFilter ClassifyImages." << std::endl;
    std::cout << e << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cout << "Error in RidgeSeedFilter ClassifyImages." << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "Update & Classification done." << std::endl;

  FeatureImageWriterType::Pointer feature2ImageWriter =
    FeatureImageWriterType::New();
  feature2ImageWriter->SetFileName( argv[6] );
  feature2ImageWriter->SetUseCompression( true );
  feature2ImageWriter->SetInput( filter->GetRidgeFeatureGenerator()
    ->GetFeatureImage( 0 ) );
  try
    {
    feature2ImageWriter->Update();
    }
  catch ( itk::ExceptionObject& e )
    {
    std::cout << "Exception caught during write:" << std::endl;
    std::cout << e << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cout << "Error in write." << std::endl;
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
  catch ( itk::ExceptionObject& e )
    {
    std::cout << "Exception caught during write:" << std::endl;
    std::cout << e << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cout << "Error in labelmap write." << std::endl;
    return EXIT_FAILURE;
    }

  FeatureImageWriterType::Pointer scaleImageWriter =
    FeatureImageWriterType::New();
  scaleImageWriter->SetFileName( argv[8] );
  scaleImageWriter->SetUseCompression( true );
  scaleImageWriter->SetInput( filter->GetOutputSeedScales() );
  try
    {
    scaleImageWriter->Update();
    }
  catch ( itk::ExceptionObject& e )
    {
    std::cout << "Exception caught during write:" << std::endl;
    std::cout << e << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cout << "Error in labelmap write." << std::endl;
    return EXIT_FAILURE;
    }

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
