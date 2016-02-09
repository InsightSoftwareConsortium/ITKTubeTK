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

#include "itktubeBasisFeatureVectorGenerator.h"
#include "itktubeNJetFeatureVectorGenerator.h"

int itktubeNJetBasisFeatureVectorGeneratorTest( int argc, char * argv[] )
{
  if( argc != 7 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
      << " inputImage maskImage objId bkgId outputLDA0Image outputLDA1Image"
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


  // Declare the type for the filter
  typedef itk::tube::NJetFeatureVectorGenerator< ImageType >
    FilterType;
  typedef itk::tube::BasisFeatureVectorGenerator< ImageType, LabelMapType >
    BasisFilterType;

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
  LabelMapReaderType::Pointer maskReader = LabelMapReaderType::New();
  maskReader->SetFileName( argv[2] );
  try
    {
    maskReader->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception caught during input mask read:" << std::endl
      << e;
    return EXIT_FAILURE;
    }
  LabelMapType::Pointer maskImage = maskReader->GetOutput();

  FilterType::NJetScalesType scales( 2 );
  scales[0] = 4;
  scales[1] = 8;

  FilterType::NJetScalesType scales2( 1 );
  scales2[0] = 8;

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inputImage );
  filter->SetZeroScales( scales );
  filter->SetFirstScales( scales );
  filter->SetSecondScales( scales2 );
  filter->SetRidgeScales( scales2 );
  filter->UpdateWhitenStatistics();

  std::cout << filter << std::endl;

  BasisFilterType::Pointer basisFilter = BasisFilterType::New();
  basisFilter->SetInputFeatureVectorGenerator( filter.GetPointer() );
  basisFilter->SetInput( inputImage );
  basisFilter->SetLabelMap( maskImage );
  int objId = std::atoi( argv[3] );
  int bkgId = std::atoi( argv[4] );
  basisFilter->SetObjectId( objId );
  basisFilter->AddObjectId( bkgId );
  basisFilter->SetNumberOfLDABasisToUseAsFeatures( 1 );
  basisFilter->SetNumberOfPCABasisToUseAsFeatures( 3 );
  basisFilter->GenerateBasis();

  std::cout << basisFilter << std::endl;

  basisFilter->SetLabelMap( NULL );

  for( unsigned int i=0; i<basisFilter->GetNumberOfFeatures(); ++i )
    {
    WriterType::Pointer featureImage0Writer = WriterType::New();
    char fname[80];
    sprintf( fname, "%s.%d.mha", argv[5], i );
    featureImage0Writer->SetFileName( fname );
    featureImage0Writer->SetUseCompression( true );
    featureImage0Writer->SetInput( basisFilter->GetFeatureImage( i ) );
    try
      {
      featureImage0Writer->Update();
      }
    catch (itk::ExceptionObject& e)
      {
      std::cerr << "Exception caught during write:" << std::endl << e;
      return EXIT_FAILURE;
      }
    }

  for( unsigned int i=0; i<filter->GetNumberOfFeatures(); ++i )
    {
    WriterType::Pointer featureImage1Writer = WriterType::New();
    char fname[80];
    sprintf( fname, "%s.%d.mha", argv[6], i );
    featureImage1Writer->SetFileName( fname );
    featureImage1Writer->SetUseCompression( true );
    featureImage1Writer->SetInput( filter->GetFeatureImage( i ) );
    try
      {
      featureImage1Writer->Update();
      }
    catch (itk::ExceptionObject& e)
      {
      std::cerr << "Exception caught during write:" << std::endl << e;
      return EXIT_FAILURE;
      }
    }

  return EXIT_SUCCESS;
}
