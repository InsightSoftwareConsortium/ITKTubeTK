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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkImageToImageDiffusiveDeformableRegistrationFilter.h"

#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageFileWriter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

int itkImageToImageDiffusiveDeformableRegistrationTest(int argc, char* argv [] )
{
  if( argc < 5 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << "original motion field image, "
              << "noise variance, "
              << "smoothed motion field image, "
              << "number of iterations"
              << std::endl;
    return EXIT_FAILURE;
    }

  // Typedefs
  const unsigned int                                      Dimension = 3;
  typedef double                                          PixelType;
  typedef itk::Image< PixelType, Dimension >              FixedImageType;
  typedef itk::Image< PixelType, Dimension >              MovingImageType;
  typedef itk::Vector< PixelType, Dimension >             VectorType;
  typedef itk::Image< VectorType, Dimension >             DeformationFieldType;
  typedef itk::ImageLinearIteratorWithIndex< DeformationFieldType >
                                                          IteratorType;

  // Image parameters
  int         startValue = 0;
  int         sizeValue = 50;
  double      spacingValue = 1.0;
  double      originValue = 0.0;

  // Allocate the motion field image
  DeformationFieldType::Pointer      deformationField
                                                  = DeformationFieldType::New();
  DeformationFieldType::IndexType    start;
  start.Fill( startValue );

  DeformationFieldType::SizeType     size;
  size.Fill( sizeValue );

  DeformationFieldType::RegionType   region;
  region.SetSize( size );
  region.SetIndex( start );

  DeformationFieldType::SpacingType  spacing;
  spacing.Fill( spacingValue );

  DeformationFieldType::PointType    origin;
  origin.Fill( originValue);

  deformationField->SetRegions( region );
  deformationField->SetSpacing( spacing );
  deformationField->SetOrigin( origin );
  deformationField->Allocate();

  // Fill the motion field image:
  // Top half is vectors like \, bottom half is vectors like /,
  // plus noise
  VectorType pixel;
  IteratorType it( deformationField, region );
  it.SetDirection(2);
  int numLines = 0;
  int halfway = sizeValue * sizeValue / 2;

  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer randGenerator
    = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  randGenerator->Initialize( 137593424 );
  PixelType randX = 0;
  PixelType randY = 0;
  PixelType randZ = 0;
  double mean = 0;
  double variance = atof(argv[2]);

  for( it.GoToBegin(); ! it.IsAtEnd(); it.NextLine() )
    {
    it.GoToBeginOfLine();
    while ( ! it.IsAtEndOfLine() )
      {
      if (numLines < halfway)
        {
        pixel[0] = 0.5;
        }
      else
        {
        pixel[0] = -0.5;
        }
      pixel[1] = 0.5;
      pixel[2] = 0.0;

      // add random noise
      randX = randGenerator->GetNormalVariate( mean, variance );
      randY = randGenerator->GetNormalVariate( mean, variance );
      randZ = randGenerator->GetNormalVariate( mean, variance );
      pixel[0] += randX;
      pixel[1] += randY;
      pixel[2] += randZ;

      it.Set(pixel);
      ++it;
      }
    numLines++;
    }

  // We know the normal since we've created artificial data
  VectorType normal;
  normal[0] = 0;
  normal[1] = 1;
  normal[2] = 0;

  // Save the motion field image
  typedef itk::ImageFileWriter< DeformationFieldType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[1] );
  writer->SetInput( deformationField );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // Setup the images to be registered
  FixedImageType::Pointer fixedImage      = FixedImageType::New();
  MovingImageType::Pointer movingImage    = MovingImageType::New();

  fixedImage->SetLargestPossibleRegion( region );
  fixedImage->SetSpacing( spacing );
  fixedImage->SetOrigin( origin );
  fixedImage->Allocate();

  movingImage->SetLargestPossibleRegion( region );
  movingImage->SetSpacing( spacing );
  movingImage->SetOrigin( origin );
  movingImage->Allocate();

  // Setup the registrator object
  typedef itk::ImageToImageDiffusiveDeformableRegistrationFilter
                                                      < FixedImageType,
                                                        MovingImageType,
                                                        DeformationFieldType >
                                                        RegistrationFilterType;
  RegistrationFilterType::Pointer registrator = RegistrationFilterType::New();

  registrator->SetInitialDeformationField( deformationField );
  registrator->SetMovingImage( movingImage );
  registrator->SetFixedImage( fixedImage );
  registrator->SetNormalVectors( normal );
  int numIterations = atoi( argv[4] );
  registrator->SetNumberOfIterations( numIterations );

  // Save the smoothed deformation field
  writer->SetFileName( argv[3] );
  writer->SetInput( registrator->GetOutput() );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }














  return EXIT_SUCCESS;

}

