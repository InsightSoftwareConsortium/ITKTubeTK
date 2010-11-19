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

#include "itkImageToImageAnisotropicDiffusiveDeformableRegistrationFilter.h"
#include "itkOrientedImage.h"
#include "itkVectorCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

int itkImageToImageAnisotropicDiffusiveDeformableRegistrationExecution(
                                                      int argc, char* argv [] )
{
  if( argc < 11 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << "original fixed image, "
              << "original moving image, "
              << "surface border normal vector image, "
              << "surface border weight image, "
              << "resulting motion field image, "
              << "resulting transformed moving image, "
              << "number of iterations, "
              << "time step, "
              << "compute regularization term, "
              << "should use anisotropic regularization"
              << std::endl;
    return EXIT_FAILURE;
    }

  // Typedefs
  const unsigned int                                      ImageDimension = 3;
  typedef double                                          FixedPixelType;
  typedef double                                          MovingPixelType;
  typedef MovingPixelType                                 OutputPixelType;
  typedef double                                          VectorScalarType;
  typedef itk::Image< FixedPixelType, ImageDimension >    FixedImageType;
  typedef itk::Image< MovingPixelType, ImageDimension >   MovingImageType;
  typedef itk::Image< OutputPixelType, ImageDimension >   OutputImageType;
  typedef itk::Vector< VectorScalarType, ImageDimension > VectorType;
  typedef itk::Image< VectorType, ImageDimension >        VectorImageType;
  typedef itk::Image< double, ImageDimension >            WeightImageType;
  typedef itk::Image< VectorType, ImageDimension >        FieldType;


  //--------------------------------------------------------
  std::cout << "Read the input images" << std::endl;

  typedef itk::ImageFileReader< FixedImageType > FixedImageReaderType;
  FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
  fixedImageReader->SetFileName( argv[1] );
  fixedImageReader->Update();
  FixedImageType::Pointer fixed = fixedImageReader->GetOutput();

  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  MovingImageReaderType::Pointer movingImageReader
      = MovingImageReaderType::New();
  movingImageReader->SetFileName( argv[2] );
  movingImageReader->Update();
  MovingImageType::Pointer moving = movingImageReader->GetOutput();

  FieldType::Pointer initField = FieldType::New();
  initField->SetSpacing( fixed->GetSpacing() );
  initField->SetOrigin( fixed->GetOrigin() );
  initField->SetLargestPossibleRegion( fixed->GetLargestPossibleRegion() );
  initField->SetRequestedRegion( fixed->GetRequestedRegion() );
  initField->SetBufferedRegion( fixed->GetBufferedRegion() );
  initField->Allocate();

  // fill initial deformation with zero vectors
  VectorType zeroVec;
  zeroVec.Fill( 0.0 );
  initField->FillBuffer( zeroVec );

  typedef itk::VectorCastImageFilter<FieldType,FieldType> CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( initField );
  caster->InPlaceOff();

  //--------------------------------------------------------
  std::cout << "Read the input normal vector image and weight image" << std::endl;
  typedef itk::ImageFileReader< VectorImageType > NormalVectorImageReaderType;
  NormalVectorImageReaderType::Pointer normalVectorImageReader =
      NormalVectorImageReaderType::New();
  normalVectorImageReader->SetFileName( argv[3] );
  normalVectorImageReader->Update();
  VectorImageType::Pointer normalVectorImage
      = normalVectorImageReader->GetOutput();

  typedef itk::ImageFileReader< WeightImageType > WeightImageReaderType;
  WeightImageReaderType::Pointer weightImageReader
      = WeightImageReaderType::New();
  weightImageReader->SetFileName( argv[4] );
  weightImageReader->Update();
  WeightImageType::Pointer weightImage = weightImageReader->GetOutput();

  //-------------------------------------------------------------
  std::cout << "Run registration and warp moving" << std::endl;

  typedef itk::ImageToImageAnisotropicDiffusiveDeformableRegistrationFilter
      < FixedImageType, MovingImageType, FieldType > RegistrationType;
  RegistrationType::Pointer registrator = RegistrationType::New();

  registrator->SetInitialDeformationField( caster->GetOutput() );
  registrator->SetMovingImage( moving );
  registrator->SetFixedImage( fixed );
  registrator->SetNormalVectorImage( normalVectorImage );
  registrator->SetWeightImage( weightImage );
  //registrator->SetBorderSurface( border );
  int numberOfIterations = atoi( argv[7] );
  registrator->SetNumberOfIterations( numberOfIterations );

  int compute = atoi( argv[9] );
  if (compute)
    {
    registrator->SetComputeRegularizationTerm( true );
    }
  else
    {
    registrator->SetComputeRegularizationTerm( false );
    }

  int useAnisotropic = atoi( argv[10] );
  if ( useAnisotropic )
    {
    registrator->SetUseAnisotropicRegularization( true );
    }
  else
    {
    registrator->SetUseAnisotropicRegularization( false );
    }

  registrator->SetTimeStep( atof( argv[8] ) );
  registrator->SetLambda( -0.1 );

  // warp moving image
  typedef itk::WarpImageFilter< MovingImageType, MovingImageType, FieldType >
      WarperType;
  WarperType::Pointer warper = WarperType::New();

  typedef WarperType::CoordRepType CoordRepType;
  typedef itk::LinearInterpolateImageFunction<MovingImageType,CoordRepType>
      InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  warper->SetInput( moving );
  warper->SetDeformationField( registrator->GetOutput() );
  warper->SetInterpolator( interpolator );
  warper->SetOutputSpacing( fixed->GetSpacing() );
  warper->SetOutputOrigin( fixed->GetOrigin() );
  warper->SetOutputDirection( fixed->GetDirection() );
  warper->SetEdgePaddingValue( 0 );

  // Update triggers the registration
  warper->Update();

  // ---------------------------------------------------------
  std::cout << "Printing the deformation field and transformed moving image"
            << std::endl;

  typedef itk::ImageFileWriter< FieldType > FieldWriterType;
  FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
  fieldWriter->SetFileName( argv[5] );
  fieldWriter->SetInput( registrator->GetOutput() );
  try
    {
    fieldWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileWriter< MovingImageType > ImageWriterType;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetFileName( argv[6] );
  imageWriter->SetInput( warper->GetOutput() );
  try
    {
    imageWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // ---------------------------------------------------------
  std::cout << "Compare warped moving and fixed." << std::endl;

  // compare the warp and fixed images
  itk::ImageRegionIterator<FixedImageType> fixedIter( fixed,
      fixed->GetBufferedRegion() );
  itk::ImageRegionIterator<MovingImageType> warpedIter( warper->GetOutput(),
      fixed->GetBufferedRegion() );

  unsigned int numPixelsDifferent = 0;
  while( !fixedIter.IsAtEnd() )
    {
    if( fixedIter.Get() != warpedIter.Get() )
      {
      numPixelsDifferent++;
      }
    ++fixedIter;
    ++warpedIter;
    }

  std::cout << "Number of pixels different: " << numPixelsDifferent;
  std::cout << std::endl;

//  if( numPixelsDifferent > 10 )
//    {
//    std::cout << "Test failed - too many pixels different." << std::endl;
//    return EXIT_FAILURE;
//    }

  return EXIT_SUCCESS;

}
