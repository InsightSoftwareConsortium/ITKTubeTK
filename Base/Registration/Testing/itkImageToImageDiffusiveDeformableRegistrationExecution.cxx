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

#include "itkImageFileWriter.h"
#include "itkVectorCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

int itkImageToImageDiffusiveDeformableRegistrationExecution(
                                                      int argc, char* argv [] )
{
  if( argc < 13 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << "original fixed image, "
              << "original moving image, "
              << "surface border polydata, "
              << "resulting motion field image, "
              << "resulting transformed moving image, "
              << "weight image, "
              << "normal vector image, "
              << "surface normals border polydata, "
              << "number of iterations, "
              << "time step, "
              << "compute regularization term, "
              << "should use diffusive regularization"
              << std::endl;
    return EXIT_FAILURE;
    }

  // Typedefs
  const unsigned int                                      ImageDimension = 3;
  typedef double                                          PixelType;
  typedef double                                          VectorScalarType;
  typedef itk::Image< PixelType, ImageDimension >         ImageType;
  typedef itk::Vector< VectorScalarType, ImageDimension > VectorType;
  typedef itk::Image< VectorType, ImageDimension >        VectorImageType;
  typedef itk::Image< double, ImageDimension >            WeightImageType;
  typedef itk::Image< VectorType, ImageDimension >        FieldType;
  typedef itk::ImageFileReader<ImageType>                 ImageReaderType;
  typedef itk::ImageFileWriter<ImageType>                 ImageWriterType;
  typedef ImageType::IndexType                            IndexType;


  //--------------------------------------------------------
  std::cout << "Read the input images" << std::endl;

  ImageReaderType::Pointer fixedImageReader = ImageReaderType::New();
  fixedImageReader->SetFileName( argv[1] );
  fixedImageReader->Update();
  ImageType::Pointer fixed = fixedImageReader->GetOutput();

  ImageReaderType::Pointer movingImageReader = ImageReaderType::New();
  movingImageReader->SetFileName( argv[2] );
  movingImageReader->Update();
  ImageType::Pointer moving = movingImageReader->GetOutput();

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
  std::cout << "Read the input polydata" << std::endl;
  vtkPolyDataReader * polyReader = vtkPolyDataReader::New();
  polyReader->SetFileName( argv[3] );
  polyReader->Update();
  vtkPolyData * border = polyReader->GetOutput();

  //-------------------------------------------------------------
  std::cout << "Run registration and warp moving" << std::endl;

  typedef itk::ImageToImageDiffusiveDeformableRegistrationFilter< ImageType,
                                                                  ImageType,
                                                                  FieldType >
                                                                  RegistrationType;
  RegistrationType::Pointer registrator = RegistrationType::New();

  registrator->SetInitialDeformationField( caster->GetOutput() );
  registrator->SetMovingImage( moving );
  registrator->SetFixedImage( fixed );
  registrator->SetBorderSurface( border );
  int numberOfIterations = atoi( argv[9] );
  registrator->SetNumberOfIterations( numberOfIterations );

  int compute = atoi( argv[11] );
  if (compute)
    {
    registrator->SetComputeRegularizationTerm( true );
    }
  else
    {
    registrator->SetComputeRegularizationTerm( false );
    }

  int useDiffusive = atoi( argv[12] );
  if ( useDiffusive )
    {
    registrator->SetUseDiffusiveRegularization( true );
    }
  else
    {
    registrator->SetUseDiffusiveRegularization( false );
    }

  registrator->SetTimeStep( atof( argv[10] ) );
  registrator->SetLambda( -0.1 );

  // warp moving image
  typedef itk::WarpImageFilter<ImageType,ImageType,FieldType> WarperType;
  WarperType::Pointer warper = WarperType::New();

  typedef WarperType::CoordRepType CoordRepType;
  typedef itk::LinearInterpolateImageFunction<ImageType,CoordRepType>
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
  std::cout << "Printing the normal surface border, normal vector image "
      << "and weight image" << std::endl;

  vtkPolyData * normalPolyData = registrator->GetBorderNormalsSurface();
  vtkPolyDataWriter * polyWriter = vtkPolyDataWriter::New();
  polyWriter->SetFileName( argv[8] );
  polyWriter->SetInput( normalPolyData );
  polyWriter->Write();

  typedef itk::ImageFileWriter< VectorImageType > VectorWriterType;
  VectorWriterType::Pointer vectorWriter = VectorWriterType::New();
  vectorWriter->SetFileName( argv[7] );
  vectorWriter->SetInput( registrator->GetNormalVectorImage() );
  vectorWriter->Write();

  typedef itk::ImageFileWriter< WeightImageType > WeightWriterType;
  WeightWriterType::Pointer weightWriter = WeightWriterType::New();
  weightWriter->SetFileName( argv[6] );
  weightWriter->SetInput( registrator->GetWeightImage() );
  weightWriter->Write();

  // ---------------------------------------------------------
  std::cout << "Printing the deformation field and transformed moving image"
            << std::endl;

  typedef itk::ImageFileWriter< FieldType > FieldWriterType;
  FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
  fieldWriter->SetFileName( argv[4] );
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

  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetFileName( argv[5] );
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
  itk::ImageRegionIterator<ImageType> fixedIter( fixed,
      fixed->GetBufferedRegion() );
  itk::ImageRegionIterator<ImageType> warpedIter( warper->GetOutput(),
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

