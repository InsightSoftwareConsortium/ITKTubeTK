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

#include "vtkPlaneSource.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkPolyDataWriter.h"

int itkImageToImageDiffusiveDeformableRegistrationRegularizationTest(
                                                      int argc, char* argv [] )
{
  if( argc < 12 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << "original motion field image, "
              << "smoothed motion field image, "
              << "smoothed motion field image (tangential), "
              << "smoothed motion field image (normal), "
              << "normal vector image, "
              << "normal border surface, "
              << "noise variance, "
              << "border slope, "
              << "number of iterations, "
              << "time step, "
              << "should use diffusive regularization"
              << std::endl;
    return EXIT_FAILURE;
    }

  // Typedefs
  const unsigned int                                      Dimension = 3;
  typedef double                                          PixelType;
  typedef double                                          VectorScalarType;
  typedef itk::Image< PixelType, Dimension >              FixedImageType;
  typedef itk::Image< PixelType, Dimension >              MovingImageType;
  typedef itk::Vector< VectorScalarType, Dimension >      VectorType;
  typedef itk::Image< VectorType, Dimension >             VectorImageType;
  typedef itk::Image< VectorType, Dimension >             DeformationFieldType;
  typedef itk::ImageRegionIterator< DeformationFieldType >
                                                          IteratorType;
  typedef itk::Index< Dimension >                         IndexType;

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
  // Top half has vectors like \, bottom half has vectors like /,
  // with additional noise

  PixelType   borderSlope;
  VectorType  borderN; // normal to the border
  VectorType  perpN;   // perpendicular to the border

  borderSlope = atof( argv[8] );
  if( borderSlope == 0 )
    {
    borderN[0] = 0.0;
    borderN[1] = 1.0;
    borderN[0] = 0.0;
    perpN[0] = 1.0;
    perpN[1] = 0.0;
    perpN[2] = 0.0;
    }
  else
    {
    borderN[0] = -1.0;
    borderN[1] = 1.0 / borderSlope;
    borderN[2] = 0.0;
    perpN[0] = -1.0 / borderSlope;
    perpN[1] = -1.0;
    perpN[2] = 0.0;
    }

  // Normalize the normals
  borderN = borderN / borderN.GetNorm();
  perpN = perpN / perpN.GetNorm();

  VectorType topHalfPixel;
  VectorType bottomHalfPixel;
  topHalfPixel = borderN + perpN;
  bottomHalfPixel = borderN - perpN;

  // The index at the center of the image is on the plane
  VectorType center;
  for( unsigned int i = 0; i < Dimension; i++ )
    {
    center[i] = deformationField->GetLargestPossibleRegion().GetSize()[i] / 2.0;
    }

  // Create the polydata for the surface
  vtkSmartPointer< vtkPlaneSource > plane = vtkPlaneSource::New();
  plane->SetCenter( center[0], center[1], center[2] );
  plane->SetNormal( borderN[0], borderN[1], borderN[2] );
  plane->Update();

  VectorType    pixel;
  IteratorType  it( deformationField,
                    deformationField->GetLargestPossibleRegion() );
  IndexType     index;
  VectorType    indexAsVector;
  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer randGenerator
    = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  randGenerator->Initialize( 137593424 );
  PixelType     randX = 0;
  PixelType     randY = 0;
  PixelType     randZ = 0;
  double        mean = 0;
  double        variance = atof(argv[7]);

  for( it.GoToBegin(); ! it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();

    for( unsigned int i = 0; i < Dimension; i++ )
      {
      indexAsVector[i] = index[i];
      }

    // Use definition of a plane to decide which side we are on
    if ( borderN * ( center - indexAsVector ) < 0 )
      {
      pixel = bottomHalfPixel;
      }
    else
      {
      pixel = topHalfPixel;
      }

    // Add random noise
    randX = randGenerator->GetNormalVariate( mean, variance );
    randY = randGenerator->GetNormalVariate( mean, variance );
    randZ = randGenerator->GetNormalVariate( mean, variance );
    pixel[0] += randX;
    pixel[1] += randY;
    pixel[2] += randZ;

    it.Set(pixel);
    }

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

  // Save the border polydata
  vtkPolyDataWriter * polyDataWriter = vtkPolyDataWriter::New();
  polyDataWriter->SetFileName( argv[6] );
  polyDataWriter->SetInput( plane->GetOutput() );
  polyDataWriter->Update();

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
  registrator->SetBorderSurface( plane->GetOutput() );
  int numIterations = atoi( argv[9] );
  registrator->SetNumberOfIterations( numIterations );
  // because we are just doing motion field regularization in this test:
  registrator->SetComputeIntensityDistanceTerm( false );
  registrator->SetTimeStep( atof( argv[10] ) );
  int useDiffusive = atoi( argv[11] );
  if ( useDiffusive )
    {
    registrator->SetUseDiffusiveRegularization( true );
    }
  else
    {
    registrator->SetUseDiffusiveRegularization( false );
    }

  // Save the smoothed deformation field
  writer->SetFileName( argv[2] );
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

  // Check to make sure the border normals were calculated correctly by the
  // registrator
  vtkPolyData * normalPolyData = registrator->GetBorderNormalsSurface();
  vtkSmartPointer< vtkDataArray > normalData
                                = normalPolyData->GetPointData()->GetNormals();

  // test to make sure the extracted normals match the known normals
  double ep = 0.00005;
  double test[3];
  for( int i = 0; i < normalData->GetNumberOfTuples(); i++ )
    {
    normalData->GetTuple( i, test );
    if( fabs( test[0] - borderN[0] ) > ep
        || fabs( test[1] - borderN[1] ) > ep
        || fabs( test[2] - borderN[2] ) > ep )
      {
      std::cerr << "index i=" << i << ": extracted normal [" << test[0] << " "
          << test[1] << " " << test[2] << "]" << std::endl;
      std::cerr << "does not match known border normal [" << borderN[0] << " "
          << borderN[1] << " " << borderN[2] << "]" << std::endl;
      return EXIT_FAILURE;
      }
    }
  if( plane->GetOutput()->GetPointData()->GetNumberOfTuples() !=
      normalPolyData->GetPointData()->GetNumberOfTuples() )
    {
    std::cerr << "number of tuples in original plane does not match number of "
                                    << "tuples in border normal" << std::endl;
    return EXIT_FAILURE;
    }

  // Save the normal vector image
  typedef itk::ImageFileWriter< VectorImageType > VectorWriterType;
  VectorWriterType::Pointer vectorWriter = VectorWriterType::New();
  vectorWriter->SetFileName( argv[5] );
  vectorWriter->SetInput( registrator->GetNormalVectorImage() );
  vectorWriter->Write();

  // Save the output deformation field tangential and normal images
  writer->SetFileName( argv[3] );
  writer->SetInput( registrator->GetTangentialDeformationFieldImage() );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }
  writer->SetFileName( argv[4] );
  writer->SetInput( registrator->GetNormalDeformationFieldImage() );
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

