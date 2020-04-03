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

#include "itktubeAnisotropicDiffusiveRegistrationFilter.h"

#include <itkImageFileWriter.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

#include <vtkPlaneSource.h>

int itktubeAnisotropicDiffusiveRegistrationRegularizationTest( int argc, char * argv[] )
{
  if( argc < 7 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << "smoothed motion field image, "
              << "noise variance, "
              << "border slope, "
              << "number of iterations, "
              << "time step, "
              << "should use anisotropic regularization"
              << std::endl;
    return EXIT_FAILURE;
    }

  // Typedefs
  enum { Dimension = 3 };
  typedef double                                          PixelType;
  typedef double                                          VectorScalarType;
  typedef itk::Image< PixelType, Dimension >              FixedImageType;
  typedef itk::Image< PixelType, Dimension >              MovingImageType;
  typedef itk::Vector< VectorScalarType, Dimension >      VectorType;
  typedef itk::Image< VectorType, Dimension >             DeformationFieldType;
  typedef itk::ImageRegionIterator< DeformationFieldType >
                                                          IteratorType;
  typedef itk::Index< Dimension >                         IndexType;

  // Image parameters
  int         startValue = 0;
  int         sizeValue = 10;
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
  origin.Fill( originValue );

  deformationField->SetRegions( region );
  deformationField->SetSpacing( spacing );
  deformationField->SetOrigin( origin );
  deformationField->Allocate();

  // Fill the motion field image:
  // Top half has vectors like \, bottom half has vectors like /,
  // with additional noise

  PixelType   borderSlope;
  VectorType  borderN( 0.0 ); // normal to the border
  VectorType  perpN;   // perpendicular to the border

  borderSlope = std::atof( argv[3] );
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
  vtkSmartPointer< vtkPlaneSource > plane = vtkSmartPointer< vtkPlaneSource >::New();
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
  double        mean = 0;
  double        variance = std::atof( argv[2] );

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();

    for( unsigned int i = 0; i < Dimension; i++ )
      {
      indexAsVector[i] = index[i];
      }

    // Use definition of a plane to decide which side we are on
    if( borderN * ( center - indexAsVector ) < 0 )
      {
      pixel = bottomHalfPixel;
      }
    else
      {
      pixel = topHalfPixel;
      }

    // Add random noise
    PixelType randX = randGenerator->GetNormalVariate( mean, variance );
    PixelType randY = randGenerator->GetNormalVariate( mean, variance );
    PixelType randZ = randGenerator->GetNormalVariate( mean, variance );
    pixel[0] += randX;
    pixel[1] += randY;
    pixel[2] += randZ;

    it.Set( pixel );
    }

  // Setup the images to be registered
  FixedImageType::Pointer fixedImage      = FixedImageType::New();
  MovingImageType::Pointer movingImage    = MovingImageType::New();

  fixedImage->SetRegions( region );
  fixedImage->SetSpacing( spacing );
  fixedImage->SetOrigin( origin );
  fixedImage->Allocate();
  fixedImage->FillBuffer( 0.0 );

  movingImage->SetRegions( region );
  movingImage->SetSpacing( spacing );
  movingImage->SetOrigin( origin );
  movingImage->Allocate();
  movingImage->FillBuffer( 0.0 );

  // Setup the registrator object
  typedef itk::tube::DiffusiveRegistrationFilter
      < FixedImageType, MovingImageType, DeformationFieldType >
      DiffusiveRegistrationFilterType;
  typedef itk::tube::AnisotropicDiffusiveRegistrationFilter
      < FixedImageType, MovingImageType, DeformationFieldType >
      AnisotropicDiffusiveRegistrationFilterType;

  DiffusiveRegistrationFilterType::Pointer registrator = nullptr;
  AnisotropicDiffusiveRegistrationFilterType::Pointer anisotropicRegistrator =
    nullptr;
  int useAnisotropic = std::atoi( argv[6] );
  if( useAnisotropic )
    {
    registrator = AnisotropicDiffusiveRegistrationFilterType::New();
    anisotropicRegistrator
        = dynamic_cast < AnisotropicDiffusiveRegistrationFilterType * >(
            registrator.GetPointer() );
    }
  else
    {
    registrator = DiffusiveRegistrationFilterType::New();
    }
  registrator->SetInitialDisplacementField( deformationField );
  registrator->SetMovingImage( movingImage );
  registrator->SetFixedImage( fixedImage );
  // because we are just doing motion field regularization in this test:
  registrator->SetComputeIntensityDistanceTerm( false );
  registrator->SetTimeStep( std::atof( argv[5] ) );
  registrator->SetNumberOfIterations( std::atoi( argv[4] ) );
  if( anisotropicRegistrator )
    {
    anisotropicRegistrator->SetBorderSurface( plane->GetOutput() );
    }

  // Save the smoothed deformation field
  typedef itk::ImageFileWriter< DeformationFieldType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[1] );
  writer->SetInput( registrator->GetOutput() );
  writer->SetUseCompression( true );
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
  if( anisotropicRegistrator )
    {
    vtkPolyData * normalPolyData = anisotropicRegistrator->GetBorderSurface();
    vtkSmartPointer< vtkDataArray > normalData
                                  = normalPolyData->GetPointData()->GetNormals();

    // test to make sure the extracted normals match the known normals
    double ep = 0.00005;
    double test[3];
    for( int i = 0; i < normalData->GetNumberOfTuples(); i++ )
      {
      normalData->GetTuple( i, test );
      if( std::fabs( test[0] - borderN[0] ) > ep
          || std::fabs( test[1] - borderN[1] ) > ep
          || std::fabs( test[2] - borderN[2] ) > ep )
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
    }

  return EXIT_SUCCESS;
}
