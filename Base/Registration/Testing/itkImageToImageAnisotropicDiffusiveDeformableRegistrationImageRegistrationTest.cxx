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

#include "itkImageFileWriter.h"
#include "itkVectorCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"

#include "vtkSphereSource.h"
#include "vtkPlaneSource.h"
#include "vtkCubeSource.h"
#include "vtkPolyDataWriter.h"
#include "vtkAppendPolyData.h"
#include "vtkDensifyPolyData.h"
#include "vtkSmoothPolyDataFilter.h"

// Template function to fill in an image with a sphere.
template <class TImage>
void
FillWithSphere(
TImage * image,
double * center,
double radius,
typename TImage::PixelType foregnd,
typename TImage::PixelType backgnd )
{
  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  Iterator it( image, image->GetBufferedRegion() );
  it.Begin();

  typename TImage::IndexType index;
  double r2 = vnl_math_sqr( radius );

  for( ; !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    double distance = 0;
    for( unsigned int j = 0; j < TImage::ImageDimension; j++ )
      {
      distance += vnl_math_sqr((double) index[j] - center[j]);
      }
    if( distance <= r2 ) it.Set( foregnd );
    else it.Set( backgnd );
    }
}

// Template function to fill in an image with two boxes
template <class TImage>
void
FillWithBox(
TImage * image,
double * bottomBox,
double * topBox,
double * size,
typename TImage::PixelType backgnd,
typename TImage::PixelType bottomStart,
typename TImage::PixelType bottomEnd,
typename TImage::PixelType topStart,
typename TImage::PixelType topEnd
)
{
  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  Iterator it( image, image->GetBufferedRegion() );
  it.Begin();

  typename TImage::IndexType index;
  typename TImage::PixelType bottomRange = bottomEnd - bottomStart;
  typename TImage::PixelType topRange = topEnd - topStart;
  double intensity;

  for( ; !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    if( index[0] <= bottomBox[0]
        && index[1] >= bottomBox[1] && index[1] < bottomBox[1] + size[1]
        && index[2] > bottomBox[2] && index[2] < bottomBox[2] + size[2] )
      {
      it.Set( bottomStart );
      }
    else if( index[0] > bottomBox[0] && index[0] < bottomBox[0] + size[0]
        && index[1] >= bottomBox[1] && index[1] < bottomBox[1] + size[1]
        && index[2] > bottomBox[2] && index[2] < bottomBox[2] + size[2] )
      {
      intensity = ( (index[0] - bottomBox[0] ) / size[0] ) * bottomRange + bottomStart;
      it.Set( intensity ); // or bottomStart for solid blocks
      }
    else if( index[0] >= bottomBox[0] + size[0]
        && index[1] >= bottomBox[1] && index[1] < bottomBox[1] + size[1]
        && index[2] > bottomBox[2] && index[2] < bottomBox[2] + size[2] )
      {
      it.Set( bottomEnd );
      }
    else if( index[0] <= topBox[0]
        && index[1] >= topBox[1] && index[1] < topBox[1] + size[1]
        && index[2] > topBox[2] && index[2] < topBox[2] + size[2] )
      {
      it.Set( topStart );
      }
    else if( index[0] > topBox[0] && index[0] < topBox[0] + size[0]
        && index[1] >= topBox[1] && index[1] < topBox[1] + size[1]
        && index[2] > topBox[2] && index[2] < topBox[2] + size[2] )
      {
      intensity = ( (index[0] - topBox[0] ) / size[0] ) * topRange + topStart;
      it.Set( intensity ); // or topEnd for solid blocks
      }
    else if( index[0] >= topBox[0] + size[0]
      && index[1] >= topBox[1] && index[1] < topBox[1] + size[1]
      && index[2] > topBox[2] && index[2] < topBox[2] + size[2] )
      {
      it.Set( topEnd );
      }
    else
      {
      it.Set( backgnd );
      }
    }
}

// Function to create the spherical polydata
vtkPolyData* CreateSpherePolydata( double * center, double radius )
{
  vtkSphereSource * sphere = vtkSphereSource::New();
  sphere->SetRadius( radius );
  sphere->SetCenter( center );
  sphere->SetThetaResolution( 18 );
  sphere->SetPhiResolution( 18 );
  sphere->Update();
  return sphere->GetOutput();
}

// Function to create the planar polydata
vtkPolyData* CreatePlanePolydata( double * origin,
                                  double * point1, double * point2,
                                  int resolution )
{
  vtkPlaneSource * plane = vtkPlaneSource::New();
  plane->SetOrigin( origin );
  plane->SetPoint1( point1 );
  plane->SetPoint2( point2 );
  plane->SetXResolution( resolution );
  plane->SetYResolution( resolution );
  plane->Update();
  return plane->GetOutput();
}

// Function to create the cube polydata
vtkPolyData* CreateCubePolydata( double * bottomBox, double * topBox,
                                 double * size )
{
  vtkCubeSource * topCube = vtkCubeSource::New();
  topCube->SetBounds( topBox[0], topBox[0] + size[0],
                      topBox[1] - 0.5, topBox[1] + size[1] - 0.5,
                      topBox[2], topBox[2] + size[2]);
  topCube->Update();

  vtkCubeSource * bottomCube = vtkCubeSource::New();
  bottomCube->SetBounds( bottomBox[0], bottomBox[0] + size[0],
                         bottomBox[1] - 0.5, bottomBox[1] + size[1] - 0.5,
                         bottomBox[2], bottomBox[2] + size[2] );
  bottomCube->Update();

  vtkAppendPolyData * append = vtkAppendPolyData::New();
  append->AddInput( topCube->GetOutput() );
  append->AddInput( bottomCube->GetOutput() );
  append->Update();

  vtkDensifyPolyData * densify = vtkDensifyPolyData::New();
  densify->AddInput( append->GetOutput() );
  densify->SetNumberOfSubdivisions( 6 );
  densify->Update();

  return densify->GetOutput();
}

int itkImageToImageAnisotropicDiffusiveDeformableRegistrationImageRegistrationTest(
                                                      int argc, char* argv [] )
{
  if( argc < 13 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << "original fixed image, "
              << "original moving image, "
              << "resulting motion field image, "
              << "resulting transformed moving image, "
              << "weight image, "
              << "normal vector image, "
              << "normal surface border polydata, "
              << "test type (0 for circles, 1 for squares), "
              << "number of iterations, "
              << "time step, "
              << "compute regularization term, "
              << "should use anisotropic regularization"
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
  typedef ImageType::IndexType                            IndexType;
  typedef ImageType::SizeType                             SizeType;
  typedef ImageType::SpacingType                          SpacingType;
  typedef ImageType::RegionType                           RegionType;

  //--------------------------------------------------------
  std::cout << "Generate input images and initial deformation field";
  std::cout << std::endl;

  bool circles = atoi( argv[8] ) == 0;

  // Image parameters
  double      sizeValue;
  if( circles )
    {
    sizeValue = 30;
    }
  else
    {
    sizeValue = 30;
    }

  double      originValue = 0.0;
  double      spacingValue = 1.0;

  ImageType::SizeValueType sizeArray[ImageDimension];
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    sizeArray[i] = sizeValue;
    }

  SizeType size;
  size.SetSize( sizeArray );

  IndexType index;
  index.Fill( originValue );

  SpacingType spacing;
  spacing.Fill( spacingValue );

  RegionType region;
  region.SetSize( size );
  region.SetIndex( index );

  ImageType::Pointer moving = ImageType::New();
  ImageType::Pointer fixed = ImageType::New();
  FieldType::Pointer initField = FieldType::New();

  moving->SetSpacing( spacing );
  moving->SetLargestPossibleRegion( region );
  moving->SetBufferedRegion( region );
  moving->Allocate();

  fixed->SetSpacing( spacing );
  fixed->SetLargestPossibleRegion( region );
  fixed->SetBufferedRegion( region );
  fixed->Allocate();

  initField->SetSpacing( spacing );
  initField->SetLargestPossibleRegion( region );
  initField->SetBufferedRegion( region );
  initField->Allocate();

  PixelType bgnd = 15;
  vtkPolyData * border;

  if ( circles )
    {
    double movingCenter[ImageDimension];
    double fixedCenter[ImageDimension];
    PixelType fgnd = 250;

    // fill moving with sphere
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      movingCenter[i] = sizeValue / 2.0;
      }
    double movingRadius = ( sizeValue / 3.0 ) - 2.0;
    FillWithSphere<ImageType>( moving, movingCenter, movingRadius, fgnd, bgnd );

    // fill fixed with sphere
    fixedCenter[0] = (sizeValue / 2.0) - 2.0;
    for ( unsigned int i = 1; i < ImageDimension; i++ )
      {
      fixedCenter[i] = sizeValue / 2.0;
      }
    double fixedRadius = sizeValue / 3.0;;
    FillWithSphere<ImageType>( fixed, fixedCenter, fixedRadius, fgnd, bgnd );

    // setup the normals
    border = CreateSpherePolydata( fixedCenter, fixedRadius );
    if( !border )
      {
      std::cerr << "Could not generate sphere surface" << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    double boxSize[3] = { sizeValue / (8.0/3.0),
                          sizeValue / 5.0,
                          sizeValue / 5.0 };
    double center[3] = {sizeValue / 2.0, sizeValue / 2.0, sizeValue / 2.0 };
    double offset = sizeValue / 8.0;
    PixelType bottomStart = 30;
    PixelType bottomEnd = 120;
    PixelType topStart = 130;
    PixelType topEnd = 220;

    // Create the two boxes on the fixed image
    double fixedBottomBox[3] = { center[0] - boxSize[0] + offset,
                                 center[1] - boxSize[1],
                                 center[2] - boxSize[2] / 2 };
    double fixedTopBox[3] = { center[0] - offset,
                              center[1],
                              center[2] - boxSize[2] / 2 };

    FillWithBox<ImageType>( fixed, fixedBottomBox, fixedTopBox, boxSize,
                            bgnd, bottomStart, bottomEnd, topStart, topEnd );

    // Create the two boxes on the moving image
    double shift = sizeValue / 20.0; // <-- this is the transformation!!!
    double movingBottomBox[3] = { fixedBottomBox[0] - shift,
                                  fixedBottomBox[1],
                                  fixedBottomBox[2] };
    double movingTopBox[3] = { fixedTopBox[0] + shift,
                                fixedTopBox[1],
                                fixedTopBox[2] };
    FillWithBox<ImageType>( moving, movingBottomBox, movingTopBox, boxSize,
                            bgnd, bottomStart, bottomEnd, topStart, topEnd );

    // setup the normals for the cubes
    double spacer = sizeValue / 8.0;
    double bottomCubeBorder[3] = { 0.0 - spacer,
                                   fixedBottomBox[1],
                                   fixedBottomBox[2] };
    double topCubeBorder[3] = { 0.0 - spacer,
                                fixedTopBox[1],
                                fixedTopBox[2] };
    double cubeSize[3] = { sizeValue + ( spacer * 2.0 ),
                           boxSize[1],
                           boxSize[2] };

    border = CreateCubePolydata( bottomCubeBorder, topCubeBorder, cubeSize );

//    // setup the normals for a plane
//    double origin[3] = { 0.0, center[1] - 0.5, 0.0 };
//    double point1[3] = {sizeValue, center[1] - 0.5, 0.0 };
//    double point2[3] = {0.0, center[1] - 0.5, sizeValue };
//    int resolution = 30;
//    border = CreatePlanePolydata( origin, point1, point2, resolution );

    if( !border )
      {
      std::cerr << "Could not generate planar surface" << std::endl;
      return EXIT_FAILURE;
      }
    }

  // fill initial deformation with zero vectors
  VectorType zeroVec;
  zeroVec.Fill( 0.0 );
  initField->FillBuffer( zeroVec );

  // ---------------------------------------------------------
  std::cout << "Printing the initial fixed and moving images" << std::endl;

  // Save the initial fixed and moving images
  typedef itk::ImageFileWriter< ImageType > ImageWriterType;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetFileName( argv[1] );
  imageWriter->SetInput( fixed );
  try
    {
    imageWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }
  imageWriter->SetFileName( argv[2] );
  imageWriter->SetInput( moving );
  try
    {
    imageWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::VectorCastImageFilter<FieldType,FieldType> CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( initField );
  caster->InPlaceOff();

  //-------------------------------------------------------------
  std::cout << "Run registration and warp moving" << std::endl;

  typedef itk::ImageToImageAnisotropicDiffusiveDeformableRegistrationFilter
      < ImageType, ImageType, FieldType > RegistrationType;
  RegistrationType::Pointer registrator = RegistrationType::New();

  registrator->SetInitialDeformationField( caster->GetOutput() );
  registrator->SetMovingImage( moving );
  registrator->SetFixedImage( fixed );
  registrator->SetBorderSurface( border );
  int numberOfIterations = atoi( argv[9] );
  registrator->SetNumberOfIterations( numberOfIterations );

  int compute = atoi( argv[11] );
  if ( compute )
    {
    registrator->SetComputeRegularizationTerm( true );
    }
  else
    {
    registrator->SetComputeRegularizationTerm( false );
    }

  int useAnisotropic = atoi( argv[12] );
  if ( useAnisotropic )
    {
    registrator->SetUseAnisotropicRegularization( true );
    }
  else
    {
    registrator->SetUseAnisotropicRegularization( false );
    }

  registrator->SetTimeStep( atof( argv[10] ) );
  registrator->SetLambda( -0.05 );

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
  warper->SetEdgePaddingValue( bgnd );

  // Update triggers the registration
  warper->Update();

  // ---------------------------------------------------------
  std::cout << "Printing the normal surface border, normal vector image "
      << "and weight image" << std::endl;

  if( compute && useAnisotropic )
    {
    vtkPolyData * normalPolyData = registrator->GetBorderNormalsSurface();
    vtkPolyDataWriter * polyWriter = vtkPolyDataWriter::New();
    polyWriter->SetFileName( argv[7] );
    polyWriter->SetInput( normalPolyData );
    polyWriter->Write();

    typedef itk::ImageFileWriter< VectorImageType > VectorWriterType;
    VectorWriterType::Pointer vectorWriter = VectorWriterType::New();
    vectorWriter->SetFileName( argv[6] );
    vectorWriter->SetInput( registrator->GetNormalVectorImage() );
    vectorWriter->Write();

    typedef itk::ImageFileWriter< WeightImageType > WeightWriterType;
    WeightWriterType::Pointer weightWriter = WeightWriterType::New();
    weightWriter->SetFileName( argv[5] );
    weightWriter->SetInput( registrator->GetWeightImage() );
    weightWriter->Write();
    }

  // ---------------------------------------------------------
  std::cout << "Printing the deformation field and transformed moving image"
            << std::endl;

  typedef itk::ImageFileWriter< FieldType > FieldWriterType;
  FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
  fieldWriter->SetFileName( argv[3] );
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

  imageWriter->SetFileName( argv[4] );
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
