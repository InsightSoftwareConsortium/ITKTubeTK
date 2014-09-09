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

#include <itkImageFileWriter.h>
#include <itkSpatialObjectWriter.h>

#include <vtkAppendPolyData.h>
#include <vtkCubeSource.h>
#include <vtkDensifyPolyData.h>
#include <vtkPlaneSource.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <vtkVersion.h>

// Template function to fill in an image with a sphere.
template< class TImage >
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
  it.GoToBegin();

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
template< class TImage >
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
typename TImage::PixelType topEnd )
{
  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  Iterator it( image, image->GetBufferedRegion() );
  it.GoToBegin();

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

// Determines whether a point is in a tube
template< class TIndex >
bool PointInTube( TIndex index, double * tubeLeftPoint, double radius)
{
  TIndex centerPoint;
  centerPoint[0] = index[0];
  centerPoint[1] = tubeLeftPoint[1];
  centerPoint[2] = tubeLeftPoint[2];

  double distance = 0;
  for(int i = 0; i < 3; i++)
    {
    distance += vnl_math_sqr(index[i] - centerPoint[i]);
    }

  return vcl_sqrt(distance) <= radius;
}

// Template function to fill in an image with two tubes
template< class TImage >
void
FillWithTubes(
TImage * image,
double * bottomTubeLeftPoint,
double * topTubeLeftPoint,
double length,
double radius,
typename TImage::PixelType backgnd,
typename TImage::PixelType bottomIntensity,
typename TImage::PixelType topIntensity )
{
  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  Iterator it( image, image->GetBufferedRegion() );
  it.GoToBegin();

  typename TImage::IndexType index;

  for( ; !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    bool inTube = false;

    if( index[0] > bottomTubeLeftPoint[0] && index[0] < bottomTubeLeftPoint[0] + length )
      {
      if( PointInTube(index, bottomTubeLeftPoint, radius) )
        {
        it.Set( bottomIntensity );
        inTube = true;
        }
      }

    if( index[0] > topTubeLeftPoint[0] && index[0] < topTubeLeftPoint[0] + length )
      {
      if( PointInTube(index, topTubeLeftPoint, radius) )
        {
        it.Set( topIntensity );
        inTube = true;
        }
      }

    if( !inTube )
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
#if VTK_MAJOR_VERSION > 5
  append->AddInputData( topCube->GetOutput() );
  append->AddInputData( bottomCube->GetOutput() );
#else
  append->AddInput( topCube->GetOutput() );
  append->AddInput( bottomCube->GetOutput() );
#endif
  append->Update();

  vtkDensifyPolyData * densify = vtkDensifyPolyData::New();
#if VTK_MAJOR_VERSION > 5
  densify->AddInputData( append->GetOutput() );
#else
  densify->AddInput( append->GetOutput() );
#endif
  densify->SetNumberOfSubdivisions( 6 );
  densify->Update();

  return densify->GetOutput();
}


// Intensity windowing from 0..255 to 0..1
template< class TImage >
void IntensityWindow(
TImage * image )
{
  float valMin = 0;
  float valMax = 255;
  float outMin = 0;
  float outMax = 1;
  itk::ImageRegionIterator< TImage > it2( image,
                                          image->GetLargestPossibleRegion() );
  it2.GoToBegin();
  while( !it2.IsAtEnd() )
    {
    double tf = it2.Get();
    tf = ( tf-valMin )/( valMax-valMin );
    if( tf<0 )
      {
      tf = 0;
      }
    if( tf>1 )
      {
      tf = 1;
      }
    tf = ( tf * ( outMax-outMin ) ) + outMin;
    it2.Set( ( typename TImage::PixelType )tf );
    ++it2;
    }
}


int itkAnisotropicDiffusiveRegistrationGenerateTestingImages( int argc, char * argv[] )
{
  if( argc < 7 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << "output fixed image, "
              << "output moving image, "
              << "output surface border polydata or tube spatial objects, "
              << "test type (0 for circles, 1 for boxes, 2 for tubes), "
              << "intensity window to [0..1] (0 = no, 1 = yes), "
              << "image size (creates square images)"
              << std::endl;
    return EXIT_FAILURE;
    }

  // Typedefs
  const unsigned int                                          ImageDimension = 3;
  typedef double                                              PixelType;
  typedef itk::Image< PixelType, ImageDimension >             ImageType;
  typedef ImageType::IndexType                                IndexType;
  typedef ImageType::SizeType                                 SizeType;
  typedef ImageType::SpacingType                              SpacingType;
  typedef ImageType::RegionType                               RegionType;
  typedef itk::GroupSpatialObject< ImageDimension >           GroupType;
  typedef itk::VesselTubeSpatialObjectPoint< ImageDimension > VectorTubePointType;
  typedef itk::VesselTubeSpatialObject< ImageDimension >      VesselTubeType;

  //--------------------------------------------------------
  std::cout << "Generate registration input images" << std::endl;

  // Image parameters
  ImageType::SizeValueType sizeValue = std::atoi( argv[6] );
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

  moving->SetSpacing( spacing );
  moving->SetLargestPossibleRegion( region );
  moving->SetBufferedRegion( region );
  moving->Allocate();

  fixed->SetSpacing( spacing );
  fixed->SetLargestPossibleRegion( region );
  fixed->SetBufferedRegion( region );
  fixed->Allocate();

  PixelType bgnd = 15;
  vtkPolyData * border = 0;
  GroupType::Pointer group = 0;

  enum geometryTypes { circles, boxes, tubes };

  int geometry = std::atoi( argv[4] );
  if( geometry == circles )
    {
    double movingCenter[ImageDimension];
    double fixedCenter[ImageDimension];
    PixelType fgnd = 250;

    // fill moving with sphere
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      movingCenter[i] = sizeValue / 2.0;
      }
    double movingRadius = ( sizeValue / 3.0 ) - 2.0;
    FillWithSphere<ImageType>( moving, movingCenter, movingRadius, fgnd, bgnd );

    // fill fixed with sphere
    fixedCenter[0] = (sizeValue / 2.0) - 2.0;
    for( unsigned int i = 1; i < ImageDimension; i++ )
      {
      fixedCenter[i] = sizeValue / 2.0;
      }
    double fixedRadius = sizeValue / 3.0;
    FillWithSphere<ImageType>( fixed, fixedCenter, fixedRadius, fgnd, bgnd );

    // setup the normals
    border = CreateSpherePolydata( fixedCenter, fixedRadius );
    if( !border )
      {
      std::cerr << "Could not generate sphere surface" << std::endl;
      return EXIT_FAILURE;
      }
    }
  else if( geometry == boxes )
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

    if( !border )
      {
      std::cerr << "Could not generate planar surface" << std::endl;
      return EXIT_FAILURE;
      }
    }
  else if( geometry == tubes )
    {
    double length = sizeValue / (8.0/3.0);
    double radius = sizeValue / 10.0;
    double center[3] = {sizeValue / 2.0, sizeValue / 2.0, sizeValue / 2.0 };
    double offset = sizeValue / 8.0;
    PixelType bottomIntensity = 30;
    PixelType topIntensity = 220;

    // Create the two tubes on the fixed image
    double fixedBottomTubeLeftPoint[3] = { center[0] - ( length / 2 ) - offset,
                                           center[1] - radius,
                                           center[2] };
    double fixedTopTubeLeftPoint[3] = { center[0] - ( length / 2 ) + offset,
                                        center[1] + radius,
                                        center[2] };
    FillWithTubes<ImageType>( fixed, fixedBottomTubeLeftPoint, fixedTopTubeLeftPoint,
                              length, radius, bgnd, bottomIntensity, topIntensity );

    // Create the two boxes on the moving image
    double shift = sizeValue / 20.0; // <-- this is the transformation!!!
    double movingBottomTubeLeftPoint[3] = { fixedBottomTubeLeftPoint[0] - shift,
                                            fixedBottomTubeLeftPoint[1],
                                            fixedBottomTubeLeftPoint[2] };
    double movingTopTubeLeftPoint[3] = { fixedTopTubeLeftPoint[0] + shift,
                                         fixedTopTubeLeftPoint[1],
                                         fixedTopTubeLeftPoint[2] };
    FillWithTubes<ImageType>( moving, movingBottomTubeLeftPoint, movingTopTubeLeftPoint,
                              length, radius, bgnd, bottomIntensity, topIntensity );

    // Create tube spatial objects
    int numPoints = 0;
    double bottomX = 0;
    double topX = 0;

    bool defineShortTube = false;
    if( defineShortTube )
      {
      numPoints = static_cast<int>(length) + 1;
      bottomX = vnl_math_ceil( fixedBottomTubeLeftPoint[0] );
      topX = vnl_math_ceil( fixedTopTubeLeftPoint[0] );
      }
    else
      {
      numPoints = sizeValue + 1;
      bottomX = 0;
      topX = 0;
      }

    VesselTubeType::PointListType bottomTubePoints;
    bottomTubePoints.resize(numPoints);
    VesselTubeType::PointListType topTubePoints;
    topTubePoints.resize(numPoints);

    for( int i = 0; i < numPoints; i++ )
      {
      VectorTubePointType bottomPoint;
      bottomPoint.SetPosition(bottomX, fixedBottomTubeLeftPoint[1], fixedBottomTubeLeftPoint[2]);
      bottomPoint.SetNormal1(0, 1, 0);
      bottomPoint.SetNormal2(0, 0, 1);
      bottomPoint.SetTangent(1, 0, 0);
      bottomPoint.SetRadius(radius);
      bottomTubePoints[i] = bottomPoint;

      VectorTubePointType topPoint;
      topPoint.SetPosition(topX, fixedTopTubeLeftPoint[1], fixedTopTubeLeftPoint[2]);
      topPoint.SetNormal1(0, 1, 0);
      topPoint.SetNormal2(0, 0, 1);
      topPoint.SetTangent(1, 0, 0);
      topPoint.SetRadius(radius);
      topTubePoints[i] = topPoint;

      bottomX += spacingValue;
      topX += spacingValue;
      }

    VesselTubeType::Pointer bottomTube = VesselTubeType::New();
    bottomTube->SetPoints(bottomTubePoints);
    VesselTubeType::Pointer topTube = VesselTubeType::New();
    topTube->SetPoints(topTubePoints);

    group = GroupType::New();
    group->AddSpatialObject(bottomTube);
    group->AddSpatialObject(topTube);
    }

  // Scale the images to 0..1
  bool intensityWindow = std::atoi( argv[5] ) == 1;
  if(intensityWindow)
    {
    IntensityWindow<ImageType>( fixed );
    IntensityWindow<ImageType>( moving );
    }

  // ---------------------------------------------------------
  std::cout << "Saving the initial fixed and moving images" << std::endl;

  // Save the initial fixed and moving images
  typedef itk::ImageFileWriter< ImageType > ImageWriterType;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetFileName( argv[1] );
  imageWriter->SetUseCompression( true );
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
  imageWriter->SetUseCompression( true );
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
  if( geometry == circles || geometry == boxes )
    {
    std::cout << "Saving the normal surface border" << std::endl;

    // Compute normals on the polydata
    vtkPolyDataNormals * normalsFilter = vtkPolyDataNormals::New();
    normalsFilter->ComputePointNormalsOn();
    normalsFilter->ComputeCellNormalsOff();
#if VTK_MAJOR_VERSION > 5
    normalsFilter->SetInputData( border );
#else
    normalsFilter->SetInput( border );
#endif
    normalsFilter->Update();

    vtkPolyDataWriter * polyWriter = vtkPolyDataWriter::New();
    polyWriter->SetFileName( argv[3] );
#if VTK_MAJOR_VERSION > 5
    polyWriter->SetInputData( normalsFilter->GetOutput() );
#else
    polyWriter->SetInput( normalsFilter->GetOutput() );
#endif
    polyWriter->Write();
    }
  else if( geometry == tubes )
    {
    std::cout << "Saving the tube spatial objects" << std::endl;

    typedef itk::SpatialObjectWriter< ImageDimension > TubeWriterType;
    TubeWriterType::Pointer tubeWriter = TubeWriterType::New();
    tubeWriter->SetFileName( argv[3] );
    tubeWriter->SetInput( group );
    tubeWriter->Update();
    }

  return EXIT_SUCCESS;
}
