/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#include "tubeTubeMathFilters.h"

#include <itkImageFileWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectToImageFilter.h>

int itktubePointsToImageTest( int argc, char * argv[] )
{
  if( argc != 3 )
    {
    std::cout
      << "itktubePointsToImageTest <vessel.tre> <outputImage>"
      << std::endl;
    return EXIT_FAILURE;
    }

  enum { Dimension = 3 };

  typedef itk::SpatialObjectReader< Dimension >      ReaderType;
  typedef itk::SpatialObject< Dimension >::ChildrenListType
      ObjectListType;
  typedef itk::SpatialObject< Dimension >            SpatialObjectType;
  typedef itk::GroupSpatialObject< Dimension >       GroupType;
  typedef itk::TubeSpatialObject< Dimension >        TubeType;
  typedef TubeType::TubePointListType                TubePointListType;
  typedef TubeType::TubePointType                    TubePointType;
  typedef itk::Image< float, Dimension >             ImageType;

  // Read the tube
  ReaderType::Pointer tubeReader = ReaderType::New();
  tubeReader->SetFileName( argv[1] );
  tubeReader->Update();

  // Get the group
  GroupType::Pointer group;
  SpatialObjectType::Pointer soScene = tubeReader->GetGroup();
  if( soScene->GetTypeName() == "GroupSpatialObject" )
    {
    group = static_cast< GroupType * >( soScene.GetPointer() );
    }
  else
    {
    group = tubeReader->GetGroup();
    }

  const unsigned int numberOfChildren = group->GetNumberOfChildren(-1);
  std::cout << "Number of children = " << numberOfChildren
    << std::endl;
  if( numberOfChildren != 2 )
    {
    std::cerr << "Wrong number of children" << std::endl;
    return EXIT_FAILURE;
    }

  // Get the list of tubes from the group
  char tubeName[17];
  std::strcpy( tubeName, "Tube" );
  ObjectListType * tubeList = group->GetChildren( -1, tubeName );
  const unsigned int numTubes = tubeList->size();
  std::cout << "Number of tubes: " << numTubes << std::endl;
  if( numTubes != 2 )
    {
    std::cerr << "Wrong number of tubes" << std::endl;
    delete tubeList;
    return EXIT_FAILURE;
    }

  // Iterate through the tubes
  ObjectListType::iterator            tubeIt = tubeList->begin();
  TubeType::Pointer                   tube;
  TubePointListType                   tubePointList;
  TubePointListType::iterator         pointIt;
  TubePointType *                     point;
  TubePointType::PointType            position;
  float                               radius;
  TubePointType::CovariantVectorType  normal1;
  TubePointType::CovariantVectorType  normal2;

  for( unsigned int i = 0; i < numTubes; i++ )
    {
    std::cout << "----tube " << i << std::endl;
    tube = static_cast< TubeType * >( tubeIt->GetPointer() );
    tube->ComputeTangentAndNormals();

    // Get the list of points in the tube
    tubePointList = tube->GetPoints();
    unsigned int numPoints = tubePointList.size();

    // Iterate through the points
    pointIt = tubePointList.begin();
    for( unsigned int j = 0; j < numPoints; j++ )
      {
      point = static_cast< TubePointType * >( &( *pointIt ) );
      position = point->GetPositionInObjectSpace();
      radius = point->GetRadiusInObjectSpace();
      normal1 = point->GetNormal1InObjectSpace();
      normal2 = point->GetNormal2InObjectSpace();

      std::cout << "point " << j << "; "
                << "position " << position
                << "radius " << radius
                << "normal1 " << normal1
                << "normal2 " << normal2 << std::endl;
      ++pointIt;
      }
    ++tubeIt;
    }
  delete tubeList;

  typedef itk::SpatialObjectToImageFilter<GroupType, ImageType>
                                              SpatialObjectToImageFilterType;
  SpatialObjectToImageFilterType::Pointer tubeToImageFilter =
    SpatialObjectToImageFilterType::New();

  tubeToImageFilter->SetInput( group );

  ImageType::SizeType   size;
  size[0] = 100;
  size[1] = 100;
  size[2] = 100;
  tubeToImageFilter->SetSize( size );

  tubeToImageFilter->SetInsideValue( 1.0 );
  tubeToImageFilter->SetOutsideValue( 0.0 );
  tubeToImageFilter->Update();
  typedef itk::ImageFileWriter< ImageType > ImageWriterType;

  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetFileName( argv[2] );
  imageWriter->SetUseCompression( true );
  imageWriter->SetInput( tubeToImageFilter->GetOutput() );
  imageWriter->Update();
  return EXIT_SUCCESS;
}
