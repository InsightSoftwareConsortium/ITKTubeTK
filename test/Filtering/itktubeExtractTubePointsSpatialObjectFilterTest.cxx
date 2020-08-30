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

#include "itktubeExtractTubePointsSpatialObjectFilter.h"

#include <itkSpatialObjectReader.h>

int itktubeExtractTubePointsSpatialObjectFilterTest( int argc, char* argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: "
      << "inputTubeTree "
      << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputTubeTree = argv[1];

  static const unsigned int Dimension = 3;
  typedef itk::TubeSpatialObject< Dimension >   TubeSpatialObjectType;
  typedef itk::GroupSpatialObject< Dimension >  GroupSpatialObjectType;

  // Read input tube tree.
  typedef itk::SpatialObjectReader< Dimension >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputTubeTree );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Read Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cerr << "Read failed." << std::endl;
    return EXIT_FAILURE;
    }

  GroupSpatialObjectType::Pointer groupSpatialObject = reader->GetGroup();
  std::cout << "Number of children = "
    << groupSpatialObject->GetNumberOfChildren() << std::endl;

  typedef itk::tube::ExtractTubePointsSpatialObjectFilter<
    TubeSpatialObjectType > ExtractTubePointsSpatialObjectFilterType;
  ExtractTubePointsSpatialObjectFilterType::Pointer extractTubePointsFilter =
    ExtractTubePointsSpatialObjectFilterType::New();
  extractTubePointsFilter->SetInput( reader->GetGroup() );
  try
    {
    extractTubePointsFilter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Update Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cerr << "Update failed." << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "extractTubePointsFilter->Update() SUCCESS" << std::endl;

  typedef ExtractTubePointsSpatialObjectFilterType::PointsContainerType
    PointsContainerType;

  const PointsContainerType * pointsContainer =
    extractTubePointsFilter->GetPointsContainer();
  if( pointsContainer == NULL )
    {
    std::cerr << "Point container is null" << std::endl;
    return EXIT_FAILURE;
    }

  const itk::SizeValueType numberOfPoints = pointsContainer->Size();
  std::cout << "Output points container size: " << numberOfPoints
    << std::endl;
  const itk::SizeValueType expectedNumberOfPoints = 110658;
  if( numberOfPoints != expectedNumberOfPoints )
    {
    std::cerr << "Did not get expected size." << std::endl;
    return EXIT_FAILURE;
    }
  for( unsigned int ii = 0; ii < 3; ++ii )
    {
    const ExtractTubePointsSpatialObjectFilterType::TubePointType &
      tubePoint = pointsContainer->ElementAt( ii );
    std::cout << "Point index: " << ii << '\n'
              << "  Point Id:       " << tubePoint.GetId() << '\n'
              << "  Point Position: " << tubePoint.GetPositionInObjectSpace() << '\n'
              << "  Point Radius:   " << tubePoint.GetRadiusInObjectSpace() << '\n'
              << "  Point Normal1:  " << tubePoint.GetNormal1InObjectSpace() << '\n'
              << "  Point Normal2:  " << tubePoint.GetNormal2InObjectSpace() << '\n'
              << "  Point Tangent:  " << tubePoint.GetTangentInObjectSpace() << '\n'
              << std::endl;
    }

  //delete pointsContainer;

  return EXIT_SUCCESS;
}
