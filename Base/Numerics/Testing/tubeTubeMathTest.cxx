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

#include "tubeTubeMath.h"

#include "itkTubeSpatialObject.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

double tubeLength( itk::TubeSpatialObject<3> * tube )
{
  typedef itk::TubeSpatialObject<3>   TubeType;

  TubeType::PointListType pointList = tube->GetPoints();
  TubeType::PointListType::iterator pointItr;

  TubeType::PointListType::iterator beginItr = pointList.begin();
  TubeType::PointListType::iterator endItr = pointList.end();
  pointItr = beginItr;

  TubeType::PointType lastX = pointItr->GetPosition();
  ++pointItr;

  double length = 0;
  while( pointItr != endItr )
    {
    double dist = 0;
    for( unsigned int d=0; d<3; ++d )
      {
      dist += ( pointItr->GetPosition()[d] - lastX[d] ) *
        ( pointItr->GetPosition()[d] - lastX[d] );
      }
    length += std::sqrt( dist );
    ++pointItr;
    }

  return length;
}


int tubeTubeMathTest( int tubeNotUsed( argc ),
  char * tubeNotUsed( argv )[] )
{
  bool returnStatus = EXIT_SUCCESS;

  typedef itk::TubeSpatialObject<3>   TubeType;

  TubeType::Pointer tube = TubeType::New();
  TubeType::Pointer tube0 = TubeType::New();

  TubeType::PointListType pointList;
  TubeType::TubePointType point;

  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator RandomType;
  RandomType::Pointer rndGen = RandomType::New();
  rndGen->Initialize( 1 );

  // Build the vessel
  double x = 0;
  double y = 0;
  double z = 0;
  for( unsigned int i = 0; i < 100; ++i )
    {
    double rndX = rndGen->GetNormalVariate( 0, 1 );
    double rndY = rndGen->GetNormalVariate( 0, 1 );
    double rndZ = rndGen->GetNormalVariate( 0, 1 );
    x += rndX;
    y += rndY;
    z += rndZ;
    point.SetPosition( x, y, z );
    pointList.push_back( point );
    }

  std::cout << "Tangents and normals..." << std::endl;
  tube0->SetPoints( pointList );
  ::tube::ComputeTubeTangentsAndNormals< TubeType >( tube0 );

  // Run tests
  std::cout << std::endl << "*********** Smoothing Tests *************" << std::endl;

  std::cout << std::endl << "-----> Smoothing by average on index" << std::endl;
  tube = tube0;
  double tLength = tubeLength( tube );
  std::cout << "Length = " << tLength << std::endl;

  for( unsigned int i = 0; i < 20; ++i )
    {
    std::cout << "Smoothing..." << std::endl;
    TubeType::Pointer tube2 = ::tube::SmoothTube< TubeType >( tube, 2,
      ::tube::SMOOTH_TUBE_USING_INDEX_AVERAGE );

    double t2Length = tubeLength( tube2 );
    std::cout << "Length = " << t2Length << std::endl;
    if( tLength <= t2Length )
      {
      std::cerr << "ERROR: Raw length = " << tLength
        << " <= Smooth length = " << t2Length << std::endl;
      returnStatus = EXIT_FAILURE;
      }

    tube = tube2;
    tLength = t2Length;
    }

  std::cout << std::endl << "-----> Smoothing by gaussian on index" << std::endl;
  tube = tube0;
  tLength = tubeLength( tube );
  std::cout << "Length = " << tLength << std::endl;

  for( unsigned int i = 0; i < 20; ++i )
    {
    std::cout << "Smoothing..." << std::endl;
    TubeType::Pointer tube2 = ::tube::SmoothTube< TubeType >( tube, 4,
      ::tube::SMOOTH_TUBE_USING_INDEX_GAUSSIAN );

    double t2Length = tubeLength( tube2 );
    std::cout << "Length = " << t2Length << std::endl;
    if( tLength <= t2Length )
      {
      std::cerr << "ERROR: Raw length = " << tLength
        << " <= Smooth length = " << t2Length << std::endl;
      returnStatus = EXIT_FAILURE;
      }

    tube = tube2;
    tLength = t2Length;
    }


  std::cout << std::endl << "*********** Subsampling Tests *************" << std::endl;
  tube = tube0;
  tLength = tubeLength( tube );
  int tNumPoints = tube->GetNumberOfPoints();
  std::cout << "Length = " << tLength << std::endl;
  std::cout<< "Number of Points = " << tNumPoints <<std::endl;

  for( unsigned int i = 0; i < 5; ++i )
    {
    std::cout << "Subsampling..." << std::endl;
    TubeType::Pointer tube3 = ::tube::SubsampleTube< TubeType >( tube, 2 );

    double t3Length = tubeLength( tube3 );
    int t3NumPoints = tube3->GetNumberOfPoints();
    std::cout << "Length = " << t3Length << std::endl;
    std::cout<< "Number of Points = " << t3NumPoints <<std::endl;
    if( tLength <= t3Length )
      {
      std::cerr << "ERROR: Raw length = " << tLength
        << " <= Smooth length = " << t3Length << std::endl;
      returnStatus = EXIT_FAILURE;
      }
    if( t3NumPoints >= tNumPoints )
      {
      std::cerr << "ERROR: Raw: NumberOfPoints = " << tube->GetNumberOfPoints()
        << " <= Subsampled: NumberOfPoints = " << tube3->GetNumberOfPoints()
        << std::endl;
      returnStatus = EXIT_FAILURE;
      }

    tube = tube3;
    tLength = t3Length;
    tNumPoints = t3NumPoints;
    }

  return returnStatus;
}
