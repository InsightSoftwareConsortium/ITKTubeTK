/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include "itkTubeSpatialObject.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

#include "tubeMacro.h"
#include "tubeTubeMathFilters.h"

double
tubeLength(itk::TubeSpatialObject<3> * tube)
{
  typedef itk::TubeSpatialObject<3> TubeType;

  TubeType::TubePointListType           tubePointList = tube->GetPoints();
  TubeType::TubePointListType::iterator tubePointItr;

  TubeType::TubePointListType::iterator beginItr = tubePointList.begin();
  TubeType::TubePointListType::iterator endItr = tubePointList.end();
  tubePointItr = beginItr;

  TubeType::PointType lastX = tubePointItr->GetPositionInObjectSpace();
  ++tubePointItr;

  double length = 0;
  while (tubePointItr != endItr)
  {
    double dist = 0;
    for (unsigned int d = 0; d < 3; ++d)
    {
      dist += (tubePointItr->GetPositionInObjectSpace()[d] - lastX[d]) *
              (tubePointItr->GetPositionInObjectSpace()[d] - lastX[d]);
    }
    length += std::sqrt(dist);
    ++tubePointItr;
  }

  return length;
}


int
tubeTubeMathFiltersTest(int tubeNotUsed(argc), char * tubeNotUsed(argv)[])
{
  bool returnStatus = EXIT_SUCCESS;

  typedef itk::TubeSpatialObject<3> TubeType;

  TubeType::Pointer tube = TubeType::New();
  TubeType::Pointer tube0 = TubeType::New();

  TubeType::TubePointListType tubePointList;
  TubeType::TubePointType     tubePoint;

  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator RandomType;
  RandomType::Pointer                                            rndGen = RandomType::New();
  rndGen->Initialize(1);

  // Build the vessel
  TubeType::TubePointType::PointType pnt;
  pnt.Fill(0);
  for (unsigned int i = 0; i < 100; ++i)
  {
    double rndX = rndGen->GetNormalVariate(0, 1);
    double rndY = rndGen->GetNormalVariate(0, 1);
    double rndZ = rndGen->GetNormalVariate(0, 1);
    pnt[0] += rndX;
    pnt[1] += rndY;
    pnt[2] += rndZ;
    tubePoint.SetPositionInObjectSpace(pnt);
    tubePointList.push_back(tubePoint);
  }

  std::cout << "Tangents and normals..." << std::endl;
  tube0->SetPoints(tubePointList);
  tube0->ComputeTangentsAndNormals();

  // Run tests
  std::cout << std::endl << "*********** Smoothing Tests *************" << std::endl;

  std::cout << std::endl << "-----> Smoothing by average on index" << std::endl;
  tube = tube0->Clone();
  double tLength = tubeLength(tube);
  std::cout << "Length = " << tLength << std::endl;

  ::tube::TubeMathFilters<3> filter;
  for (unsigned int i = 0; i < 20; ++i)
  {
    std::cout << "Smoothing..." << std::endl;
    TubeType::Pointer tube2 = tube->Clone();
    filter.SetInputTube(tube2);
    filter.SmoothTube(2, ::tube::TubeMathFilters<3>::SMOOTH_TUBE_USING_INDEX_AVERAGE);
    // tube2 = filter.GetOutput();

    tLength = tubeLength(tube);
    std::cout << "Comp Length = " << tLength << std::endl;
    double t2Length = tubeLength(tube2);
    std::cout << "   with Length = " << t2Length << std::endl;
    if (tLength <= t2Length)
    {
      std::cerr << "ERROR: Raw length = " << tLength << " <= Smooth length = " << t2Length << std::endl;
      returnStatus = EXIT_FAILURE;
    }

    tube = tube2;
    tLength = t2Length;
  }

  std::cout << std::endl << "-----> Smoothing by gaussian on index" << std::endl;
  tube = tube0->Clone();
  tLength = tubeLength(tube);
  std::cout << "Length = " << tLength << std::endl;

  filter.SetInputTube(tube);
  for (unsigned int i = 0; i < 20; ++i)
  {
    std::cout << "Smoothing..." << std::endl;
    filter.SmoothTube(4, ::tube::TubeMathFilters<3>::SMOOTH_TUBE_USING_INDEX_GAUSSIAN);

    double t2Length = tubeLength(tube);
    std::cout << "Length = " << t2Length << std::endl;
    if (tLength <= t2Length)
    {
      std::cerr << "ERROR: Raw length = " << tLength << " <= Smooth length = " << t2Length << std::endl;
      returnStatus = EXIT_FAILURE;
    }

    tLength = t2Length;
  }


  std::cout << std::endl << "*********** Subsampling Tests *************" << std::endl;
  tube = tube0->Clone();
  tLength = tubeLength(tube);
  int tNumPoints = tube->GetNumberOfPoints();
  std::cout << "Length = " << tLength << std::endl;
  std::cout << "Number of Points = " << tNumPoints << std::endl;

  filter.SetInputTube(tube);
  for (unsigned int i = 0; i < 5; ++i)
  {
    std::cout << "Subsampling..." << std::endl;
    filter.SubsampleTube(2);

    double t3Length = tubeLength(tube);
    int    t3NumPoints = tube->GetNumberOfPoints();
    std::cout << "Length = " << t3Length << std::endl;
    std::cout << "Number of Points = " << t3NumPoints << std::endl;
    if (tLength <= t3Length)
    {
      std::cerr << "ERROR: Raw length = " << tLength << " <= Smooth length = " << t3Length << std::endl;
      returnStatus = EXIT_FAILURE;
    }
    if (t3NumPoints >= tNumPoints)
    {
      std::cerr << "ERROR: Raw: NumberOfPoints = " << tNumPoints << " <= Subsampled: NumberOfPoints = " << t3NumPoints
                << std::endl;
      returnStatus = EXIT_FAILURE;
    }

    tLength = t3Length;
    tNumPoints = t3NumPoints;
  }

  return returnStatus;
}
