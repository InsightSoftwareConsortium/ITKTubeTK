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

// Tortuostity includes
#include <itktubeTortuositySpatialObjectFilter.h>
#include <itkVesselTubeSpatialObject.h>

// ITK includes
#include <itkMath.h>

typedef itk::VesselTubeSpatialObject<3> VesselTubeType;
typedef VesselTubeType::VectorType VectorType;
typedef itk::tube::TortuositySpatialObjectFilter<VesselTubeType> FilterType;

//-----------------------------------------------------------------------------
VesselTubeType::Pointer GenerateStraightTube(VectorType start,
                                             VectorType increment,
                                             int numberOfPoints)
{
  VesselTubeType::Pointer vessel = VesselTubeType::New();
  VesselTubeType::PointListType pointList;

  VectorType pos = start;
  for (int i = 0; i < numberOfPoints; ++i)
    {
    VesselTubeType::TubePointType point;
    point.SetPosition(pos[0], pos[1], pos[2]);
    pointList.push_back(point);

    pos += increment;
    }

  vessel->SetPoints(pointList);
  return vessel;
}

double Modulus(double v, double mod)
{
  if (v < mod)
    {
    return v;
    }
  return Modulus(v-mod, mod);
}

//-----------------------------------------------------------------------------
VesselTubeType::Pointer GenerateCosTube(double length,
                                        double amplitude,
                                        double frequency
                                        )
{
  VesselTubeType::Pointer vessel = VesselTubeType::New();
  VesselTubeType::PointListType pointList;

  VectorType pos(0.0), lastPos(0.0);
  double curveLength = 0.0;
  bool lastPosValid = false;

  int numberOfSamples = 1e3;
  for (double l = 0.0; l < length; l += length / numberOfSamples)
    {
    double t = Modulus(l, itk::Math::pi * 2);

    VesselTubeType::TubePointType point;
    point.SetPosition(pos[0], pos[1], pos[2]);
    pointList.push_back(point);

    if (lastPosValid)
      {
      curveLength += (pos - lastPos).GetNorm();
      }
    else
      {
      lastPosValid = true;
      }
    lastPos = pos;

    pos[0] = l;
    pos[1] = amplitude * sin(frequency * t);
    pos[2] = 0.0;
    }

  VectorType origin(0.0);
  VectorType vec = lastPos - origin;
  length = (lastPos - origin).GetNorm();

  vessel->SetPoints(pointList);
  return vessel;
}

//-----------------------------------------------------------------------------
bool TestVesselMetrics(VesselTubeType::Pointer vessel, double results[3], int i)
{
  // Run the metrics
  FilterType::Pointer filter = FilterType::New();
  filter->SetMeasureFlag(
    FilterType::DISTANCE_METRIC
    | FilterType::INFLECTION_COUNT_METRIC
    | FilterType::SUM_OF_ANGLES_METRIC);
  filter->SetInput(vessel);
  filter->Update();
  if (! filter->GetOutput() )
    {
    std::cerr<<"Error while running the metrics !"<<std::endl;
    return false;
    }

  double dm = filter->GetDistanceMetric();
  double icm = filter->GetInflectionCountMetric();
  double soam = filter->GetSumOfAnglesMetric();

  bool success = true;
  for (int i = 0; i < 3; ++i)
    {
    success &= itk::Math::FloatAlmostEqual(dm, results[0]);
    success &= itk::Math::FloatAlmostEqual(icm, results[1]);
    success &= itk::Math::FloatAlmostEqual(soam, results[2]);

    if (!success)
      {
      std::cerr<<"Error at index #"<<i
        <<", an array doesn't have the expected value !"<<std::endl;
      std::cerr<<"  DM: expected "<<results[0]<<" got "<<dm<<std::endl;
      std::cerr<<"  ICM: expected "<<results[1]<<" got "<<icm<<std::endl;
      std::cerr<<"  SOAM: expected "<<results[2]<<" got "<<soam<<std::endl;
      return success;
      }
    }

  return success;
}

//-----------------------------------------------------------------------------
int itktubeTortuositySpatialObjectFilterTest( int argc, char * argv[] )
{
  //
  // Test straight spatial objects
  //
  const int NUMBER_OF_STRAIGHT_OBJECT_TESTS = 3;
  double start[NUMBER_OF_STRAIGHT_OBJECT_TESTS][3] = {
      {26.2, -9601.0, -42.01},
      {0.0, 0.0, 0.0},
      {3.20, -89.0, 1.0},
    };

  double increment[NUMBER_OF_STRAIGHT_OBJECT_TESTS][3] = {
      {0.05, 80.4, 1.0},
      {0.01, 100.0, 0.0},
      {0.0, 0.0, -0.007},
    };
  double numberOfPoints[NUMBER_OF_STRAIGHT_OBJECT_TESTS] =
    {
    25,
    4,
    10
    };
  double straightObjectResults[NUMBER_OF_STRAIGHT_OBJECT_TESTS][3] = {
    {1.0, 1.0, 0.0},
    {1.0, 1.0, 0.0},
    {1.0, 1.0, 0.0},
    };

  for (int i = 0; i < NUMBER_OF_STRAIGHT_OBJECT_TESTS; ++i)
    {
    //std::cerr<<"Straight object test #"<<i<<std::endl;
    VesselTubeType::Pointer vessel =
      GenerateStraightTube(start[i], increment[i], numberOfPoints[i]);

    if (!TestVesselMetrics(vessel, straightObjectResults[i], i))
      {
      std::cerr<<"Error in straight object test for object #"<<i<<std::endl;
      return EXIT_FAILURE;
      }
    }

  //
  // Test cosinus spatial objects
  //
  const int NUMBER_OF_COSINUS_OBJECT_TESTS = 4;
  double length[NUMBER_OF_COSINUS_OBJECT_TESTS] = {
      itk::Math::pi * 2,
      2.0 * itk::Math::pi * 2,
      itk::Math::pi * 2,
      itk::Math::pi * 2,
    };
  double amplitude[NUMBER_OF_COSINUS_OBJECT_TESTS] = {
      1.0,
      1.0,
      9.0,
      1.0,
    };
  double frequency[NUMBER_OF_COSINUS_OBJECT_TESTS] = {
      1.0,
      1.0,
      1.0,
      5.0,
    };
  double cosResults[NUMBER_OF_COSINUS_OBJECT_TESTS][3] = {
    {
      1.21581,
      1.21581 * 2.0,
      0.411187
    },
    {
      1.21581,
      1.21581 * 4.0,
      0.411187
    },
    {
      5.87042,
      5.87042 * 2.0,
      0.158497
    },
    {
      3.40308,
      3.40308 * 10.0,
      1.28584
    },
    };

  for (int i = 0; i < NUMBER_OF_COSINUS_OBJECT_TESTS; ++i)
    {
    //std::cerr<<"Cos object test #"<<i<<std::endl;
    VesselTubeType::Pointer vessel =
      GenerateCosTube(length[i], amplitude[i], frequency[i]);

    if (!TestVesselMetrics(vessel, cosResults[i], i))
      {
      std::cerr<<"Error in straight object test for object #"<<i<<std::endl;
      return EXIT_FAILURE;
      }
    }


  return EXIT_SUCCESS;
}
