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
#include <vtkSlicerTortuosityLogic.h>

// MRML includes
#include <vtkMRMLSpatialObjectsNode.h>
#include <vtkMRMLSpatialObjectsStorageNode.h>

// TubeTK includes
#include <itkVesselTubeSpatialObject.h>
#include <itkVesselTubeSpatialObjectPoint.h>

// VTK includes
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkMathUtilities.h>
#include <vtkNew.h>

#include <sstream>

typedef itk::VesselTubeSpatialObject<3> VesselTubeType;
double EPSILON = 1e-3;
double TWOPI = 2.0*vtkMath::Pi();

//-----------------------------------------------------------------------------
VesselTubeType::Pointer GenerateStraightTube(double start[3],
                                             double increment[3],
                                             int numberOfPoints)
{
  VesselTubeType::Pointer vessel = VesselTubeType::New();
  VesselTubeType::PointListType pointList;

  double pos[3];
  pos[0] = start[0];
  pos[1] = start[1];
  pos[2] = start[2];

  for (int i = 0; i < numberOfPoints; ++i)
    {
    VesselTubeType::TubePointType point;
    point.SetPosition(pos[0], pos[1], pos[2]);
    pointList.push_back(point);

    vtkMath::Add(pos, increment, pos);
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

  double pos[3] = {0.0 , 0.0, 0.0};
  double lastPos[3] = {0.0 , 0.0, 0.0};
  double curveLength = 0.0;
  bool lastPosValid = false;

  int numberOfSamples = 1e3;
  for (double l = 0.0; l < length; l += length / numberOfSamples)
    {
    double t = Modulus(l, TWOPI);

    VesselTubeType::TubePointType point;
    point.SetPosition(pos[0], pos[1], pos[2]);
    pointList.push_back(point);

    if (lastPosValid)
      {
      double lastToCurrent[3];
      vtkMath::Subtract(pos, lastPos, lastToCurrent);
      curveLength += vtkMath::Norm(lastToCurrent);
      }
    else
      {
      lastPosValid = true;
      }
    lastPos[0] = pos[0];
    lastPos[1] = pos[1];
    lastPos[2] = pos[2];

    pos[0] = l;
    pos[1] = amplitude * sin(frequency * t);
    pos[2] = 0.0;
    }

  double origin[3] = {0.0, 0.0, 0.0};
  double vec[3];
  vtkMath::Subtract(lastPos, origin, vec);
  length = vtkMath::Norm(vec);

  vessel->SetPoints(pointList);
  return vessel;
}

//-----------------------------------------------------------------------------
bool TestVesselMetrics(VesselTubeType::Pointer vessel, double results[3], int i)
{
  // 1- Create the spatial object node
  typedef vtkMRMLSpatialObjectsNode::TubeNetType GroupType;
  vtkNew<vtkMRMLSpatialObjectsNode> node;

  GroupType::Pointer group = GroupType::New();
  group->AddSpatialObject(vessel);
  node->SetSpatialObject(group);

  // DEBUG: save the node to file
  //vtkNew<vtkMRMLSpatialObjectsStorageNode> storageNode;
  //std::stringstream ss;
  //ss << "W:/Prometheus/Data/SegmentTubesTest/test"<<i<<".tre";
  //storageNode->SetFileName(ss.str().c_str());
  //assert(storageNode->WriteData(node.GetPointer()));

  // Run the metrics
  vtkNew<vtkSlicerTortuosityLogic> logic;
  int flag = vtkSlicerTortuosityLogic::DistanceMetric
    | vtkSlicerTortuosityLogic::InflectionCountMetric
    | vtkSlicerTortuosityLogic::SumOfAnglesMetric;
  if (! logic->RunMetrics(node.GetPointer(), flag) )
    {
    std::cerr<<"Error while running the metrics !"<<std::endl;
    return false;
    }

  vtkDoubleArray* dm = logic->GetDistanceMetricArray(node.GetPointer());
  vtkDoubleArray* icm = logic->GetInflectionCountMetricArray(node.GetPointer());
  vtkDoubleArray* soam = logic->GetSumOfAnglesMetricArray(node.GetPointer());
  if (!dm || !icm || !soam)
    {
    std::cerr<<"Cannot find all the metric arrays !"<<std::endl;
    return false;
    }

  if (dm->GetNumberOfTuples() != vessel->GetNumberOfPoints()
    || dm->GetNumberOfTuples() != icm->GetNumberOfTuples()
    || dm->GetNumberOfTuples() != soam->GetNumberOfTuples())
    {
    std::cerr<<"All the arrays don't have the same size !"<<std::endl;
    return false;
    }

  bool success = true;
  for (int i = 0; i < dm->GetNumberOfTuples(); ++i)
    {
    success &= vtkMathUtilities::FuzzyCompare(dm->GetValue(i), results[0], EPSILON);
    success &= vtkMathUtilities::FuzzyCompare(icm->GetValue(i), results[1], EPSILON);
    success &= vtkMathUtilities::FuzzyCompare(soam->GetValue(i), results[2], EPSILON);

    if (!success)
      {
      std::cerr<<"Error at index #"<<i
        <<", an array doesn't have the expected value !"<<std::endl;
      std::cerr<<"  DM: expected "<<results[0]
        <<" got "<<dm->GetValue(i)<<std::endl;
      std::cerr<<"  ICM: expected "<<results[1]
        <<" got "<<icm->GetValue(i)<<std::endl;
      std::cerr<<"  SOAM: expected "<<results[2]
        <<" got "<<soam->GetValue(i)<<std::endl;
      return success;
      }
    }

  return success;
}

//-----------------------------------------------------------------------------
int TestTortuosityMetrics( int argc, char * argv[] )
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
      TWOPI,
      2.0 * TWOPI,
      TWOPI,
      TWOPI,
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
