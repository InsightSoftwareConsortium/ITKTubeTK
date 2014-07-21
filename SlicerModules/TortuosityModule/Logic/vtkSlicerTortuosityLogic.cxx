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

#include "vtkSlicerTortuosityLogic.h"

// ITK includes
#include "itkVesselTubeSpatialObject.h"

// Spatial object includes
#include "vtkMRMLSpatialObjectsNode.h"

// VTK includes
#include "vtkObjectFactory.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"

vtkCxxRevisionMacro(vtkSlicerTortuosityLogic, "$Revision: 1.9.12.1 $");
vtkStandardNewMacro(vtkSlicerTortuosityLogic);

//------------------------------------------------------------------------------
vtkSlicerTortuosityLogic::vtkSlicerTortuosityLogic( void )
{
  this->FlagToArrayNames[DistanceMetric] = "DistanceMetric";
  this->FlagToArrayNames[InflectionCountMetric] = "InflectionCountMetric";
  this->FlagToArrayNames[InflectionPoints] = "InflectionPoints";
  this->FlagToArrayNames[SumOfAnglesMetric] = "SumOfAnglesMetric";
}

//------------------------------------------------------------------------------
vtkSlicerTortuosityLogic::~vtkSlicerTortuosityLogic( void )
{
}

//------------------------------------------------------------------------------
void vtkSlicerTortuosityLogic::PrintSelf(ostream& os, vtkIndent indent)
{
}

//------------------------------------------------------------------------------
bool vtkSlicerTortuosityLogic::UniqueMeasure(int flag)
{
  return flag == vtkSlicerTortuosityLogic::DistanceMetric ||
    flag == vtkSlicerTortuosityLogic::InflectionCountMetric ||
    flag == vtkSlicerTortuosityLogic::InflectionPoints ||
    flag == vtkSlicerTortuosityLogic::SumOfAnglesMetric;
}

//------------------------------------------------------------------------------
vtkDoubleArray* vtkSlicerTortuosityLogic
::GetDistanceMetricArray(vtkMRMLSpatialObjectsNode* node, int flag)
{
  return this->GetArray(node, flag & vtkSlicerTortuosityLogic::DistanceMetric);
}

//------------------------------------------------------------------------------
vtkDoubleArray* vtkSlicerTortuosityLogic
::GetInflectionCountMetricArray(vtkMRMLSpatialObjectsNode* node, int flag)
{
  return this->GetArray(
    node, flag & vtkSlicerTortuosityLogic::InflectionCountMetric);
}

//------------------------------------------------------------------------------
vtkDoubleArray* vtkSlicerTortuosityLogic
::GetInflectionPointsArray(vtkMRMLSpatialObjectsNode* node, int flag)
{
  return this->GetArray(
    node, flag & vtkSlicerTortuosityLogic::InflectionPoints);
}

//------------------------------------------------------------------------------
vtkDoubleArray* vtkSlicerTortuosityLogic
::GetSumOfAnglesMetricArray(vtkMRMLSpatialObjectsNode* node, int flag)
{
  return this->GetArray(
    node, flag & vtkSlicerTortuosityLogic::SumOfAnglesMetric);
}

//------------------------------------------------------------------------------
vtkDoubleArray* vtkSlicerTortuosityLogic
::GetArray(vtkMRMLSpatialObjectsNode* node, int flag)
{
  if (!node || !flag || !this->UniqueMeasure(flag))
    {
    return NULL;
    }

  vtkPolyData* polydata = node->GetPolyData();
  if (!polydata)
    {
    return NULL;
    }
  vtkPointData* pointData = polydata->GetPointData();
  if (!pointData)
    {
    return NULL;
    }

  std::string name = this->FlagToArrayNames[flag];
  vtkDoubleArray* metricArray =
    vtkDoubleArray::SafeDownCast(pointData->GetArray(name.c_str()));
  if (!metricArray)
    {
    vtkDoubleArray* ids =
      vtkDoubleArray::SafeDownCast(pointData->GetArray("TubeIDs"));
    assert(ids);

    vtkNew<vtkDoubleArray> newMetricArray;
    newMetricArray->SetName(name.c_str());
    newMetricArray->SetNumberOfValues(ids->GetSize());
    pointData->AddArray(newMetricArray.GetPointer());
    return newMetricArray.GetPointer();
    }
  return metricArray;
}

//------------------------------------------------------------------------------
bool vtkSlicerTortuosityLogic
::RunDistanceMetric(vtkMRMLSpatialObjectsNode* node)
{
  return this->RunMetrics(
    node, vtkSlicerTortuosityLogic::DistanceMetric);
}

//------------------------------------------------------------------------------
bool vtkSlicerTortuosityLogic
::RunInflectionCountMetric(vtkMRMLSpatialObjectsNode* node)
{
  return this->RunMetrics(
    node, vtkSlicerTortuosityLogic::InflectionCountMetric);
}

//------------------------------------------------------------------------------
bool vtkSlicerTortuosityLogic
::RunInflectionPoints(vtkMRMLSpatialObjectsNode* node)
{
  return this->RunMetrics(
    node, vtkSlicerTortuosityLogic::InflectionPoints);
}

//------------------------------------------------------------------------------
bool vtkSlicerTortuosityLogic
::RunSumOfAnglesMetric(vtkMRMLSpatialObjectsNode* node)
{
  return this->RunMetrics(
    node, vtkSlicerTortuosityLogic::SumOfAnglesMetric);
}

//------------------------------------------------------------------------------
bool vtkSlicerTortuosityLogic
::RunMetrics(vtkMRMLSpatialObjectsNode* node, int flag)
{
  if (!node)
    {
    return false;
    }

  // 1 - Get the metric arrays
  vtkDoubleArray* dm = this->GetDistanceMetricArray(node, flag);
  vtkDoubleArray* icm = this->GetInflectionCountMetricArray(node, flag);
  vtkDoubleArray* ip = this->GetInflectionPointsArray(node, flag);
  vtkDoubleArray* soam = this->GetSumOfAnglesMetricArray(node, flag);

  // 2 - (Re)Fill the metric arrays
  typedef vtkMRMLSpatialObjectsNode::TubeNetType  TubeNetType;
  typedef itk::VesselTubeSpatialObject<3>         VesselTubeType;
  typedef VesselTubeType::TubePointType           VesselTubePointType;
  typedef VesselTubePointType::PointType          PointType;

  TubeNetType* spatialObject = node->GetSpatialObject();

  char childName[] = "Tube";
  TubeNetType::ChildrenListType* tubeList =
    spatialObject->GetChildren(spatialObject->GetMaximumDepth(), childName);

  vtkIdType tubeID = 0;
  size_t totalNumberOfPointsAdded = 0;
  for(TubeNetType::ChildrenListType::iterator tubeIt = tubeList->begin();
        tubeIt != tubeList->end(); ++tubeIt, ++tubeID )
    {
    VesselTubeType* currTube =
      dynamic_cast<VesselTubeType*>((*tubeIt).GetPointer());
    if (!currTube)
      {
      continue;
      }
    assert(currTube->GetNumberOfPoints() >= 2);

    //
    // DM variables
    double start[3], end[3];
    double pathLength = 0.0;

    // ICM variables
    double previousN[3] = {0.0, 0.0, 0.0};
    int inflectionCount = 1;

    size_t numberOfPoints = currTube->GetPoints().size();
    for(size_t index = 0; index < numberOfPoints; ++index)
      {
      VesselTubePointType* point =
        dynamic_cast<VesselTubePointType*>(currTube->GetPoint(index));
      assert(point);

      double currentPoint[3];
      currentPoint[0] = currTube->GetPoint(index)->GetPosition()[0];
      currentPoint[1] = currTube->GetPoint(index)->GetPosition()[1];
      currentPoint[2] = currTube->GetPoint(index)->GetPosition()[2];

      //
      // General variables
      bool nextPointAvailable = (index < numberOfPoints - 1);
      double nextPoint[3] = {0.0, 0.0, 0.0};
      if (nextPointAvailable)
        {
        nextPoint[0] = currTube->GetPoint(index + 1)->GetPosition()[0];
        nextPoint[1] = currTube->GetPoint(index + 1)->GetPosition()[1];
        nextPoint[2] = currTube->GetPoint(index + 1)->GetPosition()[2];
        }
      bool previousPointAvailable = (index > 0);
      double previousPoint[3] = {0.0, 0.0, 0.0};
      if (previousPointAvailable)
        {
        previousPoint[0] = currTube->GetPoint(index - 1)->GetPosition()[0];
        previousPoint[1] = currTube->GetPoint(index - 1)->GetPosition()[1];
        previousPoint[2] = currTube->GetPoint(index - 1)->GetPosition()[2];
        }

      //
      // DM Computations
      if (index == 0)
        {
        start[0] = currentPoint[0];
        start[1] = currentPoint[1];
        start[2] = currentPoint[2];
        }
      if (index == numberOfPoints - 1)
        {
        end[0] = currentPoint[0];
        end[1] = currentPoint[1];
        end[2] = currentPoint[2];
        }

      if ((dm || icm || ip) && nextPointAvailable)
        {
        double currentToNext[3];
        vtkMath::Subtract(nextPoint, currentPoint, currentToNext);
        pathLength += vtkMath::Norm(currentToNext);
        }

      //
      // ICM Computations
      double inflectionValue = 0.0;
      if (previousPointAvailable && nextPointAvailable)
        {
        // Compute velocity and acceleration
        double v[3], t1[3], t2[3], a[3];
        vtkMath::Subtract(nextPoint, previousPoint, v);
        vtkMath::Subtract(currentPoint, previousPoint, t1);
        vtkMath::Subtract(nextPoint, currentPoint, t2);
        vtkMath::Subtract(t2, t1, a);

        // Compute n
        double n[3];
        vtkMath::Cross(v, a, n);
        vtkMath::Cross(v, n, n);
        vtkMath::Normalize(n);

        if (vtkMath::Norm(a) > 1e-6) // make sure acceleration isn't to small
          {
          // Check for inflection
          double deltaN[3];
          vtkMath::Subtract(n, previousN, deltaN);

          inflectionValue = vtkMath::Dot(deltaN, deltaN);
          if (inflectionValue > 1.0)
            {
            inflectionCount += 1;
            }
          }

        previousN[0] = n[0];
        previousN[1] = n[1];
        previousN[2] = n[2];
        }

      // Set the inflection value for this point
      if (ip)
        {
        ip->SetValue(index, inflectionValue);
        }
      }

    //
    // DM final calculation
    double dmResult = -1.0;
    if (dm || icm)
      {
      double startToEnd[3];
      vtkMath::Subtract(start, end, startToEnd);
      double straighLineLength = vtkMath::Norm(startToEnd);
      if (straighLineLength > 0.0)
        {
        dmResult = pathLength / straighLineLength;
        }

      if (dmResult < 1.0)
        {
        std::cerr<<"Error while computing the distance metric."
          <<"DM (="<<dmResult<<") > 1.0"<<std::endl;
        }
      }

    // ICM final calculation
    // Nothing to do

    //
    // Fill the arrays
    for(size_t index = totalNumberOfPointsAdded;
      index < numberOfPoints + totalNumberOfPointsAdded; ++index)
      {
      if (dm)
        {
        dm->SetValue(index, dmResult);
        }
      if (icm)
        {
        icm->SetValue(index, inflectionCount * dmResult);
        }
      }
    totalNumberOfPointsAdded += numberOfPoints;
    }

  return true;
}
