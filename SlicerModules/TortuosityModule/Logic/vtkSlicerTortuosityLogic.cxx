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

#include <math.h>

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
namespace
{
template<typename PointType> void CopyPoints(PointType* point, double* p)
{
  for (int i = 0; i < 3; ++i)
    {
    p[i] = static_cast<double>(point->GetPosition()[i]);
    }
}

void CopyVector3(double* from, double* to)
{
  for (int i = 0; i < 3; ++i)
    {
    to[i] = from[i];
    }
}

void InitVector3(double* v)
{
  for (int i = 0; i < 3; ++i)
    {
    v[i] = 0.0;
    }
}

double SafeAcos(double x)
{
  if (x < -1.0)
    {
    x = -1.0;
    }
  else if (x > 1.0)
    {
    x = 1.0;
    }
  return acos(x);
}

} // end namespace

//------------------------------------------------------------------------------
bool vtkSlicerTortuosityLogic
::RunMetrics(vtkMRMLSpatialObjectsNode* node, int flag)
{
  if (!node)
    {
    return false;
    }

  bool noProblem = true;

  // 1 - Get the metric arrays
  vtkDoubleArray* dm = this->GetDistanceMetricArray(node, flag);
  vtkDoubleArray* icm = this->GetInflectionCountMetricArray(node, flag);
  vtkDoubleArray* ip = this->GetInflectionPointsArray(node, flag);
  vtkDoubleArray* soam = this->GetSumOfAnglesMetricArray(node, flag);

  // 2 - (Re)Fill the metric arrays
  typedef vtkMRMLSpatialObjectsNode::TubeNetType  TubeNetType;
  typedef itk::VesselTubeSpatialObject<3>         VesselTubeType;
  typedef VesselTubeType::TubePointType           VesselTubePointType;
  typedef itk::SpatialObjectPoint<3>              PointType;

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
    double T[3], N[3], B[3]; // for the Frenet frame
    bool vectorBIsValid = false;
    int inflectionCount = 1;

    InitVector3(T);
    InitVector3(N);
    InitVector3(B);

    // SOAM variables
    double totalCurvature = 0.0;

    size_t numberOfPoints = currTube->GetPoints().size();
    for(size_t index = 0; index < numberOfPoints; ++index)
      {
      VesselTubePointType* point =
        dynamic_cast<VesselTubePointType*>(currTube->GetPoint(index));
      assert(point);

      double currentPoint[3];
      CopyPoints<PointType>(currTube->GetPoint(index), currentPoint);

      //
      // General variables
      bool nextPointAvailable = (index < numberOfPoints - 1);
      double nextPoint[3] = {0.0, 0.0, 0.0};
      if (nextPointAvailable)
        {
        CopyPoints<PointType>(currTube->GetPoint(index + 1), nextPoint);
        }
      bool previousPointAvailable = (index > 0);
      double previousPoint[3] = {0.0, 0.0, 0.0};
      if (previousPointAvailable)
        {
        CopyPoints<PointType>(currTube->GetPoint(index - 1), previousPoint);
        }
      // t1 and t2, used both in icm and soam
      double t1[3] = {0.0, 0.0, 0.0};
      double t2[3] = {0.0, 0.0, 0.0};
      if (previousPointAvailable && nextPointAvailable)
        {
        vtkMath::Subtract(currentPoint, previousPoint, t1);
        vtkMath::Subtract(nextPoint, currentPoint, t2);
        }

      bool nPlus2PointAvailable = (index < numberOfPoints - 2);
      double nPlus2Point[3] = {0.0, 0.0, 0.0};
      if (nPlus2PointAvailable)
        {
        CopyPoints<PointType>(currTube->GetPoint(index + 2), nPlus2Point);
        }

      //
      // DM Computations
      if (index == 0)
        {
        CopyVector3(currentPoint, start);
        }
      if (index == numberOfPoints - 1)
        {
        CopyVector3(currentPoint, end);
        }

      if ((dm || icm || ip || soam) && nextPointAvailable)
        {
        double currentToNext[3];
        vtkMath::Subtract(nextPoint, currentPoint, currentToNext);
        pathLength += vtkMath::Norm(currentToNext);
        }

      //
      // ICM Computations
      double inflectionValue = 0.0;
      if ((icm || ip) && previousPointAvailable && nextPointAvailable)
        {
        // Compute velocity and acceleration
        double v[3], a[3];
        vtkMath::Subtract(nextPoint, previousPoint, v);
        vtkMath::Subtract(t2, t1, a);

        // Compute the Frenet frame
        // 1 - T = v / |v|
        CopyVector3(v, T);
        vtkMath::Normalize(T);

        // 2 - N = v × a × v / |v × a × v|
        bool canCheckForinflection = vtkMath::Norm(a) > 1e-6;
        if (canCheckForinflection)
          {
          vtkMath::Cross(v, a, N);
          vtkMath::Cross(v, N, N);
          vtkMath::Normalize(N);
          vectorBIsValid = true;
          }
        else if (vectorBIsValid) // 2nd chance
          {
          // Acceleration can be null when the curve approximates a straight
          // line (sin around pi for example). Unfortunately that could happen
          // when the curve is crossing the straight line and the inflection
          // would be missed...
          // This assumes that no pure torsion along the N vector happened.
          // Note that this is only valid is B was already computed at least
          // once.
          vtkMath::Cross(B, T, N);
          vtkMath::Normalize(N);
          canCheckForinflection = true;
          }

        // 3 - B = T x N (in case of null acceleration. See above)
        vtkMath::Cross(T, N, B);
        vtkMath::Normalize(B);

        if (canCheckForinflection)
          {
          // Check for inflection
          double deltaN[3];
          vtkMath::Subtract(N, previousN, deltaN);

          inflectionValue = vtkMath::Dot(deltaN, deltaN);
          if (inflectionValue > 1.0 + 1e-6)
            {
            inflectionCount += 1;
            }
          }

        CopyVector3(N, previousN);
        }

      // Set the inflection value for this point
      if (ip)
        {
        ip->SetValue(index, inflectionValue);
        }

      //
      // SOAM Computations
      if (soam &&
        previousPointAvailable && nextPointAvailable && nPlus2PointAvailable)
        {
        double t3[3];
        vtkMath::Subtract(nPlus2Point, nextPoint, t3);

        // Compute in-plane angle
        double normT1[3], normT2[3];
        CopyVector3(t1, normT1);
        vtkMath::Normalize(normT1);
        CopyVector3(t2, normT2);
        vtkMath::Normalize(normT2);

        double inPlaneAngle = SafeAcos( vtkMath::Dot(normT1, normT2) );

        // Compute torsionnal angle
        double t1t2Cross[3];
        vtkMath::Cross(t1, t2, t1t2Cross);
        double t2t3Cross[3];
        vtkMath::Cross(t2, t3, t2t3Cross);

        double torsionAngle = 0.0;
        double t1t2Norm = vtkMath::Normalize(t1t2Cross);
        double t2t3Norm = vtkMath::Normalize(t2t3Cross);
        if (t1t2Norm > 1e-6 || t2t3Norm > 1e-6)
          {
          // We need to make sure we don't artificially scale those vectors
          // back to life. Otherwise we end up with a torsion angles where
          // there should not be one
          torsionAngle = SafeAcos( vtkMath::Dot(t1t2Cross, t2t3Cross) );
          }
        // Finally get the curvature
        totalCurvature +=
          sqrt(inPlaneAngle*inPlaneAngle + torsionAngle*torsionAngle);
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
        noProblem = false;
        }
      }

    // ICM final calculation
    double icmResult = inflectionCount * dmResult;

    // SOAM final calculation
    double soamResult = 0.0;
    if (soam && pathLength > 0.0)
      {
      soamResult = totalCurvature / pathLength;
      }
    else if (soam && pathLength <= 0.0)
      {
      std::cerr<<"Cannot compute SOAM, total tube path (="
        <<pathLength<<") <= 0.0"<<std::endl;
      noProblem = false;
      }

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
        icm->SetValue(index, icmResult);
        }
      if (soam)
        {
        soam->SetValue(index, soamResult);
        }
      }
    totalNumberOfPointsAdded += numberOfPoints;
    }

  return noProblem;
}
