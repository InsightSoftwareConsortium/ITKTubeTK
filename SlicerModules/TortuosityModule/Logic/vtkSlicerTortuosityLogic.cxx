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
#include "itktubeTortuositySpatialObjectFilter.h"

// Spatial object includes
#include "vtkMRMLSpatialObjectsNode.h"

// VTK includes
#include "vtkObjectFactory.h"
#include "vtkDelimitedTextWriter.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTable.h"

#include <math.h>

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
  return this->GetOrCreateArray(
    node, flag & vtkSlicerTortuosityLogic::DistanceMetric);
}

//------------------------------------------------------------------------------
vtkDoubleArray* vtkSlicerTortuosityLogic
::GetInflectionCountMetricArray(vtkMRMLSpatialObjectsNode* node, int flag)
{
  return this->GetOrCreateArray(
    node, flag & vtkSlicerTortuosityLogic::InflectionCountMetric);
}

//------------------------------------------------------------------------------
vtkDoubleArray* vtkSlicerTortuosityLogic
::GetInflectionPointsArray(vtkMRMLSpatialObjectsNode* node, int flag)
{
  return this->GetOrCreateArray(
    node, flag & vtkSlicerTortuosityLogic::InflectionPoints);
}

//------------------------------------------------------------------------------
vtkDoubleArray* vtkSlicerTortuosityLogic
::GetSumOfAnglesMetricArray(vtkMRMLSpatialObjectsNode* node, int flag)
{
  return this->GetOrCreateArray(
    node, flag & vtkSlicerTortuosityLogic::SumOfAnglesMetric);
}

//------------------------------------------------------------------------------
vtkDoubleArray* vtkSlicerTortuosityLogic
::GetOrCreateArray(vtkMRMLSpatialObjectsNode* node, int flag)
{
  if (!flag || !this->UniqueMeasure(flag))
    {
    return NULL;
    }

  std::string name = this->FlagToArrayNames[flag];
  vtkDoubleArray* metricArray =
    this->GetOrCreateArray<vtkDoubleArray>(node, name.c_str());
  if (!metricArray)
    {
    return NULL;
    }

  // If it's new, make it the correct size
  vtkDoubleArray* ids = this->GetArray<vtkDoubleArray>(node, "TubeIDs");
  assert(ids);
  if (metricArray->GetSize() != ids->GetSize());
    {
    metricArray->SetNumberOfValues(ids->GetSize());
    }
  return metricArray;
}

//------------------------------------------------------------------------------
template<typename T> T* vtkSlicerTortuosityLogic
::GetArray(vtkMRMLSpatialObjectsNode* node, const char* name)
{
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

  return T::SafeDownCast(pointData->GetArray(name));
}

//------------------------------------------------------------------------------
template<typename T> T* vtkSlicerTortuosityLogic
::GetOrCreateArray(vtkMRMLSpatialObjectsNode* node, const char* name)
{
  T* metricArray = this->GetArray<T>(node, name);
  if (!metricArray)
    {
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

    vtkNew<T> newMetricArray;
    newMetricArray->SetName(name);
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
  vtkIntArray* nop =
    this->GetOrCreateArray<vtkIntArray>(node, "NumberOfPoints");
  assert(nop);

  // 2 - Fill the metric arrays
  typedef vtkMRMLSpatialObjectsNode::TubeNetType                    TubeNetType;
  typedef itk::VesselTubeSpatialObject<3>                           VesselTubeType;
  typedef itk::tube::TortuositySpatialObjectFilter<VesselTubeType>  FilterType;

  TubeNetType* spatialObject = node->GetSpatialObject();

  char childName[] = "Tube";
  TubeNetType::ChildrenListType* tubeList =
    spatialObject->GetChildren(spatialObject->GetMaximumDepth(), childName);

  size_t totalNumberOfPointsAdded = 0;
  for(TubeNetType::ChildrenListType::iterator tubeIt = tubeList->begin();
        tubeIt != tubeList->end(); ++tubeIt)
    {
    VesselTubeType* currTube =
      dynamic_cast<VesselTubeType*>((*tubeIt).GetPointer());
    if (!currTube)
      {
      continue;
      }
    assert(currTube->GetNumberOfPoints() >= 2);

    // Update filter
    FilterType::Pointer filter = FilterType::New();
    filter->SetMeasureFlag(flag);
    filter->SetInput(currTube);
    filter->Update();

    if (filter->GetOutput()->GetId() != currTube->GetId())
      {
      std::cerr<<"Error while running filter on tube."<<std::endl;
      return false;
      }

    // Fill the arrays
    size_t numberOfPoints = currTube->GetPoints().size();
    for(size_t index = totalNumberOfPointsAdded;
      index < numberOfPoints + totalNumberOfPointsAdded; ++index)
      {
      if (dm)
        {
        dm->SetValue(index, filter->GetDistanceMetric());
        }
      if (icm)
        {
        icm->SetValue(index, filter->GetInflectionCountMetric());
        }
      if (soam)
        {
        soam->SetValue(index, filter->GetSumOfAnglesMetric());
        }
      if (ip)
        {
        ip->SetValue(index, filter->GetInflectionPointValue(index));
        }
      }
    nop->InsertNextValue(numberOfPoints);

    totalNumberOfPointsAdded += numberOfPoints;
    }

  return true;
}

//------------------------------------------------------------------------------
std::vector<std::string> vtkSlicerTortuosityLogic::GetNamesFromFlag(int flag)
{
  std::vector<std::string> names;
  for (int compareFlag = vtkSlicerTortuosityLogic::DistanceMetric;
    compareFlag <= vtkSlicerTortuosityLogic::SumOfAnglesMetric;
    compareFlag = compareFlag << 1)
    {
    if (flag & compareFlag)
      {
      names.push_back(this->FlagToArrayNames[compareFlag]);
      }
    }
  return names;
}

//------------------------------------------------------------------------------
bool vtkSlicerTortuosityLogic
::SaveAsCSV(vtkMRMLSpatialObjectsNode* node, const char* filename, int flag)
{
  if (!node || !filename)
    {
    return false;
    }

  // Get the metric arrays
  std::vector<vtkDoubleArray*> metricArrays;
  std::vector<std::string> names = this->GetNamesFromFlag(flag);
  for (std::vector<std::string>::iterator it = names.begin();
    it != names.end(); ++it)
    {
    vtkDoubleArray* metricArray =
      this->GetArray<vtkDoubleArray>(node, it->c_str());
    if (metricArray)
      {
      metricArrays.push_back(metricArray);
      }
    }

  if (metricArrays.size() <= 0)
    {
    return false;
    }

  // Create  the table. Each column has only one value per vessel
  // instead of one value per each point of the vessel.
  vtkNew<vtkTable> table;
  vtkIntArray* numberOfPointsArray =
    this->GetArray<vtkIntArray>(node, "NumberOfPoints");
  for(std::vector<vtkDoubleArray*>::iterator it = metricArrays.begin();
    it != metricArrays.end(); ++it)
    {
    vtkNew<vtkDoubleArray> newArray;
    newArray->SetName((*it)->GetName());
    newArray->SetNumberOfValues(numberOfPointsArray->GetNumberOfTuples());

    int index = 0;
    for (int j = 0; j < numberOfPointsArray->GetNumberOfTuples(); ++j)
      {
      newArray->SetValue(j, (*it)->GetValue(index));

      std::cout<<(*it)->GetName()<<": "<<index<<" = "<<(*it)->GetValue(index)<<std::endl;

      index += numberOfPointsArray->GetValue(j);
      }

    table->AddColumn(newArray.GetPointer());
    }

  // Write out the table to file
  vtkNew<vtkDelimitedTextWriter> writer;
  writer->SetFileName(filename);

#if (VTK_MAJOR_VERSION < 6)
  writer->SetInput(table.GetPointer());
#else
  writer->SetInputData(table.GetPointer());
#endif

  return writer->Write();
}
