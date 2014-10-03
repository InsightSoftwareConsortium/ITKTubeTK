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

// .NAME vtkSlicerTortuosityLogic -
// Interface with TubeTK Tortuosity methods and provides additionnal support
// .SECTION Description
// The tortuosity logic is in charge of interfacing with the TubeTK Tortuosity
// module. It is also used for convience methods on tube objects.

#ifndef __vtkSlicerTortuosityLogic_h
#define __vtkSlicerTortuosityLogic_h

#include <vtkSlicerModuleLogic.h>
#include <vtkSlicerTortuosityModuleLogicExport.h>

#include <map>
#include <vector>

class vtkDoubleArray;
class vtkMRMLSpatialObjectsNode;

class VTK_SLICER_TORTUOSITY_MODULE_LOGIC_EXPORT vtkSlicerTortuosityLogic
 : public vtkSlicerModuleLogic
{
public:
  static vtkSlicerTortuosityLogic *New( void );
  vtkTypeMacro(vtkSlicerTortuosityLogic,vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Different kind of metrics that can be run on a spatial object node.
  enum MeasureTypes
    {
    DistanceMetric = 0x01,
    InflectionCountMetric = 0x02,
    InflectionPoints = 0x04,
    SumOfAnglesMetric = 0x08,
    All = 0xFF,
    };

  // Return whether the flag asks for an unique measure or not.
  bool UniqueMeasure(int flag);

  // Return the array corresponding to the metric array. If a flag is
  // passed, the array is valid only if the flag uniquely corresponds
  // to a metric.
  // \sa GetArray()
  vtkDoubleArray* GetDistanceMetricArray(
    vtkMRMLSpatialObjectsNode* node, int flag = DistanceMetric);
  vtkDoubleArray* GetInflectionCountMetricArray(
    vtkMRMLSpatialObjectsNode* node, int flag = InflectionCountMetric);
  vtkDoubleArray* GetInflectionPointsArray(
    vtkMRMLSpatialObjectsNode* node, int flag = InflectionPoints);
  vtkDoubleArray* GetSumOfAnglesMetricArray(
    vtkMRMLSpatialObjectsNode* node, int flag = SumOfAnglesMetric);

  // Get the metric array on the given node. If no array corresponding to
  // the flag (or name), an empty array will be created.
  vtkDoubleArray* GetOrCreateArray(vtkMRMLSpatialObjectsNode* node, int flag);

  // Run the metric on the given spatial object node.
  // Running any metric will create a field data array on the node if none
  // already exists. A field data array called NumberOfPoints will also be
  // created automatically to help when exporting.
  // \sa RunMetrics()
  bool RunDistanceMetric(vtkMRMLSpatialObjectsNode* node);
  bool RunInflectionCountMetric(vtkMRMLSpatialObjectsNode* node);
  bool RunInflectionPoints(vtkMRMLSpatialObjectsNode* node);
  bool RunSumOfAnglesMetric(vtkMRMLSpatialObjectsNode* node);

  // Run the metric specified by the flag on the given spatial object node.
  bool RunMetrics(vtkMRMLSpatialObjectsNode* node, int flag);

  // Save the given metrics to CSV. Only value will be saved per vessel.
  // The export MUST find the "NumberOfPoints" array generated during the
  // metric run on the node's polydata point data.
  bool SaveAsCSV(
    vtkMRMLSpatialObjectsNode* node, const char* filename, int flag = All);

protected:
  vtkSlicerTortuosityLogic( void );
  ~vtkSlicerTortuosityLogic( void );
  vtkSlicerTortuosityLogic(const vtkSlicerTortuosityLogic&);
  void operator=(const vtkSlicerTortuosityLogic&);

  // Get names from the given flag
  std::vector<std::string> GetNamesFromFlag(int flag);

  // Get the array of the type T with the given name on the node's polydata
  // pointdata.
  template<typename T>
    T* GetArray(vtkMRMLSpatialObjectsNode* node, const char* name);

  // Same than GetArray() but if no array exists, one will be created.
  template<typename T>
    T* GetOrCreateArray(vtkMRMLSpatialObjectsNode* node, const char* name);

private:
  std::map<int, std::string> FlagToArrayNames;

}; // End class vtkSlicerTortuosityLogic

#endif // End !defined(__vtkSlicerTortuosityLogic_h)
