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

// .NAME vtkSlicerTortuosityLogic -
// Interface with TubeTK Tortuosity methods and provides additionnal support
// .SECTION Description
// The tortuosity logic is in charge of interfacing with the TubeTK Tortuosity
// module. It is also used for convience methods on tube objects.

#ifndef __vtkSlicerTortuosityLogic_h
#define __vtkSlicerTortuosityLogic_h

#include <vtkSlicerModuleLogic.h>
#include <vtkSlicerTortuosityModuleLogicExport.h>

// ITK includes
#include "itkVesselTubeSpatialObject.h"
#include "itktubeTortuositySpatialObjectFilter.h"
#include "itkVector.h"

// Spatial object includes
#include "vtkMRMLSpatialObjectsNode.h"

// TubeTK includes
#include "tubeTubeMath.h"

#include <map>
#include <vector>

class vtkDoubleArray;
class vtkMRMLSpatialObjectsNode;

class VTK_SLICER_TORTUOSITY_MODULE_LOGIC_EXPORT vtkSlicerTortuosityLogic
 : public vtkSlicerModuleLogic
{
public:
  static vtkSlicerTortuosityLogic *New( void );
  vtkTypeMacro( vtkSlicerTortuosityLogic,vtkSlicerModuleLogic );
  void PrintSelf( ostream& os, vtkIndent indent );

  // typdefs
  typedef vtkMRMLSpatialObjectsNode::TubeNetType                      TubeNetType;
  typedef itk::VesselTubeSpatialObject< 3 >                           VesselTubeType;
  typedef itk::tube::TortuositySpatialObjectFilter< VesselTubeType >  FilterType;
  typedef VesselTubeType::TubePointType                               VesselTubePointType;

  // Different groups of metrics that can be run on a spatial object node.
  // See FilterType::MeasureType for more info on the metrics.
  // Tau4 is an additional metric that isn't computed by the filter.
  // - Old metrics:
  //                DistanceMetric
  //                InflectionCountMetric
  //                InflectionPointsMetric
  //                SumOfAnglesMetric
  //                SumOfTorsionMetric
  // - Basic Metrics
  //                AverageRadiusMetric
  //                ChordLengthMetric
  //                PathLengthMetric
  // - Curvature Metrics
  //                CurvatureScalarMetric
  //                InflectionCount1Metric
  //                InflectionCount2Metric
  //                Percentile95Metric
  //                Tau4Metric
  //                TotalCurvatureMetric
  //                TotalSquared CurvatureMetric
  // - Histogram Metrics
  //                CurvatureHistogramMetrics
  // - All Metrics
  //                All of the above
  enum MetricGroupsFlag
    {
    BasicMetricsGroup      = 0x01,
    OldMetricsGroup        = 0x02,
    CurvatureMetricsGroup  = 0x04,
    HistogramMetricsGroup  = 0x08,
    AllMetricsGroup        = 0xFF
    };

  // Return an array of pointers to the arrays to compute using the flag
  std::vector< vtkDoubleArray* >
    GetMetricArraysToCompute( vtkMRMLSpatialObjectsNode* node, int flag );

  // Get the metric double array on the given node and name. If no array corresponds to
  // the name, an empty array will be created.
  vtkDoubleArray*
    GetOrCreateDoubleArray( vtkMRMLSpatialObjectsNode* node, const char* name );

  // Run the metric specified by the flag on the given spatial object node.
  // Before running the metrics, smoothing and subsampling is applied to the tube
  //
  // smoothingMethod: enum that specifies which smoothing method to apply
  //
  // smoothingScale: Depending on the smoothingMethod, this has different roles:
  //   smoothingMethod == tube::SMOOTH_TUBE_USING_INDEX_GAUSSIAN:
  //      -> smothingScale is the std deviation of the gaussian
  //   smoothingMethod == tube::SMOOTH_TUBE_USING_INDEX_AVERAGE:
  //      -> smothingScale is the half average-window size
  //   smoothingMethod == tube::SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN:
  //      -> smothingScale is the std deviation of the gaussian
  //
  // subSampling: The subsampling factor.
  //    1 = no subsampling
  //    2 = divide number of points by 2
  //    etc...
  bool RunMetrics( vtkMRMLSpatialObjectsNode* node,
                  int groupFlag,
                  tube::SmoothTubeFunctionEnum smoothingMethod =
                    tube::SMOOTH_TUBE_USING_INDEX_GAUSSIAN,
                  double smoothingScale = 0.0,
                  int subsampling = 1 );

  // Save the given metrics to CSV. Only value will be saved per vessel.
  // The export MUST find the "NumberOfPoints" array generated during the
  // metric run on the node's polydata point data.
  bool SaveAsCSV( vtkMRMLSpatialObjectsNode* node,
                 const char* filename,
                 int groupFlag = AllMetricsGroup );

  // Load a CSV file with two columns : ID and Value, and assign the values
  // to the passed node tubes with corresponding IDs, as a point data.
  // If there is more ( ID, value ) pairs in the file than tubes in the passed
  // node, they will be ignored. If there is less, they will be all assigned,
  // and the tubes that are missing a value will be assigned a default value.
  bool LoadColorsFromCSV( vtkMRMLSpatialObjectsNode* node, const char *filename );


protected:
  vtkSlicerTortuosityLogic( void );
  ~vtkSlicerTortuosityLogic( void );
  vtkSlicerTortuosityLogic( const vtkSlicerTortuosityLogic& );
  void operator=( const vtkSlicerTortuosityLogic& );

  // Get names from the given flag
  std::vector< std::string > GetPrintableNamesFromMetricFlag( int metricFlag );

  // Convert the group flag ( vtkSlicerTortuosityLogic::MetricGroupsFlag ) into a
  // metric flag ( FilterType::MeasureType )
  int GetMetricFlagFromGroupFlag( int groupFlag );

  // Get the array of the type T with the given name on the node's polydata
  // pointdata.
  template<typename T>
    T* GetArray( vtkMRMLSpatialObjectsNode* node, const char* name );

  // Same than GetArray() but if no array exists, one will be created.
  template<typename T>
    T* GetOrCreateArray( vtkMRMLSpatialObjectsNode* node, const char* name );

private:
  std::map<int, std::string> m_MetricFlagToArrayNames;
  std::map<int, int>         m_GroupFlagToMetricFlag;

  std::vector< vtkSmartPointer< vtkIntArray > > m_HistogramArrays;

}; // End class vtkSlicerTortuosityLogic

#endif // End !defined( __vtkSlicerTortuosityLogic_h )
