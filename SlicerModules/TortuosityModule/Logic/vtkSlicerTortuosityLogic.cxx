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

#include "vtkSlicerTortuosityLogic.h"

// TubeTK includes
#include "tubeTubeMath.h"

// VTK includes
#include "vtkObjectFactory.h"
#include "vtkDelimitedTextWriter.h"
#include "vtkDelimitedTextReader.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTable.h"

#include <math.h>

vtkStandardNewMacro( vtkSlicerTortuosityLogic );

//------------------------------------------------------------------------------
vtkSlicerTortuosityLogic::vtkSlicerTortuosityLogic( void )
{
  // Vessel-wise metrics
  this->m_MetricFlagToArrayNames
      [ FilterType::AVERAGE_RADIUS_METRIC ] = "AverageRadiusMetric";
  this->m_MetricFlagToArrayNames
      [ FilterType::CHORD_LENGTH_METRIC ] = "ChordLengthMetric";
  this->m_MetricFlagToArrayNames
      [ FilterType::DISTANCE_METRIC ] = "DistanceMetric";
  this->m_MetricFlagToArrayNames
      [ FilterType::INFLECTION_COUNT_METRIC ] = "InflectionCountMetric";
  this->m_MetricFlagToArrayNames
      [ FilterType::INFLECTION_COUNT_1_METRIC ] = "InflectionCount1Metric";
  this->m_MetricFlagToArrayNames
      [ FilterType::INFLECTION_COUNT_2_METRIC ] = "InflectionCount2Metric";
  this->m_MetricFlagToArrayNames
      [ FilterType::PATH_LENGTH_METRIC ] = "PathLengthMetric";
  this->m_MetricFlagToArrayNames
      [ FilterType::PERCENTILE_95_METRIC ] = "Percentile95Metric";
  this->m_MetricFlagToArrayNames
      [ FilterType::SUM_OF_ANGLES_METRIC ] = "SumOfAnglesMetric";
  this->m_MetricFlagToArrayNames
      [ FilterType::SUM_OF_TORSION_METRIC ] = "SumOfTorsionMetric";
  this->m_MetricFlagToArrayNames
      [ FilterType::TOTAL_CURVATURE_METRIC ] = "TotalCurvatureMetric";
  this->m_MetricFlagToArrayNames
      [ FilterType::TOTAL_SQUARED_CURVATURE_METRIC ] = "TotalSquaredCurvatureMetric";

  // Point-wise metrics
  this->m_MetricFlagToArrayNames
      [ FilterType::CURVATURE_SCALAR_METRIC ] = "CurvatureScalarMetric";
//  this->m_MetricFlagToArrayNames
//      [ FilterType::CURVATURE_VECTOR_METRIC ] = "CurvatureVectorMetric";
  this->m_MetricFlagToArrayNames
      [ FilterType::INFLECTION_POINTS_METRIC ] = "InflectionPointsMetric";

  // Fill map of GroupFlag to metricFlag ( TubeTK flag )
  this->m_GroupFlagToMetricFlag[ BasicMetricsGroup ] =
                                    FilterType::AVERAGE_RADIUS_METRIC
                                  | FilterType::CHORD_LENGTH_METRIC
                                  | FilterType::PATH_LENGTH_METRIC;

  this->m_GroupFlagToMetricFlag[ OldMetricsGroup ] =
                                    FilterType::DISTANCE_METRIC
                                  | FilterType::INFLECTION_COUNT_METRIC
                                  | FilterType::INFLECTION_POINTS_METRIC
                                  | FilterType::SUM_OF_ANGLES_METRIC
                                  | FilterType::SUM_OF_TORSION_METRIC;

  this->m_GroupFlagToMetricFlag[ CurvatureMetricsGroup ] =
                                    FilterType::INFLECTION_COUNT_1_METRIC
                                  | FilterType::INFLECTION_COUNT_2_METRIC
                                  | FilterType::PERCENTILE_95_METRIC
                                  | FilterType::TOTAL_CURVATURE_METRIC
                                  | FilterType::TOTAL_SQUARED_CURVATURE_METRIC
                                  | FilterType::CURVATURE_SCALAR_METRIC
                                  | FilterType::CURVATURE_VECTOR_METRIC;

  this->m_GroupFlagToMetricFlag[ HistogramMetricsGroup ] =
                                    FilterType::CURVATURE_HISTOGRAM_METRICS;

}

//------------------------------------------------------------------------------
vtkSlicerTortuosityLogic::~vtkSlicerTortuosityLogic( void )
{
}

//------------------------------------------------------------------------------
void vtkSlicerTortuosityLogic::PrintSelf( ostream& os, vtkIndent indent )
{
  Superclass::PrintSelf( os,indent );
}

//------------------------------------------------------------------------------
std::vector< vtkDoubleArray* > vtkSlicerTortuosityLogic
::GetMetricArraysToCompute( vtkMRMLSpatialObjectsNode* node, int metricFlag )
{
  std::vector< vtkDoubleArray* > metricArraysVector;

  // Iterate through all the metric flags
  for ( long int compareFlag = 0x01;
    compareFlag <= FilterType::BITMASK_VESSEL_WISE_METRICS;
    compareFlag = compareFlag << 1 )
    {

    // If single metric flag is on, create array and add it to the vector
    if( ( metricFlag & compareFlag ) > 0 )
      {
      std::string metricName = m_MetricFlagToArrayNames[ compareFlag ];
      if( metricName != "" )
        {
        vtkDoubleArray* metricArray =
          this->GetOrCreateDoubleArray( node, metricName.c_str() );
        if( metricArray )
          {
          metricArraysVector.push_back( metricArray );
          }
        }
      }
    }

  return metricArraysVector;
}

//------------------------------------------------------------------------------
template<typename T> T* vtkSlicerTortuosityLogic
::GetArray( vtkMRMLSpatialObjectsNode* node, const char* name )
{
  vtkPolyData* polydata = node->GetPolyData();
  if ( !polydata )
    {
    return NULL;
    }
  vtkPointData* pointData = polydata->GetPointData();
  if ( !pointData )
    {
    return NULL;
    }

  return T::SafeDownCast( pointData->GetArray( name ) );
}

//------------------------------------------------------------------------------
template<typename T> T* vtkSlicerTortuosityLogic
::GetOrCreateArray( vtkMRMLSpatialObjectsNode* node, const char* name )
{
  T* metricArray = this->GetArray<T>( node, name );
  if ( !metricArray )
    {
    vtkPolyData* polydata = node->GetPolyData();
    if ( !polydata )
      {
      return NULL;
      }
    vtkPointData* pointData = polydata->GetPointData();
    if ( !pointData )
      {
      return NULL;
      }

    vtkNew<T> newMetricArray;
    newMetricArray->SetName( name );
    pointData->AddArray( newMetricArray.GetPointer() );
    return newMetricArray.GetPointer();
    }
  return metricArray;
}

//------------------------------------------------------------------------------
vtkDoubleArray* vtkSlicerTortuosityLogic
::GetOrCreateDoubleArray( vtkMRMLSpatialObjectsNode* node, const char* name )
{
  vtkDoubleArray* metricArray =
    this->GetOrCreateArray< vtkDoubleArray >( node, name );
  if ( !metricArray )
    {
    std::cerr<<"The array "<<name<<" couldn't be created or fetched"
             <<std::endl;
    return NULL;
    }

  // If it's new, make it the correct size
  vtkDoubleArray* ids = this->GetArray< vtkDoubleArray >( node, "TubeIDs" );
  assert( ids );
  if ( metricArray->GetSize() != ids->GetSize() )
    {
    metricArray->Initialize();
    metricArray->SetNumberOfValues( ids->GetSize() );
    }
  return metricArray;
}

//------------------------------------------------------------------------------
int vtkSlicerTortuosityLogic
::GetMetricFlagFromGroupFlag( int groupFlag )
{
  // Iterate through all the group flags
  int metricFlag = 0x00;
  for ( long int compareFlag = 0x01;
    compareFlag <= AllMetricsGroup;
    compareFlag = compareFlag << 1 )
    {
    // If group metric is on, add the corresponding metric flags
    if( ( groupFlag & compareFlag ) > 0 )
      {
      metricFlag |= m_GroupFlagToMetricFlag[ compareFlag ];
      }
    }
  return metricFlag;
}

//------------------------------------------------------------------------------
bool vtkSlicerTortuosityLogic
::RunMetrics( vtkMRMLSpatialObjectsNode* node, int groupFlag,
             tube::SmoothTubeFunctionEnum smoothingMethod, double smoothingScale,
             int subsampling )
{
  if ( !node )
    {
    return false;
    }

  if( groupFlag == 0 )
    {
    std::cout<<"No metrics to compute"<<std::endl;
    return true;
    }

  // Convert group flag to metric flag
  int metricFlag = GetMetricFlagFromGroupFlag( groupFlag );
  TubeNetType* spatialObject = node->GetSpatialObject();

  char childName[] = "Tube";
  TubeNetType::ChildrenListType* tubeList =
    spatialObject->GetChildren( spatialObject->GetMaximumDepth(), childName );

  // 1 - Get the metric arrays vector
  std::vector< vtkDoubleArray* > metricsVector;
  metricsVector = this->GetMetricArraysToCompute( node, metricFlag );

  // Additionnal metrics that don't use the tortuosity filter
  // Tau4 Metric is part of the CurvatureMetrics group
  if( ( groupFlag & CurvatureMetricsGroup ) > 0 )
    {
    // Tau4 need the path length metric, so we add it manually
    // it will be computed in the filter, but won't be added as an array
    // or printed in the CSV file.
    metricFlag |= FilterType::PATH_LENGTH_METRIC;
    vtkDoubleArray * tau4Array = this->GetOrCreateDoubleArray( node, "Tau4Metric" );
    if( tau4Array )
      {
      metricsVector.push_back( tau4Array );
      }
    }
  // Additional metrics that don't use the tortuosity filter
  // Volume Metric is part of the BasicMetrics Group
  if( ( groupFlag & BasicMetricsGroup ) > 0 )
    {
    vtkDoubleArray * volumeArray = this->GetOrCreateDoubleArray( node, "VolumeMetric" );
    if( volumeArray )
      {
      metricsVector.push_back( volumeArray );
      }
    }

  // Histogram Metric
  unsigned int numberOfBins = 20;
  double histMin = 0.0;
  double histMax = 1.0;
  double histStep = ( histMax - histMin ) / numberOfBins;
  if( ( metricFlag & FilterType::CURVATURE_HISTOGRAM_METRICS ) > 0 )
    {
    // Create the data arrays for the histogram
    for( unsigned int i = 0; i < numberOfBins; i++ )
      {
      vtkSmartPointer< vtkIntArray > histArray;
      std::ostringstream oss;
      oss << "Hist-Bin#"
          << i << ": "
          << i*histStep <<" - "
          << ( i+1 )*histStep;
      std::string binArrayName =oss.str();

      if( m_HistogramArrays.size() < numberOfBins )
        {
        histArray = vtkSmartPointer< vtkIntArray >::New();
        histArray->SetName( binArrayName.c_str() );
        histArray->SetNumberOfValues( tubeList->size() );
        m_HistogramArrays.push_back( histArray );
        }
      else
        {
        histArray = m_HistogramArrays[i];
        histArray->Initialize();
        histArray->SetName( binArrayName.c_str() );
        histArray->SetNumberOfValues( tubeList->size() );
        }
      }
    }

  // Rewrite number of points array everytime
  vtkIntArray* nop = this->GetOrCreateArray< vtkIntArray >( node, "NumberOfPoints" );
  nop->Initialize();

  vtkIntArray* ids = this->GetOrCreateArray< vtkIntArray >( node, "TubeID" );
  ids->Initialize();

  // 2 - Fill the metric arrays
  int tubeNumber = 0;
  int totalNumberOfPointsAdded = 0;
  for( TubeNetType::ChildrenListType::iterator tubeIt = tubeList->begin();
        tubeIt != tubeList->end(); ++tubeIt )
    {
    VesselTubeType::Pointer currTube =
      dynamic_cast< VesselTubeType* >( ( *tubeIt ).GetPointer() );
    if ( !currTube )
      {
      continue;
      }

    if ( currTube->GetNumberOfPoints() < 2 )
      {
      std::cerr<<"Error, vessel #"<<currTube->GetId()
        <<" has less than 2 points !"<<std::endl
        <<"Skipping the vessel."<<std::endl;
      continue;
      }

    // Set filter parameters and update it
    FilterType::Pointer filter = FilterType::New();
    filter->SetMeasureFlag( metricFlag );
    filter->SetSmoothingMethod( smoothingMethod );
    filter->SetSmoothingScale( smoothingScale );
    filter->SetSubsamplingScale( subsampling );
    filter->SetNumberOfBins( numberOfBins );
    filter->SetHistogramMin( histMin );
    filter->SetHistogramMax( histMax );
    filter->SetInput( currTube );
    filter->Update();

    if ( filter->GetOutput()->GetId() != currTube->GetId() )
      {
      std::cerr<<"Error while running filter on tube."<<std::endl;
      return false;
      }

    //Update tube to get the preprocessed tube.
    currTube = filter->GetOutput();

    // Fill the arrays

    // tubeIndex: index of points in the whole tube data set ( same across all tubes )
    // filterIndex: index of points in a specific tube ( reset to 0 when changing tube )
    int numberOfPoints = currTube->GetPoints().size();
    for( int filterIndex = 0, tubeIndex = totalNumberOfPointsAdded;
      filterIndex < numberOfPoints; ++filterIndex, ++tubeIndex )
      {
      for( unsigned int i = 0; i < metricsVector.size(); i++ )
        {
        if( metricsVector[i] != NULL )
          {
          std::string arrayName = metricsVector[i]->GetName();

          // Printable metrics
          if( arrayName == "AverageRadiusMetric" )
            {
            metricsVector[i]->SetValue( tubeIndex, filter->GetAverageRadiusMetric() );
            }
          if( arrayName == "ChordLengthMetric" )
            {
            metricsVector[i]->SetValue( tubeIndex, filter->GetChordLengthMetric() );
            }
          if( arrayName == "DistanceMetric" )
            {
            metricsVector[i]->SetValue( tubeIndex, filter->GetDistanceMetric() );
            }
          if( arrayName == "InflectionCountMetric" )
            {
            metricsVector[i]->SetValue( tubeIndex, filter->GetInflectionCountMetric() );
            }
          if( arrayName == "InflectionCount1Metric" )
            {
            metricsVector[i]->SetValue( tubeIndex, filter->GetInflectionCount1Metric() );
            }
          if( arrayName == "InflectionCount2Metric" )
            {
            metricsVector[i]->SetValue( tubeIndex, filter->GetInflectionCount2Metric() );
            }
          if( arrayName == "PathLengthMetric" )
            {
            metricsVector[i]->SetValue( tubeIndex, filter->GetPathLengthMetric() );
            }
          if( arrayName == "Percentile95Metric" )
            {
            metricsVector[i]->SetValue( tubeIndex, filter->GetPercentile95Metric() );
            }
          if( arrayName == "SumOfAnglesMetric" )
            {
            metricsVector[i]->SetValue( tubeIndex, filter->GetSumOfAnglesMetric() );
            }
          if( arrayName == "SumOfTorsionMetric" )
            {
            metricsVector[i]->SetValue( tubeIndex, filter->GetSumOfTorsionMetric() );
            }
          if( arrayName == "Tau4Metric" )
            {
            metricsVector[i]->SetValue( tubeIndex,
              filter->GetTotalCurvatureMetric() / filter->GetPathLengthMetric() );
            }
          if( arrayName == "TotalCurvatureMetric" )
            {
            metricsVector[i]->SetValue( tubeIndex, filter->GetTotalCurvatureMetric() );
            }
          if( arrayName == "TotalSquaredCurvatureMetric" )
            {
            metricsVector[i]->SetValue( tubeIndex,
              filter->GetTotalSquaredCurvatureMetric() );
            }

          // Not printable metrics
          if( arrayName == "CurvatureScalarMetric" )
            {
            metricsVector[i]->SetValue( tubeIndex,
              filter->GetCurvatureScalarMetric( filterIndex ) );
            }
          if( arrayName == "InflectionPointsMetric" )
            {
            metricsVector[i]->SetValue( tubeIndex,
              filter->GetInflectionPointValue( filterIndex ) );
            }
          if( arrayName == "VolumeMetric" )
            {
            if( filterIndex == 0 )
              {
              metricsVector[i]->SetValue( tubeIndex, 0 );
              break;
              }
            VesselTubePointType* tubePointCurrent = dynamic_cast< VesselTubePointType* >
              ( currTube->GetPoint( filterIndex ) );
            VesselTubePointType* tubePointBefore = dynamic_cast< VesselTubePointType* >
              ( currTube->GetPoint( filterIndex - 1 ) );
            double height = tubePointCurrent->GetPosition().EuclideanDistanceTo
              ( tubePointBefore->GetPosition() );
            double radius =
              ( tubePointCurrent->GetRadius() + tubePointBefore->GetRadius() ) / 2;
            double volume = itk::Math::pi * radius * radius * height;
            metricsVector[i]->SetValue( tubeIndex, volume );
            double totalVolume =
              metricsVector[i]->GetValue( tubeIndex - filterIndex ) + volume;
            metricsVector[i]->SetValue( tubeIndex - filterIndex, totalVolume );
            }
          }
        }
      nop->InsertNextValue( numberOfPoints );
      ids->InsertNextValue( currTube->GetId() );
      }

    // Set Histogram metrics values
    if( ( metricFlag & FilterType::CURVATURE_HISTOGRAM_METRICS ) > 0 )
      {
      // Get the histogram features
      for( unsigned int i = 0; i < numberOfBins; i++ )
        {
        m_HistogramArrays[i]->SetValue( tubeNumber,
          filter->GetCurvatureHistogramMetric( i ) );
        }
      }

    tubeNumber++;
    totalNumberOfPointsAdded += numberOfPoints;
    }

  // Update the arrays to recompute the range

  for( unsigned int i = 0; i < metricsVector.size(); i++ )
    {
    metricsVector[i]->Modified();
    }
  if ( ( metricFlag & FilterType::CURVATURE_HISTOGRAM_METRICS ) > 0 )
    {
    for ( unsigned int i = 0; i < numberOfBins; i++ )
      {
      m_HistogramArrays[i]->Modified();
      }
    }

  return true;
}

//------------------------------------------------------------------------------
std::vector<std::string>
vtkSlicerTortuosityLogic::GetPrintableNamesFromMetricFlag( int metricFlag )
{
  std::vector< std::string > names;
  names.push_back( "TubeID" );
  names.push_back( "NumberOfPoints" );
  for ( long int compareFlag = 0x01;
    compareFlag <= FilterType::BITMASK_VESSEL_WISE_METRICS;
    compareFlag = compareFlag << 1 )
    {

    // If metric is asked to print and is printable
    if( ( metricFlag &
        compareFlag &
        FilterType::BITMASK_VESSEL_WISE_METRICS ) > 0 )
      {
      names.push_back( this->m_MetricFlagToArrayNames[ compareFlag ] );
      }
    }
  return names;
}

//------------------------------------------------------------------------------
bool vtkSlicerTortuosityLogic
::SaveAsCSV( vtkMRMLSpatialObjectsNode* node, const char* filename, int groupFlag )
{
  if ( !node || !filename )
    {
    return false;
    }

  // Convert group flag to metric flag
  int metricFlag = GetMetricFlagFromGroupFlag( groupFlag );

  // Get the metric arrays
  std::vector< vtkDataArray* > metricArrays;
  std::vector< std::string > names = this->GetPrintableNamesFromMetricFlag( metricFlag );
  // Additionnal metrics that don't use the tortuosity filter
  if( ( groupFlag & CurvatureMetricsGroup ) > 0 )
    {
    names.push_back( "Tau4Metric" );
    }
  if( ( groupFlag & BasicMetricsGroup ) > 0 )
    {
    names.push_back( "VolumeMetric" );
    }
  for ( std::vector<std::string>::iterator it = names.begin();
    it != names.end(); ++it )
    {
    vtkDataArray* metricArray =
      this->GetArray< vtkDataArray >( node, it->c_str() );
    if ( metricArray )
      {
      metricArrays.push_back( metricArray );
      }
    }

  // Make sure we have everything we need for export
  if ( metricArrays.size() <= 0 )
    {
    std::cout<<"No array found for given flag: "<<metricFlag<<std::endl;
    return false;
    }

  vtkIntArray* numberOfPointsArray =
    this->GetArray< vtkIntArray >( node, "NumberOfPoints" );
  if ( !numberOfPointsArray )
    {
    std::cerr<<"Expected ''NumberOfPoints'' array on the node point data."
      <<std::endl<<"Cannot proceed."<<std::endl;
    return false;
    }

  // Create  the table. Each column has only one value per vessel
  // instead of one value per each point of the vessel.
  vtkNew< vtkTable > table;
  for( std::vector< vtkDataArray* >::iterator it = metricArrays.begin();
    it != metricArrays.end(); ++it )
    {
    vtkNew< vtkDoubleArray > newArray;
    newArray->SetName( ( *it )->GetName() );

    for ( int j = 0; j < numberOfPointsArray->GetNumberOfTuples();
      j += numberOfPointsArray->GetValue( j ) )
      {
      newArray->InsertNextTuple( ( *it )->GetTuple( j ) );
      }

    table->AddColumn( newArray.GetPointer() );
    }


  // Add the histogram features to the table
  if( ( metricFlag & FilterType::CURVATURE_HISTOGRAM_METRICS ) > 0 )
    {
    for( std::vector< vtkSmartPointer< vtkIntArray > >::iterator it =
      m_HistogramArrays.begin();
      it != m_HistogramArrays.end(); ++it )
      {
      table->AddColumn( ( *it ) );
      }
    }

  // Write out the table to file
  vtkNew< vtkDelimitedTextWriter > writer;
  writer->SetFileName( filename );

#if ( VTK_MAJOR_VERSION < 6 )
  writer->SetInput( table.GetPointer() );
#else
  writer->SetInputData( table.GetPointer() );
#endif

  return writer->Write();
}

//------------------------------------------------------------------------------
bool vtkSlicerTortuosityLogic::LoadColorsFromCSV( 
  vtkMRMLSpatialObjectsNode *node, const char* filename )
{
  if ( !node || !filename )
    {
    return false;
    }

  typedef vtkMRMLSpatialObjectsNode::TubeNetType  TubeNetType;
  typedef itk::VesselTubeSpatialObject<3>         VesselTubeType;

  // Load the table from file
  vtkNew< vtkDelimitedTextReader > reader;
  reader->SetFileName( filename );
  reader->SetFieldDelimiterCharacters( "," );
  reader->SetHaveHeaders( true );
  reader->SetDetectNumericColumns( true );
  reader->Update();
  vtkTable* colorTable = reader->GetOutput();
  if ( !colorTable )
    {
    std::cerr<<"Error in reading CSV file"<<std::endl;
    return false;
    }

  // Check if table is valid
  if ( colorTable->GetNumberOfColumns() != 2 )
    {
    std::cerr<<"Expected 2 columns in CSV file."
      <<std::endl<<"Cannot proceed."<<std::endl;
    return false;
    }

  // Get the tube list of the spatial object
  TubeNetType* spatialObject = node->GetSpatialObject();
  char childName[] = "Tube";
  TubeNetType::ChildrenListType* tubeList =
    spatialObject->GetChildren( spatialObject->GetMaximumDepth(), childName );

  // Create a new data array in the node
  double defaultValue = 0.0;
  vtkDoubleArray* customColorScaleArray =
    this->GetOrCreateArray< vtkDoubleArray >( node, "CustomColorScale" );

  // Set the size of the array
  vtkDoubleArray* ids = this->GetArray< vtkDoubleArray >( node, "TubeIDs" );
  assert( ids );
  if ( customColorScaleArray->GetNumberOfTuples() != ids->GetNumberOfTuples() )
    {
    customColorScaleArray->SetNumberOfTuples( ids->GetNumberOfTuples() );
    }

  // Initialize the array with the default value
  customColorScaleArray->FillComponent( 0, defaultValue );

  // Iterate through tubeList
  size_t totalNumberOfPoints = 0;
  TubeNetType::ChildrenListType::iterator tubeIt;
  for ( tubeIt = tubeList->begin(); tubeIt != tubeList->end(); ++tubeIt )
    {
    VesselTubeType* currTube =
      dynamic_cast< VesselTubeType* >( ( *tubeIt ).GetPointer() );
    if ( !currTube )
      {
      continue;
      }
    if ( currTube->GetNumberOfPoints() < 2 )
      {
      std::cerr<<"Error, vessel #"<<currTube->GetId()
        <<" has less than 2 points !"<<std::endl
        <<"Skipping the vessel."<<std::endl;
      continue;
      }

    // Get the current tube ID
    int tubeId = currTube->GetId();
    vtkDebugMacro( <<"Tube ID "<<tubeId );

    // Look for the ID in the table and get the corresponding value
    double valueToAssign = 0.0; //Default value for not specified tubes
    int tubeIndex = -1;
    for ( int i = 0; i < colorTable->GetNumberOfRows(); i++ )
      {
      if ( colorTable->GetValue( i, 0 ).ToInt() == tubeId )
        {
        tubeIndex = i;
        valueToAssign = colorTable->GetValue( tubeIndex, 1 ).ToDouble();
        vtkDebugMacro( <<" found in CSV : value = "<<valueToAssign );
        break;
        }
      }
    if ( tubeIndex == -1 )
      {
      vtkDebugMacro( <<" not found in the CSV file" );
      totalNumberOfPoints += currTube->GetPoints().size();
      continue;
      }

    // Fill the array of that tube
    size_t numberOfPoints = currTube->GetPoints().size();
    for ( size_t j = totalNumberOfPoints; j < totalNumberOfPoints + numberOfPoints; j++ )
      {
      customColorScaleArray->SetValue( j, valueToAssign );
      }
    totalNumberOfPoints += numberOfPoints;
    }

  // Notify the array of the changes
  customColorScaleArray->Modified();

  return true;
}
