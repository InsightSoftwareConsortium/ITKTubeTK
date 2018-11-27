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

// TubeTK INCLUDES
#include "tubeMacro.h"
#include "tubeMessage.h"
#include "../CLI/tubeCLIProgressReporter.h"
#include "itktubeTortuositySpatialObjectFilter.h"
#include "tubeTubeMath.h"

// ITK INCLUDES
#include "itkTimeProbesCollectorBase.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkGroupSpatialObject.h"
#include "itkVesselTubeSpatialObject.h"
#include "metaScene.h"

// VTK INCLUDES
#include "vtkNew.h"
#include "vtkTable.h"
#include "vtkVersion.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkSmartPointer.h"
#include "vtkDelimitedTextWriter.h"

// std includes
#include <sstream>
#include <string.h>
#include <vector>
#include <map>
#include <string>

#include "ComputeTubeTortuosityMeasuresCLP.h"

template< unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // Ensure that the input image dimension is valid
  // We only support 2D and 3D Images due to the
  // limitation of itkTubeSpatialObject
  if( VDimension != 2 && VDimension != 3 )
    {
    tube::ErrorMessage( "Error: Only 2D and 3D data is currently supported." );
    return EXIT_FAILURE;
    }

  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;

  // typedefs
  typedef itk::VesselTubeSpatialObject< VDimension >  TubeType;
  typedef itk::SpatialObjectReader< VDimension >      TubesReaderType;
  typedef itk::GroupSpatialObject< VDimension >       TubeGroupType;
  typedef typename TubeGroupType::ChildrenListPointer TubeListPointerType;

  typedef itk::tube::TortuositySpatialObjectFilter< TubeType >
    TortuosityFilterType;

  // Load TRE File
  tubeStandardOutputMacro( << "\n>> Loading TRE File" );

  timeCollector.Start( "Loading Input TRE File" );

  typename TubesReaderType::Pointer tubeFileReader = TubesReaderType::New();

  try
    {
    tubeFileReader->SetFileName( inputTREFile.c_str() );
    tubeFileReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error loading TRE File: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  typename TubeGroupType::Pointer pTubeGroup = tubeFileReader->GetGroup();

  timeCollector.Stop( "Loading Input TRE File" );

  // Get specified smoothing method
  tube::SmoothTubeFunctionEnum smoothingMethodEnum;

  if( smoothingMethod == "SMOOTH_TUBE_USING_INDEX_AVERAGE" )
    {
    smoothingMethodEnum =
      tube::SMOOTH_TUBE_USING_INDEX_AVERAGE;
    }
  else if( smoothingMethod == "SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN" )
    {
    smoothingMethodEnum =
      tube::SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN;
    }
  else
    {
    smoothingMethodEnum =
      tube::SMOOTH_TUBE_USING_INDEX_GAUSSIAN;
    }

  // Prepare tortuosity measure flag
  int metricFlag = 0;

  if( basicMetrics )
    {
    metricFlag = metricFlag
                 | TortuosityFilterType::AVERAGE_RADIUS_METRIC
                 | TortuosityFilterType::CHORD_LENGTH_METRIC
                 | TortuosityFilterType::PATH_LENGTH_METRIC;
    }

  if( oldMetrics )
    {
    metricFlag = metricFlag
                 | TortuosityFilterType::DISTANCE_METRIC
                 | TortuosityFilterType::INFLECTION_COUNT_METRIC
                 | TortuosityFilterType::INFLECTION_POINTS_METRIC
                 | TortuosityFilterType::SUM_OF_ANGLES_METRIC;
    }

  if( curvatureMetrics )
    {
    metricFlag = metricFlag
                 | TortuosityFilterType::INFLECTION_COUNT_1_METRIC
                 | TortuosityFilterType::INFLECTION_COUNT_2_METRIC
                 | TortuosityFilterType::PERCENTILE_95_METRIC
                 | TortuosityFilterType::TOTAL_CURVATURE_METRIC
                 | TortuosityFilterType::TOTAL_SQUARED_CURVATURE_METRIC
                 | TortuosityFilterType::CURVATURE_SCALAR_METRIC
                 | TortuosityFilterType::CURVATURE_VECTOR_METRIC;
    }

  if( histogramMetrics )
    {
    metricFlag |= TortuosityFilterType::CURVATURE_HISTOGRAM_METRICS;
    }

  std::map<int, std::string> MetricFlagToNameMap;
  MetricFlagToNameMap
    [TortuosityFilterType::AVERAGE_RADIUS_METRIC] = "AverageRadiusMetric";
  MetricFlagToNameMap
    [TortuosityFilterType::CHORD_LENGTH_METRIC] = "ChordLengthMetric";
  MetricFlagToNameMap
    [TortuosityFilterType::DISTANCE_METRIC] = "DistanceMetric";
  MetricFlagToNameMap
    [TortuosityFilterType::INFLECTION_COUNT_METRIC] = "InflectionCountMetric";
  MetricFlagToNameMap
    [TortuosityFilterType::INFLECTION_COUNT_1_METRIC]
      = "InflectionCount1Metric";
  MetricFlagToNameMap
    [TortuosityFilterType::INFLECTION_COUNT_2_METRIC]
      = "InflectionCount2Metric";
  MetricFlagToNameMap
    [TortuosityFilterType::PATH_LENGTH_METRIC] = "PathLengthMetric";
  MetricFlagToNameMap
    [TortuosityFilterType::PERCENTILE_95_METRIC] = "Percentile95Metric";
  MetricFlagToNameMap
    [TortuosityFilterType::SUM_OF_ANGLES_METRIC] = "SumOfAnglesMetric";
  MetricFlagToNameMap
    [TortuosityFilterType::TOTAL_CURVATURE_METRIC] = "TotalCurvatureMetric";
  MetricFlagToNameMap
    [TortuosityFilterType::TOTAL_SQUARED_CURVATURE_METRIC]
      = "TotalSquaredCurvatureMetric";

  // Run tortuosity filter
  tubeStandardOutputMacro( << "\n>> Computing tortuosity measures" );

  timeCollector.Start( "Computing tortuosity measures" );

  char childName[] = "Tube";
  TubeListPointerType tubeList = pTubeGroup->GetChildren(
    pTubeGroup->GetMaximumDepth(), childName );

  vtkSmartPointer< vtkIntArray > tubeIdArray =
    vtkSmartPointer<vtkIntArray>::New();
  tubeIdArray->Initialize();
  tubeIdArray->SetName( "TubeIDs" );
  tubeIdArray->SetNumberOfValues( tubeList->size() );

  vtkSmartPointer< vtkIntArray > numPointsArray =
    vtkSmartPointer<vtkIntArray>::New();
  numPointsArray->Initialize();
  numPointsArray->SetName( "NumberOfPoints" );
  numPointsArray->SetNumberOfValues( tubeList->size() );

  std::vector< vtkSmartPointer< vtkDoubleArray >  > metricArrayVec;
  for( int compareFlag = 0x01; compareFlag <=
    static_cast< int >( TortuosityFilterType::BITMASK_ALL_METRICS );
    compareFlag = compareFlag << 1 )
    {
    // If metric is asked to print and is printable
    if( ( metricFlag & compareFlag &
        TortuosityFilterType::BITMASK_VESSEL_WISE_METRICS ) > 0 )
      {
      vtkSmartPointer< vtkDoubleArray > metricArray =
        vtkSmartPointer< vtkDoubleArray >::New();
      metricArray->Initialize();
      metricArray->SetName( MetricFlagToNameMap[compareFlag].c_str() );
      metricArray->SetNumberOfValues( tubeList->size() );
      metricArrayVec.push_back( metricArray );
      }
    }

  if( curvatureMetrics )
    {
    vtkSmartPointer< vtkDoubleArray > tau4Array =
      vtkSmartPointer< vtkDoubleArray >::New();
    tau4Array->Initialize();
    tau4Array->SetName( "Tau4Metric" );
    tau4Array->SetNumberOfValues( tubeList->size() );
    metricArrayVec.push_back( tau4Array );
    }

  std::vector< vtkSmartPointer<vtkIntArray> > histogramArrays;
  if( ( metricFlag & TortuosityFilterType::CURVATURE_HISTOGRAM_METRICS ) > 0 )
    {
    double histStep = ( histogramMax - histogramMin ) / numberOfHistogramBins;

    for( int i = 0; i < numberOfHistogramBins; i++ )
      {
      std::ostringstream oss;
      oss << "Hist-Bin#" << i << ": " << i*histStep <<" - "
        << ( i+1 ) * histStep;
      std::string binArrayName = oss.str();

      vtkSmartPointer< vtkIntArray > histArray =
        vtkSmartPointer< vtkIntArray >::New();
      histArray->Initialize();
      histArray->SetName( binArrayName.c_str() );
      histArray->SetNumberOfValues( tubeList->size() );
      histogramArrays.push_back( histArray );
      }
    }

  int tubeIndex = 0;
  for( typename TubeGroupType::ChildrenListType::iterator
    itTubes = tubeList->begin(); itTubes != tubeList->end(); ++itTubes )
    {
    TubeType* curTube = dynamic_cast<TubeType*>( ( *itTubes ).GetPointer() );
    if( !curTube )
      {
      continue;
      }

    typename TortuosityFilterType::Pointer tortuosityFilter =
      TortuosityFilterType::New();
    tortuosityFilter->SetMeasureFlag( metricFlag );
    tortuosityFilter->SetSmoothingScale( smoothingScale );
    tortuosityFilter->SetSmoothingMethod( smoothingMethodEnum );
    tortuosityFilter->SetNumberOfBins( numberOfHistogramBins );
    tortuosityFilter->SetHistogramMin( histogramMin );
    tortuosityFilter->SetHistogramMax( histogramMax );
    tortuosityFilter->SetInput( curTube );
    tortuosityFilter->Update();

    tubeIdArray->SetValue( tubeIndex, tubeIndex );
    numPointsArray->SetValue( tubeIndex, curTube->GetNumberOfPoints() );

    std::cout << "vess = " << curTube->GetId() << std::endl;

    for( unsigned int i = 0; i < metricArrayVec.size(); i++ )
      {
      std::string metricName = metricArrayVec[i]->GetName();

      if( metricName == "AverageRadiusMetric" )
        {
        metricArrayVec[i]->SetValue( tubeIndex,
          tortuosityFilter->GetAverageRadiusMetric() );
        }
      if( metricName == "ChordLengthMetric" )
        {
        metricArrayVec[i]->SetValue( tubeIndex,
          tortuosityFilter->GetChordLengthMetric() );
        }
      if( metricName == "DistanceMetric" )
        {
        metricArrayVec[i]->SetValue( tubeIndex,
          tortuosityFilter->GetDistanceMetric() );
        }
      if( metricName == "InflectionCountMetric" )
        {
        metricArrayVec[i]->SetValue( tubeIndex,
          tortuosityFilter->GetInflectionCountMetric() );
        }
      if( metricName == "InflectionCount1Metric" )
        {
        metricArrayVec[i]->SetValue( tubeIndex,
          tortuosityFilter->GetInflectionCount1Metric() );
        }
      if( metricName == "InflectionCount2Metric" )
        {
        metricArrayVec[i]->SetValue( tubeIndex,
          tortuosityFilter->GetInflectionCount2Metric() );
        }
      if( metricName == "PathLengthMetric" )
        {
        metricArrayVec[i]->SetValue( tubeIndex,
          tortuosityFilter->GetPathLengthMetric() );
        }
      if( metricName == "Percentile95Metric" )
        {
        metricArrayVec[i]->SetValue( tubeIndex,
          tortuosityFilter->GetPercentile95Metric() );
        }
      if( metricName == "SumOfAnglesMetric" )
        {
        metricArrayVec[i]->SetValue( tubeIndex,
          tortuosityFilter->GetSumOfAnglesMetric() );
        }
      if( metricName == "Tau4Metric" )
        {
        metricArrayVec[i]->SetValue( tubeIndex,
          tortuosityFilter->GetTotalCurvatureMetric()
            / tortuosityFilter->GetPathLengthMetric() );
        }
      if( metricName == "TotalCurvatureMetric" )
        {
        metricArrayVec[i]->SetValue( tubeIndex,
          tortuosityFilter->GetTotalCurvatureMetric() );
        }
      if( metricName == "TotalSquaredCurvatureMetric" )
        {
        metricArrayVec[i]->SetValue( tubeIndex,
          tortuosityFilter->GetTotalSquaredCurvatureMetric() );
        }
      }

    if( ( metricFlag & TortuosityFilterType::CURVATURE_HISTOGRAM_METRICS )
      > 0 )
      {
      // Get the histogram features
      for( int i = 0; i < numberOfHistogramBins; i++ )
        {
        histogramArrays[i]->SetValue( tubeIndex,
          tortuosityFilter->GetCurvatureHistogramMetric( i ) );
        }
      }

    tubeIndex++;
    }

  timeCollector.Stop( "Computing tortuosity measures" );

  // Write tortuosity measures to a CSV file
  tubeStandardOutputMacro( << "\n>> Writing tortuosity measures to CSV" );

  timeCollector.Start( "Writing tortuosity measures to CSV" );

  vtkNew< vtkTable > table;

  table->AddColumn( tubeIdArray.GetPointer() );
  table->AddColumn( numPointsArray.GetPointer() );

  for( unsigned int i = 0; i < metricArrayVec.size(); i++ )
    {
    table->AddColumn( metricArrayVec[i].GetPointer() );
    }

  if( ( metricFlag & TortuosityFilterType::CURVATURE_HISTOGRAM_METRICS )
    > 0 )
    {
    // Get the histogram features
    for( int i = 0; i < numberOfHistogramBins; i++ )
      {
      table->AddColumn( histogramArrays[i].GetPointer() );
      }
    }

  vtkNew< vtkDelimitedTextWriter > writer;
  writer->SetFileName( outputCSVFile.c_str() );
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput( table.GetPointer() );
#else
  writer->SetInputData( table.GetPointer() );
#endif
  writer->Write();

  timeCollector.Stop( "Writing tortuosity measures to CSV" );

  // All done
  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  try
    {
    PARSE_ARGS;
    }
  catch( const std::exception & err )
    {
    tube::ErrorMessage( err.what() );
    return EXIT_FAILURE;
    }
  PARSE_ARGS;

  MetaScene *mScene = new MetaScene;
  mScene->Read( inputTREFile.c_str() );
  if( mScene->GetObjectList()->empty() )
    {
    tubeWarningMacro( << "Input TRE file has no spatial objects" );
    return EXIT_SUCCESS;
    }

  switch( mScene->GetObjectList()->front()->NDims() )
    {
    case 3:
      {
      bool result = DoIt<3>( argc, argv );
      delete mScene;
      return result;
      break;
      }
    default:
      {
      tubeErrorMacro( << "Error: Only 3D data is currently supported." );
      delete mScene;
      return EXIT_FAILURE;
      }
    }

  delete mScene;
  return EXIT_FAILURE;
}
