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

// TubeTK includes
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeMacro.h"
#include "itktubeTortuositySpatialObjectFilter.h"
#include "tubeTubeMath.h"

// ITK includes
#include "itkTimeProbesCollectorBase.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkGroupSpatialObject.h"
#include "itkVesselTubeSpatialObject.h"
#include "metaScene.h"

// std includes
#include <sstream>
#include <string.h>

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
    tube::ErrorMessage("Error: Only 2D and 3D data is currently supported.");
    return EXIT_FAILURE;
    }

  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;

  // Load TRE File
  tubeStandardOutputMacro( << "\n>> Loading TRE File" );

  typedef itk::SpatialObjectReader< VDimension > TubesReaderType;
  typedef itk::GroupSpatialObject< VDimension >  TubeGroupType;

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

  // Prepare tortuosity measure flag
  int measureFlag = 0;

  if( basicMetrics )
    {
    measureFlag = measureFlag
                  | TortuosityFilterType::AVERAGE_RADIUS_METRIC
                  | TortuosityFilterType::CHORD_LENGTH_METRIC
                  | TortuosityFilterType::PATH_LENGTH_METRIC;
    }

  if( oldMetrics )
    {
    measureFlag = measureFlag
                  | TortuosityFilterType::DISTANCE_METRIC
                  | TortuosityFilterType::INFLECTION_COUNT_METRIC
                  | TortuosityFilterType::INFLECTION_POINTS_METRIC
                  | TortuosityFilterType::SUM_OF_ANGLES_METRIC;
    }

  if( curvatureMetrics )
    {
    measureFlag = measureFlag
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
    measureFlag |= TortuosityFilterType::CURVATURE_HISTOGRAM_METRICS;
    }

  // Get specified smoothing method
  ::tube::SmoothTubeFunctionEnum smoothingMethodEnum;

  if( strcmp( smoothingMethod, "SMOOTH_TUBE_USING_INDEX_AVERAGE" ) == 0 )
    {
    smoothingMethodEnum =
      ::tube::SmoothTubeFunctionEnum::SMOOTH_TUBE_USING_INDEX_AVERAGE;
    }
  else if( strcmp( smoothingMethod, "SMOOTH_TUBE_USING_INDEX_GAUSSIAN" ) == 0 )
    {
    smoothingMethodEnum =
      ::tube::SmoothTubeFunctionEnum::SMOOTH_TUBE_USING_INDEX_GAUSSIAN;
    }
  else if( strcmp( smoothingMethod, "SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN" ) == 0 )
    {
    smoothingMethodEnum =
      ::tube::SmoothTubeFunctionEnum::SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN;
    }

  // Run tortuosity filter
  tubeStandardOutputMacro( << "\n>> Computing tortuosity measures" );

  timeCollector.Start( "Computing tortuosity measures" );

  typedef itk::VesselTubeSpatialObject< VDimension >  TubeType;
  typedef typename TubeType::Pointer                  TubePointerType;
  typedef typename TubeGroupType::ChildrenListPointer TubeListPointerType;
  typedef itk::tube::TortuositySpatialObjectFilter< TubeType >
  TortuosityFilterType;

  char childName[] = "Tube";
  TubeListPointerType pTubeList = pTubeGroup->GetChildren(
    pTubeGroup->GetMaximumDepth(), childName );

  for( typename TubeGroupType::ChildrenListType::iterator
       itTubes = pTubeList->begin();
       itTubes != pTubeList->end(); ++itTubes )
    {
    TubeType* pCurTube = dynamic_cast<TubeType*>((*itTubes).GetPointer());
    if (!pCurTube)
      {
      continue;
      }

    typename TortuosityFilterType::Pointer pTortuosityFilter =
    TortuosityFilterType::New();

    pTortuosityFilter->SetInput( pCurTube );
    pTortuosityFilter->SetMeasureFlag( measureFlag );
    pTortuosityFilter->SetSmoothingScale( smoothingScale );
    pTortuosityFilter->SetSmoothingMethod( smoothingMethodEnum );
    pTortuosityFilter->SetNumberOfBins( numberOfHistogramBins );
    pTortuosityFilter->SetHistogramMin( histogramMin );
    pTortuosityFilter->SetHistogramMax( histogramMax );
    pTortuosityFilter->Update();
    }

  timeCollector.Stop( "Computing tortuosity measures" );

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
    case 2:
      return DoIt<2>( argc, argv );
      break;

    case 3:
      return DoIt<3>( argc, argv );
      break;

    default:
      tubeErrorMacro(<< "Error: Only 2D and 3D data is currently supported.");
      return EXIT_FAILURE;
    }
}
