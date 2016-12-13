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

#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include "itkTimeProbesCollectorBase.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMinimumMaximumImageFilter.h"

#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkGroupSpatialObject.h"
#include "itkTubeSpatialObject.h"
#include "itkTubeSpatialObjectPoint.h"
#include "itkVesselTubeSpatialObject.h"
#include "itkVesselTubeSpatialObjectPoint.h"
#include "tubeMacro.h"
#include "metaScene.h"
#include "tubeMacro.h"

#include "itktubeMinimumSpanningTreeVesselConnectivityFilter.h"

#include <sstream>

#include "ConvertTubesToTubeTreeCLP.h"

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

  // Run vessel connecitivity filter
  tubeStandardOutputMacro( << "\n>> Running vessel connectivity filter" );

  timeCollector.Start( "Running vessel connectivity filter" );

  typedef itk::tube::MinimumSpanningTreeVesselConnectivityFilter< VDimension >
  VesselConnectivityFilterType;
  typename VesselConnectivityFilterType::Pointer vesselConnectivityFilter =
  VesselConnectivityFilterType::New();

  vesselConnectivityFilter->SetInput( pTubeGroup );
  vesselConnectivityFilter->SetMaxTubeDistanceToRadiusRatio( 
    maxTubeDistanceToRadiusRatio );
  vesselConnectivityFilter->SetMaxContinuityAngleError( 
    maxContinuityAngleError );
  vesselConnectivityFilter->SetRemoveOrphanTubes( removeOrphanTubes );

  if( !rootTubeIdList.empty() )
    {
    typename VesselConnectivityFilterType::TubeIdListType
      IdList( rootTubeIdList.begin(), rootTubeIdList.end() );
    vesselConnectivityFilter->SetRootTubeIdList( IdList );
    }

  vesselConnectivityFilter->Update();

  timeCollector.Stop( "Running vessel connectivity filter" );

  // Write tube group with connectivity information
  tubeStandardOutputMacro( 
    << "\n>> Writing TRE file with connectivity information" );

  timeCollector.Start( "Writing TRE file with connectivity information" );

  typedef itk::SpatialObjectWriter< VDimension > TubeWriterType;
  typename TubeWriterType::Pointer tubeWriter = TubeWriterType::New();

  try
    {
    tubeWriter->SetFileName( outputTREFile.c_str() );
    tubeWriter->SetInput( vesselConnectivityFilter->GetOutput() );
    tubeWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error writing TRE file: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Writing TRE file with connectivity information" );

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
      {
      bool result = DoIt<2>( argc, argv );
      delete mScene;
      return result;
      }
    case 3:
      {
      bool result = DoIt<3>( argc, argv );
      delete mScene;
      return result;
      }
    default:
      {
      tubeErrorMacro( 
        << "Error: Only 2D and 3D data is currently supported." );
      delete mScene;
      return EXIT_FAILURE;
      }
    }
}
