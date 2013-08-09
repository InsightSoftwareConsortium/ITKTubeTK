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

#include "itktubeSubSampleTubeTreeSpatialObjectFilter.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkGroupSpatialObject.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkVesselTubeSpatialObject.h>

#include "SubSampleTubesCLP.h"

template< unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  const unsigned int Dimension = VDimension;

  PARSE_ARGS;

  typedef itk::VesselTubeSpatialObject< Dimension > TubeSpatialObjectType;
  typedef itk::GroupSpatialObject< Dimension >      GroupSpatialObjectType;

  double progress = 0.0;
  itk::TimeProbesCollectorBase timeCollector;

  tube::CLIProgressReporter progressReporter( "SubSampleTubes",
                                        CLPProcessInformation );

  progressReporter.Start();
  progressReporter.Report( progress );

  // Read in the tubes
  timeCollector.Start( "Read tubes" );
  typedef itk::SpatialObjectReader< Dimension >  ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputTubeFile );
  try
    {
    reader->Update();
    }
  catch( const std::exception & error )
    {
    tube::ErrorMessage( error.what() );
    return EXIT_FAILURE;
    }
  typename GroupSpatialObjectType::Pointer groupSpatialObject = reader->GetGroup();
  std::ostringstream ostrm;
  ostrm << "Number of children = "
    << groupSpatialObject->GetNumberOfChildren();
  tube::InformationMessage( ostrm.str() );
  timeCollector.Stop( "Read tubes" );

  progress = 0.3;
  progressReporter.Report( progress );

  timeCollector.Start( "Sub-sample tubes" );
  typedef itk::tube::SubSampleTubeTreeSpatialObjectFilter< GroupSpatialObjectType,
    TubeSpatialObjectType >
      SubSampleTubeTreeFilterType;
  typename SubSampleTubeTreeFilterType::Pointer subSampleTubeTreeFilter =
    SubSampleTubeTreeFilterType::New();
  subSampleTubeTreeFilter->SetInput( reader->GetGroup() );
  subSampleTubeTreeFilter->SetSampling( samplingFactor );
  try
    {
    subSampleTubeTreeFilter->Update();
    }
  catch( const std::exception &e )
    {
    tube::ErrorMessage( e.what() );
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Sub-sample tubes" );

  progress = 0.9;
  progressReporter.Report( progress );

  timeCollector.Start( "Write output");
  typedef itk::SpatialObjectWriter< Dimension > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputTubeFile );
  writer->SetInput( subSampleTubeTreeFilter->GetOutput() );
  try
    {
    writer->Update();
    }
  catch( const std::exception &e )
    {
    tube::ErrorMessage( e.what() );
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Write output" );

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  //PARSE_ARGS;

  // TODO: work on 2D, too
  return DoIt< 3 >( argc, argv );
}
