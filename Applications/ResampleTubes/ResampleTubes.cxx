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
#include "itktubeTubeToTubeTransformFilter.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkVesselTubeSpatialObject.h>

#include <itkImageFileReader.h>
#include <itkDisplacementFieldTransform.h>
#include <itkTransformFileReader.h>
#include <itktubeResampleTubesFilter.h>

#include "ResampleTubesCLP.h"

template< unsigned int Dimension >
void WriteOutput( typename itk::GroupSpatialObject<Dimension>::Pointer
  tubesGroup, const char * fileName )
{
  typedef itk::SpatialObjectWriter< Dimension > SpatialObjectWriterType;

  typename SpatialObjectWriterType::Pointer writer =
    SpatialObjectWriterType::New();
  writer->SetInput( tubesGroup );
  writer->SetFileName( fileName );
  writer->Update();
}

template< unsigned int Dimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef itk::tube::ResampleTubesFilter< Dimension > FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  typedef FilterType::TubeGroupType   GroupSpatialObjectType;

  double progress = 0.0;
  itk::TimeProbesCollectorBase timeCollector;

  tube::CLIProgressReporter progressReporter( "TubeTransform",
                                        CLPProcessInformation );
  progressReporter.Start();
  progressReporter.Report( progress );

  // Read in the tubes
  timeCollector.Start( "Read tubes" );
  typename GroupSpatialObjectType::Pointer tubesGroup =
    GroupSpatialObjectType::New();

  typedef itk::SpatialObjectReader< Dimension > SpatialObjectReaderType;
  typename SpatialObjectReaderType::Pointer reader =
    SpatialObjectReaderType::New();
  reader->SetFileName( inputTubeFile.c_str() );
  reader->Update();

  tubesGroup = reader->GetGroup();
  if( tubesGroup.IsNotNull() )
    {
    tubesGroup->ComputeObjectToWorldTransform();
    filter->SetInput( tubesGroup );
    filter->SetInputSpatialObject( tubesGroup );
    }
  else
    {
    std::cerr << "Cannot read tubes from file: " << inputTubeFile
      << std::endl;
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Read tubes" );

  progress = 0.3;
  progressReporter.Report( progress );

  timeCollector.Start( "Read Parameters" );

  if( !matchImage.empty() )
    {
    typedef itk::Image< char, Dimension >    ImageType;
    typedef itk::ImageFileReader< ImageType> ImageReaderType;

    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName( matchImage.c_str() );
    reader->Update();
    filter->SetMatchImage( reader->GetOutput() );
    }

  if( !loadDisplacementField.empty() )
      {
      typedef typename FilterType::DisplacementFieldType DisplacementFieldType;
      typedef typename itk::ImageFileReader< DisplacementFieldType >
        DisplacementFieldReaderType;

      // Read the displacement field
      typename DisplacementFieldReaderType::Pointer dfReader =
        DisplacementFieldReaderType::New();
      dfReader->SetFileName( loadDisplacementField.c_str() );
      dfReader->Update();
      filter->SetDisplacementField( dfReader->GetOutput() );
      }

  if( !loadTransform.empty() )
    {
     // Read transform from file
    itk::TransformFileReader::Pointer reader =
      itk::TransformFileReader::New();
    reader->SetFileName( loadTransform.c_str() );
    reader->Update();
    filter->SetReadTransformList( reader->GetTransformList() );
    filter->SetUseInverseTransform( useInverseTransform );
    }
  filter->SetSamplingFactor( samplingFactor );

  timeCollector.Stop( "Read Parameters" );
  progress = 0.4;
  progressReporter.Report( progress );

  timeCollector.Start( "Run Filter");
  filter->Update();
  timeCollector.Stop( "Run Filter" );

  progress = 0.9;
  progressReporter.Report( progress );

  timeCollector.Start( "Write output");
  WriteOutput< Dimension >( filter->GetOutput(), outputTubeFile.c_str() );
  timeCollector.Stop( "Write output" );

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

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
  mScene->Read( inputTubeFile.c_str() );

  if( mScene->GetObjectList()->empty() )
    {
    tubeWarningMacro( << "Input TRE file has no spatial objects" );
    delete mScene;
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
      tubeErrorMacro(
        << "Error: Only 3D data is currently supported." );
      delete mScene;
      return EXIT_FAILURE;
      break;
      }
    }
  return EXIT_FAILURE;
}
