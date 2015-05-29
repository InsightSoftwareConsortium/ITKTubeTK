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

#include <itkGroupSpatialObject.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkVesselTubeSpatialObject.h>

#include <itkImageFileReader.h>
#include <itkDisplacementFieldTransform.h>
#include <itkTransformFactory.h>
#include <itkTransformFileReader.h>

#include "ResampleTubesCLP.h"

template< unsigned int Dimension >
typename itk::GroupSpatialObject< Dimension >::Pointer
ProcessTubes( itk::TransformFileReader::TransformPointer inputTransform,
  typename itk::GroupSpatialObject< Dimension >::Pointer inputTubes,
  typename itk::GroupSpatialObject< Dimension >::TransformType::Pointer
  outputIndexToObjectTransform, bool useInverseTransform )
{
  typedef typename itk::GroupSpatialObject< Dimension >::TransformType
    TransformType;
  typedef itk::tube::TubeToTubeTransformFilter< TransformType, Dimension >
    TransformFilterType;

  typename TransformType::Pointer transform =
    dynamic_cast<TransformType *>( inputTransform.GetPointer() );
  typename TransformFilterType::Pointer filter = TransformFilterType::New();

  if( useInverseTransform )
    {
    typename TransformType::InverseTransformBaseType::Pointer ivT =
      transform->GetInverseTransform();
    transform = ( TransformType * )ivT.GetPointer();
    }

  filter->SetInput( inputTubes );
  filter->SetTransform( transform );
  filter->SetOutputIndexToObjectTransform( outputIndexToObjectTransform );
  filter->Update();
  return filter->GetOutput();
}

template< unsigned int Dimension >
typename itk::GroupSpatialObject< Dimension >::Pointer
ApplyTransform( typename itk::GroupSpatialObject< Dimension >::Pointer
  inputTubes, const std::string &transformFile,
  typename itk::GroupSpatialObject< Dimension >::TransformType::Pointer
  outputIndexToObjectTransform, bool useInverseTransform )
{
  // Read transform from file
  itk::TransformFileReader::Pointer reader =
    itk::TransformFileReader::New();
  reader->SetFileName( transformFile.c_str() );
  reader->Update();

  // Transformed tubes
  typename itk::GroupSpatialObject< Dimension >::Pointer outputTubes;

  // Get list of transforms
  itk::TransformFileReader::TransformListType *tList =
    reader->GetTransformList();

  if( tList->size() != 1 )
    {
    tube::ErrorMessage( "#Transforms > 1!" );
    }

  itk::TransformFileReader::TransformListType::const_iterator tListIt;
  for( tListIt = tList->begin(); tListIt != tList->end(); ++tListIt )
    {
    outputTubes = ProcessTubes< Dimension >( *tListIt, inputTubes,
      outputIndexToObjectTransform, useInverseTransform );
    }
  return outputTubes;
}

template< unsigned int Dimension >
typename itk::GroupSpatialObject< Dimension >::Pointer
ApplyDisplacementField(
  typename itk::GroupSpatialObject< Dimension >::Pointer inputTubes,
  typename itk::GroupSpatialObject< Dimension >::TransformType::Pointer
  outputIndexToObjectTransform, const std::string & displacementFieldFile )
{
  typedef itk::DisplacementFieldTransform< double, Dimension >
    DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType
    DisplacementFieldType;
  typedef typename itk::ImageFileReader< DisplacementFieldType >
    DisplacementFieldReaderType;

  // Read the displacement field
  typename DisplacementFieldReaderType::Pointer dfReader =
    DisplacementFieldReaderType::New();
  dfReader->SetFileName( displacementFieldFile.c_str() );
  dfReader->Update();

  // Create new transform
  typename DisplacementFieldTransformType::Pointer dft =
    DisplacementFieldTransformType::New();
  dft->SetDisplacementField( dfReader->GetOutput() );

  // Define the transform filter type
  typedef itk::tube::TubeToTubeTransformFilter<
    DisplacementFieldTransformType,
    Dimension> DisplacementFieldTransformFilterType;

  // Create the filter and apply
  typename DisplacementFieldTransformFilterType::Pointer filter
    = DisplacementFieldTransformFilterType::New();
  filter->SetInput( inputTubes );
  filter->SetTransform( dft );
  filter->SetOutputIndexToObjectTransform( outputIndexToObjectTransform
    .GetPointer() );
  filter->Update();

  return filter->GetOutput();
}

template< unsigned int Dimension >
typename itk::GroupSpatialObject< Dimension >::Pointer
ReadTubeFile( const char * fileName )
{
  typedef itk::SpatialObjectReader< Dimension > SpatialObjectReaderType;

  typename SpatialObjectReaderType::Pointer reader =
    SpatialObjectReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  typename SpatialObjectReaderType::GroupType::Pointer group =
    reader->GetGroup();
  group->ComputeObjectToWorldTransform();
  return group;
}

template< unsigned int Dimension, class TransformType >
void
ReadImageTransform( const char * fileName,
  typename TransformType::Pointer outputTransform )
{
  typedef itk::Image< char, Dimension >    ImageType;
  typedef itk::ImageFileReader< ImageType> ImageReaderType;

  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  typename ImageType::SpacingType spacing =
    reader->GetOutput()->GetSpacing();
  typename ImageType::PointType origin = reader->GetOutput()->GetOrigin();
  typename ImageType::DirectionType directions =
    reader->GetOutput()->GetDirection();

  outputTransform->SetScale( spacing );
  outputTransform->SetCenter( origin );
  outputTransform->SetMatrix( directions );
}

template< unsigned int Dimension >
void WriteOutput( typename itk::GroupSpatialObject<Dimension>::Pointer
  object, const char * fileName )
{
  typedef itk::SpatialObjectWriter< Dimension > SpatialObjectWriterType;

  typename SpatialObjectWriterType::Pointer writer =
    SpatialObjectWriterType::New();
  writer->SetInput( object );
  writer->SetFileName( fileName );
  writer->Update();
}

template< unsigned int Dimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef itk::MatrixOffsetTransformBase< double, Dimension, Dimension >
    MatrixOffsetTransformType;

  typedef itk::GroupSpatialObject< Dimension > GroupSpatialObjectType;

  typedef itk::VesselTubeSpatialObject< Dimension > TubeSpatialObjectType;

  // Register transform type (this is the one produced by ANTS)
  itk::TransformFactory< MatrixOffsetTransformType >::RegisterTransform();

  double progress = 0.0;
  itk::TimeProbesCollectorBase timeCollector;

  tube::CLIProgressReporter progressReporter( "TubeTransform",
                                        CLPProcessInformation );

  progressReporter.Start();
  progressReporter.Report( progress );

  // Read in the tubes
  timeCollector.Start( "Read tubes" );
  typename GroupSpatialObjectType::Pointer tubesGroup =
    ReadTubeFile<Dimension>( inputTubeFile.c_str() );
  timeCollector.Stop( "Read tubes" );

  progress = 0.3;
  progressReporter.Report( progress );

  typename GroupSpatialObjectType::Pointer outputTubes;
  typename GroupSpatialObjectType::TransformType::Pointer
    outputIndexToObjectTransform;
  if( !matchImage.empty() )
    {
    ReadImageTransform< Dimension, typename
      GroupSpatialObjectType::TransformType >( matchImage.c_str(),
      outputIndexToObjectTransform );
    }
  else
    {
    char soTypeName[80];
    strcpy( soTypeName, "VesselTubeSpatialObject" );
    typename TubeSpatialObjectType::ChildrenListPointer tubeList =
      tubesGroup->GetChildren( tubesGroup->GetMaximumDepth(), soTypeName );
    outputIndexToObjectTransform = (*(tubeList->begin()))
      ->GetIndexToObjectTransform();
    }
  try
    {
    if( !loadDisplacementField.empty() )
      {
      timeCollector.Start(" Apply displacement field ");
      outputTubes = ApplyDisplacementField< Dimension >( tubesGroup,
        outputIndexToObjectTransform, loadDisplacementField );
      timeCollector.Stop(" Apply displacement field ");
      }
    else if( !loadTransform.empty() )
      {
      timeCollector.Start( "Apply transform" );
      outputTubes = ApplyTransform< Dimension >( tubesGroup, loadTransform,
        outputIndexToObjectTransform, useInverseTransform );
      timeCollector.Stop( "Apply transform" );
      }
    }
  catch( const std::exception &e )
    {
    tube::ErrorMessage( e.what() );
    return EXIT_FAILURE;
    }

  if( samplingFactor != 1 )
    {
    timeCollector.Start( "Sub-sample tubes" );
    typedef itk::tube::SubSampleTubeTreeSpatialObjectFilter<
      GroupSpatialObjectType, TubeSpatialObjectType >
      SubSampleTubeTreeFilterType;
    typename SubSampleTubeTreeFilterType::Pointer subSampleTubeTreeFilter =
      SubSampleTubeTreeFilterType::New();
    subSampleTubeTreeFilter->SetInput( outputTubes );
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
    outputTubes = subSampleTubeTreeFilter->GetOutput();
    timeCollector.Stop( "Sub-sample tubes" );
    }

  progress = 0.9;
  progressReporter.Report( progress );

  timeCollector.Start( "Write output");
  WriteOutput< Dimension >( outputTubes, outputTubeFile.c_str() );
  timeCollector.Stop( "Write output" );

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  // To-Do: Detect spatial object dimensionality
  return DoIt< 3 >( argc, argv );
}
