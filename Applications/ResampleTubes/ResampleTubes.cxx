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
#include <itkTransformFactory.h>
#include <itkTransformFileReader.h>

#include "ResampleTubesCLP.h"

template< unsigned int Dimension >
typename itk::GroupSpatialObject< Dimension >::Pointer
ProcessTubes( itk::TransformFileReader::TransformType::Pointer
  inputTransform, typename itk::GroupSpatialObject< Dimension >::Pointer
  inputTubes, typename itk::GroupSpatialObject< Dimension >::TransformType::
  Pointer outputTransform, bool useInverseTransform )
{
  typedef itk::MatrixOffsetTransformBase< double, Dimension >
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
  filter->SetOutputIndexToObjectTransform( outputTransform );
  filter->Update();
  return filter->GetOutput();
}

template< unsigned int Dimension >
typename itk::GroupSpatialObject< Dimension >::Pointer
ApplyTransform( typename itk::GroupSpatialObject< Dimension >::Pointer
  inputTubes, const std::string &transformFile,
  typename itk::GroupSpatialObject< Dimension >::TransformType::
  Pointer outputTransform, bool useInverseTransform )
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

  typename itk::GroupSpatialObject< Dimension >::Pointer tmpTubes;
  tmpTubes = inputTubes;

  itk::TransformFileReader::TransformListType::const_iterator tListIt;
  for( tListIt = tList->begin(); tListIt != tList->end(); ++tListIt )
    {
    outputTubes = ProcessTubes< Dimension >( *tListIt, tmpTubes,
      outputTransform, useInverseTransform );
    tmpTubes = outputTubes;
    }

  return outputTubes;
}

template< unsigned int Dimension >
typename itk::GroupSpatialObject< Dimension >::Pointer
ApplyDisplacementField( typename itk::GroupSpatialObject< Dimension >::
  Pointer inputTubes, typename itk::GroupSpatialObject< Dimension >::
  TransformType::Pointer outputTransform, const std::string &
  displacementFieldFile )
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
  filter->SetOutputIndexToObjectTransform( outputTransform
    .GetPointer() );
  filter->Update();

  return filter->GetOutput();
}

template< unsigned int Dimension >
bool
ReadTubeFile( const char * fileName, typename itk::GroupSpatialObject<
  Dimension >::Pointer & tubesGroup )
{
  typedef itk::SpatialObjectReader< Dimension > SpatialObjectReaderType;

  typename SpatialObjectReaderType::Pointer reader =
    SpatialObjectReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  tubesGroup = reader->GetGroup();
  if( tubesGroup.IsNotNull() )
    {
    tubesGroup->ComputeObjectToWorldTransform();
    return true;
    }
  else
    {
    return false;
    }
}

template< unsigned int Dimension, class TransformType >
void
ReadImageTransform( const char * fileName,
  typename TransformType::Pointer & outputTransform )
{
  typedef itk::Image< char, Dimension >    ImageType;
  typedef itk::ImageFileReader< ImageType> ImageReaderType;

  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  typename ImageType::SpacingType spacing =
    reader->GetOutput()->GetSpacing();
  typename ImageType::PointType origin =
    reader->GetOutput()->GetOrigin();
  typename ImageType::DirectionType directions =
    reader->GetOutput()->GetDirection();

  outputTransform = TransformType::New();
  outputTransform->SetIdentity();
  itk::Vector< double, Dimension > offset;
  for( unsigned int i=0; i<Dimension; ++i )
    {
    offset[i] = origin[i];
    }
  outputTransform->SetScale( spacing );
  outputTransform->SetMatrix( directions );
  outputTransform->SetOffset( offset );
}

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
    GroupSpatialObjectType::New();
  if( !ReadTubeFile<Dimension>( inputTubeFile.c_str(), tubesGroup ) )
    {
    std::cerr << "Cannot read tubes from file: " << inputTubeFile
      << std::endl;
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Read tubes" );

  progress = 0.3;
  progressReporter.Report( progress );

  std::cout << "Here" << std::endl;
  typename GroupSpatialObjectType::Pointer outputTubes;
  typename GroupSpatialObjectType::TransformType::Pointer
    outputTransform;
  if( !matchImage.empty() )
    {
    ReadImageTransform< Dimension, typename
      GroupSpatialObjectType::TransformType >( matchImage.c_str(),
      outputTransform );
    }
  else
    {
    char soTypeName[80];
    strcpy( soTypeName, "VesselTubeSpatialObject" );
    typename TubeSpatialObjectType::ChildrenListPointer tubeList =
      tubesGroup->GetChildren( tubesGroup->GetMaximumDepth(), soTypeName );
    (*(tubeList->begin()))->ComputeObjectToWorldTransform();
    outputTransform = (*(tubeList->begin()))->GetIndexToWorldTransform();
    tubeList->clear();
    delete tubeList;
    }
  try
    {
    if( !loadDisplacementField.empty() )
      {
      timeCollector.Start(" Apply displacement field ");
      outputTubes = ApplyDisplacementField< Dimension >( tubesGroup,
        outputTransform, loadDisplacementField );
      timeCollector.Stop(" Apply displacement field ");
      }
    else if( !loadTransform.empty() )
      {
      timeCollector.Start( "Apply transform" );
      outputTubes = ApplyTransform< Dimension >( tubesGroup, loadTransform,
        outputTransform, useInverseTransform );
      timeCollector.Stop( "Apply transform" );
      }
    else if( !matchImage.empty() )
      {
      typedef itk::AffineTransform< double, Dimension >
        TransformType;
      typename TransformType::Pointer identityTransform =
        TransformType::New();
      identityTransform->SetIdentity();

      typedef itk::tube::TubeToTubeTransformFilter< TransformType,
        Dimension>
        TransformFilterType;

      typename TransformFilterType::Pointer filter =
        TransformFilterType::New();

      filter->SetInput( tubesGroup );
      filter->SetTransform( identityTransform );
      filter->SetOutputIndexToObjectTransform( outputTransform );
      filter->Update();
      outputTubes = filter->GetOutput();
      }
    else
      {
      outputTubes = tubesGroup;
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
