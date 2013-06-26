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

#include "itktubeTubeToTubeTransformFilter.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkDisplacementFieldTransform.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkInverseDisplacementFieldImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkTransformFactory.h>
#include <itkTransformFactoryBase.h>
#include <itkTransformFileReader.h>
#include <itkVector.h>
#include <itkVectorImage.h>
#include <itkVesselTubeSpatialObject.h>

#include <vnl/vnl_inverse.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "TubeTransformCLP.h"

using namespace tube;

enum { Dimension = 3 };

typedef itk::GroupSpatialObject< Dimension > GroupSpatialObjectType;

/** ProcessTubes handles the actual transformation of the tubes; if required,
  * the inverse transform is computed and applied.
  */
GroupSpatialObjectType::Pointer ProcessTubes(
  itk::TransformFileReader::TransformPointer genericInputTransform,
  GroupSpatialObjectType::Pointer inputTubes, bool useInverseTransform = false);

GroupSpatialObjectType::Pointer ApplyTransform(
  GroupSpatialObjectType::Pointer inputTubes,
  const std::string &transformFile,
  bool useInverseTransform = false );

/** Displacement fields are handled differently than transforms. The user is
  * responsible for providing the inverse displacement field (if this is what
  * is needed). ApplyDisplacementField reads in a displacement field from a
  * file, creates the correct filter type for transformation and executes the
  * filter.
  */
GroupSpatialObjectType::Pointer ApplyDisplacementField(
  GroupSpatialObjectType::Pointer inputTubes,
  const std::string &displacementFieldFile );

/** ReadTubeFile is responsible for reading in the spatial objects from
  * a file.
  */
GroupSpatialObjectType::Pointer ReadTubeFile( const char *tubeFile );

/** WriteOutput writes transformed tubes to an output file. */
void WriteOutput( GroupSpatialObjectType::Pointer, const char * );

/* Main routine */
int DoIt( int argc, char * argv[] );

GroupSpatialObjectType::Pointer
ProcessTubes( itk::TransformFileReader::TransformPointer genericInputTransform,
              GroupSpatialObjectType::Pointer inputTubes,
              bool useInverseTransform )
{
  typedef itk::Transform< double, Dimension > TransformType;
  typedef itk::tube::TubeToTubeTransformFilter< TransformType, Dimension >
    TransformFilterType;

  TransformType::Pointer transform =
    dynamic_cast<TransformType *>( genericInputTransform.GetPointer() );
  TransformFilterType::Pointer filter = TransformFilterType::New();

  if( useInverseTransform )
    {
      TransformType::InverseTransformBaseType::Pointer ivT =
      transform->GetInverseTransform();
    transform = ( TransformType * )ivT.GetPointer();
    }

  filter->SetInput( inputTubes );
  filter->SetTransform( transform );
  filter->Update();
  return filter->GetOutput();
}

GroupSpatialObjectType::Pointer
ApplyTransform( GroupSpatialObjectType::Pointer inputTubes,
                const std::string &transformFile,
                bool useInverseTransform )
{
  // Read transform from file
  itk::TransformFileReader::Pointer reader =
    itk::TransformFileReader::New();
  reader->SetFileName( transformFile.c_str() );
  reader->Update();

  // Transformed tubes
  GroupSpatialObjectType::Pointer outputTubes;

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
    outputTubes = ProcessTubes( *tListIt, inputTubes, useInverseTransform );
    }
  return outputTubes;
}

GroupSpatialObjectType::Pointer
ApplyDisplacementField( GroupSpatialObjectType::Pointer inputTubes,
                        const std::string &displacementFieldFile )
{
  typedef itk::DisplacementFieldTransform< double, Dimension > DisplacementFieldTransformType;
  typedef DisplacementFieldTransformType::DisplacementFieldType
    DisplacementFieldType;
  typedef itk::ImageFileReader< DisplacementFieldType >
    DisplacementFieldReaderType;

  // Read the displacement field
  DisplacementFieldReaderType::Pointer dfReader =
    DisplacementFieldReaderType::New();
  dfReader->SetFileName( displacementFieldFile.c_str() );
  dfReader->Update();

  // Create new transform
  DisplacementFieldTransformType::Pointer dft =
    DisplacementFieldTransformType::New();
  dft->SetDisplacementField( dfReader->GetOutput() );

  // Define the transform filter type
  typedef itk::tube::TubeToTubeTransformFilter<
    DisplacementFieldTransformType,
    Dimension> DisplacementFieldTransformFilterType;

  // Create the filter and apply
  DisplacementFieldTransformFilterType::Pointer filter
    = DisplacementFieldTransformFilterType::New();
  filter->SetInput( inputTubes );
  filter->SetTransform( dft );
  filter->Update();

  return filter->GetOutput();
}

int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef itk::MatrixOffsetTransformBase< double, Dimension, Dimension >
    MatrixOffsetTransformType;

  // Register transform type (this is the one produced by ANTS)
  itk::TransformFactory< MatrixOffsetTransformType >::RegisterTransform();

  double progress = 0.0;
  itk::TimeProbesCollectorBase timeCollector;

  CLIProgressReporter progressReporter( "TubeTransform",
                                        CLPProcessInformation );

  progressReporter.Start();
  progressReporter.Report( progress );

  // Read in the tubes
  timeCollector.Start( "Read tubes" );
  GroupSpatialObjectType::Pointer tubes =
    ReadTubeFile( inputTubeFile.c_str() );
  timeCollector.Stop( "Read tubes" );

  progress = 0.3;
  progressReporter.Report( progress );

  GroupSpatialObjectType::Pointer outputTubes;
  try
    {
    /** NOTE: In case a displacement field file is given by the user, we
      * ignore any other given transform and just apply the displacement field
      * to the tubes.
      */
    if( !displacementField.empty() )
      {
      timeCollector.Start(" Apply displacement field ");
      outputTubes = ApplyDisplacementField( tubes,
                                            displacementField );
      timeCollector.Stop(" Apply displacement field ");
      }
    else
      {
      timeCollector.Start( "Apply transform" );
      outputTubes = ApplyTransform( tubes,
                                    transformFile,
                                    useInverseTransform );
      timeCollector.Stop( "Apply transform" );
      }
    }
  catch( const std::exception &e )
    {
    tube::ErrorMessage( e.what() );
    return EXIT_FAILURE;
    }

  progress = 0.9;
  progressReporter.Report( progress );

  timeCollector.Start( "Write output");
  WriteOutput( outputTubes, outputTubeFile.c_str() );
  timeCollector.Stop( "Write output" );

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

GroupSpatialObjectType::Pointer ReadTubeFile( const char * fileName )
{
  typedef itk::SpatialObjectReader< Dimension > SpatialObjectReaderType;

  SpatialObjectReaderType::Pointer reader =
    SpatialObjectReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  SpatialObjectReaderType::GroupType::Pointer group =
    reader->GetGroup();
  group->ComputeObjectToWorldTransform();
  return group;
}

void WriteOutput( GroupSpatialObjectType::Pointer object,
                  const char * fileName )
{
  typedef itk::SpatialObjectWriter< Dimension > SpatialObjectWriterType;

  SpatialObjectWriterType::Pointer writer =
    SpatialObjectWriterType::New();
  writer->SetInput( object );
  writer->SetFileName( fileName );
  writer->Update();
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return DoIt( argc, argv );
}
