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

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include <itkTimeProbesCollectorBase.h>

#include <itkVesselTubeSpatialObject.h>
#include <itkSpatialObjectWriter.h>
#include <itkSpatialObjectReader.h>

#include <itkVector.h>
#include <itkImage.h>
#include <itkDisplacementFieldTransform.h>
#include <itkInverseDisplacementFieldImageFilter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFactoryBase.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkVectorImage.h>

#include <itkTransformFileReader.h>
#include <itkResampleImageFilter.h>
#include <itkTransformFactory.h>
#include <itkVersorRigid3DTransform.h>

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_inverse.h>

#include "tubeMessage.h"
#include "TubeTransformCLP.h"
#include "itkTubeToTubeTransformFilter.h"


using namespace tube;

enum { TDimension = 3 };

typedef itk::GroupSpatialObject< TDimension >       GroupSpatialObjectType;
typedef itk::VesselTubeSpatialObject< TDimension >  TubeSpatialObjectType;
typedef itk::SpatialObjectReader< TDimension >      SpatialObjectReaderType;
typedef itk::SpatialObjectWriter< TDimension >      SpatialObjectWriterType;
//typedef itk::AffineTransform< double, TDimension >  AffineTransformType;
//typedef itk::VersorRigid3DTransform< double >       VersorRigidTransformType;

itk::Transform<double, TDimension>            TransformType;


typedef itk::DisplacementFieldTransform<double,
                                        TDimension> DisplacementFieldTransformType;
typedef DisplacementFieldTransformType
  ::DisplacementFieldType                           DisplacementFieldType;
typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFieldReaderType;
typedef itk::MatrixOffsetTransformBase<double,
                                       TDimension,
                                       TDimension>  MatrixOffsetTransformType;

typedef itk::tube::TubeToTubeTransformFilter< TransformType, TDimension >
  TransformFilterType;


/** ProcessTubes handles the actual transformation of the tubes; if required,
  * the inverse transform is computed and applied.
  */
//template < typename TransformType >
GroupSpatialObjectType::Pointer ProcessTubes(
  itk::TransformFileReader::TransformPointer genericInputTransform,
  GroupSpatialObjectType::Pointer inputTubes,
  bool useInverseTransform = false);

/** ApplyTransform 1) reads in a transform from a file, 2) chooses the right
  * filter type (according to the transform information) and 3) eventually
  * passes the right transform and the tubes to ProcessTubes.
  */
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
int DoIt( int, char *[] );


//template < typename TransformType >
GroupSpatialObjectType::Pointer
ProcessTubes( itk::TransformFileReader::TransformPointer genericInputTransform,
              GroupSpatialObjectType::Pointer inputTubes,
              bool useInverseTransform )
{
  // Get pointer to transform
  //typename TransformType::Pointer transform =
  TransformType::Pointer transform =
    dynamic_cast<TransformType *>( genericInputTransform.GetPointer() );

  // Define filter type and create filter
  //typedef itk::tube::TubeToTubeTransformFilter<TransformType, TDimension>
  //  TransformFilterType;
  //typename TransformFilterType::Pointer filter = TransformFilterType::New();
  TransformFilterType::Pointer filter = TransformFilterType::New();

  if( useInverseTransform )
    {
	  //typename TransformType::InverseTransformBaseType::Pointer ivT =
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
  for( tListIt=tList->begin(); tListIt != tList->end(); tListIt++)
    {
    outputTubes = ProcessTubes( *tListIt, inputTubes, useInverseTransform ); 
    }  
  return outputTubes;
}


GroupSpatialObjectType::Pointer
ApplyDisplacementField( GroupSpatialObjectType::Pointer inputTubes,
                        const std::string &displacementFieldFile )
{
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
    TDimension> TransformFilterType;

  // Create the filter and apply
  TransformFilterType::Pointer filter = TransformFilterType::New();
  filter->SetInput( inputTubes );
  filter->SetTransform( dft );
  filter->Update();

  return filter->GetOutput();
}


int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // Register transform type (this is the one produced by ANTS)
  itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();

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
  catch( std::exception &e )
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
