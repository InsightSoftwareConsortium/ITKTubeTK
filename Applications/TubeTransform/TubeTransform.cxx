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

#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

#include <itkVesselTubeSpatialObject.h>
#include <itkSpatialObjectWriter.h>
#include <itkSpatialObjectReader.h>

#include <OptionList.h>
#include "itkTubeToTubeTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactoryBase.h"

#include "TubeTransformCLP.h"


using namespace tube;

enum { TDimension = 3 };

/** Data typedefs **/
typedef itk::GroupSpatialObject<TDimension>       GroupSpatialObjectType;
typedef itk::VesselTubeSpatialObject<TDimension>  TubeSpatialObjectType;

/** IO typedefs */
typedef itk::SpatialObjectReader<TDimension>  SpatialObjectReaderType;
typedef itk::SpatialObjectWriter<TDimension>  SpatialObjectWriterType;
typedef TubeSpatialObjectType::TransformType  TransformType;
typedef itk::tube::TubeToTubeTransformFilter< TransformType, TDimension > TransformFilterType;

typedef itk::TransformFileReader                    TransformFileReaderType;
typedef itk::TransformFileReader::TransformType     BaseTransformType;
typedef itk::TransformFileReader::TransformListType TransformListType;
typedef itk::ScalableAffineTransform<double, 3>     ScalableAffineTransformType;
typedef itk::AffineTransform<double, 3>             AffineTransformType;

/** Forward decl. */
GroupSpatialObjectType::Pointer ProcessTube( const char *, const char * );
GroupSpatialObjectType::Pointer ReadTubeFile( const char * );
TransformType::Pointer          ReadTransformFile( const char * );
void WriteOutput( GroupSpatialObjectType::Pointer, const char * );
int DoIt( int, char *[] );


int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return DoIt( argc, argv );
}


/** The work happens here */
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::TimeProbesCollectorBase timeCollector;

  tube::CLIProgressReporter progressReporter(
    "tubeTransform",
    CLPProcessInformation );

  progressReporter.Start();

  timeCollector.Start( "Load data" );
  GroupSpatialObjectType::Pointer tube = ProcessTube(
      inputTubeFile.c_str(),
      transformFile.c_str()
      );

  timeCollector.Stop( "Load data" );

  double progress = 0.5;
  progressReporter.Report( progress );

  timeCollector.Start( "Save data" );
  WriteOutput( tube, outputTubeFile.c_str() );
  timeCollector.Stop( "Save data" );

  progress = 1.0;
  progressReporter.Report( progress );

  progressReporter.End();
  timeCollector.Report();

  return EXIT_SUCCESS;
}


/** Pipeline = Read tubes; Read transform; Perform transform */
GroupSpatialObjectType::Pointer ProcessTube(
  const char * tubeFileName, const char * transFileName )
{
  GroupSpatialObjectType::Pointer tubes = ReadTubeFile( tubeFileName );
  TransformType::Pointer trans = ReadTransformFile( transFileName );

  TransformFilterType::Pointer filter = TransformFilterType::New();
  filter->SetInput( tubes );
  filter->SetTransform( trans );
  filter->Update();

  return filter->GetOutput();
}


/** Reads tubes from file */
GroupSpatialObjectType::Pointer ReadTubeFile( const char * fileName )
{
  SpatialObjectReaderType::Pointer reader = SpatialObjectReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  SpatialObjectReaderType::GroupType::Pointer group = reader->GetGroup();

  group->ComputeObjectToWorldTransform();
  return group;
}


/** Reads a 'ScalableAffineTransform' from file */
TransformType::Pointer ReadTransformFile( const char * fileName )
{
  TransformFileReaderType::Pointer reader = itk::TransformFileReader::New();
  reader->SetFileName( fileName );
  reader->Update();

  TransformListType *tList =  reader->GetTransformList();
  TransformListType::const_iterator it = tList->begin();

  /**
   * CAUTION: We assume an ITK 'AffineTransform' to be read in - which
   * is usually the output of an image registration algorithm, configured
   * to use affine registration.
   */
  AffineTransformType::Pointer atp =
    static_cast< AffineTransformType * >( (*it).GetPointer() );

  // Build an ITK 'ScalableAffineTransform' based on the affine one
  ScalableAffineTransformType::Pointer satp = ScalableAffineTransformType::New();
  satp->SetTranslation( atp->GetTranslation() );
  satp->SetCenter( atp->GetCenter() );
  satp->SetMatrix( atp->GetMatrix() );

  double scale[3] = { 1, 1, 1 };
  itk::Vector<double,3> v ( scale );
  satp->SetScale( v ) ;

  std::cout << "Transl " << satp->GetTranslation() << std::endl;
  std::cout << "Center " << satp->GetCenter() << std::endl;
  std::cout << "Matrix " << satp->GetMatrix() << std::endl;
  std::cout << "Offset " << satp->GetOffset() << std::endl;
  return satp;
}


void WriteOutput(
  GroupSpatialObjectType::Pointer object,
  const char * fileName )
{
  SpatialObjectWriterType::Pointer writer = SpatialObjectWriterType::New();
  writer->SetInput( object );
  writer->SetFileName( fileName );
  writer->Update();
}
