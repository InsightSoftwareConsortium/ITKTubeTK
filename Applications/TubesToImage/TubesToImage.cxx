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

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGroupSpatialObject.h"
#include "itkSpatialObjectReader.h"
#include "itkAffineTransform.h"
#include "itkTubeSpatialObjectToImageFilter.h"

#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

#include "itkTubeToTubeTransformFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactoryBase.h"

// CLI parsing
#include "TubesToImageCLP.h"

using namespace tube;

const int Dimensions = 3;

typedef unsigned char                                   OutputPixelType;

typedef itk::Image< unsigned char, Dimensions >         TemplateImageType;
typedef itk::ImageFileReader< TemplateImageType >       TemplateImageReaderType;
typedef itk::Image< OutputPixelType, Dimensions >       OutputImageType;

typedef itk::tube::TubeSpatialObjectToImageFilter<
    Dimensions,
    OutputImageType >                                   TubetoImageFilterType;

typedef itk::GroupSpatialObject< Dimensions >           TubeGroupType;
typedef TubeGroupType::Pointer                          TubeGroupPointer;
typedef TubeGroupType                                   TubesType;
typedef itk::SpatialObjectReader< Dimensions >          TubesReaderType;

/** Read tube file from disk */
TubesType::Pointer ReadTubes( const char * file );

/** Write tube image to disk */
void WriteImage( const char * file, OutputImageType::Pointer image );

/** Main routine */
int DoIt( int, char *[] );


int main(int argc, char *argv[])
{
  PARSE_ARGS;
  return DoIt( argc, argv );
}


/** Main work happens here */
int DoIt( int argc, char *argv[] )
{
  PARSE_ARGS;

  double progress = 0.0;
  itk::TimeProbesCollectorBase timeCollector;

  CLIProgressReporter progressReporter(
    "TubesToImage",
    CLPProcessInformation );

  progressReporter.Start();
  progressReporter.Report( progress );

  // read tubes
  TubesType::Pointer tubes = ReadTubes( inputTubeFile.c_str() );
  progress = 0.1; // At about 10% done
  progressReporter.Report( progress );

  // configure template image reader
  TemplateImageReaderType::Pointer imTemplateReader;
  imTemplateReader = TemplateImageReaderType::New();
  imTemplateReader->SetFileName(inputTemplateImage.c_str());

  timeCollector.Start( "Loading template image" );
  imTemplateReader->Update();
  timeCollector.Stop( "Loading template image" );

  progress = 0.2; // At about 20% done
  progressReporter.Report( progress );

  // get image's spacing and size
  TemplateImageType::Pointer imT = imTemplateReader->GetOutput();
  TubetoImageFilterType::SizeType size;
  double spacing[Dimensions];
  for(int i = 0; i < Dimensions; i++ )
    {
    size[i] = imT->GetLargestPossibleRegion().GetSize()[i];
    spacing[i] = imT->GetSpacing()[i];
    }

  // configure filter
  TubetoImageFilterType::Pointer
    tubeFilter = TubetoImageFilterType::New();
  tubeFilter->SetBuildRadiusImage( false );
  tubeFilter->SetBuildTangentImage( false );
  tubeFilter->SetUseRadius( false );
  tubeFilter->SetSize( size );
  tubeFilter->SetSpacing( spacing );
  tubeFilter->SetInput( tubes );

  // update filter
  timeCollector.Start( "Update filter" );
  tubeFilter->Update();
  timeCollector.Stop( "Update filter" );

  progress = 0.8; // At about 80% done after filter
  progressReporter.Report( progress );

  timeCollector.Start( "Save data" );
  WriteImage( outputImageFile.c_str(), tubeFilter->GetOutput() );
  timeCollector.Stop( "Save data" );

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();
  timeCollector.Report();
  return EXIT_SUCCESS;
}


TubesType::Pointer ReadTubes( const char * file )
{
  TubesReaderType::Pointer reader = TubesReaderType::New();
  try
    {
    reader->SetFileName( file );
    reader->Update();
    }
  catch( ... )
    {
    tube::ErrorMessage( "No readable tubes found!" );
    return NULL;
    }

  return reader->GetGroup();
}


void WriteImage( const char * file, OutputImageType::Pointer image )
{
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer  writer = WriterType::New();

  writer->SetFileName( file );
  writer->SetInput( image );
  writer->SetUseCompression( true );
  writer->Update();
}
