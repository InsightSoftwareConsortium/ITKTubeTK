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

#include "itktubeTubeXIO.h"

#include "itkGroupSpatialObject.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkImageFileReader.h"

#include "ConvertTRECLP.h"

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef itk::tube::TubeXIO< 3 >       TubeXIOType;
  typedef itk::SpatialObjectReader< 3 > SOReaderType;
  typedef itk::SpatialObjectWriter< 3 > SOWriterType;

  // The image type does not matter. We're only interested in the image size.
  typedef itk::Image< float, 3 >            ImageType;
  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  typedef TubeXIOType::SizeType             SizeType;


  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter progressReporter( "ConvertTRE",
    CLPProcessInformation );
  progressReporter.Start();
  float progress = 0;

  if( !reverse )
    {
    timeCollector.Start( "Load data" );
    TubeXIOType::Pointer reader = TubeXIOType::New();
    try
      {
      reader->Read( inputTREFileName.c_str() );
      }
    catch( ... )
      {
      tube::ErrorMessage( "Error reading TubeX file. " );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Load data" );

    progress = 0.5;
    progressReporter.Report( progress );

    timeCollector.Start( "Save data" );
    SOWriterType::Pointer writer = SOWriterType::New();
    writer->SetFileName( outputTREFileName.c_str() );
    writer->SetInput( reader->GetTubeGroup() );
    try
      {
      writer->Update();
      }
    catch( ... )
      {
      tube::ErrorMessage( "Error writing spatial objects file." );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Save data" );

    progress = 1.0;
    progressReporter.Report( progress );
    progressReporter.End();

    timeCollector.Report();
    return EXIT_SUCCESS;
    }
  else
    {
    timeCollector.Start( "Load data" );
    SOReaderType::Pointer reader = SOReaderType::New();
    reader->SetFileName( inputTREFileName.c_str() );
    try
      {
      reader->Update();
      }
    catch( ... )
      {
      tube::ErrorMessage( "Error reading spatial objects file." );
      timeCollector.Report();
      return EXIT_FAILURE;
      }

    ImageReaderType::Pointer inputImageReader = ImageReaderType::New();
    inputImageReader->SetFileName( inputImageFileName.c_str() );
    ImageType* inputImage = 0;
    try
      {
      inputImageReader->Update();
      inputImage = inputImageReader->GetOutput();
      }
    catch( ... )
      {
      tube::WarningMessage(
        "No input image found. Defaulting to a ( 1, 1, 1 ) size." );
      timeCollector.Report();
      }
    timeCollector.Stop( "Load data" );

    progress = 0.5;
    progressReporter.Report( progress );

    timeCollector.Start( "Save data" );
    TubeXIOType::Pointer writer = TubeXIOType::New();
    writer->SetTubeGroup( reader->GetGroup() );

    SizeType size;
    size.Fill( 1 );
    if( inputImage )
      {
      size = inputImageReader->GetOutput()->GetLargestPossibleRegion().
        GetSize();
      }
    writer->SetDimensions( size );

    try
      {
      writer->Write( outputTREFileName.c_str() );
      }
    catch( ... )
      {
      tube::ErrorMessage( "Error writing TubeX file. " );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Save data" );

    progress = 1.0;
    progressReporter.Report( progress );
    progressReporter.End();

    timeCollector.Report();
    return EXIT_SUCCESS;
    }
}
