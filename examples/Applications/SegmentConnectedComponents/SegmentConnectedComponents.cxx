/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#include "../CLI/tubeCLIFilterWatcher.h"
#include "../CLI/tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeSegmentConnectedComponents.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>


#include "SegmentConnectedComponentsCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#define PARSE_ARGS_INT_ONLY
#include "../CLI/tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "ConnectedComponents",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef itk::Image< TPixel, VDimension >         MaskType;
  typedef itk::ImageFileReader< MaskType >         MaskReaderType;

  //
  //
  //
  timeCollector.Start( "Load mask" );
  typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
  maskReader->SetFileName( inputMask.c_str() );
  try
    {
    maskReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading mask: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Load mask" );
  double progress = 0.1;
  progressReporter.Report( progress );

  typename MaskType::Pointer curMask = maskReader->GetOutput();

  //
  //
  //
  timeCollector.Start( "Connected Components" );

  typedef tube::SegmentConnectedComponents< MaskType, MaskType >
    FilterType;
  typename FilterType::Pointer filter;

  filter = FilterType::New();
  filter->SetInput( curMask );

  filter->SetMinimumVolume( minSize );

  if( seedMask.size() > 0 )
    {
    timeCollector.Start( "Load seed mask" );
    typename MaskReaderType::Pointer seedMaskReader = MaskReaderType::New();
    seedMaskReader->SetFileName( seedMask.c_str() );
    try
      {
      seedMaskReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading seed mask: Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Load seed mask" );

    filter->SetSeedMask( seedMaskReader->GetOutput() );
    }

  filter->Update();

  typedef itk::ImageFileWriter< MaskType  >   ImageWriterType;

  timeCollector.Start( "Save data" );
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputMask.c_str() );
  writer->SetInput( filter->GetOutput() );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume: Exception caught: "
      + std::string( err.GetDescription() ) );
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

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputMask, argc, argv );
}
