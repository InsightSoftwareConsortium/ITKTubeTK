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

#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkTimeProbesCollectorBase.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>

// Must include CLP before including tubeCLIHelperFunctions
#include "SampleCLIApplicationCLP.h"

// Must do a forward declaration of DoIt before including
// tubeCLIHelperFunctions
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter progressReporter( "SampleCLIApplication",
    CLPProcessInformation );
  progressReporter.Start();

  typedef TPixel                                    InputPixelType;
  typedef itk::Image< InputPixelType, VDimension >  InputImageType;
  typedef itk::ImageFileReader< InputImageType >    ReaderType;

  typedef float                                     OutputPixelType;
  typedef itk::Image< OutputPixelType, VDimension > OutputImageType;
  typedef itk::ImageFileWriter< OutputImageType  >  WriterType;

  typedef itk::DiscreteGaussianImageFilter< InputImageType,
    OutputImageType > FilterType;

  timeCollector.Start( "Load data" );
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFileName.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Load data" );

  double progress = 0.1;
  progressReporter.Report( progress );

  if( gaussianBlurStdDev > 0 )
    {
    timeCollector.Start( "Gaussian Blur" );

    typename FilterType::Pointer filter;
    filter = FilterType::New();
    filter->SetInput( reader->GetOutput() );
    filter->SetVariance( gaussianBlurStdDev * gaussianBlurStdDev );

    double progressFraction = 0.8;
    tube::CLIFilterWatcher watcher( filter, "Blur Filter 1D",
      CLPProcessInformation, progressFraction, progress, true );

    try
      {
      filter->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Processing volume: Exception caught: "
        + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Gaussian Blur" );

    timeCollector.Start( "Save data" );
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputImageFileName.c_str() );
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
  else
    {
    tubeErrorMacro( << "Specified scale <= 0 is not supported." );
    progress = 1.0;
    progressReporter.Report( progress );
    progressReporter.End();
    timeCollector.Report();
    return EXIT_FAILURE;
    }
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputImageFileName.
  return tube::ParseArgsAndCallDoIt( inputImageFileName, argc, argv );
}
