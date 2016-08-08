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

#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkTimeProbesCollectorBase.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itktubeComputeTubeMeasuresFilter.h>

// Must include CLP before including tubeCLIHelperFunctions
#include "ComputeTubeMeasuresCLP.h"

// Must do a forward declaration of DoIt before including
// tubeCLIHelperFunctions
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
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
  tube::CLIProgressReporter progressReporter( "Compute Tube Measures",
    CLPProcessInformation );
  progressReporter.Start();

  typedef TPixel                                       InputPixelType;
  typedef itk::tube::ComputeTubeMeasuresFilter
  < InputPixelType, VDimension > FilterType;
  FilterType::Pointer filter = FilterType::New();

  typedef typename FilterType::InputImageType       InputImageType;
  typedef itk::ImageFileReader< InputImageType >    ReaderType;
  typedef typename FilterType::OutputImageType      OutputImageType;
  typedef itk::ImageFileWriter< OutputImageType  >  WriterType;

  timeCollector.Start("Load data");
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFileName.c_str() );
  try
    {
    reader->Update();
    filter->SetInputImage(  reader->GetOutput() );
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Load data");

  double progress = 0.1;
  progressReporter.Report( progress );

  if( scale > 0 )
    {
    timeCollector.Start( "Gaussian Blur" );

    filter->SetScale( scale );
    try
      {
      filter->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Processing volume: Exception caught: "
        + std::string( err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    timeCollector.Stop("Gaussian Blur");

    timeCollector.Start("Save data");
    typename WriterType::Pointer writer = WriterType::New();

    writer->SetFileName( ridgeFilename );
    writer->SetInput( filter->GetRidgeness() );
    writer->SetUseCompression( true );
    writer->Update();
  
    writer->SetFileName( roundFilename );
    writer->SetInput( filter->GetRoundness() );
    writer->SetUseCompression( true );
    writer->Update();
  
    writer->SetFileName( curveFilename );
    writer->SetInput( filter->GetCurvature() );
    writer->SetUseCompression( true );
    writer->Update();
  
    writer->SetFileName( levelFilename );
    writer->SetInput( filter->GetLevelness() );
    writer->SetUseCompression( true );
    writer->Update();

    timeCollector.Stop("Save data");

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
