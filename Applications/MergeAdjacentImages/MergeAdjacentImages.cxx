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

// ITK includes
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>

// TubeTKITK includes
#include "tubeMergeAdjacentImages.h"

// TubeTK includes
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include "MergeAdjacentImagesCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

#define PARSE_ARGS_FLOAT_ONLY 1

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // setup progress reporting
  double progress = 0.0;

  tube::CLIProgressReporter progressReporter( 
    "MergeAdjacentImages", CLPProcessInformation );
  progressReporter.Start();
  progressReporter.Report( progress );

  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;

  // define types
  typedef TPixel                                        PixelType;
  typedef itk::Image< PixelType, VDimension >           ImageType;
  typedef itk::ImageFileReader< ImageType >             ReaderType;
  typedef itk::ImageFileWriter< ImageType  >            WriterType;

  // Load input image 1
  timeCollector.Start( "Loading input image 1" );

  typename ReaderType::Pointer reader1 = ReaderType::New();

  try
    {
    reader1->SetFileName( inputVolume1.c_str() );
    reader1->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error loading input image 1: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Loading input image 1" );
  progress = 0.05;
  progressReporter.Report( progress );

  // Load input image 2
  timeCollector.Start( "Loading input image 2" );

  typename ReaderType::Pointer reader2 = ReaderType::New();

  try
    {
    reader2->SetFileName( inputVolume2.c_str() );
    reader2->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error loading input image 2: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Loading input image 2" );
  progress = 0.1;
  progressReporter.Report( progress );

  // Perform merging
  timeCollector.Start( "Merging images" );

  typedef tube::MergeAdjacentImages< ImageType >
    MergeAdjacentImagesFilterType;

  typename MergeAdjacentImagesFilterType::Pointer filter =
    MergeAdjacentImagesFilterType::New();

  filter->SetInput1( reader1->GetOutput() );
  filter->SetInput2( reader2->GetOutput() );

  filter->SetBackground( background );
  filter->SetMaskZero( mask );
  filter->SetMaxIterations( iterations );
  filter->SetExpectedOffset( expectedOffset );
  filter->SetExpectedRotation( expectedRotation );
  filter->SetSamplingRatio( samplingRatio );
  filter->LoadTransform( loadTransform );
  filter->SaveTransform( saveTransform );
  filter->SetBlendUsingAverage( averagePixels );
  filter->SetUseFastBlending( useFastBlending );

  if( boundary.size() == VDimension )
    {
    filter->SetPadding( boundary );
    }

  filter->Update();

  timeCollector.Stop( "Merging images" );
  progress = 0.9;
  progressReporter.Report( progress );

  // Write output image
  timeCollector.Start( "Writing output image" );

  typename WriterType::Pointer writer = WriterType::New();

  try
    {
    writer->SetFileName( outputVolume.c_str() );
    writer->SetInput( filter->GetOutput() );
    writer->SetUseCompression( true );
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error writing output image: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Writing output image" );
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
  return tube::ParseArgsAndCallDoIt( inputVolume1, argc, argv );
}
