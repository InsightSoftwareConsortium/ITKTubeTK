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

#include "../CLI/tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>
#include <tubeEnhanceContrastUsingPrior.h>

#include "EnhanceContrastUsingPriorCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"


template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with Slicer GUI
  tube::CLIProgressReporter    progressReporter(
    "ContrastImage", CLPProcessInformation );
  progressReporter.Start();

  typedef float                                PixelType;

  typedef tube::EnhanceContrastUsingPrior
    < PixelType, VDimension >                  FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  typedef typename FilterType::ImageType       ImageType;

  /** Read input images */
  typename ImageType::Pointer inputImage;
  typename ImageType::Pointer inputMask;

  timeCollector.Start( "Read" );
    {
    typedef itk::ImageFileReader< ImageType >   ReaderType;

    typename ReaderType::Pointer readerInputImage = ReaderType::New();
    typename ReaderType::Pointer readerInputMask = ReaderType::New();

    //read input image
    readerInputImage->SetFileName( inputImageName.c_str() );
    readerInputMask->SetFileName( inputMaskName.c_str() );

    try
      {
      readerInputImage->Update();
      filter->SetInput( readerInputImage->GetOutput() );
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading input image. Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }

    try
      {
      readerInputMask->Update();
      filter->SetInputMaskImage( readerInputMask->GetOutput() );
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading input mask. Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    }
  progressReporter.Report( 0.1 );
  timeCollector.Stop( "Read" );

  filter->SetObjectScale( objectScale );
  filter->SetBackgroundScale( backgroundScale );
  filter->SetMaskObjectValue( maskObjectValue );
  filter->SetMaskBackgroundValue( maskBackgroundValue );
  filter->SetOptimizationIterations( iterations );
  filter->SetOptimizationSeed( seed );

  timeCollector.Start( "Run Filter" );

  filter->Update();

  timeCollector.Stop( "Run Filter" );
  progressReporter.Report( 0.8 );

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetFileName( outputImageName.c_str() );
  writer->SetInput( filter->GetOutput() );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume. Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  progressReporter.Report( 1.0 );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( inputImageName, argc, argv );
}
