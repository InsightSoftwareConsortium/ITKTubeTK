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
#include "tubeSegmentUsingOtsuThreshold.h"

// TubeTK includes
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include "SegmentUsingOtsuThresholdCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // setup progress reporting
  double progress = 0.0;

  tube::CLIProgressReporter progressReporter(
    "SegmentUsingOtsuThreshold", CLPProcessInformation );
  progressReporter.Start();
  progressReporter.Report( progress );

  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;

  // typedefs
  typedef tube::SegmentUsingOtsuThreshold< TPixel, VDimension > FilterType;

  // Load input image
  timeCollector.Start( "Load data" );

  typedef typename FilterType::InputImageType     InputImageType;
  typedef itk::ImageFileReader< InputImageType >  ImageReaderType;

  typename ImageReaderType::Pointer inputReader = ImageReaderType::New();

  try
    {
    inputReader->SetFileName( inputVolume.c_str() );
    inputReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error loading input image: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  // Load mask image if provided
  typedef typename FilterType::MaskImageType     MaskImageType;
  typedef itk::ImageFileReader< MaskImageType >  MaskReaderType;

  typename MaskReaderType::Pointer maskReader = MaskReaderType::New();

  if( maskVolume.size() > 0 )
    {
    try
      {
      maskReader->SetFileName( maskVolume.c_str() );
      maskReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Error reading input mask: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    }

  timeCollector.Stop( "Load data" );
  progress = 0.1;
  progressReporter.Report( progress );

  // run otsu thresholding
  timeCollector.Start( "Otsu thresholding" );

  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( inputReader->GetOutput() );

  if( maskVolume.size() > 0 )
    {
    filter->SetMaskValue( maskValue );
    filter->SetMaskImage( maskReader->GetOutput() );
    }

  filter->Update();

  std::cout << "Chosen threshold = " << filter->GetThreshold() << std::endl;

  timeCollector.Stop( "Otsu thresholding" );
  progress = 0.8; // At about 80% done
  progressReporter.Report( progress );

  // write output
  typedef typename FilterType::OutputImageType      OutputImageType;
  typedef itk::ImageFileWriter< OutputImageType >   OutputWriterType;

  timeCollector.Start( "Write segmentation mask" );

  typename OutputWriterType::Pointer writer = OutputWriterType::New();

  try
    {
    writer->SetFileName( outputVolume.c_str() );
    writer->SetInput( filter->GetOutput() );
    writer->SetUseCompression( true );
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error writing segmentation mask: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Write segmentation mask" );
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
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
