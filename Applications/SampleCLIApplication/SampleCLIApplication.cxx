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

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// The following three should be used in every CLI application
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include <itkTimeProbesCollectorBase.h>

// Includes specific to this CLI application
#include <itkRecursiveGaussianImageFilter.h>

// Must include CLP before including tubeCLIHelperFunctions
#include "SampleCLIApplicationCLP.h"

// Must do a forward declaraction of DoIt before including
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
  tube::CLIProgressReporter    progressReporter( "SampleCLIApplication",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef float                                 PixelType;
  typedef itk::Image< PixelType, VDimension >   ImageType;
  typedef itk::ImageFileReader< ImageType >     ReaderType;

  timeCollector.Start("Load data");
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Load data");
  double progress = 0.1;
  progressReporter.Report( progress );

  typename ImageType::Pointer curImage = reader->GetOutput();

  if( gaussianBlurStdDev > 0 )
    {
    timeCollector.Start("Gaussian Blur");

    typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType >
      FilterType;
    typename FilterType::Pointer filter;

    // Progress per iteration
    double progressFraction = 0.8/VDimension;

    for(unsigned int i=0; i<VDimension; i++)
      {
      filter = FilterType::New();
      filter->SetInput( curImage );
      filter->SetNormalizeAcrossScale( true );
      filter->SetSigma( gaussianBlurStdDev );

      filter->SetOrder(
               itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
      filter->SetDirection( i );
      tube::CLIFilterWatcher watcher( filter,
                                      "Blur Filter 1D",
                                      CLPProcessInformation,
                                      progressFraction,
                                      progress,
                                      true );

      filter->Update();
      curImage = filter->GetOutput();
      }

    timeCollector.Stop("Gaussian Blur");
    }

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;

  timeCollector.Start("Save data");
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( curImage );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume: Exception caught: "
      + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Save data");
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
