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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// It is important to use OrientedImages
#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Includes specific to this CLI application
#include "itkAnisotropicCoherenceEnhancingDiffusionImageFilter.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "CoherenceEnhancingAnisotropicDiffusionCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;
  
  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter("CoherenceEnhancingAnisotropicDiffusion",
                                                CLPProcessInformation );
  progressReporter.Start();

  // Define the types and dimension of the images
  const unsigned int Dimension                = 3;
  typedef double                              PixelType;
  typedef itk::Image< PixelType, Dimension >  InputImageType;
  typedef itk::Image< PixelType, Dimension >  OutputImageType;

  // Read the input volume
  timeCollector.Start("Load data");
  typedef itk::ImageFileReader< InputImageType  >      ImageReaderType;
  typename ImageReaderType::Pointer reader = ImageReaderType::New();
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

  // Perform the coherence enhancing anisotropic diffusion
  timeCollector.Start("Coherence enhancing anisotropic diffusion");

  // Declare the anisotropic diffusion coherence enhancing filter
  typedef itk::AnisotropicCoherenceEnhancingDiffusionImageFilter< InputImageType,
                                           OutputImageType>  CoherenceEnhancingFilterType;

  // Create a coherence enhancing filter
  typename CoherenceEnhancingFilterType::Pointer CoherenceEnhancingFilter =
      CoherenceEnhancingFilterType::New();

  CoherenceEnhancingFilter->SetInput( reader->GetOutput() );

  //Set/Get CED parameters
  CoherenceEnhancingFilter->SetSigma( scaleParameter );
  CoherenceEnhancingFilter->SetAlpha( alpha );
  CoherenceEnhancingFilter->SetContrastParameterLambdaC( cedContrastParameter );
  CoherenceEnhancingFilter->SetTimeStep( timeStep );
  CoherenceEnhancingFilter->SetNumberOfIterations( numberOfIterations );

  double progressFraction = 0.8;
  tube::CLIFilterWatcher watcher( CoherenceEnhancingFilter,
                                  "Coherence enhancing anisotropic diffusion",
                                  CLPProcessInformation,
                                  progressFraction,
                                  progress,
                                  true );

  try
    {
    CoherenceEnhancingFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Coherence enhancing anisotropic diffusion: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop("Coherence enhancing anisotropic diffusion");
  progress = 0.9;
  progressReporter.Report( progress );

  // Save output data
  timeCollector.Start("Save data");
  typedef itk::ImageFileWriter< OutputImageType  >      ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput ( CoherenceEnhancingFilter->GetOutput() );

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

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
