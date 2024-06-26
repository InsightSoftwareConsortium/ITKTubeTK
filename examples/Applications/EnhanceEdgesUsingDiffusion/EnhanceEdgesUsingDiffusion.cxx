/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
// TubeTKITK includes
#include "tubeEnhanceEdgesUsingDiffusion.h"

#include "../CLI/tubeCLIFilterWatcher.h"
#include "../CLI/tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkImageFileReader.h>
#include <itkOrientImageFilter.h>
#include <itkTimeProbesCollectorBase.h>

#include "EnhanceEdgesUsingDiffusionCLP.h"

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

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter(
    "EdgeEnhancingAnisotropicDiffusion", CLPProcessInformation );
  progressReporter.Start();

  // Define the types and dimension of the images
  // Use the input image to dictate the type of the image reader/writer,
  // but use double for the filter to avoid rounding off errors in the
  // filter's floating point operations
  enum { Dimension = 3 };
  typedef TPixel                                   ImagePixelType;
  typedef itk::Image< ImagePixelType, Dimension >  InputImageType;
  typedef itk::Image< ImagePixelType, Dimension >  OutputImageType;
  typedef double                                   FilterPixelType;
  typedef itk::Image< FilterPixelType, Dimension > FilterInputImageType;
  typedef itk::Image< FilterPixelType, Dimension > FilterOutputImageType;

  // Read the input volume
  timeCollector.Start( "Load data" );
  typedef itk::ImageFileReader< InputImageType  >  ImageReaderType;
  typename ImageReaderType::Pointer   reader = ImageReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
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

  // C-style cast from input image type to type 'double'
  typedef itk::CastImageFilter< InputImageType, FilterInputImageType >
    CastInputImageFilterType;
  typename CastInputImageFilterType::Pointer castInputImageFilter =
    CastInputImageFilterType::New();
  castInputImageFilter->SetInput( reader->GetOutput() );
  castInputImageFilter->Update();

  // Reorient to axial because the anisotropic diffusion tensor function does
  // not handle direction
  typedef itk::OrientImageFilter< FilterInputImageType, FilterInputImageType >
      OrientInputFilterType;
  typename OrientInputFilterType::Pointer orientInputFilter
      = OrientInputFilterType::New();
  orientInputFilter->UseImageDirectionOn();
  orientInputFilter->SetDesiredCoordinateOrientationToAxial();
  orientInputFilter->SetInput( castInputImageFilter->GetOutput() );
  orientInputFilter->Update();

  // Perform the edge enhancing anisotropic diffusion
  timeCollector.Start( "Edge enhancing anisotropic diffusion" );

  // Declare the anisotropic diffusion edge enhancement filter
  typedef tube::EnhanceEdgesUsingDiffusion<
    FilterInputImageType, FilterOutputImageType>  EdgeEnhancementFilterType;

  // Create a edge enhancement Filter
  typename EdgeEnhancementFilterType::Pointer EdgeEnhancementFilter =
    EdgeEnhancementFilterType::New();

  EdgeEnhancementFilter->SetInput( orientInputFilter->GetOutput() );

  //Set/Get EED parameters
  EdgeEnhancementFilter->SetSigma( scaleParameter );
  EdgeEnhancementFilter->SetContrastParameterLambdaE( eedContrastParameter );
  EdgeEnhancementFilter->SetTimeStep( timeStep );
  EdgeEnhancementFilter->SetNumberOfIterations( numberOfIterations );

  double progressFraction = 0.8;
  tube::CLIFilterWatcher watcher( EdgeEnhancementFilter,
    "Edge enhancement anisotropic diffusion", CLPProcessInformation,
    progressFraction, progress, true );

  try
    {
    EdgeEnhancementFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage(
      "Edge enhancing anisotropic diffusion: Exception caught: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Edge enhancing anisotropic diffusion" );
  progress = 0.9;
  progressReporter.Report( progress );

  // C-style cast from double to output image type
  typedef itk::CastImageFilter< FilterOutputImageType, OutputImageType >
    CastOutputImageFilterType;
  typename CastOutputImageFilterType::Pointer castOutputImageFilter =
    CastOutputImageFilterType::New();
  castOutputImageFilter->SetInput( EdgeEnhancementFilter->GetOutput() );
  castOutputImageFilter->Update();

  // Reorient back from axial to whatever direction we had before
  typedef itk::OrientImageFilter< OutputImageType, OutputImageType >
      OrientOutputFilterType;
  typename OrientOutputFilterType::Pointer orientOutputFilter
      = OrientOutputFilterType::New();
  orientOutputFilter->UseImageDirectionOn();
  orientOutputFilter->SetDesiredCoordinateDirection(
      reader->GetOutput()->GetDirection() );
  orientOutputFilter->SetInput( castOutputImageFilter->GetOutput() );
  orientOutputFilter->Update();

  // Save output data
  timeCollector.Start( "Save data" );
  typedef itk::ImageFileWriter< OutputImageType  >      ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput ( orientOutputFilter->GetOutput() );

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

  // All objects should be automatically destroyed at this point
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
