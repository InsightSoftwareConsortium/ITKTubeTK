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

// TubeTK includes
#include "tubeMessage.h"
#include "tubeCLIProgressReporter.h"

// TubeTKITK includes
#include "tubeComputeTubeFlyThroughImage.h"

// ITK includes
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>

#include "ComputeTubeFlyThroughImageCLP.h"

template< class TPixel, unsigned int TDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "<ModuleName>CLP.h"
//   and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // Ensure that the input image dimension is valid
  // We only support 2D and 3D Images due to the
  // limitation of itkTubeSpatialObject
  if( VDimension != 2 && VDimension != 3 )
    {
    tube::ErrorMessage(
      "Error: Only 2D and 3D Images are currently supported." );
    return EXIT_FAILURE;
    }

  // setup progress reporting
  double progress = 0.0;

  tube::CLIProgressReporter progressReporter(
    "ComputeTubeFlyThroughImage", CLPProcessInformation );
  progressReporter.Start();
  progressReporter.Report( progress );

  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;

  // Load Input Image
  //std::cout << "Loading Input Image" << std::endl;

  typedef itk::Image< TPixel, VDimension >              ImageType;
  typedef itk::ImageFileReader< ImageType >             ImageReaderType;

  timeCollector.Start( "Loading input image" );

  typename ImageReaderType::Pointer pImageReader = ImageReaderType::New();

  try
    {
    pImageReader->SetFileName( inputImageFile.c_str() );
    pImageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error loading input image: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Loading input image" );
  progress = 0.1; // At about 10% done
  progressReporter.Report( progress );

  // Load TRE File
  //std::cout << "Loading TRE File" << std::endl;

  typedef itk::SpatialObjectReader< VDimension >        TubesReaderType;

  timeCollector.Start( "Loading input TRE file" );

  typename TubesReaderType::Pointer pTubeFileReader = TubesReaderType::New();

  try
    {
    pTubeFileReader->SetFileName( inputTREFile.c_str() );
    pTubeFileReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error loading TRE file: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Loading input TRE file" );
  progress = 0.2; // At about 20% done
  progressReporter.Report( progress );

  // call ComputeTubeFlyThroughImage
  typedef tube::ComputeTubeFlyThroughImage< TPixel, VDimension >
    ComputeFlyThroughImageFilterType;

  timeCollector.Start( "Computing tube fly through images" );

  typename ComputeFlyThroughImageFilterType::Pointer pFlyThroughImageFilter
    = ComputeFlyThroughImageFilterType::New();

  pFlyThroughImageFilter->SetTubeId( inputTubeId );
  pFlyThroughImageFilter->SetInputImage( pImageReader->GetOutput() );
  pFlyThroughImageFilter->SetInput( pTubeFileReader->GetGroup() );
  pFlyThroughImageFilter->Update();

  timeCollector.Stop( "Computing tube fly through images" );
  progress = 0.8; // At about 80% done
  progressReporter.Report( progress );

  // Write fly through image
  typedef itk::ImageFileWriter< ImageType > ImageWriterType;

  timeCollector.Start( "Writing tube fly through image" );

  typename ImageWriterType::Pointer imageWriter = ImageWriterType::New();

  try
    {
    imageWriter->SetFileName( outputImageFile.c_str() );
    imageWriter->SetInput( pFlyThroughImageFilter->GetOutput() );
    imageWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error writing tube fly through image: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Writing tube fly through image" );
  progress = 0.9;
  progressReporter.Report( progress );

  // Write tube mask fly through image
  typedef typename ComputeFlyThroughImageFilterType::OutputMaskType
    MaskType;
  typedef itk::ImageFileWriter< MaskType > MaskWriterType;

  timeCollector.Start( "Writing tube mask fly through image" );

  typename MaskWriterType::Pointer maskWriter = MaskWriterType::New();

  try
    {
    maskWriter->SetFileName( outputTubeMaskFile.c_str() );
    maskWriter->SetInput( pFlyThroughImageFilter->GetOutputMask() );
    maskWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error writing tube mask fly through image: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Writing tube mask fly through image" );
  progress = 1.0;
  progressReporter.Report( progress );

  // All done
  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  try
    {
    PARSE_ARGS;
    }
  catch( const std::exception & err )
    {
    tube::ErrorMessage( err.what() );
    return EXIT_FAILURE;
    }

  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( inputImageFile, argc, argv );
}
