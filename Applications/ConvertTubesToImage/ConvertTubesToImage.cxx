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
#include "tubeConvertTubesToImage.h"

// ITK includes
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>

#include "ConvertTubesToImageCLP.h"

template< class TPixel, unsigned int TDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "<ModuleName>CLP.h"
//   and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int Dimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // setup progress reporting
  double progress = 0.0;
  itk::TimeProbesCollectorBase timeCollector;

  tube::CLIProgressReporter progressReporter( 
    "ConvertTubesToImage",
    CLPProcessInformation );

  progressReporter.Start();
  progressReporter.Report( progress );

  // read tubes
  typedef itk::SpatialObjectReader< Dimension > TubesReaderType;

  timeCollector.Start( "Reading tubes file" );

  typename TubesReaderType::Pointer tubeFileReader = TubesReaderType::New();

  try
    {
    tubeFileReader->SetFileName( inputTubeFile.c_str() );
    tubeFileReader->Update();
    }
  catch( ... )
    {
    tube::ErrorMessage( "No readable tubes found!" );
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Reading tubes file" );
  progress = 0.1; // At about 10% done
  progressReporter.Report( progress );

  // read template image
  typedef itk::Image< TPixel, Dimension >            TemplateImageType;
  typedef itk::ImageFileReader< TemplateImageType >  TemplateImageReaderType;

  timeCollector.Start( "Reading template image" );

  typename TemplateImageReaderType::Pointer templateImageReader =
    TemplateImageReaderType::New();

  try
    {
    templateImageReader->SetFileName( inputTemplateImage.c_str() );
    templateImageReader->Update();
    }
  catch( ... )
    {
    tube::ErrorMessage( "Could not read template image!" );
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Reading template image" );
  progress = 0.2; // At about 20% done
  progressReporter.Report( progress );

  // call TubesToImageFilter
  typedef tube::ConvertTubesToImage< Dimension, TPixel >
    TubesToImageFilterType;

  timeCollector.Start( "Converting Tubes To Image" );

  typename TubesToImageFilterType::Pointer tubesToImageFilter =
    TubesToImageFilterType::New();

  tubesToImageFilter->SetUseRadius( useRadii );
  tubesToImageFilter->SetTemplateImage( templateImageReader->GetOutput() );
  tubesToImageFilter->SetInput( tubeFileReader->GetGroup() );
  tubesToImageFilter->Update();

  timeCollector.Stop( "Converting Tubes To Image" );
  progress = 0.8; // At about 80% done after filter
  progressReporter.Report( progress );

  // write tube image to file
  timeCollector.Start( "Writing tube image to file" );

  typedef itk::ImageFileWriter< TemplateImageType > TubeImageWriterType;
  typename TubeImageWriterType::Pointer tubeImageWriter =
    TubeImageWriterType::New();

  tubeImageWriter->SetFileName( outputImageFile.c_str() );
  tubeImageWriter->SetInput( tubesToImageFilter->GetOutput() );
  tubeImageWriter->SetUseCompression( true );
  tubeImageWriter->Update();

  timeCollector.Stop( "Writing tube image to file" );

  progress = 1.0;
  progressReporter.Report( progress );

  // Finish progress reporting
  progressReporter.End();
  timeCollector.Report();

  return EXIT_SUCCESS;
}

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

  return tube::ParseArgsAndCallDoIt( inputTemplateImage, argc, argv );
}
