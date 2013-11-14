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

#include "itktubeInnerOpticToPlusImageReader.h"

#include <itkTimeProbesCollectorBase.h>
#include <itkImageFileWriter.h>
#include <itkMetaImageIO.h>
#include <itkRGBToLuminanceImageFilter.h>

// Must include CLP before including tubeCLIHelperFunctions
#include "ConvertInnerOpticToPlusCLP.h"

// Must do a forward declaration of DoIt before including
// tubeCLIHelperFunctions
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
//#include "tubeCLIHelperFunctions.h"

int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "ConvertInnerOpticToPlus",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef itk::tube::InnerOpticToPlusImageReader ReaderType;
  typedef ReaderType::OutputImageType RGBImageType;
  typedef itk::Image< ReaderType::PixelComponentType, ReaderType::ImageDimension >
    OutputImageType;

  timeCollector.Start("Load data");
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( innerOpticMetaDataFileName.c_str() );
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
  RGBImageType::Pointer inputImage = reader->GetOutput();
  timeCollector.Stop("Load data");
  double progress = 0.1;
  progressReporter.Report( progress );

  timeCollector.Start("Get luminance");
  typedef itk::RGBToLuminanceImageFilter< RGBImageType, OutputImageType >
    LuminanceFilterType;

  LuminanceFilterType::Pointer luminanceFilter =
    LuminanceFilterType::New();
  luminanceFilter->SetInput( inputImage );
  try
    {
    luminanceFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Get luminance: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Get luminance");
  progress = 0.1;
  progressReporter.Report( progress );


  timeCollector.Start("Save data");
  typedef itk::ImageFileWriter< OutputImageType > ImageWriterType;
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( metaImageFileName.c_str() );
  writer->SetInput( luminanceFilter->GetOutput() );
  writer->SetUseInputMetaDataDictionary( false );
  typedef itk::MetaImageIO ImageIOType;
  ImageIOType::Pointer metaIO = ImageIOType::New();
  metaIO->SetMetaDataDictionary( inputImage->GetMetaDataDictionary() );
  writer->SetImageIO( metaIO );
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


int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return DoIt( argc, argv );
}
