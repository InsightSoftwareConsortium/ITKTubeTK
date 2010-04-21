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
#include "itkOrImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "OrImagesCLP.h"

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
  tube::CLIProgressReporter    progressReporter( "OrImages",
                                                 CLPProcessInformation );
  progressReporter.Start();

  
  typedef bool                                               PixelType;
  typedef unsigned char                                      OutputPixelType;
  typedef itk::OrientedImage< PixelType,  dimensionT >       ImageType;
  typedef itk::OrientedImage< OutputPixelType, dimensionT >  OutputImageType;
  typedef itk::ImageFileReader< ImageType >                  ReaderType;
  typedef itk::ImageFileWriter< OutputImageType  >           OutputWriterType;
  typedef itk::OrImageFilter< ImageType, ImageType, OutputImageType >
                                                             FilterType;
  typedef itk::MultiplyByConstantImageFilter< OutputImageType, double,
    OutputImageType >                                        MultiplierType;
  
  const OutputPixelType OUTPUT_MAX = 
    itk::NumericTraits< OutputPixelType >::max();

  timeCollector.Start("Load data");
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume1.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume 1: Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  typename ImageType::Pointer curImage1 = reader->GetOutput();
  

  reader = ReaderType::New();
  reader->SetFileName( inputVolume2.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume 2: Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  typename ImageType::Pointer curImage2 = reader->GetOutput();

  timeCollector.Stop("Load data");
  double progress = 0.1;
  progressReporter.Report( progress );

  timeCollector.Start("Or Images");

  typename FilterType::Pointer filter;

  // Progress per iteration
  double progressFraction = 0.6;


  filter = FilterType::New();
  filter->SetInput1( curImage1 );
  filter->SetInput2( curImage2 );
  tube::CLIFilterWatcher watcher( filter,
                                  "OrImageFilter",
                                  CLPProcessInformation,
                                  progressFraction,
                                  progress,
                                  true );
  progress = 0.7;
  progressFraction = 0.2;  
  
  filter->Update();
  typename OutputImageType::Pointer curImage = filter->GetOutput();
  timeCollector.Stop("Or Images");


  timeCollector.Start("Multiply By a Constant");
  typename MultiplierType::Pointer multiplier = MultiplierType::New();
  multiplier->SetInput( curImage );
  multiplier->SetConstant( OUTPUT_MAX );
  tube::CLIFilterWatcher watcher2( multiplier,
                                   "Multiplier",
                                   CLPProcessInformation,
                                   progressFraction,
                                   progress,
                                   true );
  multiplier->Update();
  curImage = multiplier->GetOutput();
  timeCollector.Stop("Multiply By a Constant");

  timeCollector.Start("Save data");
  typename OutputWriterType::Pointer writer = OutputWriterType::New();
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
  progressReporter.End( );
  
  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume1, argc, argv );
}
