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
#include "itkNJetImageFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "RidgeBasedDetectionCLP.h"

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
  tube::CLIProgressReporter progressReporter( "RidgeBasedDetection",
                                              CLPProcessInformation );
  progressReporter.Start();

  // typedefs for data structures
  typedef float                                             PixelType;
  typedef itk::OrientedImage<PixelType, dimensionT>         ImageType;
  typedef itk::ImageFileReader<ImageType>                   ReaderType;
  typedef itk::ImageFileWriter<ImageType>                   WriterType;

  // typedefs for numerics
  typedef itk::NJetImageFunction<ImageType>                 CalculatorType;

  // typedefs for iterators
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> ConstIterType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType>      IterType;

  // Setup the readers to load the input data (image + prior)
  timeCollector.Start( "Load data" );
  typename ReaderType::Pointer reader = ReaderType::New();
  typename ReaderType::Pointer priorReader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  priorReader->SetFileName( inputPrior.c_str() );

  // Load the input image with exception handling
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught while reading the input image!"
              << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // Load the input prior with exception handling
  try
    {
    priorReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught while reading the input prior!"
              << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Load data" );
  double progress = 0.1;
  progressReporter.Report( progress );

  // Store the input images and prepare the target regions
  typename ImageType::Pointer curImage = reader->GetOutput();
  typename ImageType::Pointer curPrior = priorReader->GetOutput();

  // Allocate the output image and fill with 0s
  typename ImageType::Pointer outImage = ImageType::New();
  outImage->CopyInformation( curImage );
  outImage->SetRegions( curImage->GetLargestPossibleRegion() );
  outImage->Allocate();
  outImage->FillBuffer( 0 );  

  // Setup the Calculator
  typename CalculatorType::Pointer calc = CalculatorType::New();
  calc->SetInputImage( curImage );
  calc->SetInverseRidgeness( true );

  ConstIterType inputItr( curImage, curImage->GetLargestPossibleRegion() );
  IterType outputItr( outImage, outImage->GetLargestPossibleRegion() );
  double samples = 0;
  inputItr.GoToBegin();
  while( !inputItr.IsAtEnd() && !outputItr.IsAtEnd() )
    {
    ++samples;
    ++inputItr;
    }
  inputItr.GoToBegin();
  outputItr.GoToBegin();
  double portion = 0.8;
  double step = portion/samples;
  unsigned int count = 0;
  while( !inputItr.IsAtEnd() && !outputItr.IsAtEnd() )
    {
    outputItr.Set( calc->RidgenessAtIndex( inputItr.GetIndex(), scale ) );
    progress += step;
    if( count == 1000 )
      {
      progressReporter.Report( progress );
      count = 0;
      }
    else
      {
      ++count;
      }
    ++inputItr;
    ++outputItr;
    }

  timeCollector.Start("Save data");
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( outImage );
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

  // You May Need To update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
