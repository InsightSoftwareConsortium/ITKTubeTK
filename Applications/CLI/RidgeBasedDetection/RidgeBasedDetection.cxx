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
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "vnl/vnl_math.h"
#include "math.h"

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
  typedef float                                                 PixelType;
  typedef itk::OrientedImage<PixelType, dimensionT>             ImageType;
  typedef itk::ImageFileReader<ImageType>                       ReaderType;
  typedef itk::ImageFileWriter<ImageType>                       WriterType;

  // typedefs for numerics
  typedef itk::NJetImageFunction<ImageType>                     CalculatorType;

  // typedefs for filters
  typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> RescaleType;

  // typedefs for iterators
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType>     IterType;

  // Setup the readers to load the input data (image + prior)
  timeCollector.Start( "Load data" );
  typename ReaderType::Pointer reader = ReaderType::New();
  typename ReaderType::Pointer priorReader = ReaderType::New();
  typename ReaderType::Pointer defectReader = ReaderType::New();
  typename ReaderType::Pointer centerlinesReader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  priorReader->SetFileName( inputPrior.c_str() );
  defectReader->SetFileName( inputDefects.c_str() );
  centerlinesReader->SetFileName( inputCenterlines.c_str() );

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

  // Load the input prior with exception handling
  try
    {
    defectReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught while reading the input defects!"
              << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // Load the input prior with exception handling
  try
    {
    centerlinesReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught while reading the input centerlines!"
              << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Load data" );
  double progress = 0.1;
  progressReporter.Report( progress );

  // Setup the sigmas based on the wire scale
  PixelType sigmaMedium = (scale/2)*0.6667;
  PixelType sigmaSmall = 0.6667*sigmaMedium;
  PixelType sigmaLarge = 1.3333*sigmaMedium;

  // Get the centerlines
  typename ImageType::Pointer centerlines = centerlinesReader->GetOutput();
  typename ImageType::Pointer defects = defectReader->GetOutput();

  // Rescale the input
  typename RescaleType::Pointer rescale = RescaleType::New();
  rescale->SetInput( reader->GetOutput() );
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( 1 );
  rescale->Update();
  typename ImageType::Pointer curImage = rescale->GetOutput();

  // Rescale the prior
  rescale = RescaleType::New();
  rescale->SetInput( priorReader->GetOutput() );
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( 1 );
  rescale->Update();
  typename ImageType::Pointer curPrior = rescale->GetOutput();

  // Setup the Calculators
  typename CalculatorType::Pointer inputCalc = CalculatorType::New();
  inputCalc->SetInputImage( curImage );
  inputCalc->SetInverseRidgeness( true );

  typename CalculatorType::Pointer priorCalc = CalculatorType::New();
  priorCalc->SetInputImage( curPrior );
  priorCalc->SetInverseRidgeness( true );

  typename CalculatorType::Pointer defectCalc = CalculatorType::New();
  defectCalc->SetInputImage( defects );

  centerlines->GetLargestPossibleRegion().Print( std::cout );
  defects->GetLargestPossibleRegion().Print( std::cout );
  curImage->GetLargestPossibleRegion().Print( std::cout );
  curPrior->GetLargestPossibleRegion().Print( std::cout );
  IterType centerlineItr( centerlines, centerlines->GetLargestPossibleRegion() );
  double samples = 0;
  centerlineItr.GoToBegin();
  while( !centerlineItr.IsAtEnd() )
    {
    if( centerlineItr.Get() > 0 )
      {
      ++samples;
      }
    ++centerlineItr;
    }

  // Get the outfile started
  std::ofstream output( outputFile.c_str() );
  output << "% Generated by RidgeBasedDetection in TubeTK\n";
  output << "@RELATION errors\n";
  output << "\n";
  output << "@ATTRIBUTE v1g NUMERIC\n";
  output << "@ATTRIBUTE v2g NUMERIC\n";
  output << "@ATTRIBUTE v3g NUMERIC\n";
  output << "@ATTRIBUTE v1e NUMERIC\n";
  output << "@ATTRIBUTE v2e NUMERIC\n";
  output << "@ATTRIBUTE v3e NUMERIC\n";
  output << "@ATTRIBUTE v1g-v1e NUMERIC\n";
  output << "@ATTRIBUTE v2g-v2e NUMERIC\n";
  output << "@ATTRIBUTE v3g-v3e NUMERIC\n";
  output << "@ATTRIBUTE log(abs(v1g-v1e)) NUMERIC\n";
  output << "@ATTRIBUTE log(abs(v2g-v2e)) NUMERIC\n";
  output << "@ATTRIBUTE log(abs(v3g-v3e)) NUMERIC\n";
  output << "@ATTRIBUTE class {good,bad}\n";
  output << "\n";
  output << "@DATA\n";


  centerlineItr.GoToBegin();
  double portion = 0.9;
  double step = portion/samples;
  unsigned int count = 0;
  unsigned int goodCount = 0;
  while( !centerlineItr.IsAtEnd() )
    {
    if( centerlineItr.Get() > 0 )
      {
      typename CalculatorType::PointType curPoint;
      centerlines->TransformIndexToPhysicalPoint( centerlineItr.GetIndex(),
                                                  curPoint );

      // Set label value
      std::string label;
      if( defectCalc->Evaluate( curPoint, sigmaMedium ) )
        {
        label = "bad";
        }
      else
        {
        label = "good";
        ++goodCount;
        if( goodCount != 29 )
          {
          progress += step;
          ++centerlineItr;
          continue;
          }
        else
          {
          goodCount = 0;
          }
        }

      PixelType v1g, v2g, v3g, v1e, v2e, v3e;

      v1g = priorCalc->Ridgeness( curPoint, sigmaMedium );
      v1e = inputCalc->Ridgeness( curPoint, sigmaMedium );

      v2g = ( priorCalc->Evaluate( curPoint, sigmaSmall ) -
              priorCalc->Evaluate( curPoint, sigmaLarge ) ) /
        priorCalc->Evaluate( curPoint, sigmaLarge );
      v2e = ( inputCalc->Evaluate( curPoint, sigmaSmall ) -
              inputCalc->Evaluate( curPoint, sigmaLarge ) ) /
        inputCalc->Evaluate( curPoint, sigmaLarge );

      v3g = priorCalc->Evaluate( curPoint );
      v3e = inputCalc->Evaluate( curPoint );

      output << v1g << "," << v2g << "," << v3g << ","
             << v1e << "," << v2e << "," << v3e << ","
             << v1g-v1e << "," << v2g-v2e << "," << v3g-v3e << ","
             << log( vnl_math_abs( v1g-v1e ) ) << ","
             << log( vnl_math_abs( v2g-v2e ) )  << "," 
             << log( vnl_math_abs( v3g-v3e ) ) << ","
             << label << "\n";

      if( count == 1000 )
        {
        progressReporter.Report( progress );
        count = 0;
        }
      else
        {
        ++count;
        }

      }

    progress += step;
    ++centerlineItr;
    }

  output.close();

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
