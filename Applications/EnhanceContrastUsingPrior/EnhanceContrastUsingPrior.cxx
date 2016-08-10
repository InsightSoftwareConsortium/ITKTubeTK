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

#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkFRPROptimizer.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNormalVariateGenerator.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkTimeProbesCollectorBase.h>
#include <itktubeContrastCostFunction.h>

#include "EnhanceContrastUsingPriorCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"


template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with Slicer GUI
  tube::CLIProgressReporter    progressReporter(
    "ContrastImage", CLPProcessInformation );
  progressReporter.Start();

  typedef float                                PixelType;
  typedef itk::Image< PixelType, VDimension >  ImageType;

  /** Read input images */
  typename ImageType::Pointer inputImage;
  typename ImageType::Pointer inputMask;

  timeCollector.Start("Read");
    {
    typedef itk::ImageFileReader< ImageType >   ReaderType;

    typename ReaderType::Pointer readerInputImage = ReaderType::New();
    typename ReaderType::Pointer readerInputMask = ReaderType::New();

    //read input image
    readerInputImage->SetFileName( inputImageName.c_str() );
    readerInputMask->SetFileName( inputMaskName.c_str() );

    try
      {
      readerInputImage->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading input image. Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }

    try
      {
      readerInputMask->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading input mask. Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }

    inputImage = readerInputImage->GetOutput();
    inputMask = readerInputMask->GetOutput();
    }
  progressReporter.Report( 0.1 );
  timeCollector.Stop("Read");

  //
  // Generate output image
  //
  typename ImageType::Pointer outputImage = ImageType::New();
  outputImage->CopyInformation( inputImage );
  outputImage->SetRegions( inputImage->GetLargestPossibleRegion() );
  outputImage->Allocate();
  progressReporter.Report( 0.2 );

  typedef itk::tube::ContrastCostFunction< PixelType, VDimension >
                                                ContrastCostFunctionType;
  typedef itk::OnePlusOneEvolutionaryOptimizer  InitialOptimizerType;
  typedef itk::FRPROptimizer                    OptimizerType;
  typedef itk::ImageRegionIterator< ImageType > ImageIteratorType;

  ImageIteratorType iter( inputImage,
    inputImage->GetLargestPossibleRegion() );
  double inputMin = iter.Get();
  double inputMax = inputMin;
  while( !iter.IsAtEnd() )
    {
    double tf = iter.Get();
    if( tf < inputMin )
      {
      inputMin = tf;
      }
    else if( tf > inputMax )
      {
      inputMax = tf;
      }
    ++iter;
    }

  itk::Array<double> params(3);
  params[0] = objectScale;
  params[1] = backgroundScale;
  params[2] = 20*(inputMax - inputMin);

  typename ContrastCostFunctionType::Pointer costFunc =
    ContrastCostFunctionType::New();
  costFunc->SetInputImage( inputImage );
  costFunc->SetInputMask( inputMask );
  costFunc->SetOutputImage( outputImage );
  costFunc->SetMaskObjectValue( maskObjectValue );
  costFunc->SetMaskBackgroundValue( maskBackgroundValue );

  InitialOptimizerType::Pointer initOptimizer =
    InitialOptimizerType::New();
  itk::Statistics::NormalVariateGenerator::Pointer normGen =
    itk::Statistics::NormalVariateGenerator::New();
  if( seed > 0 )
    {
    normGen->Initialize( seed );
    }
  initOptimizer->SetNormalVariateGenerator( normGen );
  initOptimizer->Initialize( 1.0 );
  initOptimizer->SetMetricWorstPossibleValue( 101 );
  initOptimizer->SetMaximumIteration( iterations*0.5 );
  initOptimizer->SetMaximize( true );

  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetUseUnitLengthGradient( true );
  optimizer->SetMaximumIteration( iterations*0.4 );
  optimizer->SetMaximumLineIteration( iterations*0.2 );
  optimizer->SetStepLength( 0.1 );
  optimizer->SetStepTolerance( 0.001 );
  optimizer->SetValueTolerance( 0.01 );
  optimizer->SetMaximize( true );

  OptimizerType::ScalesType scales( 3 );
  scales[0] = 1.0 / 0.1;
  scales[1] = 1.0 / 2.0;
  scales[2] = 1.0 / (params[2]/10);

  typename ContrastCostFunctionType::ParametersType costFunctionScales( 3 );
  costFunctionScales[0] = scales[0];
  costFunctionScales[1] = scales[1];
  costFunctionScales[2] = scales[2];

  OptimizerType::ScalesType scales2( 3 );
  scales2[0] = scales[0] * scales[0];
  scales2[1] = scales[1] * scales[1];
  scales2[2] = scales[2] * scales[2];

  // OnePlusOne should be passed squared-scales
  initOptimizer->SetScales( scales2 );
  optimizer->SetScales( scales );
  costFunc->SetScales( costFunctionScales );

  initOptimizer->SetCostFunction( costFunc );
  optimizer->SetCostFunction( costFunc );

  progressReporter.Report( 0.25 );

  costFunc->SetOutputImage( outputImage );
  costFunc->Initialize();

  initOptimizer->SetInitialPosition( params );

  progressReporter.Report( 0.3 );

  initOptimizer->StartOptimization();

  progressReporter.Report( 0.5 );

  params = initOptimizer->GetCurrentPosition();
  double result = costFunc->GetValue( params );
  std::cout << "Intermediate params = " << params
            << " Result = " << result << std::endl;

  optimizer->SetInitialPosition( params );
  optimizer->StartOptimization();

  progressReporter.Report( 0.8 );

  params = optimizer->GetCurrentPosition();
  result = costFunc->GetValue( params );
  std::cout << "Winning params = " << params
            << " Result = " << result << std::endl;

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetFileName( outputImageName.c_str() );
  writer->SetInput( outputImage );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume. Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  progressReporter.Report( 1.0 );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( inputImageName, argc, argv );
}
