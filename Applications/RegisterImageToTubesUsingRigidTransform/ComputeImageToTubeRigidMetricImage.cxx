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

#include "itktubeImageToTubeRigidMetric.h"
#include "itktubeImageToTubeRigidRegistration.h"
#include "itktubeRecordOptimizationParameterProgressionCommand.h"
#include "itktubeSingleValuedCostFunctionImageSource.h"
#include "itktubeTubeToTubeTransformFilter.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"


#include <itkCastImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>

#include "ComputeImageToTubeRigidMetricImageCLP.h"
#include "PreProcessRegistrationInputs.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Does not currently use TPixel
#define PARSE_ARGS_FLOAT_ONLY 1

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter progressReporter(
    "RegisterImageToTubesUsingRigidTransform",
    CLPProcessInformation );
  progressReporter.Start();

  const unsigned int Dimension = 3;
  typedef double     FloatType;

  typedef itk::VesselTubeSpatialObject< Dimension >      TubeType;
  typedef itk::GroupSpatialObject< Dimension >           TubeNetType;
  typedef itk::Image< FloatType, Dimension >             ImageType;
  typedef itk::tube::ImageToTubeRigidRegistration< ImageType, TubeNetType, TubeType >
                                                         RegistrationMethodType;
  typedef typename RegistrationMethodType::TransformType TransformType;
  typedef itk::tube::TubeToTubeTransformFilter< TransformType, Dimension >
                                                         TubeTransformFilterType;
  typedef RegistrationMethodType::FeatureWeightsType     PointWeightsType;

  typename ImageType::Pointer currentImage;
  typename TubeNetType::Pointer tubeNet;
  PointWeightsType pointWeights;

  if( PreProcessRegistrationInputs< Dimension,
        FloatType,
        TubeType,
        TubeNetType,
        ImageType,
        RegistrationMethodType >( argc,
          argv,
          timeCollector,
          progressReporter,
          currentImage,
          tubeNet,
          pointWeights ) == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }


  timeCollector.Start("Sample parameter space");

  typedef itk::tube::ImageToTubeRigidMetric< ImageType,
    TubeNetType,
    TubeType > CostFunctionType;
  CostFunctionType::Pointer costFunction = CostFunctionType::New();
  costFunction->SetFixedImage( currentImage );
  costFunction->SetMovingSpatialObject( tubeNet );
  costFunction->SetFeatureWeights( pointWeights );
  TransformType::Pointer transform = TransformType::New();
  costFunction->SetTransform( transform );
  costFunction->Initialize();

  const unsigned int NumberOfParameters = 6;
  typedef itk::tube::SingleValuedCostFunctionImageSource< CostFunctionType,
    NumberOfParameters > CostFunctionImageSourceType;
  CostFunctionImageSourceType::Pointer costFunctionImageSource =
    CostFunctionImageSourceType::New();
  costFunctionImageSource->SetCostFunction( costFunction );

  typedef CostFunctionImageSourceType::ParametersType ParametersType;

  ParametersType parametersLowerBound( NumberOfParameters );
  parametersLowerBound[0] = -0.3;
  parametersLowerBound[1] = -0.3;
  parametersLowerBound[2] = -0.3;
  parametersLowerBound[3] = -5.0;
  parametersLowerBound[4] = -5.0;
  parametersLowerBound[5] = -5.0;

  ParametersType parametersUpperBound( NumberOfParameters );
  parametersUpperBound[0] = 0.3;
  parametersUpperBound[1] = 0.3;
  parametersUpperBound[2] = 0.3;
  parametersUpperBound[3] = 5.0;
  parametersUpperBound[4] = 5.0;
  parametersUpperBound[5] = 5.0;

  ParametersType parametersStep( NumberOfParameters );
  parametersStep[0] = 0.3;
  parametersStep[1] = 0.3;
  parametersStep[2] = 0.3;
  parametersStep[3] = 5.0;
  parametersStep[4] = 5.0;
  parametersStep[5] = 5.0;

  costFunctionImageSource->SetParametersLowerBound( parametersLowerBound );
  costFunctionImageSource->SetParametersUpperBound( parametersUpperBound );
  costFunctionImageSource->SetParametersStep( parametersStep );

  try
    {
    costFunctionImageSource->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Sampling parameter space: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  double progress = 0.9;
  progressReporter.Report( progress );

  timeCollector.Stop("Sample parameter space");


  timeCollector.Start("Save data");

  typedef itk::Image< float, NumberOfParameters > OutputImageType;
  typedef itk::CastImageFilter< CostFunctionImageSourceType::OutputImageType,
    OutputImageType > CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( costFunctionImageSource->GetOutput() );

  typedef itk::ImageFileWriter< OutputImageType >
    WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputMetricImage.c_str() );
  writer->SetInput( caster->GetOutput() );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing output image: Exception caught: "
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
