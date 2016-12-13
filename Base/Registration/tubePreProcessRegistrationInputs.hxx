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

#ifndef __tubePreProcessRegistrationInputs_hxx
#define __tubePreProcessRegistrationInputs_hxx

#include "tubePreProcessRegistrationInputs.h"

#include "itktubeSubSampleTubeTreeSpatialObjectFilter.h"
#include "itktubeTubeAngleOfIncidenceWeightFunction.h"
#include "itktubeTubeExponentialResolutionWeightFunction.h"
#include "itktubeTubePointWeightsCalculator.h"
#include "itkUltrasoundProbeGeometryCalculator.h"

#include <itkImageFileReader.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>

#include "tubeCLIProgressReporter.h"

namespace tube {

template<
  unsigned int VDimension,
  typename TFloat,
  typename TTube,
  typename TTubeNet,
  typename TImage,
  typename TRegistrationMethod >
int
PreProcessRegistrationInputs( int argc,
  char * argv[],
  itk::TimeProbesCollectorBase & timeCollector,
  tube::CLIProgressReporter & progressReporter,
  typename TImage::Pointer & currentImage,
  typename TTubeNet::Pointer & tubeNet,
  typename TRegistrationMethod::FeatureWeightsType & pointWeights )
{
  const unsigned int Dimension = VDimension;
  typedef TFloat                                     FloatType;
  typedef TTube                                      TubeType;
  typedef TTubeNet                                   TubeNetType;
  typedef TImage                                     ImageType;
  typedef TRegistrationMethod                        RegistrationMethodType;

  PARSE_ARGS;

  timeCollector.Start( "Load data" );
  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  typename ImageReaderType::Pointer reader = ImageReaderType::New();
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

  typedef itk::SpatialObjectReader< Dimension > TubeNetReaderType;
  typename TubeNetReaderType::Pointer vesselReader = TubeNetReaderType::New();
  vesselReader->SetFileName( inputVessel );
  try
    {
    vesselReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading vessel: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Load data" );
  double progress = 0.1;
  progressReporter.Report( progress );

  timeCollector.Start( "Sub-sample data" );
  typedef itk::tube::SubSampleTubeTreeSpatialObjectFilter< TubeNetType, TubeType >
                                                         SubSampleTubeTreeFilterType;
  typename SubSampleTubeTreeFilterType::Pointer subSampleTubeTreeFilter =
    SubSampleTubeTreeFilterType::New();
  subSampleTubeTreeFilter->SetInput( vesselReader->GetGroup() );
  subSampleTubeTreeFilter->SetSampling( 100 );
  try
    {
    subSampleTubeTreeFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Sub-sampling vessel: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  tubeNet = subSampleTubeTreeFilter->GetOutput();

  timeCollector.Stop( "Sub-sample data" );
  progress = 0.2;
  progressReporter.Report( progress );


  currentImage = reader->GetOutput();
  if( gaussianBlurStdDev > 0.0 )
    {
    timeCollector.Start( "Gaussian Blur" );

    typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType >
      GaussianFilterType;
    typename GaussianFilterType::Pointer gaussianFilter;

    // Progress per iteration
    const double progressFraction = 0.1/Dimension;
    for( unsigned int ii = 0; ii < Dimension; ++ii )
      {
      gaussianFilter = GaussianFilterType::New();
      gaussianFilter->SetInput( currentImage );
      gaussianFilter->SetSigma( gaussianBlurStdDev );

      gaussianFilter->SetOrder(
               itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
      gaussianFilter->SetDirection( ii );
      tube::CLIFilterWatcher watcher( gaussianFilter,
                                      "Blur Filter 1D",
                                      CLPProcessInformation,
                                      progressFraction,
                                      progress,
                                      true );

      gaussianFilter->Update();
      currentImage = gaussianFilter->GetOutput();
      }

    timeCollector.Stop( "Gaussian Blur" );
    }

  timeCollector.Start( "Compute Model Feature Weights" );

  typedef itk::tube::Function::TubeExponentialResolutionWeightFunction<
    typename TubeType::TubePointType, FloatType >
    WeightFunctionType;
  typedef typename RegistrationMethodType::FeatureWeightsType
    PointWeightsType;
  typedef itk::tube::TubePointWeightsCalculator< Dimension,
    TubeType, WeightFunctionType, PointWeightsType >
    PointWeightsCalculatorType;

  typename WeightFunctionType::Pointer weightFunction =
    WeightFunctionType::New();
  typename PointWeightsCalculatorType::Pointer resolutionWeightsCalculator =
    PointWeightsCalculatorType::New();

  resolutionWeightsCalculator->SetTubeTreeSpatialObject(
    subSampleTubeTreeFilter->GetOutput() );
  resolutionWeightsCalculator->SetPointWeightFunction( weightFunction );
  resolutionWeightsCalculator->Compute();

  pointWeights = resolutionWeightsCalculator->GetPointWeights();

  timeCollector.Stop( "Compute Model Feature Weights" );

  return EXIT_SUCCESS;
}

} //namespace

#endif // End !defined( __tubePreProcessRegistrationInputs_hxx )
