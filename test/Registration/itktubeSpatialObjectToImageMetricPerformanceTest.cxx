/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#include "itktubePointBasedSpatialObjectToImageMetric.h"

#include "itktubeSubSampleSpatialObjectFilter.h"

#include <itkImageFileReader.h>
#include <itkMemoryProbesCollectorBase.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkComposeScaleSkewVersor3DTransform.h>

/**
 *  This test exercised the metric evaluation methods in the
 *  itktubeImageToTubeRigidMetric class. The distance between
 *  a 3D binary images ( 32x32x32 ) and a .tre image is computed.
 */

int itktubeSpatialObjectToImageMetricPerformanceTest( int argc, char * argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << " Input_FixedImage "
              << "Input_SpatialObject "
              << "Output_Results"
              << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<double, 3>                             Image3DType;

  typedef itk::ImageFileReader<Image3DType>                 ImageReaderType;
  typedef itk::SpatialObjectReader<3>                       SOReaderType;

  typedef itk::tube::PointBasedSpatialObjectToImageMetric<3, Image3DType >
                                                            MetricType;
  typedef itk::ComposeScaleSkewVersor3DTransform<double>    TransformType;

  // read image ( fixedImage )
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( argv[1] );
  try
    {
    imageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // read tube ( spatialObject )
  SOReaderType::Pointer tubeReader = SOReaderType::New();
  tubeReader->SetFileName( argv[2] );
  try
    {
    tubeReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // subsample points in vessel
  typedef itk::tube::SubSampleSpatialObjectFilter<>
    SubSampleFilterType;
  SubSampleFilterType::Pointer subSampleFilter =
    SubSampleFilterType::New();
  subSampleFilter->SetInput( tubeReader->GetGroup() );
  subSampleFilter->SetSampling( 30 );
  try
    {
    subSampleFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  //------------------------------------------------------------------
  // Compute the metric for a 3D image susampled at 1/30
  //------------------------------------------------------------------
  MetricType::Pointer metric = MetricType::New();
  metric->SetExtent( 3 );

  TransformType::Pointer transform = TransformType::New();
  TransformType::ParametersType parameters = transform->GetParameters();

  metric->SetFixedImage( imageReader->GetOutput() );
  metric->SetMovingSpatialObject ( subSampleFilter->GetOutput() );
  metric->SetTransform( transform.GetPointer() );

  // Add a time probe
  itk::TimeProbesCollectorBase chronometer;
  itk::MemoryProbesCollectorBase memorymeter;

  // Create stream to record the measure
  std::ofstream measuresFile;
  measuresFile.open( argv[3] );
  if( !measuresFile.is_open() )
    {
    std::cerr << "Unable to open: " << argv[3] << std::endl;
    return EXIT_FAILURE;
    }

  try
    {
    memorymeter.Start( "MetricComputation" );
    chronometer.Start( "MetricComputation" );

    metric->Initialize();
    MetricType::MeasureType value = metric->GetValue( parameters );

    memorymeter.Stop( "MetricComputation" );
    chronometer.Stop( "MetricComputation" );

    // Report the time and memory taken by the registration
    chronometer.Report( measuresFile );
    memorymeter.Report( measuresFile );

    std::cout << "Metric value: " << value << std::endl;
    }
  catch( itk::ExceptionObject &excp )
    {
    std::cerr << "Exception caught while initializing metric." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  measuresFile.close();
  return EXIT_SUCCESS;
}
