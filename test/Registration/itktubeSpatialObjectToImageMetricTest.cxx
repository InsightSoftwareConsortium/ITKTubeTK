/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#include "itktubePointBasedSpatialObjectToImageMetric.h"
#include "itktubeSubSampleSpatialObjectFilter.h"

#include <itkImageFileReader.h>
#include <itkSpatialObjectReader.h>
#include <itkComposeScaleSkewVersor3DTransform.h>

/**
 *  This test exercised the metric evaluation methods in the
 *  itktubeSpatialObjectToImageMetric class. The distance between
 *  a 3D binary images ( 32x32x32 ) and a .tre image is computed and check with
 *  the reference for the metric.
 */

int itktubeSpatialObjectToImageMetricTest( int argc, char * argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << " Input_FixedImage "
              << "Input_SpatialObject "
              << "Input_ExpectedValue."
              << std::endl;
    return EXIT_FAILURE;
    }

  typedef double            FloatType;
  static const unsigned int ImageDimension = 3;
  static const unsigned int ObjectDimension = 3;

  typedef itk::Image< FloatType, ImageDimension >         ImageType;
  typedef itk::GroupSpatialObject< ObjectDimension >      GroupType;

  typedef itk::ImageFileReader< ImageType >               ImageReaderType;
  typedef itk::SpatialObjectReader< ObjectDimension >     GroupReaderType;

  typedef itk::tube::PointBasedSpatialObjectToImageMetric< 3, ImageType >
                                                          MetricType;
  typedef itk::ComposeScaleSkewVersor3DTransform< double >       TransformType;

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

  // read group ( spatialObject )
  GroupReaderType::Pointer groupReader = GroupReaderType::New();
  groupReader->SetFileName( argv[2] );
  try
    {
    groupReader->Update();
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
  subSampleFilter->SetInput( groupReader->GetGroup() );
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
  metric->SetTransform( transform );
  try
    {
    metric->Initialize();
    }
  catch( itk::ExceptionObject &excp )
    {
    std::cerr << "Exception caught while initializing metric." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  const double epsilonReg = 0.05; // Delta threshold on the measure checking.
  MetricType::MeasureType value = metric->GetValue( parameters );
  if( value < ( std::atof( argv[3] ) - epsilonReg ) ||
      value > ( std::atof( argv[3] ) + epsilonReg ) )
    {
    std::cerr << "Distance value different than expected: "
              << value
              << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
