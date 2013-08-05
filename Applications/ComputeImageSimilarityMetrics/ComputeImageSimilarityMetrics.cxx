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

#include <itkIdentityTransform.h>
#include <itkImageFileReader.h>
#include <itkMutualInformationImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkNormalizeImageFilter.h>

#include "ComputeImageSimilarityMetricsCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef TPixel                                              PixelType;
  typedef itk::Image< PixelType, VDimension >                 ImageType;
  typedef itk::ImageFileReader< ImageType >                   ReaderType;

  typename ReaderType::Pointer reader1 = ReaderType::New();
  typename ReaderType::Pointer reader2 = ReaderType::New();

  //read input image
  reader1->SetFileName( inputVolume1.c_str() );
  reader2->SetFileName( inputVolume2.c_str() );

  try
    {
    reader1->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught, image reader 1 !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  try
    {
    reader2->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught, image reader2 !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  typename ImageType::Pointer image1 = reader1->GetOutput();
  typename ImageType::Pointer image2 = reader2->GetOutput();

  typedef itk::NormalizeImageFilter< ImageType, ImageType > NormFilterType;

  typename NormFilterType::Pointer norm1 = NormFilterType::New();
  norm1->SetInput( image1 );
  norm1->Update();
  image1 = norm1->GetOutput();

  typename NormFilterType::Pointer norm2 = NormFilterType::New();
  norm2->SetInput( image2 );
  norm2->Update();
  image2 = norm2->GetOutput();

  typedef itk::IdentityTransform< double, VDimension >
    TransformType;
  typename TransformType::Pointer transform = TransformType::New();

  typedef itk::LinearInterpolateImageFunction< ImageType, double >
    InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( image2 );

  typedef itk::ImageToImageMetric< ImageType, ImageType >
    MetricType;
  typename MetricType::Pointer metric;

  if( !correlation )
    {
    typedef itk::MutualInformationImageToImageMetric< ImageType, ImageType >
                                                           MIMetricType;
    metric = MIMetricType::New();
    }
  else
    {
    typedef itk::NormalizedCorrelationImageToImageMetric< ImageType,
                                                           ImageType >
                                                           CorMetricType;
    metric = CorMetricType::New();
    }

  typename ImageType::SizeType size = image1->GetLargestPossibleRegion().
                                              GetSize();

  metric->SetFixedImage( image1 );
  metric->SetMovingImage( image2 );
  metric->SetFixedImageRegion( image1->GetLargestPossibleRegion() );
  metric->SetTransform( transform );
  metric->SetInterpolator( interpolator );
  metric->SetNumberOfSpatialSamples( size[0]*size[1]*samplingRate );
  metric->Initialize();
  metric->MultiThreadingInitialize();

  if( !correlation )
    {
    std::cout << metric->GetValue( transform->GetParameters() )
              << std::endl;
    }
  else
    {
    std::cout << -metric->GetValue( transform->GetParameters() )
              << std::endl;
    }

  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( inputVolume1, argc, argv );
}
