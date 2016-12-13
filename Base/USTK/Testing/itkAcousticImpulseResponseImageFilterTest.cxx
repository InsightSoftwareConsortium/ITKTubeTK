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

#include "itkAcousticImpulseResponseImageFilter.h"
#include "itkGradientBasedAngleOfIncidenceImageFilter.h"

#include <itkAbsImageAdaptor.h>
#include <itkAddImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkIntensityWindowingImageFilter.h>
#include <itkLog10ImageAdaptor.h>

// Very rough B-Mode that does not take into account system properties or a
// a realistic image formation process.
template< class TInputImage >
void
RoughBMode( TInputImage * inputImage, const char * fileName )
{
  typedef TInputImage InputImageType;
  typedef itk::AbsImageAdaptor< InputImageType,
      typename InputImageType::PixelType >
    AbsAdaptorType;
  typename AbsAdaptorType::Pointer absAdaptor = AbsAdaptorType::New();
  absAdaptor->SetImage( inputImage );

  typedef itk::AddImageFilter<
      AbsAdaptorType,
      InputImageType,
      InputImageType >
    AddConstantFilterType;
  typename AddConstantFilterType::Pointer addConstantFilter =
    AddConstantFilterType::New();
  addConstantFilter->SetInput1( absAdaptor );
  addConstantFilter->SetConstant2( 1.0e-12 );
  addConstantFilter->Update();

  typedef itk::Log10ImageAdaptor< InputImageType,
      typename InputImageType::PixelType >
    Log10AdaptorType;
  typename Log10AdaptorType::Pointer log10Adaptor =
    Log10AdaptorType::New();
  log10Adaptor->SetImage( addConstantFilter->GetOutput() );

  typedef unsigned char OutputPixelType;
  typedef itk::Image< OutputPixelType, InputImageType::ImageDimension >
    OutputImageType;
  typedef itk::IntensityWindowingImageFilter< Log10AdaptorType, OutputImageType >
    IntensityWindowingFilterType;
  typename IntensityWindowingFilterType::Pointer
    intensityWindowingFilter = IntensityWindowingFilterType::New();
  intensityWindowingFilter->SetWindowMinimum( -5.0 );
  intensityWindowingFilter->SetWindowMaximum( 2.0 );
  intensityWindowingFilter->SetInput( log10Adaptor );

  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( intensityWindowingFilter->GetOutput() );
  writer->SetFileName( fileName );
  writer->Update();
}

int itkAcousticImpulseResponseImageFilterTest( int argc, char * argv[] )
{
  // Argument parsing.
  if( argc < 7 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: "
              << argv[0]
              << " inputImage"
              << " outputImage"
              << " outputImageWithGradientRecursiveGaussian"
              << " originX"
              << " originY"
              << " angleDependence"
              << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImage = argv[1];
  const char * outputImage = argv[2];
  const char * outputImageWithGradientRecursiveGaussian = argv[3];
  const char * originX = argv[4];
  const char * originY = argv[5];
  const char * angleDependence = argv[6];

  // Types
  enum { Dimension = 2 };

  typedef float                              PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;

  // Reader
  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImage );

  // Calculate the angle of incidence
  typedef itk::GradientBasedAngleOfIncidenceImageFilter< ImageType, ImageType >
    AngleOfIncidenceFilterType;
  AngleOfIncidenceFilterType::Pointer angleOfIncidenceFilter =
    AngleOfIncidenceFilterType::New();
  angleOfIncidenceFilter->SetInput( reader->GetOutput() );
  AngleOfIncidenceFilterType::OriginType probeOrigin;
  std::istringstream istrm;
  istrm.str( originX );
  istrm >> probeOrigin[0];
  istrm.clear();
  istrm.str( originY );
  istrm >> probeOrigin[1];
  angleOfIncidenceFilter->SetUltrasoundProbeOrigin( probeOrigin );

  // Calculate the acoustic impulse response
  typedef itk::AcousticImpulseResponseImageFilter< ImageType, ImageType >
    AcousticImpulseResponseFilterType;
  AcousticImpulseResponseFilterType::Pointer acousticImpulseResponseFilter =
    AcousticImpulseResponseFilterType::New();
  acousticImpulseResponseFilter->SetInput( 0,
    reader->GetOutput() );
  acousticImpulseResponseFilter->SetInput( 1,
    angleOfIncidenceFilter->GetOutput() );
  istrm.clear();
  istrm.str( angleDependence );
  double angleDependenceDouble;
  istrm >> angleDependenceDouble;
  acousticImpulseResponseFilter->SetAngleDependence( angleDependenceDouble );

  // Write image.
  try
    {
    acousticImpulseResponseFilter->Update();
    RoughBMode( acousticImpulseResponseFilter->GetOutput(), outputImage );
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  // Use a different gradient filter.
  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<
    AcousticImpulseResponseFilterType::OperatorImageType,
    AcousticImpulseResponseFilterType::OperatorImageType
      >
      GradientMagnitudeRecursiveGaussianFilterType;
  GradientMagnitudeRecursiveGaussianFilterType::Pointer
    gradientMagnitudeRecursiveGaussianFilter =
      GradientMagnitudeRecursiveGaussianFilterType::New();
  gradientMagnitudeRecursiveGaussianFilter->SetSigma( 1.0 );
  acousticImpulseResponseFilter->SetGradientMagnitudeFilter( 
    gradientMagnitudeRecursiveGaussianFilter );

  try
    {
    acousticImpulseResponseFilter->Update();
    RoughBMode( acousticImpulseResponseFilter->GetOutput(),
      outputImageWithGradientRecursiveGaussian );
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
