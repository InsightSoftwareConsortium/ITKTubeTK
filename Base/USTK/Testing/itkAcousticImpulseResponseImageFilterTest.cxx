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

#include "itkAcousticImpulseResponseImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkGradientBasedAngleOfIncidenceImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"

#include <sstream>

int itkAcousticImpulseResponseImageFilterTest( int argc, char* argv [] )
{
  // Argument parsing.
  if( argc < 6 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: "
              << argv[0]
              << " inputImage"
              << " outputImage"
              << " outputImageWithGradientRecursiveGaussian"
              << " originX"
              << " originY"
              << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImage = argv[1];
  const char * outputImage = argv[2];
  const char * outputImageWithGradientRecursiveGaussian = argv[3];
  const char * originX = argv[4];
  const char * originY = argv[5];

  // Types
  static const unsigned int Dimension = 2;

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
  acousticImpulseResponseFilter->SetInput( 0, reader->GetOutput() );
  acousticImpulseResponseFilter->SetInput( 1, angleOfIncidenceFilter->GetOutput() );

  // Writer
  typedef itk::ImageFileWriter< ImageType > ImageWriterType;
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputImage );
  writer->SetInput( acousticImpulseResponseFilter->GetOutput() );
  try
    {
    writer->Update();
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
  gradientMagnitudeRecursiveGaussianFilter->SetSigma( 0.5 );
  acousticImpulseResponseFilter->SetGradientMagnitudeFilter( gradientMagnitudeRecursiveGaussianFilter );

  writer->SetFileName( outputImageWithGradientRecursiveGaussian );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
