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

#include "itkUltrasoundProbeGeometryCalculator.h"

#include <itkImageFileReader.h>
#include <itkThresholdImageFilter.h>

int itkUltrasoundProbeGeometryCalculatorTest( int argc, char * argv[] )
{
  // Argument parsing.
  if( argc < 2 )
    {
    std::cerr << "Usage: "
              << argv[0]
              << " inputImage"
              << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImage = argv[1];

  // Types
  enum { Dimension = 2 };

  typedef unsigned char                      PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;

  // Reader
  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImage );

  // The test input image was captured with a screen grabber, and some analog
  // noise was introduced, just using the input image with a BackgroundValue
  // will not work.  So, we threshold the input.
  typedef itk::ThresholdImageFilter< ImageType > ThresholdFilterType;
  ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
  thresholdFilter->SetInput( reader->GetOutput() );
  thresholdFilter->ThresholdBelow( 5 );

  // Calculate the probe's geometry
  typedef itk::tube::UltrasoundProbeGeometryCalculator< ImageType >
    GeometryCalculatorType;
  GeometryCalculatorType::Pointer geometryCalculator = GeometryCalculatorType::New();
  geometryCalculator->SetInput( thresholdFilter->GetOutput() );
  // The probe is oriented in the "vertical" direction in the image
  geometryCalculator->SetGeneralBeamDirection( 1 );

  try
    {
    geometryCalculator->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  const GeometryCalculatorType::OriginType ultrasoundProbeOrigin =
    geometryCalculator->GetUltrasoundProbeOrigin();
  const GeometryCalculatorType::RadiusType startOfAcquisitionRadius =
    geometryCalculator->GetStartOfAcquisitionRadius();

  std::cout << "Probe Origin: " << ultrasoundProbeOrigin << std::endl;
  std::cout << "Start of Acquisition Radius: " << startOfAcquisitionRadius << std::endl;

  const double tolerance = 1.0;
  if( std::fabs( ultrasoundProbeOrigin[0] - 331.5 ) > tolerance ||
      std::fabs( ultrasoundProbeOrigin[1] - -528.4 ) >  tolerance )
    {
    std::cerr << "Did not find the correct probe origin!" << std::endl;
    return EXIT_FAILURE;
    }

  if( std::fabs( startOfAcquisitionRadius - 555.6 ) > tolerance )
    {
    std::cerr << "Did not find the correct start of acquisition radius!" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
