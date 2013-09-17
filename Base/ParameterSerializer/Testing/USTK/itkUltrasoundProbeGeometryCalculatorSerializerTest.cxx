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

#include "itkUltrasoundProbeGeometryCalculator.h"
#include "itkUltrasoundProbeGeometryCalculatorSerializer.h"

#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkJsonCppArchiver.h>

int itkUltrasoundProbeGeometryCalculatorSerializerTest( int argc, char * argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: "
              << argv[0]
              << " archiveFileName"
              << std::endl;
    return EXIT_FAILURE;
    }
  const char * archiveFileName = argv[1];
  // Types
  enum { Dimension = 2 };

  typedef float                              PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;

  // Create an UltrasoundProbeGeometryCalculator
  typedef itk::tube::UltrasoundProbeGeometryCalculator< ImageType >
    UltrasoundProbeGeometryCalculatorType;
  UltrasoundProbeGeometryCalculatorType::Pointer geometryCalculator =
    UltrasoundProbeGeometryCalculatorType::New();

  typedef itk::tube::UltrasoundProbeGeometryCalculatorSerializer<
    UltrasoundProbeGeometryCalculatorType >
      SerializerType;
  SerializerType::Pointer serializer = SerializerType::New();

  unsigned int beamDirection = 0;
  geometryCalculator->SetGeneralBeamDirection( beamDirection );
  PixelType backgroundValue = 22;
  geometryCalculator->SetBackgroundValue( backgroundValue );
  serializer->SetTargetObject( geometryCalculator );
  serializer->Serialize();

  std::cout << serializer << std::endl;

  itk::JsonCppArchiver::Pointer archiver =
    dynamic_cast< itk::JsonCppArchiver * >( serializer->GetArchiver() );
  archiver->WriteToFile( archiveFileName );

  archiver->ReadFromFile( archiveFileName );
  geometryCalculator->SetGeneralBeamDirection( 1 );
  geometryCalculator->SetBackgroundValue( 5.0 );
  serializer->DeSerialize();
  if( geometryCalculator->GetGeneralBeamDirection() != 0 ||
      geometryCalculator->GetBackgroundValue() != 22 )
    {
    std::cerr << "DeSerialization did not occur correctly: "
      << geometryCalculator << std::endl;
    return EXIT_FAILURE;
    }

  std::ostringstream oStream;
  archiver->WriteToStdStream( oStream );

  std::istringstream iStream( oStream.str() );
  archiver->ReadFromStdStream( iStream );

  std::cout << "Writing to std::cout: " << std::endl;
  archiver->WriteToStdStream( std::cout );

  return EXIT_SUCCESS;
}
