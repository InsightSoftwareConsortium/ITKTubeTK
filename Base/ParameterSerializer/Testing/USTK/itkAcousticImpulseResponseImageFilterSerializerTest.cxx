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
#include "itkAcousticImpulseResponseImageFilterSerializer.h"

#include "itkJsonCppArchiver.h"

#include "itkAcousticImpulseResponseImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"

int itkAcousticImpulseResponseImageFilterSerializerTest( int argc, char * argv[] )
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
  static const unsigned int Dimension = 2;

  typedef float                              PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;

  // Create an AcousticImpulseResponseImageFilter
  typedef itk::AcousticImpulseResponseImageFilter< ImageType, ImageType >
    AcousticImpulseResponseFilterType;
  AcousticImpulseResponseFilterType::Pointer acousticImpulseResponseFilter =
    AcousticImpulseResponseFilterType::New();

  typedef itk::tube::AcousticImpulseResponseImageFilterSerializer<
    AcousticImpulseResponseFilterType >
      SerializerType;
  SerializerType::Pointer serializer = SerializerType::New();

  double angleDependence = 0.8;
  acousticImpulseResponseFilter->SetAngleDependence( angleDependence );
  serializer->SetTargetObject( acousticImpulseResponseFilter );
  serializer->Serialize();

  std::cout << serializer << std::endl;

  itk::JsonCppArchiver::Pointer archiver =
    dynamic_cast< itk::JsonCppArchiver * >( serializer->GetArchiver() );
  archiver->WriteToFile( archiveFileName );

  archiver->ReadFromFile( archiveFileName );
  acousticImpulseResponseFilter->SetAngleDependence( 0.2 );
  serializer->DeSerialize();
  if( acousticImpulseResponseFilter->GetAngleDependence() != 0.8 )
    {
    std::cerr << "DeSerialization did not occur correctly: "
      << acousticImpulseResponseFilter->GetAngleDependence() << std::endl;
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
  gradientMagnitudeRecursiveGaussianFilter->SetSigma( 2.3 );
  acousticImpulseResponseFilter->SetGradientMagnitudeFilter(
    gradientMagnitudeRecursiveGaussianFilter );

  serializer->Serialize();
  std::cout << "Serializer with recursive Gaussian: "
            << serializer << std::endl;

  archiver->WriteToFile( archiveFileName );
  archiver->ReadFromFile( archiveFileName );
  serializer->DeSerialize();

  // Switch the gradient filter again.
  typedef itk::GradientMagnitudeImageFilter<
    AcousticImpulseResponseFilterType::OperatorImageType,
    AcousticImpulseResponseFilterType::OperatorImageType
      >
      GradientMagnitudeImageFilterType;
  GradientMagnitudeImageFilterType::Pointer
    gradientMagnitudeImageFilter =
      GradientMagnitudeImageFilterType::New();
  acousticImpulseResponseFilter->SetGradientMagnitudeFilter(
    gradientMagnitudeImageFilter );

  serializer->Serialize();
  std::cout << "Serializer with GradientMagnitudeImageFilter: "
             << serializer << std::endl;

  acousticImpulseResponseFilter->SetGradientMagnitudeFilter(
    gradientMagnitudeRecursiveGaussianFilter );
  archiver->ReadFromFile( archiveFileName );
  serializer->DeSerialize();
  std::cout << "Serializer after revert: "
             << serializer << std::endl;

  std::ostringstream oStream;
  archiver->WriteToStdStream( oStream );

  std::istringstream iStream( oStream.str() );
  archiver->ReadFromStdStream( iStream );

  std::cout << "Writing to std::cout: " << std::endl;
  archiver->WriteToStdStream( std::cout );

  return EXIT_SUCCESS;
}
