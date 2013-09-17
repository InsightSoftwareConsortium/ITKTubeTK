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

#include "itktubeTubeAngleOfIncidenceWeightFunction.h"
#include "itktubeTubeAngleOfIncidenceWeightFunctionSerializer.h"

#include <itkJsonCppArchiver.h>
#include <itkGroupSpatialObject.h>
#include <itkVesselTubeSpatialObject.h>

int itktubeTubeAngleOfIncidenceWeightFunctionSerializerTest( int argc, char * argv[] )
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

  typedef itk::VesselTubeSpatialObject< Dimension > TubeSpatialObjectType;
  typedef float                                     WeightType;

  typedef itk::tube::Function::TubeAngleOfIncidenceWeightFunction<
    TubeSpatialObjectType::TubePointType, WeightType > WeightFunctionType;
  WeightFunctionType::Pointer weightFunction = WeightFunctionType::New();

  typedef itk::tube::TubeAngleOfIncidenceWeightFunctionSerializer<
    WeightFunctionType >
      SerializerType;
  SerializerType::Pointer serializer = SerializerType::New();

  weightFunction->SetFractionalImportance( 0.8 );
  weightFunction->SetAngleDependence( 2.5 );
  WeightFunctionType::PointType probeOrigin;
  probeOrigin[0] = 354.8;
  probeOrigin[1] = -214.8;
  weightFunction->SetUltrasoundProbeOrigin( probeOrigin );
  serializer->SetTargetObject( weightFunction );
  serializer->Serialize();

  std::cout << serializer << std::endl;

  itk::JsonCppArchiver::Pointer archiver =
    dynamic_cast< itk::JsonCppArchiver * >( serializer->GetArchiver() );
  archiver->WriteToFile( archiveFileName );

  archiver->ReadFromFile( archiveFileName );
  weightFunction->SetFractionalImportance( 0.5 );
  weightFunction->SetAngleDependence( 0.5 );
  probeOrigin[0] = 1.0;
  probeOrigin[1] = 1.0;
  probeOrigin[2] = 1.0;
  serializer->DeSerialize();
  if( weightFunction->GetFractionalImportance() != 0.8 ||
      weightFunction->GetAngleDependence() != 2.5 ||
      weightFunction->GetUltrasoundProbeOrigin()[0] != 354.8 ||
      weightFunction->GetUltrasoundProbeOrigin()[1] != -214.8 )
    {
    std::cerr << "DeSerialization did not occur correctly: "
      << weightFunction << std::endl;
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
