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

#include "tubeMessage.h"
#include "../CLI/tubeCLIProgressReporter.h"

#include <itkSpatialObjectReader.h>
#include <itkGroupSpatialObject.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkTubeSpatialObject.h>

#include <tubeWriteTubesAsPolyData.h>

#include "ConvertTubesToSurfaceCLP.h"

/** Forward declarations */
int DoIt( int argc, char * argv[] );

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return DoIt( argc, argv );
}

/** Main work happens here */
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  const unsigned int Dimension = 3;

  typedef itk::TubeSpatialObject< Dimension >       TubeSpatialObjectType;
  typedef itk::GroupSpatialObject< Dimension >      GroupSpatialObjectType;

  typedef itk::SpatialObjectReader< Dimension >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputTubeFile );
  try
    {
    reader->Update();
    }
  catch( const std::exception & error )
    {
    tube::ErrorMessage( error.what() );
    return EXIT_FAILURE;
    }
  GroupSpatialObjectType::Pointer groupSpatialObject = reader->GetGroup();
  groupSpatialObject->Update();

  tube::WriteTubesAsPolyData::Pointer writer =
    tube::WriteTubesAsPolyData::New();

  writer->SetInput( groupSpatialObject );

  if( !outputCenterline.empty() )
    {
    writer->SetCenterlineFileName( outputCenterline );
    }

  writer->SetFileName( outputSurfaceFile );

  writer->Write();

  return EXIT_SUCCESS;
}
