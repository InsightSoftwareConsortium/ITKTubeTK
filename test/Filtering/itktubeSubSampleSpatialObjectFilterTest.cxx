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

#include "itktubeSubSampleSpatialObjectFilter.h"

#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>

int itktubeSubSampleSpatialObjectFilterTest( int argc, char * argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: "
      << "inputTubeNetwork "
      << "outputTubeNetwork "
      << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputTubeNetwork = argv[1];
  const char * outputTubeNetwork = argv[2];

  enum { Dimension = 3 };
  typedef itk::TubeSpatialObject< Dimension >   TubeSpatialObjectType;
  typedef itk::GroupSpatialObject< Dimension >  GroupSpatialObjectType;

  // Read input tube tree.
  typedef itk::SpatialObjectReader< Dimension >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputTubeNetwork );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
  GroupSpatialObjectType::Pointer groupSpatialObject = reader->GetGroup();
  std::cout << "Number of children = "
    << groupSpatialObject->GetNumberOfChildren()
    << std::endl;

  // Sub-sample the tube tree.
  typedef itk::tube::SubSampleSpatialObjectFilter<>
      SubSampleFilterType;
  SubSampleFilterType::Pointer subSampleFilter =
    SubSampleFilterType::New();
  subSampleFilter->SetInput( reader->GetGroup() );

  const unsigned int sampling = 100;
  subSampleFilter->SetSampling( sampling );
  if( subSampleFilter->GetSampling() != sampling )
    {
    std::cerr << "Sampling did not get set correctly." << std::endl;
    return EXIT_FAILURE;
    }

  // Write output tube tree.
  typedef itk::SpatialObjectWriter< Dimension > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputTubeNetwork );
  writer->SetInput( subSampleFilter->GetOutput() );
  try
    {
    // Currently, there is a bug in the SpatialObjectWriter that it does not do
    // a pipeline update on its inputs.
    subSampleFilter->Update();
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
