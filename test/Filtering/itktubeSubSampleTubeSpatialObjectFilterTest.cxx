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

#include "itktubeSubSampleTubeSpatialObjectFilter.h"

#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>

int
itktubeSubSampleTubeSpatialObjectFilterTest(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cerr << "Usage: "
              << "inputTubeNetwork "
              << "outputTubeNetwork " << std::endl;
    return EXIT_FAILURE;
  }
  const char * inputTubeNetwork = argv[1];
  const char * outputTubeNetwork = argv[2];

  enum
  {
    ObjectDimension = 3
  };
  typedef itk::TubeSpatialObject<ObjectDimension>  TubeSpatialObjectType;
  typedef itk::GroupSpatialObject<ObjectDimension> GroupSpatialObjectType;

  typedef itk::SpatialObjectReader<ObjectDimension> ReaderType;
  ReaderType::Pointer                               reader = ReaderType::New();
  reader->SetFileName(inputTubeNetwork);
  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject & error)
  {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
  }
  GroupSpatialObjectType::Pointer groupSpatialObject = reader->GetGroup();
  std::cout << "Number of children = " << groupSpatialObject->GetNumberOfChildren() << std::endl;

  GroupSpatialObjectType::Pointer            output = GroupSpatialObjectType::New();
  char                                       childName[] = "Tube";
  GroupSpatialObjectType::ChildrenListType * children = groupSpatialObject->GetChildren(1000, childName);
  typedef GroupSpatialObjectType::ChildrenListType::iterator GroupChildrenIteratorType;
  for (GroupChildrenIteratorType it = children->begin(); it != children->end(); ++it)
  {
    TubeSpatialObjectType * pointBasedSpatialObject = dynamic_cast<TubeSpatialObjectType *>((*it).GetPointer());
    if (pointBasedSpatialObject)
    {
      typedef itk::tube::SubSampleTubeSpatialObjectFilter<ObjectDimension> SubSampleFilterType;
      SubSampleFilterType::Pointer                                         subSampleFilter = SubSampleFilterType::New();
      const unsigned int                                                   samplingFactor = 5;
      subSampleFilter->SetSampling(samplingFactor);
      subSampleFilter->SetInput(pointBasedSpatialObject);
      try
      {
        subSampleFilter->Update();
      }
      catch (itk::ExceptionObject & error)
      {
        std::cerr << "Error: " << error << std::endl;
        delete children;
        return EXIT_FAILURE;
      }
      output->AddChild(subSampleFilter->GetOutput());
    }
  }


  typedef itk::SpatialObjectWriter<ObjectDimension> WriterType;
  WriterType::Pointer                               writer = WriterType::New();
  writer->SetFileName(outputTubeNetwork);
  writer->SetInput(output);
  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & error)
  {
    std::cerr << "Error: " << error << std::endl;
    delete children;
    return EXIT_FAILURE;
  }

  delete children;
  return EXIT_SUCCESS;
}
