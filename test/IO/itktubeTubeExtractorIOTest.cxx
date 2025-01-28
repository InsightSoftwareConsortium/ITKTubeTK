/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#include "itktubeTubeExtractorIO.h"

int
itktubeTubeExtractorIOTest(int argc, char * argv[])
{
  if (argc != 4)
  {
    std::cout << "Missing arguments." << std::endl;
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0] << " input.mtp output1.mtp output2.mtp" << std::endl;
    return EXIT_FAILURE;
  }

  // Prep
  typedef itk::Image<float, 3> ImageType;

  // Declare the type for the Filter
  typedef itk::tube::TubeExtractorIO<ImageType> IOMethodType;

  // TubeExtractor must have been assigned an input image prior to
  //   reading parameters into it.
  ImageType::Pointer    image = ImageType::New();
  ImageType::RegionType region;
  ImageType::IndexType  indx;
  indx.Fill(0);
  ImageType::RegionType::SizeType size;
  size.Fill(100);
  ImageType::PointType origin;
  origin.Fill(0);
  ImageType::SpacingType spacing;
  spacing.Fill(1);
  region.SetIndex(indx);
  region.SetSize(size);
  image->SetOrigin(origin);
  image->SetRegions(region);
  image->SetSpacing(spacing);
  image->Allocate();
  image->FillBuffer(0);

  IOMethodType::TubeExtractorType::Pointer tubeExtractor = IOMethodType::TubeExtractorType::New();
  tubeExtractor->SetInputImage(image);

  IOMethodType ioMethod;
  ioMethod.SetTubeExtractor(tubeExtractor);

  if (!ioMethod.CanRead(argv[1]))
  {
    std::cout << "ERROR: CanRead returned false." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "CanRead." << std::endl;

  if (!ioMethod.Read(argv[1]))
  {
    std::cout << "ERROR: Read failed." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Read." << std::endl;

  IOMethodType ioMethod2;

  ioMethod2.SetTubeExtractor(tubeExtractor);
  std::cout << "SetTubeExtractor." << std::endl;

  if (!ioMethod2.Write(argv[2]))
  {
    std::cout << "ERROR: Write failed." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Write." << std::endl;

  IOMethodType ioMethod3;

  ioMethod3.CopyInfo(ioMethod);
  std::cout << "CopyInfo." << std::endl;

  if (!ioMethod3.Write(argv[3]))
  {
    std::cout << "Write file = " << argv[3] << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Write." << std::endl;

  ioMethod.PrintInfo();
  std::cout << "Print." << std::endl;

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
