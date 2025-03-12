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

#include <cstdlib>

#include "itktubeMetaTubeExtractor.h"

int
itktubeMetaTubeExtractorTest(int argc, char * argv[])
{
  if (argc != 2)
  {
    std::cout << "Usage: testname <tempfilename>" << std::endl;
    return EXIT_FAILURE;
  }

  vnl_vector<double> scales(3);
  scales[0] = 1;
  scales[1] = 2;
  scales[2] = 3;

  itk::tube::MetaTubeExtractor mtp1;
  /*
  mtp1.SetGeneralProperties( scales, 1024, 2048, 0 );
  if( mtp1.GetSeedScales() != scales
    || mtp1.GetSeedIntensityMin() != 1024
    || mtp1.GetSeedIntensityMax() != 2048
    || mtp1.GetSeedIntensityPercentile() != 0 )
    {
    std::cout << "Tube parameter seed values do not match after set"
      << std::endl;
    return EXIT_FAILURE;
    }
    */

  itk::tube::MetaTubeExtractor mtp2(mtp1);
  /*
  if( mtp2.GetSeedScales() != scales
    || mtp2.GetSeedIntensityMin() != 1024
    || mtp2.GetSeedIntensityMax() != 2048
    || mtp2.GetSeedIntensityPercentile() != 0 )
    {
    std::cout << "Tube param seed values don't match after copy constructor"
      << std::endl;
    mtp2.PrintInfo();
    return EXIT_FAILURE;
    }
    */

  mtp1.Write(argv[1]);

  itk::tube::MetaTubeExtractor mtp3(argv[1]);
  /*
  if( mtp3.GetSeedScales() != scales
    || mtp3.GetSeedIntensityMin() != 1024
    || mtp3.GetSeedIntensityMax() != 2048
    || mtp3.GetSeedIntensityPercentile() != 0 )
    {
    std::cout
      << "Tube param seed values don't match after write/read constructor"
      << std::endl;
    return EXIT_FAILURE;
    }
    */

  mtp1.Clear();
  /*
  if( mtp1.GetSeedScales().size() != 1 )
    {
    std::cout << "Tube params seed.size not 1 after clear." << std::endl;
    return EXIT_FAILURE;
    }
    */

  mtp1.Read(argv[1]);
  /*
  if( mtp1.GetSeedScales() != scales
    || mtp1.GetSeedIntensityMin() != 1024
    || mtp1.GetSeedIntensityMax() != 2048
    || mtp1.GetSeedIntensityPercentile() != 0 )
    {
    std::cout << "Tube params values don't match after read"
      << std::endl;
    return EXIT_FAILURE;
    }
    */

  return EXIT_SUCCESS;
}
