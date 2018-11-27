/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

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

#include "tubeTestMain.h"
#include "tubetklibConfigure.h"

#include <iostream>

void RegisterTests( void )
{
  REGISTER_TEST( tubeBaseSegmentationPrintTest );
  REGISTER_TEST( itktubePDFSegmenterParzenTest );
#ifdef TubeTKLib_USE_LIBSVM
  REGISTER_TEST( itktubePDFSegmenterSVMTest );
#endif
#ifdef TubeTKLib_USE_RANDOMFOREST
  REGISTER_TEST( itktubePDFSegmenterRandomForestTest );
#endif
  REGISTER_TEST( itktubeRidgeExtractorTest );
  REGISTER_TEST( itktubeRidgeExtractorTest2 );
  REGISTER_TEST( itktubeRidgeSeedFilterTest );
  REGISTER_TEST( itktubeRadiusExtractor2Test );
  REGISTER_TEST( itktubeRadiusExtractor2Test2 );
  REGISTER_TEST( itktubeTubeExtractorTest );
}
