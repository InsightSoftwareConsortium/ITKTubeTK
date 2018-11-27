/*=========================================================================

Library:   TubeTKLib

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

#include "tubetklibConfigure.h"

#include "tubeTestMain.h"

void RegisterTests( void )
{
  REGISTER_TEST( tubeBaseIOPrintTest );
  REGISTER_TEST( itktubePDFSegmenterParzenIOTest );
#ifdef TubeTKLib_USE_LIBSVM
  REGISTER_TEST( itktubePDFSegmenterSVMIOTest );
#endif
#ifdef TubeTKLib_USE_RANDOMFOREST
  REGISTER_TEST( itktubePDFSegmenterRandomForestIOTest );
#endif
  REGISTER_TEST( itktubeRidgeSeedFilterIOTest );
  REGISTER_TEST( itktubeTubeExtractorIOTest );
  REGISTER_TEST( itktubeTubeXIOTest );
}
