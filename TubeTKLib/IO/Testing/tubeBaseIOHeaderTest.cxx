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

#include "tubetkConfigure.h"

#include "itktubePDFSegmenterParzenIO.h"
#include "itktubeRidgeSeedFilterIO.h"
#ifdef TubeTK_USE_LIBSVM
#include "itktubePDFSegmenterSVMIO.h"
#endif
#ifdef TubeTK_USE_RANDOMFOREST
#  include "itktubePDFSegmenterRandomForestIO.h"
#endif
#include "itktubeTubeExtractorIO.h"
#include "itktubeTubeXIO.h"

int main( int tubeNotUsed( argc ), char * tubeNotUsed( argv )[] )
{
  return EXIT_SUCCESS;
}
