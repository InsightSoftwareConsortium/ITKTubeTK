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

#include <cstdlib>
#include <iostream>

#include "SyncRecordManager.h"

int SyncRecordTest( int argc, char * argv [] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: "
              << argv[0]
              << " metadata.txt" // :-)
              << std::endl;
    return EXIT_FAILURE;
    }
  const char * metadata = argv[1];

  SyncRecordManager syncRecordManager;
  bool loaded = syncRecordManager.load( metadata );
  if( !loaded )
    {
    std::cerr << "Error during load metadata file." << std::endl;
    return EXIT_FAILURE;
    }
  // metadata is now in memory; check whether return value is true
  SyncRecord * syncRecord;

  return EXIT_SUCCESS;
}
