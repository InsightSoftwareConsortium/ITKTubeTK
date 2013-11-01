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

  // sequential access to records
  while( ( syncRecord = syncRecordManager.getNextRecord() ) )
    {
    const int time = syncRecord->getTimestamp();
    (void) time;
    // duration in msec since the start of InnerOptic's acquisition program (at acquisition time);
    // could be used to determine how much time had passed from one tracked image to the next.

    double xInRUF[MAX_US_SCAN_CROP_POLYGON_PTS];
    double yInRUF[MAX_US_SCAN_CROP_POLYGON_PTS];

    for( int ii = 0 ; ii < MAX_US_SCAN_CROP_POLYGON_PTS ; ++ii )
      {
      syncRecord->getScanCropVertex_in_ruf( ii, xInRUF[ii], yInRUF[ii] ); // check whether return value is true
      }
    // polygonal outline of actual ultrasound pixel scan data, in drift-corrected raw ultrasound frame pixels (x going right, y down)
    // only the data within this polygon should be reconstructed and used for registration with pre-op imagery
    // currently this polygon is a trapeze whose parallel edges are aligned with the x axis of the ultrasound image

    double transformationMatrix[16];
    syncRecord->getTrackerFromRufMatrix( transformationMatrix );
    // OpenGL-style orthogonal transformation matrix (translations are in elements 12, 13, 14) which premultiplies a vector in drift-corrected
    // raw ultrasound frame pixels (z_in_ruf set to 0, add homogeneous component 1) to a vector in tracker space (tabletop field generator)

    // InnerOptic has performed a manual, approximate, constant-time-offset synchronization of these tracking matrices with the ultrasound images

    unsigned char *rgbRUFPixels = syncRecord->loadRawRgbPixels();
    (void) rgbRUFPixels;
    // points to 24-bit raw ultrasound frame pixels as described above; check for non-NULL return value
    // for now all ultrasound scan files are at 960x768 resolution, so this pointer currently points to 3*960*768 Bytes

    syncRecord->unloadRawRgbPixels();
    // in a 64-bit program you can keep thousands of images in memory, otherwise you should deallocate (which invalidates the pointer)
    }


  // in case you want to start another sequential access pass
  syncRecordManager.rewind();

  return EXIT_SUCCESS;
}
