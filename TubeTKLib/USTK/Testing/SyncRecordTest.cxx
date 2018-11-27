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

#include <cstdlib>
#include <iostream>

#include "SyncRecordManager.h"

int SyncRecordTest( int argc, char * argv [] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: "
              << argv[0]
              << " metadata.txt" // :- )
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

  const size_t numberOfRecords = syncRecordManager.getNbRecords();
  std::cout << "Number of Records: " << numberOfRecords << std::endl;
  if( numberOfRecords != 3 )
    {
    std::cerr << "Did not get the expected number of records." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Volume image path: " << syncRecordManager.getVolumeImagePath() << std::endl;
  double transformationMatrix[16];
  // exercise methods
  syncRecordManager.getTrackerFromVolumeImageMatrix( transformationMatrix );
  std::cout << "\nTrackerFromVolumeImageMatrix:";
  for( unsigned int ii = 0; ii < 4; ++ii )
    {
    std::cout << std::endl;
    for( unsigned int jj = 0; jj < 4; ++jj )
      {
      std::cout << transformationMatrix[ii * 4 + jj] << " ";
      }
    }
  syncRecordManager.getTrackerFromNavelMatrix( transformationMatrix );
  std::cout << "\n\nTrackerFromNavelMatrix:";
  for( unsigned int ii = 0; ii < 4; ++ii )
    {
    std::cout << std::endl;
    for( unsigned int jj = 0; jj < 4; ++jj )
      {
      std::cout << transformationMatrix[ii * 4 + jj] << " ";
      }
    }
  std::cout << "\n" << std::endl;

  // metadata is now in memory; check whether return value is true
  SyncRecord * syncRecord;

  // sequential access to records
  while( ( syncRecord = syncRecordManager.getNextRecord() ) )
    {
    std::cout << "Record filepath: " << syncRecord->getRufImageFilePath() << std::endl;;
    const int time = syncRecord->getTimestamp();
    ( void ) time;
    // duration in msec since the start of InnerOptic's acquisition program ( at acquisition time );
    // could be used to determine how much time had passed from one tracked image to the next.

    double xInRUF[MAX_US_SCAN_CROP_POLYGON_PTS];
    double yInRUF[MAX_US_SCAN_CROP_POLYGON_PTS];

    for( int ii = 0 ; ii < MAX_US_SCAN_CROP_POLYGON_PTS ; ++ii )
      {
      syncRecord->getScanCropVertex_in_ruf( ii, xInRUF[ii], yInRUF[ii] ); // check whether return value is true
      }
    // polygonal outline of actual ultrasound pixel scan data, in drift-corrected raw ultrasound frame pixels ( x going right, y down )
    // only the data within this polygon should be reconstructed and used for registration with pre-op imagery
    // currently this polygon is a trapeze whose parallel edges are aligned with the x axis of the ultrasound image

    syncRecord->getTrackerFromRufMatrix( transformationMatrix );
    // OpenGL-style orthogonal transformation matrix ( translations are in elements 12, 13, 14 ) which premultiplies a vector in drift-corrected
    // raw ultrasound frame pixels ( z_in_ruf set to 0, add homogeneous component 1 ) to a vector in tracker space ( tabletop field generator )

    // InnerOptic has performed a manual, approximate, constant-time-offset synchronization of these tracking matrices with the ultrasound images
    //

    syncRecord->getTransducerFromRufMatrix( transformationMatrix );

    unsigned char *rgbRUFPixels = syncRecord->loadRawRgbPixels();
    ( void ) rgbRUFPixels;
    // points to 24-bit raw ultrasound frame pixels as described above; check for non-NULL return value
    // for now all ultrasound scan files are at 960x768 resolution, so this pointer currently points to 3*960*768 Bytes

    syncRecord->unloadRawRgbPixels();
    // in a 64-bit program you can keep thousands of images in memory, otherwise you should deallocate ( which invalidates the pointer )
    }


  // in case you want to start another sequential access pass
  syncRecordManager.rewind();

  syncRecordManager.printRecords();

  // These are "INNEROPTIC_INTERNAL_ONLY" ( ? )
  // Just exercise for now.
  syncRecordManager.setVolumeImagePath( "newVolumeImagePath" );
  syncRecordManager.setTrackerFromVolumeImageMatrix( transformationMatrix );
  syncRecordManager.setTrackerFromNavelMatrix( transformationMatrix );
  syncRecord = syncRecordManager.newSequentialAppendedRecord();
  syncRecord->setTimestamp( 15 );
  syncRecord->setRufImageFilePath( "SmoothPath" );
  syncRecord->setRufImageFileIndex( 10 );
  syncRecord->setScanCropVertex_in_ruf( 0, 3.0, 4.0 );
  syncRecord->setTrackerFromRufMatrix( transformationMatrix );
  syncRecord->setTransducerFromRufMatrix( transformationMatrix );
  syncRecord->setEndoImageGeometry_in_ruf( 0, 0, 0, 0 );

  return EXIT_SUCCESS;
}
