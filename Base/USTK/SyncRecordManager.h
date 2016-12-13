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
/*! @file

@brief class for management of synchronized AIM->Kitware/Spotlight
  tracking + ultrasound data single-slice record lists

  Supports disk reading, writing and playing back of synchronous tracking +
  ultrasound data record lists initially saved asynchronously by AIM's
  Save_CMC branch, then synchronized and streamlined by the
  Replay_Reexport_CMC branch.
  Shared between AIM and Spotlight, and with Kitware.

  @note  Copyright ( c ) InnerOptic Technology, Inc.  All Rights Reserved.

  @author  Andrei State

*/
//=========================================================================
#ifndef __SyncRecordManager_h
#define __SyncRecordManager_h

#include <vector>
#include "SyncRecord.h"
using namespace std;

//=========================================================================
class SyncRecordManager{
public:
  SyncRecordManager( void );
  ~SyncRecordManager( void );

  bool load( const char *disk_load_path );

  //#   ifdef INNEROPTIC_INTERNAL_ONLY
  void setVolumeImagePath( const char *vol_path );
  void setTrackerFromVolumeImageMatrix( const double m[16] );
  void setTrackerFromNavelMatrix(       const double m[16] );
  SyncRecord *newSequentialAppendedRecord( void );
  bool dump( const char *disk_write_path, char *ident );
  //#   endif
  const char *getVolumeImagePath( void );
  void getTrackerFromVolumeImageMatrix( double m[16] );
  void getTrackerFromNavelMatrix( double m[16] );

  void printRecords( void );

  SyncRecord *getNextRecord( void );
  SyncRecord *getRecord( const int which );

  void rewind( void );
  size_t getNbRecords( void );

protected:
private:
  string                       m_Records_path;
  string                       m_Volume_path;
  vector<SyncRecord>           m_Recs;
  vector<SyncRecord>::iterator m_Curr_rec;
  double                       tracker_f_volume[16];
  double                       tracker_f_navel[16];
  unsigned int                 m_Count;

  void ostreamMatrix( ostream &os, const char *header, const double m[16] );
  bool istreamSkipPastEOL( istream &is, const int howMany );
  bool istreamMatrix( istream &is, const char *header, double m[16],
    const bool verbose );
};
#endif
//=========================================================================
// vim: set noexpandtab: REQUIRED BY ANDREI's EDITOR
