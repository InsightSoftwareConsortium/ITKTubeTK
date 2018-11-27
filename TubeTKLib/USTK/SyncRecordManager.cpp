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
/*! @file

  @brief class for management of synchronized AIM->Kitware/Spotlight
  tracking + ultrasound data single-slice record lists

   Supports disk reading, writing and playing back of synchronous tracking
   + ultrasound data record lists initially saved asynchronously by AIM's
   Save_CMC branch, then synchronized and streamlined by the
   Replay_Reexport_CMC branch.
   Shared between AIM and Spotlight, and with Kitware.

    @note  Copyright (c) InnerOptic Technology, Inc.  All Rights Reserved.

    @author  Andrei State

*/
//=========================================================================
#define _CRT_SECURE_NO_WARNINGS

/*#include <assert.h>*/
#include <stdio.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <cctype>

#include "SyncRecordManager.h"
//=========================================================================
///
SyncRecordManager::SyncRecordManager( void )
{
  m_Volume_path.assign( "VOID_VOLUME_PATH" );

  for( int i = 0 ; i < 16 ; ++i )
    {
    tracker_f_volume[i] = tracker_f_navel[i] = 0.;
    }
  m_Count = 0; rewind();
}
//-------------------------------------------------------------------------
///
SyncRecordManager::~SyncRecordManager( void )
{
  m_Recs.resize( 0 );
}
//-------------------------------------------------------------------------
///
void SyncRecordManager::ostreamMatrix( ostream &os, const char *header,
  const double m[16] )
{
  os << header << endl;
  for( int i = 0 ; i < 4 ; ++i )
    {
    for( int j = 0 ; j < 4 ; ++j )
      {
      os << m[4*i+j]; if( j < 3 ) os << " "; else os << endl;
      }
    }
}
//=========================================================================
//#ifdef INNEROPTIC_INTERNAL_ONLY
//-------------------------------------------------------------------------
///
void SyncRecordManager::setVolumeImagePath( const char *vol_path )
{
  m_Volume_path.assign( vol_path );
}
//-------------------------------------------------------------------------
///
void SyncRecordManager::setTrackerFromVolumeImageMatrix( const double m[16] )
{
  for( int i = 0 ; i < 16 ; ++i )
    {
    tracker_f_volume[i] = m[i];
    }
}
//-------------------------------------------------------------------------
///
void SyncRecordManager::setTrackerFromNavelMatrix(  const double m[16] )
{
  for( int i = 0 ; i < 16 ; ++i ) tracker_f_navel[i] = m[i];
}
//-------------------------------------------------------------------------
///
SyncRecord *SyncRecordManager::newSequentialAppendedRecord( void )
{
  ++m_Count;
  if( m_Count >= m_Recs.size() )
    {
    m_Recs.resize( m_Count );
    }
  return &m_Recs[m_Count-1];
}
//-------------------------------------------------------------------------
///
bool SyncRecordManager::dump( const char *disk_write_path, char *ident )
{
  if( m_Recs.empty() )
    {
    cerr << "SyncRecordManager::dump(): attempt to dump empty \
      SyncRecordManager object" << endl;
    return false;
    }

  ofstream ofs( disk_write_path );
  if( !ofs || !ofs.is_open() )
    {
    cerr << "SyncRecordManager::dump(): file " << disk_write_path <<
      " is not writable" << endl;
    return false;
    }

  // write a few leading comment lines starting with #, and empty
  ofs << "# Synchronized ultrasound scan / tracking data acquired and \
    processed by InnerOptic" << endl;
  ofs << "# Additional info: " << ident << endl;
  ofs << endl;
  ofs << "# Lines at the beginning of this file that are either empty or \
    start with # are considered comments" << endl;
  ofs << endl;
  ofs << "# The next 11 non-empty lines represent the volume data path and \
    the two navel registration matrices" << endl;
  ofs << "# Then follow the column title line and the actual bulk data \
    (1 line per tracked ultrasound slice), total " << m_Recs.size() <<
    " single-slice records" << endl;
  ofs << endl;

  ofs << m_Volume_path << endl;

  // navel registration matrices
  ostreamMatrix( ofs, "tracker_f_volume", tracker_f_volume );
  ostreamMatrix( ofs, "tracker_f_navel" , tracker_f_navel  );

  // column titles line
  ofs << SYNCRECORD_DUMP_ORDER << endl;

  // dump all records line by line
  for( vector<SyncRecord>::iterator i = m_Recs.begin(); i != m_Recs.end(); ++i )
    {
    (*i).dump( ofs );
    }

  ofs.close();

  return true;
}
//-------------------------------------------------------------------------
//#endif // INNEROPTIC_INTERNAL_ONLY
//=========================================================================
const char *SyncRecordManager::getVolumeImagePath( void )
{
  return m_Volume_path.c_str();
}
//-------------------------------------------------------------------------
///
void SyncRecordManager::getTrackerFromVolumeImageMatrix( double m[16] )
{
  for( int i = 0 ; i < 16 ; ++i )
    {
    m[i] = tracker_f_volume[i];
    }
}
//-------------------------------------------------------------------------
///
void SyncRecordManager::getTrackerFromNavelMatrix( double m[16] )
{
  for( int i = 0 ; i < 16 ; ++i )
    {
    m[i] = tracker_f_navel[i];
    }
}
//-------------------------------------------------------------------------
///
void SyncRecordManager::rewind( void )
{
  m_Curr_rec = m_Recs.begin();
}
//-------------------------------------------------------------------------
///
size_t SyncRecordManager::getNbRecords( void )
{
  return m_Recs.size();
}
//-------------------------------------------------------------------------
///
SyncRecord *SyncRecordManager::getNextRecord( void )
{
  if( m_Recs.empty() )
    { 
    cerr << "SyncRecordManager::getNextRecord(): manager object is empty" << endl;
    return NULL;
    }
  if( m_Curr_rec == m_Recs.end() )
    {
    return NULL;
    }
  SyncRecord *ret_rec = &*(m_Curr_rec++);

  return ret_rec;
}
//-------------------------------------------------------------------------
///
void SyncRecordManager::printRecords( void )
{
  cout << endl << "total " << m_Recs.size() << " record(s):" << endl;
  for( unsigned int i = 0 ; i < m_Recs.size() ; i++ )
    {
    cout << "---------------------------------------------- SyncRecord # " <<
      i << " ----------------------------------------------\n";
    m_Recs[i].print();
    }
  cout << endl;
}
//-------------------------------------------------------------------------
///
bool SyncRecordManager::istreamSkipPastEOL( istream &is, int howMany )
{
  int err = 0;
  for( int i = 0 ; i < howMany ; ++i )
    {
    if( !is.eof() && is.good() )
      {
      is.ignore( numeric_limits<streamsize>::max(), '\n' );
      }
    else
      {
      ++err;
      }
    }
  if( err )
    {
    cerr << "SyncRecordManager::istreamSkipPastEOL(): tried to skip past "
      << howMany << " line ends but encountered " << err << " errors" << endl;
    return false;
    }
  return true;
}
//-------------------------------------------------------------------------
///
bool SyncRecordManager::istreamMatrix( istream &is, const char *header, double m[16], const bool verbose )
{
  if( !istreamSkipPastEOL( is, 1 ) )
    {
    cerr << "SyncRecordManager::istreamMatrix(): cannot skip header line '" <<
      header << "'" << endl;
    return false;
    }
  // later we can scan in the header and actually verify it
  int err = 0;
  for( int i = 0 ; i < 4 ; ++i )
    {
    for( int j = 0 ; j < 4 ; ++j )
      {
      if( !is.eof() && is.good() )
        {
        err += !( is >> m[4*i+j] );
        }
      }
    }
  if( err )
    {
    cerr << "SyncRecordManager::istreamMatrix(): could not read " << header
      << " matrix" << endl;
    return false;
    }
  if( verbose )
    {
    cout << "SyncRecordManager::istreamMatrix(): read following matrix:" << endl;
    ostreamMatrix( cout, header, m );
    }
  if( !istreamSkipPastEOL( is, 1 ) )
    {
    cerr << "SyncRecordManager::istreamMatrix(): cannot skip to end of matrix" << endl;
    return false;
    }
  return true;
}
//-------------------------------------------------------------------------
///
bool SyncRecordManager::load( const char *disk_load_path )
{
  if( !m_Recs.empty() )
    {
    cerr << "SyncRecordManager::load(): attempt to initialize non-empty \
      SyncRecordManager object" << endl;
    return false;
    }

  m_Records_path = disk_load_path;
  ifstream ifs( m_Records_path.c_str() );
  if( !ifs || !ifs.is_open() )
    {
    cerr << "SyncRecordManager::load(): file " << m_Records_path <<
      " is not readable" << endl;
    return false;
    }

  // skip leading comment lines starting with # or empty
  while( !ifs.eof() && ifs.good() )
    {
    char c = ifs.peek();
    if( c == '#' || std::isspace(c) )
      {
      if( !istreamSkipPastEOL( ifs, 1 ) )
        {
        cerr << "SyncRecordManager::load(): problem skipping leading \
          blank/comment line(s)" << endl;
        ifs.close();
        return false;
        }
      }
    else 
      {
      break;
      }
    }

  // volume dataset path
  if( !ifs.eof() && ifs.good() )
    {
    if( !( ifs >> m_Volume_path ) )
      {
      cerr << "SyncRecordManager::load(): missing volume dataset path"
        << endl;
      ifs.close();
      return false;
      }
    }
  if( !istreamSkipPastEOL( ifs, 1 ) )
    {
    cerr << "SyncRecordManager::load(): problem skipping newline after \
      volume_path" << endl;
    ifs.close();
    return false;
    }
  // navel registration matrices
  if( !istreamMatrix( ifs, "tracker_f_volume", tracker_f_volume, false ) )
    {
    ifs.close();
    return false;
    }
  if( !istreamMatrix( ifs, "tracker_f_navel" , tracker_f_navel , false ) )
    {
    ifs.close();
    return false;
    }

  // skip column titles line
  if( !istreamSkipPastEOL( ifs, 1 ) )
    {
    cerr << "SyncRecordManager::load(): cannot skip column titles line" << endl;
    ifs.close(); return false;
    }

  // load all records line by line
  m_Count = 0;
  m_Recs.resize( 2222 );
  while( !ifs.eof() && ifs.good() )
    {
    if( !m_Recs[m_Count].load( ifs ) ) 
      {
      break;
      }
    if( ++m_Count >= m_Recs.size() )
      {
      m_Recs.resize( m_Recs.size() + 1111 );
      }
    }
  m_Recs.resize( m_Count );
  ifs.close();

  // here we may have to check and enforce some limit on the number of records (based on main or video memory sizes, texture subsampling, etc.

  /*printRecords();*/

  rewind();
  return true;
}
//=========================================================================
// vim: set noexpandtab: REQUIRED BY ANDREI's EDITOR
