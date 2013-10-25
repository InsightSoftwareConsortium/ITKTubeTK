/*! @file

    @brief class for management of synchronized AIM->Kitware/Spotlight tracking + ultrasound data single-slice record lists

    Supports disk reading, writing and playing back of synchronous tracking + ultrasound data record lists
	initially saved asynchronously by AIM's Save_CMC branch, then synchronized and streamlined by the Replay_Reexport_CMC branch.
	Shared between AIM and Spotlight, and with Kitware.

    @note  Copyright (c) InnerOptic Technology, Inc.  All Rights Reserved.

    @author  Andrei State

*/
//==========================================================================================================================================
#define _CRT_SECURE_NO_WARNINGS

/*#include <assert.h>*/
#include <stdio.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <cctype>

#include "SyncRecordManager.h"
//==========================================================================================================================================
///
SyncRecordManager::SyncRecordManager( void )
{
	volume_path.assign( "VOID_VOLUME_PATH" );

	for( int i = 0 ; i < 16 ; ++i ) tracker_f_volume[i] = tracker_f_navel[i] = 0.;
	count = 0; rewind();
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
SyncRecordManager::~SyncRecordManager( void )
{
	recs.resize( 0 );
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecordManager::ostreamMatrix( ostream &os, const char *header, const double m[16] )
{
	os << header << endl;
	for( int i = 0 ; i < 4 ; ++i ) for( int j = 0 ; j < 4 ; ++j ) { os << m[4*i+j]; if( j < 3 ) os << " "; else os << endl; }
}
//==========================================================================================================================================
#ifdef INNEROPTIC_INTERNAL_ONLY
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecordManager::setVolumeImagePath( const char *vol_path )
{
	volume_path.assign( vol_path );
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecordManager::setTrackerFromVolumeImageMatrix( const double m[16] )
{
	for( int i = 0 ; i < 16 ; ++i ) tracker_f_volume[i] = m[i];
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecordManager::setTrackerFromNavelMatrix(  const double m[16] )
{
	for( int i = 0 ; i < 16 ; ++i ) tracker_f_navel[i] = m[i];
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
SyncRecord *SyncRecordManager::newSequentialAppendedRecord( void )
{
	++count;
	if( count >= recs.size() ) recs.resize( count );
	return &recs[count-1];
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
bool SyncRecordManager::dump( const char *disk_write_path, char *ident )
{
	if( !recs.size() ){ cerr << "SyncRecordManager::dump(): attempt to dump empty SyncRecordManager object" << endl; return false; }

	ofstream ofs( disk_write_path );
	if( !ofs || !ofs.is_open() ){ cerr << "SyncRecordManager::dump(): file " << disk_write_path << " is not writable" << endl; return false; }
	cout << "SyncRecordManager::dump(): dumping to playback file " << disk_write_path << endl;

	// write a few leading comment lines starting with #, and empty
	ofs << "# Synchronized ultrasound scan / tracking data acquired and processed by InnerOptic" << endl;
	ofs << "# Additional info: " << ident << endl;
	ofs << endl;
	ofs << "# Lines at the beginning of this file that are either empty or start with # are considered comments" << endl;
	ofs << endl;
	ofs << "# The next 11 non-empty lines represent the volume data path and the two navel registration matrices" << endl;
	ofs << "# Then follow the column title line and the actual bulk data (1 line per tracked ultrasound slice), total " << recs.size() << " single-slice records" << endl;
	ofs << endl;

	ofs << volume_path << endl;

	// navel registration matrices
	ostreamMatrix( ofs, "tracker_f_volume", tracker_f_volume );
	ostreamMatrix( ofs, "tracker_f_navel" , tracker_f_navel  );

	// column titles line
	ofs << SYNCRECORD_DUMP_ORDER << endl;

	// dump all records line by line
	for( vector<SyncRecord>::iterator i = recs.begin() ; i != recs.end() ; i++ ) (*i).dump( ofs );

	ofs.close();
	cout << "SyncRecordManager::dump(): dumped " << recs.size() << " synchronized tracked ultrasound frame records" << endl;

	return true;
}
//------------------------------------------------------------------------------------------------------------------------------------------
#endif
//==========================================================================================================================================
const char *SyncRecordManager::getVolumeImagePath( void )
{
	return volume_path.c_str();
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecordManager::getTrackerFromVolumeImageMatrix( double m[16] )
{
	for( int i = 0 ; i < 16 ; ++i ) m[i] = tracker_f_volume[i];
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecordManager::getTrackerFromNavelMatrix( double m[16] )
{
	for( int i = 0 ; i < 16 ; ++i ) m[i] = tracker_f_navel[i];
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecordManager::rewind( void ){ curr_rec = recs.begin(); }
//------------------------------------------------------------------------------------------------------------------------------------------
///
int SyncRecordManager::getNbRecords( void ){ return recs.size(); }
//------------------------------------------------------------------------------------------------------------------------------------------
///
SyncRecord *SyncRecordManager::getNextRecord( void )
{
	if( !recs.size() ){ cerr << "SyncRecordManager::getNextRecord(): manager object is empty" << endl; return NULL; }
	if( curr_rec == recs.end() ) return NULL;
	SyncRecord *ret_rec = &*(curr_rec++);

	return ret_rec;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecordManager::printRecords( void )
{
	cout << endl << "total " << recs.size() << " record(s):" << endl;
	for( unsigned int i = 0 ; i < recs.size() ; i++ ){
		cout << "---------------------------------------------- SyncRecord # " << i << " ----------------------------------------------\n";
		recs[i].print();
	}
	cout << endl;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
bool SyncRecordManager::istreamSkipPastEOL( istream &is, int howMany )
{
	int err = 0;
	for( int i = 0 ; i < howMany ; ++i ) if( !is.eof() && is.good() ) is.ignore( numeric_limits<streamsize>::max(), '\n' ); else ++err;
	if( err ){
		cerr << "SyncRecordManager::istreamSkipPastEOL(): tried to skip past " << howMany << " line ends but encountered " << err << " errors" << endl;
		return false;
	}
	return true;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
bool SyncRecordManager::istreamMatrix( istream &is, const char *header, double m[16], const bool verbose )
{
	if( !istreamSkipPastEOL( is, 1 ) ){
		cerr << "SyncRecordManager::istreamMatrix(): cannot skip header line '" << header << "'" << endl; return false;
	}
	// later we can scan in the header and actually verify it
	int err = 0;
	for( int i = 0 ; i < 4 ; ++i ) for( int j = 0 ; j < 4 ; ++j ) if( !is.eof() && is.good() ) err += !( is >> m[4*i+j] );
	if( err ){     cerr << "SyncRecordManager::istreamMatrix(): could not read " << header << " matrix" << endl; return false; }
	if( verbose ){ cout <<       "SyncRecordManager::istreamMatrix(): read following matrix:" << endl; ostreamMatrix( cout, header, m ); }
	if( !istreamSkipPastEOL( is, 1 ) ){
		cerr            << "SyncRecordManager::istreamMatrix(): cannot skip to end of matrix" << endl; return false;
	}
	return true;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
bool SyncRecordManager::load( const char *disk_load_path )
{
	if( recs.size() ){ cerr << "SyncRecordManager::load(): attempt to initialize non-empty SyncRecordManager object" << endl; return false; }

	records_path = disk_load_path;
	ifstream ifs( records_path.c_str() );
	if( !ifs || !ifs.is_open() ){ cerr << "SyncRecordManager::load(): file " << records_path << " is not readable" << endl; return false; }
	cout << "SyncRecordManager::load(): reading from synchronized-records file " << records_path << endl;

	// skip leading comment lines starting with # or empty
	while( !ifs.eof() && ifs.good() ){
		char c = ifs.peek();
		if( c == '#' || std::isspace(c) ){
			if( !istreamSkipPastEOL( ifs, 1 ) ){
				cerr <<   "SyncRecordManager::load(): problem skipping leading blank/comment line(s)" << endl; ifs.close(); return false;
			}
		}else break;
	}

	// volume dataset path
	if( !ifs.eof() && ifs.good() ) if( !( ifs >> volume_path ) ){ cerr << "SyncRecordManager::load(): missing volume dataset path" << endl; ifs.close(); return false; }
	cout << "SyncRecordManager::load(): read volume dataset path " << volume_path << endl;
	if( !istreamSkipPastEOL( ifs, 1 ) ){
		cerr <<   "SyncRecordManager::load(): problem skipping newline after volume_path" << endl; ifs.close(); return false;
	}

	// navel registration matrices
	if( !istreamMatrix( ifs, "tracker_f_volume", tracker_f_volume, true ) ){ ifs.close(); return false; }
	if( !istreamMatrix( ifs, "tracker_f_navel" , tracker_f_navel , true ) ){ ifs.close(); return false; }

	// skip column titles line
	if( !istreamSkipPastEOL( ifs, 1 ) ){ cerr << "SyncRecordManager::load(): cannot skip column titles line" << endl; ifs.close(); return false; }

	// load all records line by line
	count = 0;
	recs.resize( 2222 );
	while( !ifs.eof() && ifs.good() ){
		/*cout << " " << count;*/
		if( !recs[count].load( ifs ) ) break;

#if 0
		if( count > 1150 ){
			cout << "==============> " << count << endl;
			recs[count].print();
		}
#endif

		if( ++count >= recs.size() ) recs.resize( recs.size() + 1111 );
	}
	recs.resize( count );
	ifs.close();
	cout << "SyncRecordManager::load(): loaded " << count << " synchronized tracked ultrasound frame records" << endl;

	// here we may have to check and enforce some limit on the number of records (based on main or video memory sizes, texture subsampling, etc.

	/*printRecords();*/

	rewind();
	return true;
}
//==========================================================================================================================================
// vim: set noexpandtab: REQUIRED BY ANDREI's EDITOR
