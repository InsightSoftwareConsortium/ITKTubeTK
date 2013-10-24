/*! @file

    @brief class for management of synchronized AIM->Kitware/Spotlight tracking + ultrasound data single-slice record lists

    Supports disk reading, writing and playing back of synchronous tracking + ultrasound data record lists
	initially saved asynchronously by AIM's Save_CMC branch, then synchronized and streamlined by the Replay_Reexport_CMC branch.
	Shared between AIM and Spotlight, and with Kitware.

    @note  Copyright (c) InnerOptic Technology, Inc.  All Rights Reserved.

    @author  Andrei State

*/
//==========================================================================================================================================
#ifndef INNEROPTIC_SYNCRECORDMANAGER
#define INNEROPTIC_SYNCRECORDMANAGER

#include <vector>
#include "SyncRecord.h"
using namespace std;

//==========================================================================================================================================
class SyncRecordManager{
public:
	SyncRecordManager( void );
   ~SyncRecordManager( void );

	bool load( const char *disk_load_path );

#   ifdef INNEROPTIC_INTERNAL_ONLY
		void setVolumeImagePath( const char *vol_path );
		void setTrackerFromVolumeImageMatrix( const double m[16] );
		void setTrackerFromNavelMatrix(       const double m[16] );
		SyncRecord *newSequentialAppendedRecord( void );
		bool dump( const char *disk_write_path, char *ident );
#   endif
	const char *getVolumeImagePath( void );
	void getTrackerFromVolumeImageMatrix( double m[16] );
	void getTrackerFromNavelMatrix(       double m[16] );

	void printRecords( void );

	SyncRecord *getNextRecord( void );
	SyncRecord *getRecord( const int which );

	void rewind( void );
	int getNbRecords( void );

protected:
private:
	string records_path, volume_path;
	vector<SyncRecord> recs;
	vector<SyncRecord>::iterator curr_rec;
	double tracker_f_volume[16], tracker_f_navel[16];
	unsigned int count;

	void ostreamMatrix( ostream &os, const char *header, const double m[16] );
	bool istreamSkipPastEOL( istream &is, const int howMany );
	bool istreamMatrix( istream &is, const char *header, double m[16], const bool verbose );
};
#endif
//==========================================================================================================================================
// vim: set ts=4 sw=4 tw=0 noexpandtab lines=77 columns=195: REQUIRED BY ANDREI's EDITOR
