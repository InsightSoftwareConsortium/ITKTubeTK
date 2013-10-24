/*! @file

    @brief class for management of a single-slice record of synchronized AIM->Kitware/Spotlight tracking + ultrasound data

    Supports disk reading and writing of, as well as internal access to a synchronous tracking + ultrasound data record.
	Shared between AIM and Spotlight, and with Kitware.

    @note  Copyright (c) InnerOptic Technology, Inc.  All Rights Reserved.

    @author  Andrei State

*/
//==========================================================================================================================================
#define _CRT_SECURE_NO_WARNINGS

#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "SyncRecord.h"
//==========================================================================================================================================
///
SyncRecord::SyncRecord( void )
{
	allRgbPixels = NULL;

	// hardcoded (and thus not drift-corrected) laparo image corners; top left is 0,0
	endo_x0_in_ruf =   3; // upper left corner
	endo_y0_in_ruf = 291;
	endo_x1_in_ruf = 262; // bottom right corner
	endo_y1_in_ruf = 485;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
SyncRecord::~SyncRecord( void )
{
	if( allRgbPixels ) delete allRgbPixels;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecord::ostreamMatrix( ostream &os, const char *header, const double m[16] )
{
	os << header << endl;
	for( int i = 0 ; i < 4 ; ++i ) for( int j = 0 ; j < 4 ; ++j ) { os << m[4*i+j]; if( j < 3 ) os << " "; else os << endl; }
}
//==========================================================================================================================================
#ifdef INNEROPTIC_INTERNAL_ONLY
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecord::setTimestamp( int stamp )
{
	timestamp = stamp;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecord::setRufImageFilePath( const char *ruf_path )
{
	ruf_image_path.assign( ruf_path );
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecord::setRufImageFileIndex( const int index )
{
	ostringstream ss;
	ss << "ultrasound_" << setw( 7 ) << setfill( '0' ) << index << ".ppm";
	
	ruf_image_path.assign( ss.str() );
	/*cout << ruf_image_path;*/
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
bool SyncRecord::setScanCropVertex_in_ruf( const int index, double x, double y )
{
	if( index < 0 || index >= MAX_US_SCAN_CROP_POLYGON_PTS ){
        cerr << "SyncRecord::setScanCropVertex_in_ruf(): index " << index << " outside bounds " << 0 << "-" << MAX_US_SCAN_CROP_POLYGON_PTS << endl;
		return false;
	}
	scanCropPolygon_in_ruf_x[index] = x;
	scanCropPolygon_in_ruf_y[index] = y;
	return true;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecord::setTrackerFromRufMatrix( const double m[16] )
{
	for( int i = 0 ; i < 16 ; ++i ) tracker_f_ruf[i] = m[i];
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecord::setTransducerFromRufMatrix( const double m[16] )
{
	for( int i = 0 ; i < 16 ; ++i ) transducer_f_ruf[i] = m[i];
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecord::setEndoImageGeometry_in_ruf( const int x0, const int y0, const int x1, const int y1 )
{
	endo_x0_in_ruf = x0;
	endo_y0_in_ruf = y0;
	endo_x1_in_ruf = x1;
	endo_y1_in_ruf = y1;
}
//------------------------------------------------------------------------------------------------------------------------------------------
#endif
//==========================================================================================================================================
///
int SyncRecord::getTimestamp( void )
{
	return timestamp;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
const char *SyncRecord::getRufImageFilePath( void )
{
	return ruf_image_path.c_str();
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
bool SyncRecord::getScanCropVertex_in_ruf( const int index, double &x, double &y )
{
	if( index < 0 || index >= MAX_US_SCAN_CROP_POLYGON_PTS ){
        cerr << "SyncRecord::getScanCropVertex_in_ruf(): index " << index << " outside bounds " << 0 << "-" << MAX_US_SCAN_CROP_POLYGON_PTS << endl;
		return false;
	}
	x = scanCropPolygon_in_ruf_x[index];
	y = scanCropPolygon_in_ruf_y[index];
	return true;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecord::getTrackerFromRufMatrix( double m[16] )
{
	for( int i = 0 ; i < 16 ; ++i ) m[i] = tracker_f_ruf[i];
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecord::getTransducerFromRufMatrix( double m[16] )
{
	for( int i = 0 ; i < 16 ; ++i ) m[i] = transducer_f_ruf[i];
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecord::getEndoImageGeometry_in_ruf( int &x0, int &y0, int &x1, int &y1 )
{
	x0 = endo_x0_in_ruf;
	y0 = endo_y0_in_ruf;
	x1 = endo_x1_in_ruf;
	y1 = endo_y1_in_ruf;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///                                                                                                   prints human-readable record to stdout
void SyncRecord::print( void )
{
	streamsize oldp = cout.precision(); cout.precision( 9 ); // arbitrary; the matrices come from Mat44d types, so they are double precision
	;                                                         cout << "timestamp_in_msec: " << timestamp
																   << ", u/s_image_path: "  << ruf_image_path << endl;
	;                                                         cout << "cropped_u/s_polygon_corners_x0_y0_x1_y1_x2_y2_x3_y3_in_raw_u/s_frame:";
	for( int i = 0 ; i < MAX_US_SCAN_CROP_POLYGON_PTS ; ++i ) cout << " "                   << scanCropPolygon_in_ruf_x[i] 
																   << " "                   << scanCropPolygon_in_ruf_y[i];
	;                                                         cout << endl;
	cout.setf( ios::fixed, ios::floatfield );
	ostreamMatrix( cout,    "tracker_from_raw_u/s_frame_matrix:",                                 tracker_f_ruf );
	ostreamMatrix( cout, "transducer_from_raw_u/s_frame_matrix:",                              transducer_f_ruf );
	cout.unsetf( ios::fixed );
	cout.unsetf( ios::floatfield );
	;                                                         cout << "endo_image_rectangle_x0_yo_x1_y1_in_raw_u/s_frame:";
	;                                                         cout << " "                   << endo_x0_in_ruf;
	;                                                         cout << " "                   << endo_y0_in_ruf;
	;                                                         cout << " "                   << endo_x1_in_ruf;
	;                                                         cout << " "                   << endo_y1_in_ruf;
	;                                                         cout << endl;
	cout.precision( oldp );
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
void SyncRecord::dump( ofstream &ofs )
{
	streamsize oldp = ofs.precision(); ofs.precision( 11 ); // 11 is arbitrary; the matrices come from Mat44d types, so they are double precision
	;                                                         ofs <<        timestamp
																  << " " << ruf_image_path;
	for( int i = 0 ; i < MAX_US_SCAN_CROP_POLYGON_PTS ; ++i ) ofs << " " << scanCropPolygon_in_ruf_x[i] 
																  << " " << scanCropPolygon_in_ruf_y[i];
	for( int i = 0 ; i < 16                           ; ++i ) ofs << " " <<    tracker_f_ruf[i];
	for( int i = 0 ; i < 16                           ; ++i ) ofs << " " << transducer_f_ruf[i];
	;                                                         ofs << " " << endo_x0_in_ruf;
	;                                                         ofs << " " << endo_y0_in_ruf;
	;                                                         ofs << " " << endo_x1_in_ruf;
	;                                                         ofs << " " << endo_y1_in_ruf;
	;                                                         ofs <<        endl;
	ofs.precision( oldp );
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
bool SyncRecord::loadImageAsciiHeader( void )
{
	// open PPM file in text mode to load its Save_CMC-specific ASCII header info
	ifstream ifs(                                                                     ruf_image_path.c_str() );
	if(     !ifs           ){ cerr << "SyncRecord::loadImageAsciiHeader(): file " << ruf_image_path << " not readable" << endl; return false; }
	if(     !ifs.is_open() ){ cerr << "SyncRecord::loadImageAsciiHeader(): file " << ruf_image_path << " not readable" << endl; return false; }
	// not loading anything here for now
	// we assume PPM (magic number P6) and resolution 960x768
	// Brian's code [UsVideo::parseUltrasoundFrameFileHeader() and UsVideo::getUltrasoundFrameNumber()] should be used here from a separate module

	ifs.close();
	return true;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///                                                                                              loads single-line-formatted record from ifs
bool SyncRecord::load( ifstream &ifs )
{
	int err = 0;
	;                                                         err += !( ifs >> timestamp
																			>> ruf_image_path );
	for( int i = 0 ; i < MAX_US_SCAN_CROP_POLYGON_PTS ; ++i ) err += !( ifs >> scanCropPolygon_in_ruf_x[i] 
												                            >> scanCropPolygon_in_ruf_y[i] );
	for( int i = 0 ; i < 16                           ; ++i ) err += !( ifs >>    tracker_f_ruf[i] );
	for( int i = 0 ; i < 16                           ; ++i ) err += !( ifs >> transducer_f_ruf[i] );
#if 1
	;                                                         err += !( ifs >> endo_x0_in_ruf );
	;                                                         err += !( ifs >> endo_y0_in_ruf );
	;                                                         err += !( ifs >> endo_x1_in_ruf );
	;                                                         err += !( ifs >> endo_y1_in_ruf );
#endif

	if( err ) return false;
	if( !loadImageAsciiHeader() ){ // load info from ultrasound image file
		cerr << "SyncRecord::load(): problem loading header of image file " << ruf_image_path << endl;
		ifs.close(); return false;
	}
	return true;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///                                                                                   not part of class since it should be moved out of here
static unsigned int readBinFileAnySizeWithOffset( const char *filename, unsigned char **bytes, long offset )
{
    FILE *fp;
    if( !( fp = fopen( filename, "rb" ) ) ){
        cerr << "readBinFileAnySizeWithOffset(): could not open " << filename << " for reading" << endl;
        return 0;
    }
    struct stat st;
    if( fstat( fileno( fp ), &st ) ){
        cerr << "readBinFileAnySizeWithOffset(): _fstat() failed on " << filename << endl;
		fclose( fp );
        return 0;
    }
	if( offset ){
		if( offset >= (long)st.st_size ){
			cerr << "readBinFileAnySizeWithOffset(): offset " << offset << " is too large for the size of file " << filename
				<< ", which is only " << st.st_size << " Bytes" << endl;
			fclose( fp );
			return 0;
		}
		if( fseek( fp, offset, SEEK_SET ) ){
			cerr << "readBinFileAnySizeWithOffset(): seek error w/ offset " << offset << " on file " << filename << endl;
			fclose( fp );
			return 0;
		}
	}
    *bytes = new unsigned char[st.st_size-offset];
    int totalRead = 0;
    while( !feof( fp ) ){ totalRead += fread( *bytes, sizeof( unsigned char ), st.st_size-offset, fp ); }
    fclose( fp );
    assert( totalRead == st.st_size-offset );
    return totalRead;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///                                                                                   not part of class since it should be moved out of here
static void fastWriteFile( char *filename, unsigned char *bytes, unsigned int size, bool writeThrough, bool noBuffering )
{
    FILE *fp = fopen( filename, "wb" );
    if( fp ){
        fwrite( bytes, 1, size * sizeof(unsigned char), fp );
        /*fwrite( block, size * sizeof(unsigned char), 1, fp );*/
        //fflush( fp ); // doesn't seem to make a difference
        fclose( fp );
    }
    /*slog.Exit( -1, "Wrote image to file %s\n", filename );*/
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
unsigned char *SyncRecord::loadRawRgbPixels( void )
{
	if( allRgbPixels ){
        cerr << "SyncRecord::loadRawRgbPixels(): raw RGB pixels already loaded" << endl;
		return NULL;
	}
	// open and read PPM file and retrieve its binary data, ready for glTexSubImage() (but no OpenGL here)
#if 0
	ifstream ifb(                                                  ret_rec->ruf_image_path, ios::binary );
	if(     !ifb ){ cerr << "SyncRecordManager::getNextRecord(): image file " << ret_rec->ruf_image_path << " is not readable" << endl; return NULL; }
#else
	// use fread() for max speed
	// should calculate an offset to read only from top of highest bounding box in file to bottom of lowest one
	unsigned int bytes_read;
	if( !( bytes_read = readBinFileAnySizeWithOffset( ruf_image_path.c_str(), &allRgbPixels, (long)4096 ) ) ){
        cerr << "SyncRecord::loadRawRgbPixels(): could not load raw RGB pixels" << endl;
		if( allRgbPixels ){ delete allRgbPixels; allRgbPixels = NULL; }
		return NULL;
	}
	assert( allRgbPixels );
	assert( bytes_read = 3 * 960 * 768 ); // later calculate this based on PPM image resolution
#endif
	// cout << "SyncRecordManager::loadRawRgbPixels(): read " << bytes_read << " Bytes from " << ruf_image_path  << endl;

#if 0
	// test code to write out pixel data fr verification
	static int first = 1; if( first ){ first = 0;
		fastWriteFile( "_alpha.raw", allRgbPixels, 3*960*768, true, true );
	}
#endif

	return allRgbPixels;
}
//------------------------------------------------------------------------------------------------------------------------------------------
///
bool SyncRecord::unloadRawRgbPixels( void )
{
	if( !allRgbPixels ){
        cerr << "SyncRecord::loadRawRgbPixels(): no raw RGB pixels to unload" << endl;
		return false;
	}
	delete allRgbPixels;
	allRgbPixels = NULL;
	return true;
}
//==========================================================================================================================================
// vim: set noexpandtab: REQUIRED BY ANDREI's EDITOR
