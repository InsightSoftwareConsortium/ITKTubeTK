/*! @file

    @brief class for management of a single-slice record of synchronized AIM->Kitware/Spotlight tracking + ultrasound data

    Supports disk reading and writing of, as well as internal access to a synchronous tracking + ultrasound data record.
	Shared between AIM and Spotlight, and with Kitware.

    @note  Copyright (c) InnerOptic Technology, Inc.  All Rights Reserved.

    @author  Andrei State

*/
//==========================================================================================================================================
#ifndef INNEROPTIC_SYNCRECORD
#define INNEROPTIC_SYNCRECORD

#define MAX_US_SCAN_CROP_POLYGON_PTS 4

// string representing order of values output by the dump function
#define SYNCRECORD_DUMP_ORDER "timestamp_in_msec   u/s_image_path   cropped_u/s_polygon_corners_x0_y0_x1_y1_x2_y2_x3_y3_in_raw_u/s_frame   tracker_from_RUF_matrix_by_columns   transducer_from_RUF_matrix_by_columns   endo_image_rectangle_x0_y0_x1_y1_in_raw_u/s_frame"

#include <string>
using namespace std;

//==========================================================================================================================================
class SyncRecord{
public:
	SyncRecord( void );
   ~SyncRecord( void );

#   ifdef INNEROPTIC_INTERNAL_ONLY
		void setTimestamp( int stamp );
		void setRufImageFilePath( const char *ruf_path );
		void setRufImageFileIndex( const int index );
		bool setScanCropVertex_in_ruf( const int index, double x, double y );
		void setTrackerFromRufMatrix(    const double m[16] );
		void setTransducerFromRufMatrix( const double m[16] );
		void setEndoImageGeometry_in_ruf( const int x0, const int y0, const int x1, const int y1 );
#   endif
	int getTimestamp( void );
	const char *getRufImageFilePath( void );
	bool getScanCropVertex_in_ruf( const int index, double &x, double &y );
	void getTrackerFromRufMatrix(    double m[16] );
	void getTransducerFromRufMatrix( double m[16] );
	void getEndoImageGeometry_in_ruf( int &x0, int &y0, int &x1, int &y1 );
	bool loadRufImageRgbPixels( unsigned char **rgb_pixels );

	void print( void );

	void dump( ofstream &ofs );
	bool load( ifstream &ifs );

	unsigned char *loadRawRgbPixels( void );
	bool         unloadRawRgbPixels( void );
protected:
private:
	void ostreamMatrix( ostream &os, const char *header, const double m[16] );
	bool loadImageAsciiHeader( void );

	int timestamp;
	string ruf_image_path;
	double scanCropPolygon_in_ruf_x[MAX_US_SCAN_CROP_POLYGON_PTS];
	double scanCropPolygon_in_ruf_y[MAX_US_SCAN_CROP_POLYGON_PTS];
	double tracker_f_ruf[   16]; // orthogonal matrices (contain scaling), in column-major order to match OpenGL manpage insanity
	double transducer_f_ruf[16];
	int endo_x0_in_ruf, endo_y0_in_ruf, endo_x1_in_ruf, endo_y1_in_ruf; // laparo image rectangle (0 is top left point, 1 is bottom right one)

	/*int index;*/

	// additional fields for data retrieved from PPM image header (resolution, timestamps, frame numbers, etc.)

	unsigned char *allRgbPixels; // points to 24-bit RGB pixels data from ultrasound PPM file
	unsigned char *usGreyPixels; // points to 8-bit pixels containing only the scanCropPolygon

#   ifdef INNEROPTIC_INTERNAL_ONLY
		// additional fields for two OpenGL texture IDs?
#   endif
};
#endif
//==========================================================================================================================================
// vim: set noexpandtab: REQUIRED BY ANDREI's EDITOR
