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
#include "itkTubeMetaTubeParams.h"

#include <stdio.h>
#include <ctype.h>
#include <string>
#include <string.h> // for memcpy
#include <math.h>

namespace itk {

namespace tube {

//
// MetaTubeParams Constructors
//
MetaTubeParams::
MetaTubeParams()
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaTubeParams()" << METAIO_STREAM::endl;
    }

  Clear();
}

//
MetaTubeParams::
MetaTubeParams( const char *_headerName )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaTubeParams()" << METAIO_STREAM::endl;
    }

  Clear();
  m_ReadStream = NULL;

  MetaTubeParams::Read( _headerName );
}

//
MetaTubeParams::
MetaTubeParams( const MetaTubeParams & _metaTubeParams ):
  MetaForm()
{
  if( META_DEBUG )
   {
   METAIO_STREAM::cout << "MetaTubeParams()" << METAIO_STREAM::endl;
   }

  Clear();

  CopyInfo( _metaTubeParams );
}

//
MetaTubeParams::
~MetaTubeParams()
{
  M_Destroy();
}

//
void MetaTubeParams::
PrintInfo() const
{
  MetaForm::PrintInfo();

  METAIO_STREAM::cout << "SeedScales = " << m_SeedScales
    << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "SeedIntensityMin = " << m_SeedIntensityMin
    << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "SeedIntensityMax = " << m_SeedIntensityMax
    << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "SeedIntensityPercentile = "
    << m_SeedIntensityPercentile << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeIntensityMin = "
    << m_TubeIntensityMin << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeIntensityMax = "
    << m_TubeIntensityMax << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeBright = "
    << m_TubeBright << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeColor = "
    << m_TubeColor << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeScale = "
    << m_TubeRidgeScale << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeScaleExtent = "
    << m_TubeRidgeScaleExtent << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeDynamicScale = "
    << m_TubeRidgeDynamicScale << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeStepX = "
    << m_TubeRidgeStepX << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeThresholdTangentChange = "
    << m_TubeRidgeThresholdTangentChange << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeThresholdXChange = "
    << m_TubeRidgeThresholdXChange << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeThresholdRidgeness = "
    << m_TubeRidgeThresholdRidgeness << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeThresholdRidgenessStart = "
    << m_TubeRidgeThresholdRidgenessStart << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeThresholdRoundness = "
    << m_TubeRidgeThresholdRoundness << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeThresholdRoundnessStart = "
    << m_TubeRidgeThresholdRoundnessStart << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeCurvatureMax = "
    << m_TubeRidgeCurvatureMax << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeThresholdCurvature = "
    << m_TubeRidgeThresholdCurvature << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeThresholdCurvatureStart = "
    << m_TubeRidgeThresholdCurvatureStart << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeThresholdLinearity = "
    << m_TubeRidgeThresholdLinearity << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeThresholdLinearityStart = "
    << m_TubeRidgeThresholdLinearityStart << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRidgeRecoveryMax = "
    << m_TubeRidgeRecoveryMax << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRadiusStart = "
    << m_TubeRadiusStart << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRadiusMin = "
    << m_TubeRadiusMin << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRadiusMax = "
    << m_TubeRadiusMax << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRadiusThresholdMedialness = "
    << m_TubeRadiusThresholdMedialness << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "TubeRadiusThresholdMedialnessStart = "
    << m_TubeRadiusThresholdMedialnessStart << METAIO_STREAM::endl;

}

void MetaTubeParams::
SetSeedParams( const VectorType & _seedScales,
  double _seedIntensityMin, double _seedIntensityMax,
  double _seedIntensityPercentile )
{
  m_SeedScales = _seedScales;
  m_SeedIntensityMin = _seedIntensityMin;
  m_SeedIntensityMax = _seedIntensityMax;
  m_SeedIntensityPercentile = _seedIntensityPercentile;
}

void MetaTubeParams::
SetTubeParams( double _tubeIntensityMin, double _tubeIntensityMax,
  bool _tubeBright, const VectorType & _tubeColor )
{
  m_TubeIntensityMin = _tubeIntensityMin;
  m_TubeIntensityMax = _tubeIntensityMax;
  m_TubeBright = _tubeBright;
  m_TubeColor = _tubeColor;
}

void MetaTubeParams::
SetTubeRidgeParams( double _tubeRidgeScale, double _tubeRidgeScaleExtent,
  bool _tubeRidgeDynamicScale, double _tubeRidgeStepX,
  double _tubeRidgeThresholdTangentChange,
  double _tubeRidgeThresholdXChange,
  double _tubeRidgeThresholdRidgeness,
  double _tubeRidgeThresholdRidgenessStart,
  double _tubeRidgeThresholdRoundness,
  double _tubeRidgeThresholdRoundnessStart,
  double _tubeRidgeCurvatureMax,
  double _tubeRidgeThresholdCurvature,
  double _tubeRidgeThresholdCurvatureStart,
  double _tubeRidgeThresholdLinearity,
  double _tubeRidgeThresholdLinearityStart,
  int _tubeRidgeRecoveryMax )
{
  m_TubeRidgeScale = _tubeRidgeScale;
  m_TubeRidgeScaleExtent = _tubeRidgeScaleExtent;
  m_TubeRidgeDynamicScale = _tubeRidgeDynamicScale;
  m_TubeRidgeStepX = _tubeRidgeStepX;
  m_TubeRidgeThresholdTangentChange = _tubeRidgeThresholdTangentChange;
  m_TubeRidgeThresholdXChange = _tubeRidgeThresholdXChange;
  m_TubeRidgeThresholdRidgeness = _tubeRidgeThresholdRidgeness;
  m_TubeRidgeThresholdRidgenessStart = _tubeRidgeThresholdRidgenessStart;
  m_TubeRidgeThresholdRoundness = _tubeRidgeThresholdRoundness;
  m_TubeRidgeThresholdRoundnessStart = _tubeRidgeThresholdRoundnessStart;
  m_TubeRidgeCurvatureMax = _tubeRidgeCurvatureMax;
  m_TubeRidgeThresholdCurvature = _tubeRidgeThresholdCurvature;
  m_TubeRidgeThresholdCurvatureStart = _tubeRidgeThresholdCurvatureStart;
  m_TubeRidgeThresholdLinearity = _tubeRidgeThresholdLinearity;
  m_TubeRidgeThresholdLinearityStart = _tubeRidgeThresholdLinearityStart;
  m_TubeRidgeRecoveryMax = _tubeRidgeRecoveryMax;
}

void MetaTubeParams::
SetTubeRadiusParams( double _tubeRadiusStart,
  double _tubeRadiusMin, double _tubeRadiusMax,
  double _tubeRadiusThresholdMedialness,
  double _tubeRadiusThresholdMedialnessStart )
{
  m_TubeRadiusStart = _tubeRadiusStart;
  m_TubeRadiusMin = _tubeRadiusMin;
  m_TubeRadiusMax = _tubeRadiusMax;
  m_TubeRadiusThresholdMedialness = _tubeRadiusThresholdMedialness;
  m_TubeRadiusThresholdMedialnessStart = _tubeRadiusThresholdMedialnessStart;
}

vnl_vector<double> MetaTubeParams::
GetSeedScales( void ) const
{
  return m_SeedScales;
}

double MetaTubeParams::
GetSeedIntensityMin( void ) const
{
  return m_SeedIntensityMin;
}

double MetaTubeParams::
GetSeedIntensityMax( void ) const
{
  return m_SeedIntensityMax;
}

double MetaTubeParams::
GetSeedIntensityPercentile( void ) const
{
  return m_SeedIntensityPercentile;
}

double MetaTubeParams::
GetTubeIntensityMin( void ) const
{
  return m_TubeIntensityMin;
}

double MetaTubeParams::
GetTubeIntensityMax( void ) const
{
  return m_TubeIntensityMax;
}

bool MetaTubeParams::
GetTubeBright( void ) const
{
  return m_TubeBright;
}

vnl_vector<double> MetaTubeParams::
GetTubeColor( void ) const
{
  return m_TubeColor;
}

double MetaTubeParams::
GetTubeRidgeScale( void ) const
{
  return m_TubeRidgeScale;
}

double MetaTubeParams::
GetTubeRidgeScaleExtent( void ) const
{
  return m_TubeRidgeScaleExtent;
}

bool MetaTubeParams::
GetTubeRidgeDynamicScale( void ) const
{
  return m_TubeRidgeDynamicScale;
}

double MetaTubeParams::
GetTubeRidgeStepX( void ) const
{
  return m_TubeRidgeStepX;
}

double MetaTubeParams::
GetTubeRidgeThresholdTangentChange( void ) const
{
  return m_TubeRidgeThresholdTangentChange;
}

double MetaTubeParams::
GetTubeRidgeThresholdXChange( void ) const
{
  return m_TubeRidgeThresholdXChange;
}

double MetaTubeParams::
GetTubeRidgeThresholdRidgeness( void ) const
{
  return m_TubeRidgeThresholdRidgeness;
}

double MetaTubeParams::
GetTubeRidgeThresholdRidgenessStart( void ) const
{
  return m_TubeRidgeThresholdRidgenessStart;
}

double MetaTubeParams::
GetTubeRidgeThresholdRoundness( void ) const
{
  return m_TubeRidgeThresholdRoundness;
}

double MetaTubeParams::
GetTubeRidgeThresholdRoundnessStart( void ) const
{
  return m_TubeRidgeThresholdRoundnessStart;
}

double MetaTubeParams::
GetTubeRidgeCurvatureMax( void ) const
{
  return m_TubeRidgeCurvatureMax;
}

double MetaTubeParams::
GetTubeRidgeThresholdCurvature( void ) const
{
  return m_TubeRidgeThresholdCurvature;
}

double MetaTubeParams::
GetTubeRidgeThresholdCurvatureStart( void ) const
{
  return m_TubeRidgeThresholdCurvatureStart;
}

double MetaTubeParams::
GetTubeRidgeThresholdLinearity( void ) const
{
  return m_TubeRidgeThresholdLinearity;
}

double MetaTubeParams::
GetTubeRidgeThresholdLinearityStart( void ) const
{
  return m_TubeRidgeThresholdLinearityStart;
}

int MetaTubeParams::
GetTubeRidgeRecoveryMax( void ) const
{
  return m_TubeRidgeRecoveryMax;
}

double MetaTubeParams::
GetTubeRadiusStart( void ) const
{
  return m_TubeRadiusStart;
}

double MetaTubeParams::
GetTubeRadiusMin( void ) const
{
  return m_TubeRadiusMin;
}

double MetaTubeParams::
GetTubeRadiusMax( void ) const
{
  return m_TubeRadiusMax;
}

double MetaTubeParams::
GetTubeRadiusThresholdMedialness( void ) const
{
  return m_TubeRadiusThresholdMedialness;
}

double MetaTubeParams::
GetTubeRadiusThresholdMedialnessStart( void ) const
{
  return m_TubeRadiusThresholdMedialnessStart;
}

void MetaTubeParams::
CopyInfo( const MetaTubeParams & _tubeParams )
{
  MetaForm::CopyInfo( dynamic_cast< const MetaForm * >( & _tubeParams ) );

  SetSeedParams( _tubeParams.GetSeedScales(),
    _tubeParams.GetSeedIntensityMin(),
    _tubeParams.GetSeedIntensityMax(),
    _tubeParams.GetSeedIntensityPercentile() );

  SetTubeParams( _tubeParams.GetTubeIntensityMin(),
    _tubeParams.GetTubeIntensityMax(),
    _tubeParams.GetTubeBright(),
    _tubeParams.GetTubeColor() );

  SetTubeRidgeParams( _tubeParams.GetTubeRidgeScale(),
    _tubeParams.GetTubeRidgeScaleExtent(),
    _tubeParams.GetTubeRidgeDynamicScale(),
    _tubeParams.GetTubeRidgeStepX(),
    _tubeParams.GetTubeRidgeThresholdTangentChange(),
    _tubeParams.GetTubeRidgeThresholdXChange(),
    _tubeParams.GetTubeRidgeThresholdRidgeness(),
    _tubeParams.GetTubeRidgeThresholdRidgenessStart(),
    _tubeParams.GetTubeRidgeThresholdRoundness(),
    _tubeParams.GetTubeRidgeThresholdRoundnessStart(),
    _tubeParams.GetTubeRidgeCurvatureMax(),
    _tubeParams.GetTubeRidgeThresholdCurvature(),
    _tubeParams.GetTubeRidgeThresholdCurvatureStart(),
    _tubeParams.GetTubeRidgeThresholdLinearity(),
    _tubeParams.GetTubeRidgeThresholdLinearityStart(),
    _tubeParams.GetTubeRidgeRecoveryMax() );

  SetTubeRadiusParams( _tubeParams.GetTubeRadiusStart(),
    _tubeParams.GetTubeRadiusMin(),
    _tubeParams.GetTubeRadiusMax(),
    _tubeParams.GetTubeRadiusThresholdMedialness(),
    _tubeParams.GetTubeRadiusThresholdMedialnessStart() );
}

void MetaTubeParams::
Clear( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaTubeParams: Clear" << METAIO_STREAM::endl;
    }

  m_SeedScales.set_size( 1 );
  m_SeedScales[0] = 1;
  m_SeedIntensityMin = 0;
  m_SeedIntensityMax = 0;
  m_SeedIntensityPercentile = 0;

  m_TubeIntensityMin = 0;
  m_TubeIntensityMax = 0;
  m_TubeBright = true;
  m_TubeColor.set_size( 3 );
  m_TubeColor[0] = 1;
  m_TubeColor[1] = 0;
  m_TubeColor[2] = 0;

  m_TubeRidgeScale = 1;
  m_TubeRidgeScaleExtent = 2;
  m_TubeRidgeDynamicScale = true;
  m_TubeRidgeStepX = 0.2;
  m_TubeRidgeThresholdTangentChange = 0.8;
  m_TubeRidgeThresholdXChange = 0.8;
  m_TubeRidgeThresholdRidgeness = 0.8;
  m_TubeRidgeThresholdRidgenessStart = 0.8;
  m_TubeRidgeThresholdRoundness = 0.8;
  m_TubeRidgeThresholdRoundnessStart = 0.8;
  m_TubeRidgeCurvatureMax = 0.8;
  m_TubeRidgeThresholdCurvature = 0.8;
  m_TubeRidgeThresholdCurvatureStart = 0.8;
  m_TubeRidgeThresholdLinearity = 0.8;
  m_TubeRidgeThresholdLinearityStart = 0.8;
  m_TubeRidgeRecoveryMax = 3;

  m_TubeRadiusStart = 1;
  m_TubeRadiusMin = 0.5;
  m_TubeRadiusMax = 8.0;
  m_TubeRadiusThresholdMedialness = 0.8;
  m_TubeRadiusThresholdMedialnessStart = 0.8;

  MetaForm::Clear();
}

bool MetaTubeParams::
InitializeEssential( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaTubeParams: Initialize"
      << METAIO_STREAM::endl;
    }

  MetaForm::InitializeEssential();

  Clear();

  return true;
}

//
//
//
bool MetaTubeParams::
CanRead( const char *_headerName ) const
{
  // First check the extension
  METAIO_STL::string fname = _headerName;
  if( fname == "" )
    {
    return false;
    }

  bool extensionFound = false;

  METAIO_STL::string::size_type stringPos = fname.rfind( ".mtp" );
  if ( ( stringPos != METAIO_STL::string::npos )
      && ( stringPos == fname.length() - 5 ) )
    {
    extensionFound = true;
    }

  if( !extensionFound )
    {
    return false;
    }

  // Now check the file content
  METAIO_STREAM::ifstream inputStream;

  inputStream.open( _headerName, METAIO_STREAM::ios::in |
                                 METAIO_STREAM::ios::binary );

  if( !inputStream.rdbuf()->is_open() )
    {
    return false;
    }

  bool result = !strncmp( MET_ReadForm( inputStream ).c_str(), "TubeParams", 10 );

  inputStream.close();

  return result;
}


bool MetaTubeParams::
Read( const char *_headerName )
{
  if( _headerName != NULL && strlen( _headerName ) > 1 )
    {
    FileName( _headerName );
    }

  METAIO_STREAM::ifstream * tmpStream = new METAIO_STREAM::ifstream;

  tmpStream->open( m_FileName, METAIO_STREAM::ios::in |
                              METAIO_STREAM::ios::binary );

  if( !tmpStream->rdbuf()->is_open() )
    {
    METAIO_STREAM::cout << "MetaTubeParams: Read: Cannot open file _"
                        << m_FileName << "_" << METAIO_STREAM::endl;
    delete tmpStream;
    return false;
    }

  bool result = ReadStream( tmpStream );

  tmpStream->close();

  delete tmpStream;

  return result;
}


bool MetaTubeParams::
CanReadStream( METAIO_STREAM::ifstream * _stream ) const
{
  if( !strncmp( MET_ReadForm( *_stream ).c_str(), "TubeParams", 10 ) )
    {
    return true;
    }

  return false;
}

bool MetaTubeParams::
ReadStream( METAIO_STREAM::ifstream * _stream )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaTubeParams: ReadStream"
      << METAIO_STREAM::endl;
    }

  M_Destroy();

  Clear();

  M_SetupReadFields();

  if( m_ReadStream )
    {
    METAIO_STREAM::cout << "MetaTubeParams: ReadStream: two files open?"
                        << METAIO_STREAM::endl;
    delete m_ReadStream;
    }

  m_ReadStream = _stream;

  if( !M_Read() )
    {
    METAIO_STREAM::cout << "MetaTubeParams: Read: Cannot parse file"
                        << METAIO_STREAM::endl;
    m_ReadStream = NULL;
    return false;
    }

  m_ReadStream = NULL;

  return true;
}

//
//
//
//

bool MetaTubeParams::
Write( const char *_headName )
{
  if( _headName != NULL && strlen( _headName ) > 1 )
    {
    FileName( _headName );
    }

  MET_SetFileSuffix( m_FileName, "mtp" );

  METAIO_STREAM::ofstream * tmpWriteStream = new METAIO_STREAM::ofstream;

  tmpWriteStream->open( m_FileName, METAIO_STREAM::ios::binary |
                                   METAIO_STREAM::ios::out );

  if( !tmpWriteStream->rdbuf()->is_open() )
    {
    delete tmpWriteStream;
    return false;
    }

  bool result = WriteStream( tmpWriteStream );

  tmpWriteStream->close();

  delete tmpWriteStream;

  return result;
}

bool MetaTubeParams::
WriteStream( METAIO_STREAM::ofstream * _stream )
{
  if( m_WriteStream != NULL )
    {
    METAIO_STREAM::cout << "MetaTubeParams: WriteStream: two files open?"
                        << METAIO_STREAM::endl;
    delete m_WriteStream;
    }

  m_WriteStream = _stream;

  M_SetupWriteFields();

  M_Write();

  m_WriteStream->flush();

  m_WriteStream = NULL;

  return true;
}

void MetaTubeParams::
M_Destroy( void )
{
  MetaForm::M_Destroy();
}

void MetaTubeParams::
M_SetupReadFields( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaTubeParams: M_SetupReadFields"
                        << METAIO_STREAM::endl;
    }

  MetaForm::M_SetupReadFields();

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "NSeedScales", MET_INT, true );
  m_Fields.push_back( mF );

  int nSeedScales = MET_GetFieldRecordNumber( "NSeedScales", &m_Fields );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "SeedScales", MET_FLOAT_ARRAY, true, nSeedScales );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "SeedIntensityMin", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "SeedIntensityMax", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "SeedIntensityPercentile", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeIntensityMin", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeIntensityMax", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeBright", MET_STRING, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeColor", MET_FLOAT_ARRAY, true, 3 );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeScale", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeScaleExtent", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeDynamicScale", MET_STRING, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeStepX", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeThresholdTangentChange", MET_FLOAT,
    true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeThresholdXChange", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeThresholdRidgeness", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeThresholdRidgenessStart", MET_FLOAT,
    true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeThresholdRoundness", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeThresholdRoundnessStart", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeCurvatureMax", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeThresholdCurvature", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeThresholdCurvatureStart", MET_FLOAT,
    true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeThresholdLinearity", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeThresholdLinearityStart", MET_FLOAT,
    true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRidgeRecoveryMax", MET_INT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRadiusStart", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRadiusMin", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRadiusMax", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRadiusThresholdMedialness", MET_FLOAT,
    true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeRadiusThresholdMedialnessStart", MET_FLOAT,
    true );
  m_Fields.push_back( mF );

}

void MetaTubeParams::
M_SetupWriteFields( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaTubeParams: M_SetupWriteFields"
                        << METAIO_STREAM::endl;
    }

  strcpy( m_FormTypeName, "TubeParams" );
  MetaForm::M_SetupWriteFields();

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "NSeedScales", MET_INT, m_SeedScales.size() );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "SeedScales", MET_FLOAT_ARRAY,
    m_SeedScales.size(), m_SeedScales.data_block() );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "SeedIntensityMin", MET_FLOAT,
    m_SeedIntensityMin );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "SeedIntensityMax", MET_FLOAT,
    m_SeedIntensityMax );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "SeedIntensityPercentile", MET_FLOAT,
    m_SeedIntensityPercentile );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeIntensityMin", MET_FLOAT,
    m_TubeIntensityMin );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeIntensityMax", MET_FLOAT,
    m_TubeIntensityMax );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  if( m_TubeBright )
    {
    MET_InitWriteField( mF, "TubeBright", MET_STRING, 5, "True" );
    }
  else
    {
    MET_InitWriteField( mF, "TubeBright", MET_STRING, 6, "False" );
    }
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeColor", MET_FLOAT_ARRAY, 3,
    m_TubeColor.data_block() );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeScale", MET_FLOAT, m_TubeRidgeScale );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeScaleExtent", MET_FLOAT,
    m_TubeRidgeScaleExtent );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  if( m_TubeRidgeDynamicScale )
    {
    MET_InitWriteField( mF, "TubeRidgeDynamicScale", MET_STRING, 5, "True" );
    }
  else
    {
    MET_InitWriteField( mF, "TubeRidgeDynamicScale", MET_STRING, 6, "False" );
    }
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeStepX", MET_FLOAT, m_TubeRidgeStepX );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeThresholdTangentChange", MET_FLOAT,
    m_TubeRidgeThresholdTangentChange );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeThresholdXChange", MET_FLOAT,
    m_TubeRidgeThresholdXChange );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeThresholdRidgeness", MET_FLOAT,
    m_TubeRidgeThresholdRidgeness );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeThresholdRidgenessStart", MET_FLOAT,
    m_TubeRidgeThresholdRidgenessStart );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeThresholdRoundness", MET_FLOAT,
    m_TubeRidgeThresholdRoundness );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeThresholdRoundnessStart", MET_FLOAT,
    m_TubeRidgeThresholdRoundnessStart );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeCurvatureMax", MET_FLOAT,
    m_TubeRidgeCurvatureMax );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeThresholdCurvature", MET_FLOAT,
    m_TubeRidgeThresholdCurvature );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeThresholdCurvatureStart", MET_FLOAT,
    m_TubeRidgeThresholdCurvatureStart );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeThresholdLinearity", MET_FLOAT,
    m_TubeRidgeThresholdLinearity );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeThresholdLinearityStart", MET_FLOAT,
    m_TubeRidgeThresholdLinearityStart );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRidgeRecoveryMax", MET_INT,
    m_TubeRidgeRecoveryMax );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRadiusStart", MET_FLOAT,
    m_TubeRadiusStart );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRadiusMin", MET_FLOAT, m_TubeRadiusMin );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRadiusMax", MET_FLOAT, m_TubeRadiusMax );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRadiusThresholdMedialness", MET_FLOAT,
    m_TubeRadiusThresholdMedialness );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeRadiusThresholdMedialnessStart", MET_FLOAT,
    m_TubeRadiusThresholdMedialnessStart );
  m_Fields.push_back( mF );
}


bool MetaTubeParams::
M_Read( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaTubeParams: M_Read: Loading Header"
                        << METAIO_STREAM::endl;
    }
  if( !MetaForm::M_Read() )
    {
    METAIO_STREAM::cout << "MetaTubeParams: M_Read: Error parsing file"
                        << METAIO_STREAM::endl;
    return false;
    }

  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaTubeParams: M_Read: Parsing Header"
                        << METAIO_STREAM::endl;
    }
  MET_FieldRecordType * mF;

  unsigned int nSeedScales = 0;
  mF = MET_GetFieldRecord( "NSeedScales", &m_Fields );
  nSeedScales = ( unsigned int )mF->value[0];
  m_SeedScales.set_size( nSeedScales );

  mF = MET_GetFieldRecord( "SeedScales", &m_Fields );
  for( unsigned int i=0; i<nSeedScales; i++ )
    {
    m_SeedScales[i] = ( double )mF->value[i];
    }

  mF = MET_GetFieldRecord( "SeedIntensityMin", &m_Fields );
  m_SeedIntensityMin = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "SeedIntensityMax", &m_Fields );
  m_SeedIntensityMax = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "SeedIntensityPercentile", &m_Fields );
  m_SeedIntensityPercentile = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeIntensityMin", &m_Fields );
  m_TubeIntensityMin = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeIntensityMax", &m_Fields );
  m_TubeIntensityMax = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeBright", &m_Fields );
  if( (char)(mF->value[0]) == 't' || (char)(mF->value[0]) == 'T' )
    {
    m_TubeBright = true;
    }
  else
    {
    m_TubeBright = false;
    }

  mF = MET_GetFieldRecord( "TubeColor", &m_Fields );
  m_TubeColor[0] = ( double )mF->value[0];
  m_TubeColor[1] = ( double )mF->value[1];
  m_TubeColor[2] = ( double )mF->value[2];

  mF = MET_GetFieldRecord( "TubeRidgeScale", &m_Fields );
  m_TubeRidgeScale = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeScaleExtent", &m_Fields );
  m_TubeRidgeScaleExtent = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeDynamicScale", &m_Fields );
  if( (char)(mF->value[0]) == 't' || (char)(mF->value[0]) == 'T' )
    {
    m_TubeRidgeDynamicScale = true;
    }
  else
    {
    m_TubeRidgeDynamicScale = false;
    }


  mF = MET_GetFieldRecord( "TubeRidgeStepX", &m_Fields );
  m_TubeRidgeStepX = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeThresholdTangentChange", &m_Fields );
  m_TubeRidgeThresholdTangentChange = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeThresholdXChange", &m_Fields );
  m_TubeRidgeThresholdXChange = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeThresholdRidgeness", &m_Fields );
  m_TubeRidgeThresholdRidgeness = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeThresholdRidgenessStart", &m_Fields );
  m_TubeRidgeThresholdRidgenessStart = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeThresholdRoundness", &m_Fields );
  m_TubeRidgeThresholdRoundness = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeThresholdRoundnessStart", &m_Fields );
  m_TubeRidgeThresholdRoundnessStart = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeCurvatureMax", &m_Fields );
  m_TubeRidgeCurvatureMax = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeThresholdCurvature", &m_Fields );
  m_TubeRidgeThresholdCurvature = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeThresholdCurvatureStart", &m_Fields );
  m_TubeRidgeThresholdCurvatureStart = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeThresholdLinearity", &m_Fields );
  m_TubeRidgeThresholdLinearity = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeThresholdLinearityStart", &m_Fields );
  m_TubeRidgeThresholdLinearityStart = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRidgeRecoveryMax", &m_Fields );
  m_TubeRidgeRecoveryMax = ( int )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRadiusStart", &m_Fields );
  m_TubeRadiusStart = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRadiusMin", &m_Fields );
  m_TubeRadiusMin = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRadiusMax", &m_Fields );
  m_TubeRadiusMax = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRadiusThresholdMedialness", &m_Fields );
  m_TubeRadiusThresholdMedialness = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeRadiusThresholdMedialnessStart",
    &m_Fields );
  m_TubeRadiusThresholdMedialnessStart = ( double )mF->value[0];

  return true;
}

}

}
