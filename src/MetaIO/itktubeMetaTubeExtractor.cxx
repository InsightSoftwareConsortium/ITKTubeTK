/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include <cstring>

#include "tubeStringUtilities.h"
#include "itktubeMetaTubeExtractor.h"

namespace itk
{

namespace tube
{

MetaTubeExtractor::
MetaTubeExtractor( void )
{
  if( META_DEBUG )
    {
    std::cout << "MetaTubeExtractor()" << std::endl;
    }

  Clear();
}

//
MetaTubeExtractor::
MetaTubeExtractor( const char *_headerName )
{
  if( META_DEBUG )
    {
    std::cout << "MetaTubeExtractor()" << std::endl;
    }

  Clear();
  m_ReadStream = NULL;

  MetaTubeExtractor::Read( _headerName );
}

//
MetaTubeExtractor::
MetaTubeExtractor( const MetaTubeExtractor & _metaTubeExtractor ) : MetaForm()
{
  if( META_DEBUG )
   {
   std::cout << "MetaTubeExtractor()" << std::endl;
   }

  Clear();

  CopyInfo( _metaTubeExtractor );
}

//
MetaTubeExtractor::
~MetaTubeExtractor( void )
{
  M_Destroy();
}

//
void MetaTubeExtractor::
PrintInfo( void ) const
{
  MetaForm::PrintInfo();

  std::cout << "DataMin = "
    << m_DataMin << std::endl;

  std::cout << "DataMax = "
    << m_DataMax << std::endl;

  std::cout << "TubeColor = "
    << m_TubeColor << std::endl;

  std::cout << "RidgeScale = "
    << m_RidgeScale << std::endl;

  std::cout << "RidgeScaleKernelExtent = "
    << m_RidgeScaleKernelExtent << std::endl;

  std::cout << "RidgeDynamicScale = "
    << m_RidgeDynamicScale << std::endl;

  std::cout << "RidgeDynamicStepSize = "
    << m_RidgeDynamicStepSize << std::endl;

  std::cout << "RidgeStepX = "
    << m_RidgeStepX << std::endl;

  std::cout << "RidgeMaxTangentChange = "
    << m_RidgeMaxTangentChange << std::endl;

  std::cout << "RidgeMaxXChange = "
    << m_RidgeMaxXChange << std::endl;

  std::cout << "RidgeMinRidgeness = "
    << m_RidgeMinRidgeness << std::endl;

  std::cout << "RidgeMinRidgenessStart = "
    << m_RidgeMinRidgenessStart << std::endl;

  std::cout << "RidgeMinRoundness = "
    << m_RidgeMinRoundness << std::endl;

  std::cout << "RidgeMinRoundnessStart = "
    << m_RidgeMinRoundnessStart << std::endl;

  std::cout << "RidgeMinCurvature = "
    << m_RidgeMinCurvature << std::endl;

  std::cout << "RidgeMinCurvatureStart = "
    << m_RidgeMinCurvatureStart << std::endl;

  std::cout << "RidgeMinLevelness = "
    << m_RidgeMinLevelness << std::endl;

  std::cout << "RidgeMinLevelnessStart = "
    << m_RidgeMinLevelnessStart << std::endl;

  std::cout << "RidgeMaxRecoveryAttempts = "
    << m_RidgeMaxRecoveryAttempts << std::endl;

  std::cout << "RadiusStart = "
    << m_RadiusStart << std::endl;

  std::cout << "RadiusMin = "
    << m_RadiusMin << std::endl;

  std::cout << "RadiusMax = "
    << m_RadiusMax << std::endl;

  std::cout << "RadiusMinMedialness = "
    << m_RadiusMinMedialness << std::endl;

  std::cout << "RadiusMinMedialnessStart = "
    << m_RadiusMinMedialnessStart << std::endl;

}

void MetaTubeExtractor::
SetGeneralProperties( double _dataMin, double _dataMax,
  const VectorType & _tubeColor )
{
  m_DataMin = _dataMin;
  m_DataMax = _dataMax;
  m_TubeColor = _tubeColor;
}

void MetaTubeExtractor::
SetRidgeProperties( double _ridgeScale, double _ridgeScaleKernelExtent,
  bool _ridgeDynamicScale,
  bool _ridgeDynamicStepSize,
  double _ridgeStepX,
  double _ridgeMaxTangentChange,
  double _ridgeMaxXChange,
  double _ridgeMinRidgeness,
  double _ridgeMinRidgenessStart,
  double _ridgeMinRoundness,
  double _ridgeMinRoundnessStart,
  double _ridgeMinCurvature,
  double _ridgeMinCurvatureStart,
  double _ridgeMinLevelness,
  double _ridgeMinLevelnessStart,
  int _ridgeMaxRecoveryAttempts )
{
  m_RidgeScale = _ridgeScale;
  m_RidgeScaleKernelExtent = _ridgeScaleKernelExtent;
  m_RidgeDynamicScale = _ridgeDynamicScale;
  m_RidgeDynamicStepSize = _ridgeDynamicStepSize;
  m_RidgeStepX = _ridgeStepX;
  m_RidgeMaxTangentChange = _ridgeMaxTangentChange;
  m_RidgeMaxXChange = _ridgeMaxXChange;
  m_RidgeMinRidgeness = _ridgeMinRidgeness;
  m_RidgeMinRidgenessStart = _ridgeMinRidgenessStart;
  m_RidgeMinRoundness = _ridgeMinRoundness;
  m_RidgeMinRoundnessStart = _ridgeMinRoundnessStart;
  m_RidgeMinCurvature = _ridgeMinCurvature;
  m_RidgeMinCurvatureStart = _ridgeMinCurvatureStart;
  m_RidgeMinLevelness = _ridgeMinLevelness;
  m_RidgeMinLevelnessStart = _ridgeMinLevelnessStart;
  m_RidgeMaxRecoveryAttempts = _ridgeMaxRecoveryAttempts;
}

void MetaTubeExtractor::
SetRadiusProperties( double _radiusStart,
  double _radiusMin,
  double _radiusMax,
  double _radiusMinMedialness,
  double _radiusMinMedialnessStart )
{
  m_RadiusStart = _radiusStart;
  m_RadiusMin = _radiusMin;
  m_RadiusMax = _radiusMax;
  m_RadiusMinMedialness = _radiusMinMedialness;
  m_RadiusMinMedialnessStart = _radiusMinMedialnessStart;
}

double MetaTubeExtractor::
GetDataMin( void ) const
{
  return m_DataMin;
}

double MetaTubeExtractor::
GetDataMax( void ) const
{
  return m_DataMax;
}

vnl_vector<double> MetaTubeExtractor::
GetTubeColor( void ) const
{
  return m_TubeColor;
}

double MetaTubeExtractor::
GetRidgeScale( void ) const
{
  return m_RidgeScale;
}

double MetaTubeExtractor::
GetRidgeScaleKernelExtent( void ) const
{
  return m_RidgeScaleKernelExtent;
}

bool MetaTubeExtractor::
GetRidgeDynamicScale( void ) const
{
  return m_RidgeDynamicScale;
}

bool MetaTubeExtractor::
GetRidgeDynamicStepSize( void ) const
{
  return m_RidgeDynamicStepSize;
}

double MetaTubeExtractor::
GetRidgeStepX( void ) const
{
  return m_RidgeStepX;
}

double MetaTubeExtractor::
GetRidgeMaxTangentChange( void ) const
{
  return m_RidgeMaxTangentChange;
}

double MetaTubeExtractor::
GetRidgeMaxXChange( void ) const
{
  return m_RidgeMaxXChange;
}

double MetaTubeExtractor::
GetRidgeMinRidgeness( void ) const
{
  return m_RidgeMinRidgeness;
}

double MetaTubeExtractor::
GetRidgeMinRidgenessStart( void ) const
{
  return m_RidgeMinRidgenessStart;
}

double MetaTubeExtractor::
GetRidgeMinRoundness( void ) const
{
  return m_RidgeMinRoundness;
}

double MetaTubeExtractor::
GetRidgeMinRoundnessStart( void ) const
{
  return m_RidgeMinRoundnessStart;
}

double MetaTubeExtractor::
GetRidgeMinCurvature( void ) const
{
  return m_RidgeMinCurvature;
}

double MetaTubeExtractor::
GetRidgeMinCurvatureStart( void ) const
{
  return m_RidgeMinCurvatureStart;
}

double MetaTubeExtractor::
GetRidgeMinLevelness( void ) const
{
  return m_RidgeMinLevelness;
}

double MetaTubeExtractor::
GetRidgeMinLevelnessStart( void ) const
{
  return m_RidgeMinLevelnessStart;
}

int MetaTubeExtractor::
GetRidgeMaxRecoveryAttempts( void ) const
{
  return m_RidgeMaxRecoveryAttempts;
}

double MetaTubeExtractor::
GetRadiusStart( void ) const
{
  return m_RadiusStart;
}

double MetaTubeExtractor::
GetRadiusMin( void ) const
{
  return m_RadiusMin;
}

double MetaTubeExtractor::
GetRadiusMax( void ) const
{
  return m_RadiusMax;
}

double MetaTubeExtractor::
GetRadiusMinMedialness( void ) const
{
  return m_RadiusMinMedialness;
}

double MetaTubeExtractor::
GetRadiusMinMedialnessStart( void ) const
{
  return m_RadiusMinMedialnessStart;
}

void MetaTubeExtractor::
CopyInfo( const MetaTubeExtractor & _tubeExtractor )
{
  MetaForm::CopyInfo( dynamic_cast< const MetaForm * >( & _tubeExtractor ) );

  SetGeneralProperties( _tubeExtractor.GetDataMin(),
    _tubeExtractor.GetDataMax(),
    _tubeExtractor.GetTubeColor() );

  SetRidgeProperties( _tubeExtractor.GetRidgeScale(),
    _tubeExtractor.GetRidgeScaleKernelExtent(),
    _tubeExtractor.GetRidgeDynamicScale(),
    _tubeExtractor.GetRidgeDynamicStepSize(),
    _tubeExtractor.GetRidgeStepX(),
    _tubeExtractor.GetRidgeMaxTangentChange(),
    _tubeExtractor.GetRidgeMaxXChange(),
    _tubeExtractor.GetRidgeMinRidgeness(),
    _tubeExtractor.GetRidgeMinRidgenessStart(),
    _tubeExtractor.GetRidgeMinRoundness(),
    _tubeExtractor.GetRidgeMinRoundnessStart(),
    _tubeExtractor.GetRidgeMinCurvature(),
    _tubeExtractor.GetRidgeMinCurvatureStart(),
    _tubeExtractor.GetRidgeMinLevelness(),
    _tubeExtractor.GetRidgeMinLevelnessStart(),
    _tubeExtractor.GetRidgeMaxRecoveryAttempts() );

  SetRadiusProperties( _tubeExtractor.GetRadiusStart(),
    _tubeExtractor.GetRadiusMin(),
    _tubeExtractor.GetRadiusMax(),
    _tubeExtractor.GetRadiusMinMedialness(),
    _tubeExtractor.GetRadiusMinMedialnessStart() );
}

void MetaTubeExtractor::
Clear( void )
{
  if( META_DEBUG )
    {
    std::cout << "MetaTubeExtractor: Clear"
      << std::endl;
    }

  m_DataMin = 0;
  m_DataMax = 0;
  m_TubeColor.set_size( 4 );
  m_TubeColor[0] = 1;
  m_TubeColor[1] = 0;
  m_TubeColor[2] = 0;
  m_TubeColor[3] = 1;

  m_RidgeScale = 1;
  m_RidgeScaleKernelExtent = 2;
  m_RidgeDynamicScale = true;
  m_RidgeDynamicStepSize = true;
  m_RidgeStepX = 0.2;
  m_RidgeMaxTangentChange = 0.8;
  m_RidgeMaxXChange = 0.8;
  m_RidgeMinRidgeness = 0.8;
  m_RidgeMinRidgenessStart = 0.8;
  m_RidgeMinRoundness = 0.8;
  m_RidgeMinRoundnessStart = 0.8;
  m_RidgeMinCurvature = 0.8;
  m_RidgeMinCurvatureStart = 0.8;
  m_RidgeMinLevelness = 0.8;
  m_RidgeMinLevelnessStart = 0.8;
  m_RidgeMaxRecoveryAttempts = 3;

  m_RadiusStart = 1;
  m_RadiusMin = 0.5;
  m_RadiusMax = 8.0;
  m_RadiusMinMedialness = 0.8;
  m_RadiusMinMedialnessStart = 0.8;

  MetaForm::Clear();
}

bool MetaTubeExtractor::
InitializeEssential( void )
{
  if( META_DEBUG )
    {
    std::cout << "MetaTubeExtractor: Initialize"
      << std::endl;
    }

  MetaForm::InitializeEssential();

  Clear();

  return true;
}

//
//
//
bool MetaTubeExtractor::
CanRead( const char *_headerName ) const
{
  // First check the extension
  std::string fname = _headerName;
  if( fname == "" )
    {
    return false;
    }

  bool extensionFound = false;

  std::string::size_type stringPos = fname.rfind( ".mtp" );
  if( ( stringPos != std::string::npos )
      && ( stringPos == fname.length() - 4 ) )
    {
    extensionFound = true;
    }

  if( !extensionFound )
    {
    return false;
    }

  // Now check the file content
  std::ifstream inputStream;

  inputStream.open( _headerName, std::ios::in |
                                 std::ios::binary );

  if( !inputStream.rdbuf()->is_open() )
    {
    return false;
    }

  bool result = !std::strncmp( MET_ReadForm( inputStream ).c_str(),
    "TubeExtractor", 13 );

  inputStream.close();

  return result;
}


bool MetaTubeExtractor::
Read( const char *_headerName )
{
  if( _headerName != NULL && std::strlen( _headerName ) > 1 )
    {
    FileName( _headerName );
    }

  std::ifstream * tmpStream = new std::ifstream;

  tmpStream->open( m_FileName, std::ios::in |
                              std::ios::binary );

  if( !tmpStream->rdbuf()->is_open() )
    {
    std::cout << "MetaTubeExtractor: Read: Cannot open file _"
                        << m_FileName << "_" << std::endl;
    delete tmpStream;
    return false;
    }

  bool result = ReadStream( tmpStream );

  tmpStream->close();

  delete tmpStream;

  return result;
}


bool MetaTubeExtractor::
CanReadStream( std::ifstream * _stream ) const
{
  if( !std::strncmp( MET_ReadForm( *_stream ).c_str(),
    "TubeExtractor", 10 ) )
    {
    return true;
    }

  return false;
}

bool MetaTubeExtractor::
ReadStream( std::ifstream * _stream )
{
  if( META_DEBUG )
    {
    std::cout << "MetaTubeExtractor: ReadStream"
      << std::endl;
    }

  M_Destroy();

  Clear();

  M_SetupReadFields();

  if( m_ReadStream )
    {
    std::cout << "MetaTubeExtractor: ReadStream: two files open?"
                        << std::endl;
    delete m_ReadStream;
    }

  m_ReadStream = _stream;

  if( !M_Read() )
    {
    std::cout << "MetaTubeExtractor: Read: Cannot parse file"
                        << std::endl;
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

bool MetaTubeExtractor::
Write( const char *_headName )
{
  if( _headName != NULL && std::strlen( _headName ) > 1 )
    {
    FileName( _headName );
    }

  MET_SetFileSuffix( m_FileName, ".mtp" );

  std::ofstream * tmpWriteStream = new std::ofstream;

  tmpWriteStream->open( m_FileName, std::ios::binary |
                                   std::ios::out );

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

bool MetaTubeExtractor::
WriteStream( std::ofstream * _stream )
{
  if( m_WriteStream != NULL )
    {
    std::cout << "MetaTubeExtractor: WriteStream: two files open?"
                        << std::endl;
    delete m_WriteStream;
    }

  m_WriteStream = _stream;

  M_SetupWriteFields();

  M_Write();

  m_WriteStream->flush();

  m_WriteStream = NULL;

  return true;
}

void MetaTubeExtractor::
M_Destroy( void )
{
  MetaForm::M_Destroy();
}

void MetaTubeExtractor::
M_SetupReadFields( void )
{
  if( META_DEBUG )
    {
    std::cout << "MetaTubeExtractor: M_SetupReadFields"
                        << std::endl;
    }

  MetaForm::M_SetupReadFields();

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "DataMin", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "DataMax", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TubeColor", MET_STRING, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeScale", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeScaleKernelExtent", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeDynamicScale", MET_STRING, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeDynamicStepSize", MET_STRING, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeStepX", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeMaxTangentChange", MET_FLOAT,
    true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeMaxXChange", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeMinRidgeness", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeMinRidgenessStart", MET_FLOAT,
    true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeMinRoundness", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeMinRoundnessStart", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeMinCurvature", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeMinCurvatureStart", MET_FLOAT,
    true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeMinLevelness", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeMinLevelnessStart", MET_FLOAT,
    true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeMaxRecoveryAttempts", MET_INT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RadiusStart", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RadiusMin", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RadiusMax", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RadiusMinMedialness", MET_FLOAT,
    true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RadiusMinMedialnessStart", MET_FLOAT,
    true );
  m_Fields.push_back( mF );

}

void MetaTubeExtractor::
M_SetupWriteFields( void )
{
  if( META_DEBUG )
    {
    std::cout << "MetaTubeExtractor: M_SetupWriteFields"
                        << std::endl;
    }

  strcpy( m_FormTypeName, "TubeExtractor" );

  MetaForm::M_SetupWriteFields();

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "DataMin", MET_FLOAT, m_DataMin );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "DataMax", MET_FLOAT, m_DataMax );
  m_Fields.push_back( mF );

  char colorString[80];
  snprintf( colorString, 79, "%f %f %f %f", m_TubeColor[0], m_TubeColor[1],
    m_TubeColor[2], m_TubeColor[3] );
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "TubeColor", MET_STRING,
    std::strlen( colorString ), colorString );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeScale", MET_FLOAT, m_RidgeScale );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeScaleKernelExtent", MET_FLOAT,
    m_RidgeScaleKernelExtent );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  if( m_RidgeDynamicScale )
    {
    MET_InitWriteField( mF, "RidgeDynamicScale", MET_STRING, 4, "True" );
    }
  else
    {
    MET_InitWriteField( mF, "RidgeDynamicScale", MET_STRING, 5, "False" );
    }
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  if( m_RidgeDynamicStepSize )
    {
    MET_InitWriteField( mF, "RidgeDynamicStepSize", MET_STRING, 4,
      "True" );
    }
  else
    {
    MET_InitWriteField( mF, "RidgeDynamicStepSize", MET_STRING, 5,
      "False" );
    }
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeStepX", MET_FLOAT, m_RidgeStepX );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeMaxTangentChange", MET_FLOAT,
    m_RidgeMaxTangentChange );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeMaxXChange", MET_FLOAT,
    m_RidgeMaxXChange );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeMinRidgeness", MET_FLOAT,
    m_RidgeMinRidgeness );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeMinRidgenessStart", MET_FLOAT,
    m_RidgeMinRidgenessStart );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeMinRoundness", MET_FLOAT,
    m_RidgeMinRoundness );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeMinRoundnessStart", MET_FLOAT,
    m_RidgeMinRoundnessStart );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeMinCurvature", MET_FLOAT,
    m_RidgeMinCurvature );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeMinCurvatureStart", MET_FLOAT,
    m_RidgeMinCurvatureStart );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeMinLevelness", MET_FLOAT,
    m_RidgeMinLevelness );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeMinLevelnessStart", MET_FLOAT,
    m_RidgeMinLevelnessStart );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RidgeMaxRecoveryAttempts", MET_INT,
    m_RidgeMaxRecoveryAttempts );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RadiusStart", MET_FLOAT,
    m_RadiusStart );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RadiusMin", MET_FLOAT, m_RadiusMin );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RadiusMax", MET_FLOAT, m_RadiusMax );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RadiusMinMedialness", MET_FLOAT,
    m_RadiusMinMedialness );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "RadiusMinMedialnessStart", MET_FLOAT,
    m_RadiusMinMedialnessStart );
  m_Fields.push_back( mF );
}


bool MetaTubeExtractor::
M_Read( void )
{
  if( META_DEBUG )
    {
    std::cout << "MetaTubeExtractor: M_Read: Loading Header"
                        << std::endl;
    }
  if( !MetaForm::M_Read() )
    {
    std::cout << "MetaTubeExtractor: M_Read: Error parsing file"
                        << std::endl;
    return false;
    }

  if( META_DEBUG )
    {
    std::cout << "MetaTubeExtractor: M_Read: Parsing Header"
                        << std::endl;
    }
  MET_FieldRecordType * mF;

  mF = MET_GetFieldRecord( "DataMin", &m_Fields );
  m_DataMin = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "DataMax", &m_Fields );
  m_DataMax = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "TubeColor", &m_Fields );
  std::string str = ( char * )( mF->value );
  std::vector< double > clr( 4 );
  ::tube::StringToVector( str, clr, " " );
  m_TubeColor[0] = clr[0];
  m_TubeColor[1] = clr[1];
  m_TubeColor[2] = clr[2];
  m_TubeColor[3] = clr[3];

  mF = MET_GetFieldRecord( "RidgeScale", &m_Fields );
  m_RidgeScale = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeScaleKernelExtent", &m_Fields );
  m_RidgeScaleKernelExtent = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeDynamicScale", &m_Fields );
  if( ( char )( mF->value[0] ) == 't' || ( char )( mF->value[0] ) == 'T' )
    {
    m_RidgeDynamicScale = true;
    }
  else
    {
    m_RidgeDynamicScale = false;
    }

  mF = MET_GetFieldRecord( "RidgeDynamicStepSize", &m_Fields );
  if( ( char )( mF->value[0] ) == 't' || ( char )( mF->value[0] ) == 'T' )
    {
    m_RidgeDynamicStepSize = true;
    }
  else
    {
    m_RidgeDynamicStepSize = false;
    }

  mF = MET_GetFieldRecord( "RidgeStepX", &m_Fields );
  m_RidgeStepX = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeMaxTangentChange", &m_Fields );
  m_RidgeMaxTangentChange = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeMaxXChange", &m_Fields );
  m_RidgeMaxXChange = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeMinRidgeness", &m_Fields );
  m_RidgeMinRidgeness = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeMinRidgenessStart", &m_Fields );
  m_RidgeMinRidgenessStart = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeMinRoundness", &m_Fields );
  m_RidgeMinRoundness = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeMinRoundnessStart", &m_Fields );
  m_RidgeMinRoundnessStart = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeMinCurvature", &m_Fields );
  m_RidgeMinCurvature = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeMinCurvatureStart", &m_Fields );
  m_RidgeMinCurvatureStart = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeMinLevelness", &m_Fields );
  m_RidgeMinLevelness = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeMinLevelnessStart", &m_Fields );
  m_RidgeMinLevelnessStart = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RidgeMaxRecoveryAttempts", &m_Fields );
  m_RidgeMaxRecoveryAttempts = ( int )mF->value[0];

  mF = MET_GetFieldRecord( "RadiusStart", &m_Fields );
  m_RadiusStart = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RadiusMin", &m_Fields );
  m_RadiusMin = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RadiusMax", &m_Fields );
  m_RadiusMax = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RadiusMinMedialness", &m_Fields );
  m_RadiusMinMedialness = ( double )mF->value[0];

  mF = MET_GetFieldRecord( "RadiusMinMedialnessStart",
    &m_Fields );
  m_RadiusMinMedialnessStart = ( double )mF->value[0];

  return true;
}

}

}
