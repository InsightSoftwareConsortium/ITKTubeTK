/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#include <cstring>

#include "itktubeMetaRidgeSeed.h"

namespace itk
{

namespace tube
{

MetaRidgeSeed::
MetaRidgeSeed( void )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed()" << std::endl;
    }

  Clear();
}

MetaRidgeSeed::
MetaRidgeSeed( const char * _headerName )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed()" << std::endl;
    }

  Clear();

  MetaRidgeSeed::Read( _headerName );
}

MetaRidgeSeed::
MetaRidgeSeed( const MetaRidgeSeed & _metaRidgeSeed )
: MetaLDA()
{
  if( META_DEBUG )
   {
   std::cout << "MetaRidgeSeed()" << std::endl;
   }

  Clear();

  CopyInfo( _metaRidgeSeed );
}

MetaRidgeSeed::
MetaRidgeSeed(
  const RidgeSeedScalesType & _ridgeSeedScales,
  bool _useIntensityOnly,
  bool _useFeatureMath,
  const LDAValuesType & _ldaValues,
  const LDAMatrixType & _ldaMatrix,
  const ValueListType & _inputWhitenMeans,
  const ValueListType & _inputWhitenStdDevs,
  const ValueListType & _outputWhitenMeans,
  const ValueListType & _outputWhitenStdDevs,
  const std::string & _pdfFileName )
{
  if( META_DEBUG )
   {
   std::cout << "MetaRidgeSeed()" << std::endl;
   }

  Clear();

  InitializeEssential( _ridgeSeedScales, _useIntensityOnly, _useFeatureMath,
    _ldaValues, _ldaMatrix, _inputWhitenMeans,
    _inputWhitenStdDevs, _outputWhitenMeans, _outputWhitenStdDevs,
    _pdfFileName );
}

MetaRidgeSeed::
~MetaRidgeSeed()
{
  M_Destroy();
}

void MetaRidgeSeed::
PrintInfo() const
{
  MetaLDA::PrintInfo();

  std::cout << "RidgeSeedScales = " << m_RidgeSeedScales.size() << std::endl;

  std::cout << "PDFFileaName = " << m_PDFFileName << std::endl;

  std::cout << "UseIntensityOnly = "
    << ( m_UseIntensityOnly ? "True" : "False" ) << std::endl;

  std::cout << "UseFeatureMath = "
    << ( m_UseFeatureMath ? "True" : "False" ) << std::endl;

  std::cout << "RidgeId = " << m_RidgeId << std::endl;

  std::cout << "BackgroundId = " << m_BackgroundId << std::endl;

  std::cout << "UnknownId = " << m_UnknownId << std::endl;

  std::cout << "SeedTolerance = " << m_SeedTolerance << std::endl;

  std::cout << "Skeletonize = "
    << ( m_Skeletonize ? "True" : "False" ) << std::endl;
}

void MetaRidgeSeed::
CopyInfo( const MetaRidgeSeed & _lda )
{
  MetaLDA::CopyInfo( dynamic_cast< const MetaLDA & >( _lda ) );

  SetRidgeSeedScales( _lda.GetRidgeSeedScales() );
  SetUseIntensityOnly( _lda.GetUseIntensityOnly() );
  SetUseFeatureMath( _lda.GetUseFeatureMath() );
  SetPDFFileName( _lda.GetPDFFileName() );
  SetRidgeId( _lda.GetRidgeId() );
  SetUnknownId( _lda.GetUnknownId() );
  SetBackgroundId( _lda.GetBackgroundId() );
  SetSeedTolerance( _lda.GetSeedTolerance() );
  SetSkeletonize( _lda.GetSkeletonize() );
}

void MetaRidgeSeed::
Clear( void )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: Clear" << std::endl;
    }

  MetaLDA::Clear();

  strcpy( m_FormTypeName, "RidgeSeed" );

  m_RidgeSeedScales.clear();

  m_UseIntensityOnly = false;
  m_UseFeatureMath = false;

  m_PDFFileName.clear();

  m_RidgeId = 255;
  m_BackgroundId = 127;
  m_UnknownId = 0;
  m_SeedTolerance = 1;
  m_Skeletonize = true;

  this->m_NumberOfLDABasisToUseAsFeatures = 1;
  this->m_NumberOfPCABasisToUseAsFeatures = 3;
}

bool MetaRidgeSeed::
InitializeEssential(
  const RidgeSeedScalesType & _ridgeSeedScales,
  bool _useIntensityOnly,
  bool _useFeatureMath,
  const LDAValuesType & _ldaValues,
  const LDAMatrixType & _ldaMatrix,
  const ValueListType & _inputWhitenMeans,
  const ValueListType & _inputWhitenStdDevs,
  const ValueListType & _outputWhitenMeans,
  const ValueListType & _outputWhitenStdDevs,
  const std::string & _pdfFileName )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: Initialize" << std::endl;
    }

  MetaLDA::InitializeEssential( 3, 1,
    _ldaValues, _ldaMatrix, _inputWhitenMeans, _inputWhitenStdDevs,
    _outputWhitenMeans, _outputWhitenStdDevs );

  SetRidgeSeedScales( _ridgeSeedScales );

  SetUseIntensityOnly( _useIntensityOnly );
  SetUseFeatureMath( _useFeatureMath );

  SetPDFFileName( _pdfFileName );

  return true;
}

void MetaRidgeSeed::
SetRidgeSeedScales( const RidgeSeedScalesType & _RidgeSeedScales )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: SetRidgeSeedScales" << std::endl;
    }

  m_RidgeSeedScales = _RidgeSeedScales;
}

const MetaRidgeSeed::RidgeSeedScalesType & MetaRidgeSeed::
GetRidgeSeedScales( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: GetRidgeSeedScales" << std::endl;
    }

  return m_RidgeSeedScales;
}

void MetaRidgeSeed::
SetUseIntensityOnly( bool _UseIntensityOnly )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: SetUseIntensityOnly" << std::endl;
    }

  m_UseIntensityOnly = _UseIntensityOnly;
}

bool MetaRidgeSeed::
GetUseIntensityOnly( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: GetUseIntensityOnly" << std::endl;
    }

  return m_UseIntensityOnly;
}

void MetaRidgeSeed::
SetUseFeatureMath( bool _UseFeatureMath )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: SetUseFeatureMath" << std::endl;
    }

  m_UseFeatureMath = _UseFeatureMath;
}

bool MetaRidgeSeed::
GetUseFeatureMath( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: GetUseFeatureMath" << std::endl;
    }

  return m_UseFeatureMath;
}

void MetaRidgeSeed::
SetPDFFileName( const std::string & _pdfFileName )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: SetPDFFileName" << std::endl;
    }

  m_PDFFileName = _pdfFileName;
}

const std::string & MetaRidgeSeed::
GetPDFFileName( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: GetPDFFileName" << std::endl;
    }

  return m_PDFFileName;
}

void MetaRidgeSeed::
SetRidgeId( int _RidgeId )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: SetRidgeId" << std::endl;
    }

  m_RidgeId = _RidgeId;
}

int MetaRidgeSeed::
GetRidgeId( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: GetRidgeId" << std::endl;
    }

  return m_RidgeId;
}

void MetaRidgeSeed::
SetBackgroundId( int _BackgroundId )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: SetBackgroundId" << std::endl;
    }

  m_BackgroundId = _BackgroundId;
}

int MetaRidgeSeed::
GetBackgroundId( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: GetBackgroundId" << std::endl;
    }

  return m_BackgroundId;
}

void MetaRidgeSeed::
SetUnknownId( int _UnknownId )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: SetUnknownId" << std::endl;
    }

  m_UnknownId = _UnknownId;
}

int MetaRidgeSeed::
GetUnknownId( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: GetUnknownId" << std::endl;
    }

  return m_UnknownId;
}

void MetaRidgeSeed::
SetSkeletonize( bool _Skeletonize )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: SetSkeletonize" << std::endl;
    }

  m_Skeletonize = _Skeletonize;
}

bool MetaRidgeSeed::
GetSkeletonize( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: GetSkeletonize" << std::endl;
    }

  return m_Skeletonize;
}

void MetaRidgeSeed::
SetSeedTolerance( double _SeedTolerance )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: SetSeedTolerance" << std::endl;
    }

  m_SeedTolerance = _SeedTolerance;
}

double MetaRidgeSeed::
GetSeedTolerance( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: GetSeedTolerance" << std::endl;
    }

  return m_SeedTolerance;
}

bool MetaRidgeSeed::
CanRead( const char * _headerName ) const
{
  // First check the extension
  std::string fname = _headerName;
  if( fname == "" )
    {
    return false;
    }

  bool extensionFound = false;

  std::string::size_type stringPos = fname.rfind( ".mrs" );
  if( ( stringPos != std::string::npos )
      && ( stringPos == fname.length() - 5 ) )
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

  const bool result = !std::strncmp( MET_ReadForm( inputStream ).c_str(),
    "RidgeSeed", 9 );

  inputStream.close();

  return result;
}

bool MetaRidgeSeed::
Read( const char * _headerName )
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
    std::cout << "MetaRidgeSeed: Read: Cannot open file _"
                        << m_FileName << "_" << std::endl;
    delete tmpStream;
    return false;
    }

  bool result = ReadStream( tmpStream );

  tmpStream->close();

  delete tmpStream;

  return result;
}

bool MetaRidgeSeed::
CanReadStream( std::ifstream * _stream ) const
{
  if( !std::strncmp( MET_ReadForm( * _stream ).c_str(), "RidgeSeed", 9 ) )
    {
    return true;
    }

  return false;
}

bool MetaRidgeSeed::
ReadStream( std::ifstream * _stream )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: ReadStream" << std::endl;
    }

  M_Destroy();

  Clear();

  M_SetupReadFields();

  if( m_ReadStream )
    {
    std::cout << "MetaRidgeSeed: ReadStream: two files open?" << std::endl;
    delete m_ReadStream;
    }

  m_ReadStream = _stream;

  if( !M_Read() )
    {
    std::cout << "MetaRidgeSeed: Read: Cannot parse file" << std::endl;
    m_ReadStream = NULL;
    return false;
    }

  m_ReadStream = NULL;

  InitializeEssential( m_RidgeSeedScales, m_UseIntensityOnly, m_UseFeatureMath,
    m_LDAValues, m_LDAMatrix, m_InputWhitenMeans,
    m_InputWhitenStdDevs, m_OutputWhitenMeans, m_OutputWhitenStdDevs,
    m_PDFFileName );

  return true;
}

bool MetaRidgeSeed::
Write( const char * _headName )
{
  if( _headName != NULL && std::strlen( _headName ) > 1 )
    {
    FileName( _headName );
    }

  MET_SetFileSuffix( m_FileName, ".mrs" );

  std::ofstream * tmpWriteStream = new std::ofstream;

  tmpWriteStream->open( m_FileName, std::ios::binary |
                                   std::ios::out );

  if( !tmpWriteStream->rdbuf()->is_open() )
    {
    delete tmpWriteStream;
    return false;
    }

  tmpWriteStream->precision( 10 );

  const bool result = WriteStream( tmpWriteStream );

  tmpWriteStream->close();

  delete tmpWriteStream;

  return result;
}

bool MetaRidgeSeed::
WriteStream( std::ofstream * _stream )
{
  if( m_WriteStream != NULL )
    {
    std::cout << "MetaRidgeSeed: WriteStream: two files open?" << std::endl;
    delete m_WriteStream;
    }

  m_WriteStream = _stream;

  M_SetupWriteFields();

  M_Write();

  m_WriteStream->flush();

  m_WriteStream = NULL;

  return true;
}

void MetaRidgeSeed::
M_Destroy( void )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: M_Destroy" << std::endl;
    }

  MetaLDA::M_Destroy();
}

void MetaRidgeSeed::
M_SetupReadFields( void )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: M_SetupReadFields" << std::endl;
    }

  MetaLDA::M_SetupReadFields();

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "NRidgeSeedScales", MET_INT, true );
  m_Fields.push_back( mF );
  int nScalesRecNum = MET_GetFieldRecordNumber( "NRidgeSeedScales",
    &m_Fields );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeSeedScales", MET_FLOAT_ARRAY, true,
    nScalesRecNum );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "UseIntensityOnly", MET_STRING, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "UseFeatureMath", MET_STRING, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "PDFFileName", MET_STRING, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "RidgeId", MET_INT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "BackgroundId", MET_INT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "UnknownId", MET_INT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "SeedTolerance", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "Skeletonize", MET_STRING, true );
  m_Fields.push_back( mF );
}

void MetaRidgeSeed::
M_SetupWriteFields( void )
{
  MetaLDA::M_SetupWriteFields();

  if( !m_RidgeSeedScales.empty() )
    {
    LDAValuesType ridgeSeedScales;
    ridgeSeedScales.set_size( m_RidgeSeedScales.size() );
    for( unsigned int i = 0; i < m_RidgeSeedScales.size(); i++ )
      {
      ridgeSeedScales[i] = m_RidgeSeedScales[i];
      }
    MET_FieldRecordType * mF = new MET_FieldRecordType;
    MET_InitWriteField( mF, "NRidgeSeedScales", MET_INT,
      m_RidgeSeedScales.size() );
    m_Fields.push_back( mF );

    int nRidgeSeedScales = m_RidgeSeedScales.size();
    mF = new MET_FieldRecordType;
    MET_InitWriteField( mF, "RidgeSeedScales", MET_FLOAT_ARRAY,
      nRidgeSeedScales, ridgeSeedScales.data_block() );
    m_Fields.push_back( mF );

    mF = new MET_FieldRecordType;
    if( m_UseIntensityOnly )
      {
      MET_InitWriteField( mF, "UseIntensityOnly", MET_STRING, 4, "True" );
      }
    else
      {
      MET_InitWriteField( mF, "UseIntensityOnly", MET_STRING, 5, "False" );
      }
    m_Fields.push_back( mF );

    mF = new MET_FieldRecordType;
    if( m_UseFeatureMath )
      {
      MET_InitWriteField( mF, "UseFeatureMath", MET_STRING, 4, "True" );
      }
    else
      {
      MET_InitWriteField( mF, "UseFeatureMath", MET_STRING, 5, "False" );
      }
    m_Fields.push_back( mF );

    mF = new MET_FieldRecordType;
    MET_InitWriteField( mF, "PDFFileName", MET_STRING,
      m_PDFFileName.size(), m_PDFFileName.c_str() );
    m_Fields.push_back( mF );

    mF = new MET_FieldRecordType;
    MET_InitWriteField( mF, "RidgeId", MET_INT, m_RidgeId );
    m_Fields.push_back( mF );

    mF = new MET_FieldRecordType;
    MET_InitWriteField( mF, "BackgroundId", MET_INT, m_BackgroundId );
    m_Fields.push_back( mF );

    mF = new MET_FieldRecordType;
    MET_InitWriteField( mF, "UnknownId", MET_INT, m_UnknownId );
    m_Fields.push_back( mF );

    mF = new MET_FieldRecordType;
    MET_InitWriteField( mF, "SeedTolerance", MET_FLOAT, m_SeedTolerance );
    m_Fields.push_back( mF );

    mF = new MET_FieldRecordType;
    if( m_Skeletonize )
      {
      MET_InitWriteField( mF, "Skeletonize", MET_STRING, 4, "True" );
      }
    else
      {
      MET_InitWriteField( mF, "Skeletonize", MET_STRING, 5, "False" );
      }
    m_Fields.push_back( mF );
    }

}

bool MetaRidgeSeed::
M_Read( void )
{
  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: M_Read: Loading Header" << std::endl;
    }
  if( !MetaLDA::M_Read() )
    {
    std::cout << "MetaRidgeSeed: M_Read: Error parsing file" << std::endl;
    return false;
    }

  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: M_Read: Parsing Header" << std::endl;
    }

  if( META_DEBUG )
    {
    std::cout << "MetaRidgeSeed: M_Read: num fields = "
      << m_Fields.size() << std::endl;
    for( unsigned int i = 0; i < m_Fields.size(); i++ )
      {
      std::cout << "  Field " << i << " = " << m_Fields[i]->name << std::endl;
      }
    }

  MET_FieldRecordType * mF = MET_GetFieldRecord( "NRidgeSeedScales",
    &m_Fields );
  unsigned int nRidgeSeedScales = ( unsigned int )mF->value[0];
  m_RidgeSeedScales.resize( nRidgeSeedScales, 0 );
  mF = MET_GetFieldRecord( "RidgeSeedScales", &m_Fields );
  if( mF && mF->defined )
    {
    for( unsigned int i = 0; i < nRidgeSeedScales; i++ )
      {
      m_RidgeSeedScales[i] = ( double )mF->value[i];
      }
    }

  mF = MET_GetFieldRecord( "UseIntensityOnly", &m_Fields );
  if( ( ( char * )( mF->value ) )[0] == 'T'
    || ( ( char * )( mF->value ) ) [0] == 't' )
    {
    m_UseIntensityOnly = true;
    }
  else
    {
    m_UseIntensityOnly = false;
    }

  mF = MET_GetFieldRecord( "UseFeatureMath", &m_Fields );
  if( ( ( char * )( mF->value ) )[0] == 'T'
    || ( ( char * )( mF->value ) ) [0] == 't' )
    {
    m_UseFeatureMath = true;
    }
  else
    {
    m_UseFeatureMath = false;
    }

  mF = MET_GetFieldRecord( "PDFFileName", &m_Fields );
  m_PDFFileName = ( char * )( &( mF->value[0] ) );

  mF = MET_GetFieldRecord( "RidgeId", &m_Fields );
  m_RidgeId = ( int )( mF->value[0] );

  mF = MET_GetFieldRecord( "BackgroundId", &m_Fields );
  m_BackgroundId = ( int )( mF->value[0] );

  mF = MET_GetFieldRecord( "UnknownId", &m_Fields );
  m_UnknownId = ( int )( mF->value[0] );

  mF = MET_GetFieldRecord( "SeedTolerance", &m_Fields );
  m_SeedTolerance = ( double )( mF->value[0] );

  mF = MET_GetFieldRecord( "Skeletonize", &m_Fields );
  if( ( ( char * )( mF->value ) )[0] == 'T'
    || ( ( char * )( mF->value ) ) [0] == 't' )
    {
    m_Skeletonize = true;
    }
  else
    {
    m_Skeletonize = false;
    }

  return true;
}

} // End namespace tube

} // End namespace itk
