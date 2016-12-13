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
    METAIO_STREAM::cout << "MetaRidgeSeed()" << METAIO_STREAM::endl;
    }

  Clear();
}

MetaRidgeSeed::
MetaRidgeSeed( const char * _headerName )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed()" << METAIO_STREAM::endl;
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
   METAIO_STREAM::cout << "MetaRidgeSeed()" << METAIO_STREAM::endl;
   }

  Clear();

  CopyInfo( _metaRidgeSeed );
}

MetaRidgeSeed::
MetaRidgeSeed( 
  const RidgeSeedScalesType & _ridgeSeedScales,
  bool _useIntensityOnly,
  bool _useSVM,
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
   METAIO_STREAM::cout << "MetaRidgeSeed()" << METAIO_STREAM::endl;
   }

  Clear();

  InitializeEssential( _ridgeSeedScales, _useIntensityOnly, _useSVM,
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

  METAIO_STREAM::cout << "RidgeSeedScales = " << m_RidgeSeedScales.size()
    << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "PDFFileaName = "
    << m_PDFFileName << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "UseSVM = "
    << ( m_UseSVM ? "True" : "False" ) << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "UseIntensityOnly = "
    << ( m_UseIntensityOnly ? "True" : "False" ) << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "RidgeId = " << m_RidgeId
    << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "BackgroundId = " << m_BackgroundId
    << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "UnknownId = " << m_UnknownId
    << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "SeedTolerance = " << m_SeedTolerance
    << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "Skeletonize = "
    << ( m_Skeletonize ? "True" : "False" ) << METAIO_STREAM::endl;
}

void MetaRidgeSeed::
CopyInfo( const MetaRidgeSeed & _lda )
{
  MetaLDA::CopyInfo( dynamic_cast< const MetaLDA & >( _lda ) );

  SetRidgeSeedScales( _lda.GetRidgeSeedScales() );
  SetUseIntensityOnly( _lda.GetUseIntensityOnly() );
  SetUseSVM( _lda.GetUseSVM() );
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
    METAIO_STREAM::cout << "MetaRidgeSeed: Clear" << METAIO_STREAM::endl;
    }

  MetaLDA::Clear();

  strcpy( m_FormTypeName, "RidgeSeed" );

  m_RidgeSeedScales.clear();

  m_UseIntensityOnly = false;

  m_UseSVM = false;

  m_PDFFileName.clear();

  m_RidgeId = 255;
  m_BackgroundId = 127;
  m_UnknownId = 0;
  m_SeedTolerance = 1;
  m_Skeletonize = true;
}

bool MetaRidgeSeed::
InitializeEssential( 
  const RidgeSeedScalesType & _ridgeSeedScales,
  bool _useIntensityOnly,
  bool _useSVM,
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
    METAIO_STREAM::cout << "MetaRidgeSeed: Initialize"
      << METAIO_STREAM::endl;
    }

  MetaLDA::InitializeEssential( 3, 1,
    _ldaValues, _ldaMatrix, _inputWhitenMeans, _inputWhitenStdDevs,
    _outputWhitenMeans, _outputWhitenStdDevs );

  SetRidgeSeedScales( _ridgeSeedScales );

  SetUseIntensityOnly( _useIntensityOnly );

  SetUseSVM( _useSVM );

  SetPDFFileName( _pdfFileName );

  return true;
}

void MetaRidgeSeed::
SetRidgeSeedScales( const RidgeSeedScalesType & _RidgeSeedScales )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: SetRidgeSeedScales"
      << METAIO_STREAM::endl;
    }

  m_RidgeSeedScales = _RidgeSeedScales;
}

const MetaRidgeSeed::RidgeSeedScalesType & MetaRidgeSeed::
GetRidgeSeedScales( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: GetRidgeSeedScales"
      << METAIO_STREAM::endl;
    }

  return m_RidgeSeedScales;
}

void MetaRidgeSeed::
SetUseIntensityOnly( bool _UseIntensityOnly )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: SetUseIntensityOnly"
      << METAIO_STREAM::endl;
    }

  m_UseIntensityOnly = _UseIntensityOnly;
}

bool MetaRidgeSeed::
GetUseIntensityOnly( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: GetUseIntensityOnly"
      << METAIO_STREAM::endl;
    }

  return m_UseIntensityOnly;
}

void MetaRidgeSeed::
SetUseSVM( bool _UseSVM )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: SetUseSVM"
      << METAIO_STREAM::endl;
    }

  m_UseSVM = _UseSVM;
}

bool MetaRidgeSeed::
GetUseSVM( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: GetUseSVM"
      << METAIO_STREAM::endl;
    }

  return m_UseSVM;
}


void MetaRidgeSeed::
SetPDFFileName( const std::string & _pdfFileName )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: SetPDFFileName"
      << METAIO_STREAM::endl;
    }

  m_PDFFileName = _pdfFileName;
}

const std::string & MetaRidgeSeed::
GetPDFFileName( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: GetPDFFileName"
      << METAIO_STREAM::endl;
    }

  return m_PDFFileName;
}

void MetaRidgeSeed::
SetRidgeId( int _RidgeId )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: SetRidgeId"
      << METAIO_STREAM::endl;
    }

  m_RidgeId = _RidgeId;
}

int MetaRidgeSeed::
GetRidgeId( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: GetRidgeId"
      << METAIO_STREAM::endl;
    }

  return m_RidgeId;
}

void MetaRidgeSeed::
SetBackgroundId( int _BackgroundId )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: SetBackgroundId"
      << METAIO_STREAM::endl;
    }

  m_BackgroundId = _BackgroundId;
}

int MetaRidgeSeed::
GetBackgroundId( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: GetBackgroundId"
      << METAIO_STREAM::endl;
    }

  return m_BackgroundId;
}

void MetaRidgeSeed::
SetUnknownId( int _UnknownId )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: SetUnknownId"
      << METAIO_STREAM::endl;
    }

  m_UnknownId = _UnknownId;
}

int MetaRidgeSeed::
GetUnknownId( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: GetUnknownId"
      << METAIO_STREAM::endl;
    }

  return m_UnknownId;
}

void MetaRidgeSeed::
SetSkeletonize( bool _Skeletonize )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: SetSkeletonize"
      << METAIO_STREAM::endl;
    }

  m_Skeletonize = _Skeletonize;
}

bool MetaRidgeSeed::
GetSkeletonize( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: GetSkeletonize"
      << METAIO_STREAM::endl;
    }

  return m_Skeletonize;
}

void MetaRidgeSeed::
SetSeedTolerance( double _SeedTolerance )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: SetSeedTolerance"
      << METAIO_STREAM::endl;
    }

  m_SeedTolerance = _SeedTolerance;
}

double MetaRidgeSeed::
GetSeedTolerance( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: GetSeedTolerance"
      << METAIO_STREAM::endl;
    }

  return m_SeedTolerance;
}

bool MetaRidgeSeed::
CanRead( const char * _headerName ) const
{
  // First check the extension
  METAIO_STL::string fname = _headerName;
  if( fname == "" )
    {
    return false;
    }

  bool extensionFound = false;

  METAIO_STL::string::size_type stringPos = fname.rfind( ".mrs" );
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

  METAIO_STREAM::ifstream * tmpStream = new METAIO_STREAM::ifstream;

  tmpStream->open( m_FileName, METAIO_STREAM::ios::in |
                               METAIO_STREAM::ios::binary );

  if( !tmpStream->rdbuf()->is_open() )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: Read: Cannot open file _"
                        << m_FileName << "_" << METAIO_STREAM::endl;
    delete tmpStream;
    return false;
    }

  bool result = ReadStream( tmpStream );

  tmpStream->close();

  delete tmpStream;

  return result;
}

bool MetaRidgeSeed::
CanReadStream( METAIO_STREAM::ifstream * _stream ) const
{
  if( !std::strncmp( MET_ReadForm( * _stream ).c_str(), "RidgeSeed", 9 ) )
    {
    return true;
    }

  return false;
}

bool MetaRidgeSeed::
ReadStream( METAIO_STREAM::ifstream * _stream )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: ReadStream"
      << METAIO_STREAM::endl;
    }

  M_Destroy();

  Clear();

  M_SetupReadFields();

  if( m_ReadStream )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: ReadStream: two files open?"
                        << METAIO_STREAM::endl;
    delete m_ReadStream;
    }

  m_ReadStream = _stream;

  if( !M_Read() )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: Read: Cannot parse file"
                        << METAIO_STREAM::endl;
    m_ReadStream = NULL;
    return false;
    }

  m_ReadStream = NULL;

  InitializeEssential( m_RidgeSeedScales, m_UseIntensityOnly,
    m_UseSVM, m_LDAValues, m_LDAMatrix, m_InputWhitenMeans,
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

  METAIO_STREAM::ofstream * tmpWriteStream = new METAIO_STREAM::ofstream;

  tmpWriteStream->open( m_FileName, METAIO_STREAM::ios::binary |
                                   METAIO_STREAM::ios::out );

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
WriteStream( METAIO_STREAM::ofstream * _stream )
{
  if( m_WriteStream != NULL )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: WriteStream: two files open?"
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

void MetaRidgeSeed::
M_Destroy( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: M_Destroy"
                        << METAIO_STREAM::endl;
    }

  MetaLDA::M_Destroy();
}

void MetaRidgeSeed::
M_SetupReadFields( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: M_SetupReadFields"
                        << METAIO_STREAM::endl;
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
  MET_InitReadField( mF, "UseSVM", MET_STRING, true );
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
    if( m_UseSVM )
      {
      MET_InitWriteField( mF, "UseSVM", MET_STRING, 4, "True" );
      }
    else
      {
      MET_InitWriteField( mF, "UseSVM", MET_STRING, 5, "False" );
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
    METAIO_STREAM::cout << "MetaRidgeSeed: M_Read: Loading Header"
                        << METAIO_STREAM::endl;
    }
  if( !MetaLDA::M_Read() )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: M_Read: Error parsing file"
                        << METAIO_STREAM::endl;
    return false;
    }

  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: M_Read: Parsing Header"
                        << METAIO_STREAM::endl;
    }

  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaRidgeSeed: M_Read: num fields = "
      << m_Fields.size() << METAIO_STREAM::endl;
    for( unsigned int i = 0; i < m_Fields.size(); i++ )
      {
      METAIO_STREAM::cout << "  Field " << i << " = "
        << m_Fields[i]->name << METAIO_STREAM::endl;
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

  mF = MET_GetFieldRecord( "UseSVM", &m_Fields );
  if( ( ( char * )( mF->value ) )[0] == 'T'
    || ( ( char * )( mF->value ) ) [0] == 't' )
    {
    m_UseSVM = true;
    }
  else
    {
    m_UseSVM = false;
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
