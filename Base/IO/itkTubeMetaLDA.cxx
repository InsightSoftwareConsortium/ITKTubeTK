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
#include "itkTubeMetaLDA.h"

#include <stdio.h>
#include <ctype.h>
#include <string>
#include <string.h> // for memcpy
#include <math.h>

namespace itk {

namespace tube {

//
// MetaLDA Constructors
//
MetaLDA::
MetaLDA( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA()" << METAIO_STREAM::endl;
    }

  Clear();
}

//
MetaLDA::
MetaLDA( const char *_headerName )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA()" << METAIO_STREAM::endl;
    }

  Clear();
  m_ReadStream = NULL;

  MetaLDA::Read( _headerName );
}

//
MetaLDA::
MetaLDA( const MetaLDA & _metaLDA )
  : MetaForm()
{
  if( META_DEBUG )
   {
   METAIO_STREAM::cout << "MetaLDA()" << METAIO_STREAM::endl;
   }

  Clear();

  CopyInfo( _metaLDA );
}

//
MetaLDA::
MetaLDA( const LDAValuesType & _ldaValues,
  const LDAMatrixType & _ldaMatrix,
  const ValueListType & _whitenMeans,
  const ValueListType & _whitenStdDevs )
{
  if( META_DEBUG )
   {
   METAIO_STREAM::cout << "MetaLDA()" << METAIO_STREAM::endl;
   }

  Clear();

  InitializeEssential( _ldaValues, _ldaMatrix, _whitenMeans, _whitenStdDevs );
}

//
MetaLDA::
~MetaLDA( void )
{
  M_Destroy();
}

//
void MetaLDA::
PrintInfo( void ) const
{
  MetaForm::PrintInfo();

  METAIO_STREAM::cout << "LDAValues = " << m_LDAValues
    << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "LDAMatrix = " << m_LDAMatrix
    << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "WhitenMeans = " << METAIO_STREAM::endl;
  for( unsigned int i=0; i<m_WhitenMeans.size(); i++ )
    {
    METAIO_STREAM::cout << m_WhitenMeans[i] << " ";
    }
  METAIO_STREAM::cout << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "WhitenStdDevs = " << METAIO_STREAM::endl;
  for( unsigned int i=0; i<m_WhitenStdDevs.size(); i++ )
    {
    METAIO_STREAM::cout << m_WhitenStdDevs[i] << " ";
    }
  METAIO_STREAM::cout << METAIO_STREAM::endl;
}

void MetaLDA::
CopyInfo( const MetaLDA & _lda )
{
  MetaForm::CopyInfo( dynamic_cast< const MetaForm * >( &_lda ) );

  SetLDAValues( _lda.GetLDAValues() );
  SetLDAMatrix( _lda.GetLDAMatrix() );
  SetWhitenMeans( _lda.GetWhitenMeans() );
  SetWhitenStdDevs( _lda.GetWhitenStdDevs() );
}

void MetaLDA::
Clear( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: Clear" << METAIO_STREAM::endl;
    }

  m_LDAValues.set_size( 0 );
  m_LDAMatrix.set_size( 0, 0 );
  m_WhitenMeans.clear();
  m_WhitenStdDevs.clear();

  MetaForm::Clear();
}

bool MetaLDA::
InitializeEssential( const LDAValuesType & _ldaValues,
  const LDAMatrixType & _ldaMatrix,
  const ValueListType & _whitenMeans,
  const ValueListType & _whitenStdDevs )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: Initialize" << METAIO_STREAM::endl;
    }

  MetaForm::InitializeEssential();

  SetLDAValues( _ldaValues );
  SetLDAMatrix( _ldaMatrix );
  SetWhitenMeans( _whitenMeans );
  SetWhitenStdDevs( _whitenStdDevs );

  return true;
}

//
//
//
void MetaLDA::
SetLDAValues( const LDAValuesType & _ldaValues )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: SetLDAValues" << METAIO_STREAM::endl;
    }

  m_LDAValues = _ldaValues;
}

const MetaLDA::LDAValuesType & MetaLDA::
GetLDAValues( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: GetLDAValues" << METAIO_STREAM::endl;
    }

  return m_LDAValues;
}

void MetaLDA::
SetLDAMatrix( const LDAMatrixType & _ldaMatrix )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: SetLDAMatrix" << METAIO_STREAM::endl;
    }

  m_LDAMatrix = _ldaMatrix;
}

const MetaLDA::LDAMatrixType & MetaLDA::
GetLDAMatrix( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: GetLDAMatrix" << METAIO_STREAM::endl;
    }

  return m_LDAMatrix;
}

void MetaLDA::
SetWhitenMeans( const ValueListType & _whitenMeans )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: SetWhitenMeans" << METAIO_STREAM::endl;
    }

  m_WhitenMeans = _whitenMeans;
}

const MetaLDA::ValueListType & MetaLDA::
GetWhitenMeans( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: GetWhitenMeans" << METAIO_STREAM::endl;
    }

  return m_WhitenMeans;
}

void MetaLDA::
SetWhitenStdDevs( const ValueListType & _whitenStdDevs )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: SetWhitenStdDevs" << METAIO_STREAM::endl;
    }

  m_WhitenStdDevs = _whitenStdDevs;
}

const MetaLDA::ValueListType & MetaLDA::
GetWhitenStdDevs( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: GetWhitenStdDevs" << METAIO_STREAM::endl;
    }

  return m_WhitenStdDevs;
}

bool MetaLDA::CanRead( const char *_headerName ) const
{
  // First check the extension
  METAIO_STL::string fname = _headerName;
  if( fname == "" )
    {
    return false;
    }

  bool extensionFound = false;

  METAIO_STL::string::size_type stringPos = fname.rfind( ".mlda" );
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

  bool result = !strncmp( MET_ReadForm( inputStream ).c_str(), "LDA", 3 );

  inputStream.close();

  return result;
}


bool MetaLDA::
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
    METAIO_STREAM::cout << "MetaLDA: Read: Cannot open file _"
                        << m_FileName << "_" << METAIO_STREAM::endl;
    delete tmpStream;
    return false;
    }

  bool result = ReadStream( tmpStream );

  tmpStream->close();

  delete tmpStream;

  return result;
}

bool MetaLDA::CanReadStream( METAIO_STREAM::ifstream * _stream ) const
{
  if( !strncmp( MET_ReadForm( *_stream ).c_str(), "LDA", 3 ) )
    {
    return true;
    }

  return false;
}

bool MetaLDA::
ReadStream( METAIO_STREAM::ifstream * _stream )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: ReadStream" << METAIO_STREAM::endl;
    }

  M_Destroy();

  Clear();

  M_SetupReadFields();

  if( m_ReadStream )
    {
    METAIO_STREAM::cout << "MetaLDA: ReadStream: two files open?"
                        << METAIO_STREAM::endl;
    delete m_ReadStream;
    }

  m_ReadStream = _stream;

  if( !M_Read() )
    {
    METAIO_STREAM::cout << "MetaLDA: Read: Cannot parse file"
                        << METAIO_STREAM::endl;
    m_ReadStream = NULL;
    return false;
    }

  m_ReadStream = NULL;

  InitializeEssential( m_LDAValues, m_LDAMatrix,
    m_WhitenMeans, m_WhitenStdDevs );

  return true;
}

bool MetaLDA::Write( const char *_headName )
{
  if( _headName != NULL && strlen( _headName ) > 1 )
    {
    FileName( _headName );
    }

  MET_SetFileSuffix( m_FileName, "mlda" );

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

bool MetaLDA::
WriteStream( METAIO_STREAM::ofstream * _stream )
{
  if( m_WriteStream != NULL )
    {
    METAIO_STREAM::cout << "MetaLDA: WriteStream: two files open?"
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
        
void MetaLDA::
M_Destroy( void )
{
  MetaForm::M_Destroy();
}

void MetaLDA::
M_SetupReadFields( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: M_SetupReadFields"
                        << METAIO_STREAM::endl;
    }

  MetaForm::M_SetupReadFields();

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "NDims", MET_INT, true );
  m_Fields.push_back( mF );

  int nDimsRecNum = MET_GetFieldRecordNumber( "NDims", &m_Fields );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "Values", MET_FLOAT_ARRAY, true, nDimsRecNum );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "Matrix", MET_FLOAT_MATRIX, true, nDimsRecNum );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "WhitenMeans", MET_FLOAT_ARRAY, false,
    nDimsRecNum );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "WhitenStdDevs", MET_FLOAT_ARRAY, false,
    nDimsRecNum );
  m_Fields.push_back( mF );
}

void MetaLDA::
M_SetupWriteFields( void )
{
  strcpy( m_FormTypeName, "LDA" );
  MetaForm::M_SetupWriteFields();

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "NDims", MET_INT, m_LDAValues.size() );
  m_Fields.push_back( mF );

  int nDims = m_LDAValues.size();

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "Values", MET_FLOAT_ARRAY, nDims,
    m_LDAValues.data_block() );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "Matrix", MET_FLOAT_MATRIX, nDims,
    m_LDAMatrix.data_block() );
  m_Fields.push_back( mF );

  int tfCount = m_WhitenMeans.size();

  if( tfCount > 0 )
    {
    double tf[ 4096 ];
    mF = new MET_FieldRecordType;
    for( int i=0; i<tfCount; i++ )
      {
      tf[i] = m_WhitenMeans[i];
      }
    MET_InitWriteField( mF, "WhitenMeans", MET_FLOAT_ARRAY, tfCount, tf );
    m_Fields.push_back( mF );

    mF = new MET_FieldRecordType;
    for( int i=0; i<tfCount; i++ )
      {
      tf[i] = m_WhitenStdDevs[i];
      }
    MET_InitWriteField( mF, "WhitenStdDevs", MET_FLOAT_ARRAY, tfCount, tf);
    m_Fields.push_back( mF );
    }
}


bool MetaLDA::
M_Read( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: M_Read: Loading Header"
                        << METAIO_STREAM::endl;
    }
  if( !MetaForm::M_Read() )
    {
    METAIO_STREAM::cout << "MetaLDA: M_Read: Error parsing file"
                        << METAIO_STREAM::endl;
    return false;
    }

  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: M_Read: Parsing Header"
                        << METAIO_STREAM::endl;
    }
  MET_FieldRecordType * mF;
     
  unsigned int nDims = 0;
  mF = MET_GetFieldRecord( "NDims", &m_Fields );
  if( mF && mF->defined )
    {
    nDims = ( unsigned int )mF->value[0];
    m_LDAValues.set_size( nDims );
    m_LDAValues.fill( 0 );
    m_LDAMatrix.set_size( nDims, nDims );
    m_LDAMatrix.fill( 0 );
    m_WhitenMeans.resize( nDims );
    m_WhitenStdDevs.resize( nDims );
    }
  else
    {
    METAIO_STREAM::cout << "MetaLDA: M_Read: Error: NDims required"
      << METAIO_STREAM::endl;
    return false;
    }

  mF = MET_GetFieldRecord( "Values", &m_Fields );
  if( mF && mF->defined )
    {
    for( unsigned int i=0; i<nDims; i++ )
      {
      m_LDAValues[i] = ( double )mF->value[i];
      }
    }
  else
    {
    METAIO_STREAM::cout << "MetaLDA: M_Read: Error: Values required"
      << METAIO_STREAM::endl;
    return false;
    }

  mF = MET_GetFieldRecord( "Matrix", &m_Fields );
  if( mF && mF->defined )
    {
    for( unsigned int i=0; i<nDims; i++ )
      {
      for( unsigned int j=0; j<nDims; j++ )
        {
        m_LDAMatrix[i][j] = ( double )mF->value[i*nDims + j];
        }
      }
    }
  else
    {
    METAIO_STREAM::cout << "MetaLDA: M_Read: Error: Matrix required"
      << METAIO_STREAM::endl;
    return false;
    }

  mF = MET_GetFieldRecord( "WhitenMeans", &m_Fields );
  if( mF && mF->defined )
    {
    for( unsigned int i=0; i<nDims; i++ )
      {
      m_WhitenMeans[i] = ( double )mF->value[i];
      }
    }
  else
    {
    m_WhitenMeans.clear();
    }

  mF = MET_GetFieldRecord( "WhitenStdDevs", &m_Fields );
  if( mF && mF->defined )
    {
    for( unsigned int i=0; i<nDims; i++ )
      {
      m_WhitenStdDevs[i] = ( double )mF->value[i];
      }
    }
  else
    {
    m_WhitenStdDevs.clear();
    }

  return true;
}

}

}
