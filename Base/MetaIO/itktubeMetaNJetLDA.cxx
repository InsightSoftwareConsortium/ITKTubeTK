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

#include "itktubeMetaNJetLDA.h"

namespace itk
{

namespace tube
{

MetaNJetLDA
::MetaNJetLDA( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA()" << METAIO_STREAM::endl;
    }

  this->Clear();
}

MetaNJetLDA
::MetaNJetLDA( const char * _headerName )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA()" << METAIO_STREAM::endl;
    }

  this->Clear();

  MetaNJetLDA::Read( _headerName );
}

MetaNJetLDA
::MetaNJetLDA( const MetaNJetLDA & metaNJetLDA ) : MetaLDA()
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA()" << METAIO_STREAM::endl;
    }

  this->Clear();
  this->CopyInfo( metaNJetLDA );
}

MetaNJetLDA
::MetaNJetLDA( const NJetScalesType & zeroScales,
               const NJetScalesType & firstScales,
               const NJetScalesType & secondScales,
               const NJetScalesType & ridgeScales,
               unsigned int numberOfPCABasis,
               unsigned int numberOfLDABasis,
               const LDAValuesType & ldaValues,
               const LDAMatrixType & ldaMatrix,
               const ValueListType & inputWhitenMeans,
               const ValueListType & inputWhitenStdDevs,
               const ValueListType & outputWhitenMeans,
               const ValueListType & outputWhitenStdDevs )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA()" << METAIO_STREAM::endl;
    }

  this->Clear();
  this->InitializeEssential( zeroScales, firstScales, secondScales,
    ridgeScales, numberOfPCABasis, numberOfLDABasis, ldaValues, ldaMatrix,
    inputWhitenMeans, inputWhitenStdDevs, outputWhitenMeans,
    outputWhitenStdDevs );
}

MetaNJetLDA
::~MetaNJetLDA( void )
{
  this->M_Destroy();
}

void MetaNJetLDA
::PrintInfo( void ) const
{
  MetaLDA::PrintInfo();

  METAIO_STREAM::cout << "ZeroScales = " << m_ZeroScales.size()
                      << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "FirstScales = " << m_FirstScales.size()
                      << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "SecondScales = " << m_SecondScales.size()
                      << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "RidgeScales = " << m_RidgeScales.size()
                      << METAIO_STREAM::endl;

  METAIO_STREAM::cout << "ZeroScalesTmp = " << m_ZeroScalesTmp.size()
                      << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "FirstScalesTmp = " << m_FirstScalesTmp.size()
                      << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "SecondScalesTmp = " << m_SecondScalesTmp.size()
                      << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "RidgeScalesTmp = " << m_RidgeScalesTmp.size()
                      << METAIO_STREAM::endl;
}

void MetaNJetLDA
::CopyInfo( const MetaNJetLDA & lda )
{
  MetaLDA::CopyInfo( dynamic_cast< const MetaLDA & >( lda ) );

  this->SetZeroScales( lda.GetZeroScales() );
  this->SetFirstScales( lda.GetFirstScales() );
  this->SetSecondScales( lda.GetSecondScales() );
  this->SetRidgeScales( lda.GetRidgeScales() );
}

void MetaNJetLDA
::Clear( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: Clear" << METAIO_STREAM::endl;
    }

  MetaLDA::Clear();

  strcpy( m_FormTypeName, "NJetLDA" );

  m_ZeroScales.clear();
  m_FirstScales.clear();
  m_SecondScales.clear();
  m_RidgeScales.clear();

}

bool MetaNJetLDA
::InitializeEssential( const NJetScalesType & zeroScales,
                       const NJetScalesType & firstScales,
                       const NJetScalesType & secondScales,
                       const NJetScalesType & ridgeScales,
                       unsigned int numberOfPCABasis,
                       unsigned int numberOfLDABasis,
                       const LDAValuesType & ldaValues,
                       const LDAMatrixType & ldaMatrix,
                       const ValueListType & inputWhitenMeans,
                       const ValueListType & inputWhitenStdDevs,
                       const ValueListType & outputWhitenMeans,
                       const ValueListType & outputWhitenStdDevs )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: Initialize" << METAIO_STREAM::endl;
    }

  MetaLDA::InitializeEssential( numberOfPCABasis, numberOfLDABasis,
    ldaValues, ldaMatrix, inputWhitenMeans, inputWhitenStdDevs,
    outputWhitenMeans, outputWhitenStdDevs );

  this->SetZeroScales( zeroScales );
  this->SetFirstScales( firstScales );
  this->SetSecondScales( secondScales );
  this->SetRidgeScales( ridgeScales );

  return true;
}

void MetaNJetLDA
::SetZeroScales( const NJetScalesType & zeroScales )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: SetZeroScales" << METAIO_STREAM::endl;
    METAIO_STREAM::cout << " Size = " << zeroScales.size()
                        << METAIO_STREAM::endl;
    for( unsigned int i = 0; i < zeroScales.size(); i++ )
      {
      METAIO_STREAM::cout << " Scale " << i << " = " << zeroScales[i]
                          << METAIO_STREAM::endl;
      }
    }

  m_ZeroScales = zeroScales;
}

const MetaNJetLDA::NJetScalesType & MetaNJetLDA
::GetZeroScales( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: GetZeroScales" << METAIO_STREAM::endl;
    }

  return m_ZeroScales;
}

void MetaNJetLDA
::SetFirstScales( const NJetScalesType & firstScales )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: SetFirstScales" << METAIO_STREAM::endl;
    }

  m_FirstScales = firstScales;
}

const MetaNJetLDA::NJetScalesType & MetaNJetLDA
::GetFirstScales( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: GetFirstScales" << METAIO_STREAM::endl;
    }

  return m_FirstScales;
}

void MetaNJetLDA
::SetSecondScales( const NJetScalesType & secondScales )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: SetSecondScales"
                        << METAIO_STREAM::endl;
    }

  m_SecondScales = secondScales;
}

const MetaNJetLDA::NJetScalesType & MetaNJetLDA
::GetSecondScales( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: GetSecondScales"
                        << METAIO_STREAM::endl;
    }

  return m_SecondScales;
}

void MetaNJetLDA
::SetRidgeScales( const NJetScalesType & ridgeScales )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: SetRidgeScales" << METAIO_STREAM::endl;
    }

  m_RidgeScales = ridgeScales;
}

const MetaNJetLDA::NJetScalesType & MetaNJetLDA
::GetRidgeScales( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: GetRidgeScales" << METAIO_STREAM::endl;
    }

  return m_RidgeScales;
}

bool MetaNJetLDA
::CanRead( const char * headerName ) const
{
  // First check the extension.
  METAIO_STL::string fname = headerName;
  if( fname == "" )
    {
    return false;
    }

  bool extensionFound = false;

  METAIO_STL::string::size_type stringPos = fname.rfind( ".mnda" );
  if( ( stringPos != METAIO_STL::string::npos )
      && ( stringPos == fname.length() - 5 ) )
    {
    extensionFound = true;
    }

  if( !extensionFound )
    {
    return false;
    }

  // Now check the file content.
  METAIO_STREAM::ifstream inputStream;

  inputStream.open( headerName,
                    METAIO_STREAM::ios::in | METAIO_STREAM::ios::binary );

  if( !inputStream.rdbuf()->is_open() )
    {
    return false;
    }

  const bool result = !std::strncmp( MET_ReadForm(
                                       inputStream ).c_str(), "NJetLDA", 7 );

  inputStream.close();

  return result;
}

bool MetaNJetLDA
::Read( const char * headerName )
{
  if( headerName != NULL && std::strlen( headerName ) > 1 )
    {
    this->FileName( headerName );
    }

  METAIO_STREAM::ifstream * const tmpStream = new METAIO_STREAM::ifstream();

  tmpStream->open( m_FileName,
                   METAIO_STREAM::ios::in | METAIO_STREAM::ios::binary );

  if( !tmpStream->rdbuf()->is_open() )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: Read: Cannot open file _"
                        << m_FileName << "_" << METAIO_STREAM::endl;
    delete tmpStream;
    return false;
    }

  const bool result = ReadStream( tmpStream );

  tmpStream->close();

  delete tmpStream;

  return result;
}

bool MetaNJetLDA
::CanReadStream( METAIO_STREAM::ifstream * stream ) const
{
  if( !std::strncmp( MET_ReadForm( *stream ).c_str(), "NJetLDA", 7 ) )
    {
    return true;
    }

  return false;
}

bool MetaNJetLDA
::ReadStream( METAIO_STREAM::ifstream * stream )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: ReadStream" << METAIO_STREAM::endl;
    }

  this->M_Destroy();
  this->Clear();
  this->M_SetupReadFields();

  if( m_ReadStream )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: ReadStream: Are two files open?"
                        << METAIO_STREAM::endl;
    delete m_ReadStream;
    }

  m_ReadStream = stream;

  if( !this->M_Read() )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: Read: Cannot parse file."
                        << METAIO_STREAM::endl;
    m_ReadStream = NULL;
    return false;
    }

  m_ReadStream = NULL;

  this->InitializeEssential( m_ZeroScales, m_FirstScales, m_SecondScales,
    m_RidgeScales, m_NumberOfPCABasisToUseAsFeatures,
    m_NumberOfLDABasisToUseAsFeatures, m_LDAValues, m_LDAMatrix,
    m_InputWhitenMeans, m_InputWhitenStdDevs, m_OutputWhitenMeans,
    m_OutputWhitenStdDevs );

  return true;
}

bool MetaNJetLDA
::Write( const char * headerName )
{
  if( headerName != NULL && std::strlen( headerName ) > 1 )
    {
    this->FileName( headerName );
    }

  MET_SetFileSuffix( m_FileName, "mnda" );

  METAIO_STREAM::ofstream * const tmpWriteStream = new
    METAIO_STREAM::ofstream();

  tmpWriteStream->open( m_FileName,
    METAIO_STREAM::ios::binary | METAIO_STREAM::ios::out );

  if( !tmpWriteStream->rdbuf()->is_open() )
    {
    delete tmpWriteStream;
    return false;
    }

  tmpWriteStream->precision( 10 );

  const bool result = this->WriteStream( tmpWriteStream );

  tmpWriteStream->close();
  delete tmpWriteStream;

  return result;
}

bool MetaNJetLDA
::WriteStream( METAIO_STREAM::ofstream * stream )
{
  if( m_WriteStream != NULL )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: WriteStream: Are two files open?"
                        << METAIO_STREAM::endl;
    delete m_WriteStream;
    }

  m_WriteStream = stream;

  this->M_SetupWriteFields();
  this->M_Write();

  m_WriteStream->flush();
  m_WriteStream = NULL;

  return true;
}

void MetaNJetLDA
::M_Destroy( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: M_Destroy" << METAIO_STREAM::endl;
    }

  MetaLDA::M_Destroy();
}

void MetaNJetLDA
::M_SetupReadFields( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: M_SetupReadFields"
                        << METAIO_STREAM::endl;
    }

  MetaLDA::M_SetupReadFields();

  MET_FieldRecordType * mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "NZeroScales", MET_INT, false );
  m_Fields.push_back( mF );
  int nScalesRecNum = MET_GetFieldRecordNumber( "NZeroScales", &m_Fields );

  if( META_DEBUG )
    {
    METAIO_STREAM::cout << " nZeroScalesRecNum = " << nScalesRecNum
                        << METAIO_STREAM::endl;
    }

  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "ZeroScales", MET_FLOAT_ARRAY, false, nScalesRecNum );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "NFirstScales", MET_INT, false );
  m_Fields.push_back( mF );
  nScalesRecNum = MET_GetFieldRecordNumber( "NFirstScales", &m_Fields );
  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "FirstScales", MET_FLOAT_ARRAY, false, nScalesRecNum );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "NSecondScales", MET_INT, false );
  m_Fields.push_back( mF );
  nScalesRecNum = MET_GetFieldRecordNumber( "NSecondScales", &m_Fields );
  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "SecondScales", MET_FLOAT_ARRAY, false,
                     nScalesRecNum );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "NRidgeScales", MET_INT, false );
  m_Fields.push_back( mF );
  nScalesRecNum = MET_GetFieldRecordNumber( "NRidgeScales", &m_Fields );
  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "RidgeScales", MET_FLOAT_ARRAY, false, nScalesRecNum );
  m_Fields.push_back( mF );
}

void MetaNJetLDA
::M_SetupWriteFields( void )
{
  std::strcpy( m_FormTypeName, "NJetLDA" );
  MetaLDA::M_SetupWriteFields();

  MET_FieldRecordType * mF;

  if( !m_ZeroScales.empty() )
    {
    m_ZeroScalesTmp.set_size( m_ZeroScales.size() );
    for( unsigned int i = 0; i < m_ZeroScales.size(); i++ )
      {
      m_ZeroScalesTmp[i] = m_ZeroScales[i];
      }

    mF = new MET_FieldRecordType();
    MET_InitWriteField( mF, "NZeroScales", MET_INT, m_ZeroScales.size() );
    m_Fields.push_back( mF );
    int nZeroScales = m_ZeroScales.size();
    mF = new MET_FieldRecordType();
    MET_InitWriteField( mF, "ZeroScales", MET_FLOAT_ARRAY, nZeroScales,
                        m_ZeroScalesTmp.data_block() );
    m_Fields.push_back( mF );
    }

  if( !m_FirstScales.empty() )
    {
    m_FirstScalesTmp.set_size( m_FirstScales.size() );
    for( unsigned int i = 0; i < m_FirstScales.size(); i++ )
      {
      m_FirstScalesTmp[i] = m_FirstScales[i];
      }

    mF = new MET_FieldRecordType();
    MET_InitWriteField( mF, "NFirstScales", MET_INT, m_FirstScales.size() );
    m_Fields.push_back( mF );
    int nFirstScales = m_FirstScales.size();
    mF = new MET_FieldRecordType();
    MET_InitWriteField( mF, "FirstScales", MET_FLOAT_ARRAY, nFirstScales,
                        m_FirstScalesTmp.data_block() );
    m_Fields.push_back( mF );
    }

  if( !m_SecondScales.empty() )
    {
    m_SecondScalesTmp.set_size( m_SecondScales.size() );
    for( unsigned int i = 0; i < m_SecondScales.size(); i++ )
      {
      m_SecondScalesTmp[i] = m_SecondScales[i];
      }

    mF = new MET_FieldRecordType();
    MET_InitWriteField( mF, "NSecondScales", MET_INT, m_SecondScales.size() );
    m_Fields.push_back( mF );
    int nSecondScales = m_SecondScales.size();
    mF = new MET_FieldRecordType();
    MET_InitWriteField( mF, "SecondScales", MET_FLOAT_ARRAY, nSecondScales,
                        m_SecondScalesTmp.data_block() );
    m_Fields.push_back( mF );
    }

  if( !m_RidgeScales.empty() )
    {
    m_RidgeScalesTmp.set_size( m_RidgeScales.size() );
    for( unsigned int i = 0; i < m_RidgeScales.size(); i++ )
      {
      m_RidgeScalesTmp[i] = m_RidgeScales[i];
      }

    mF = new MET_FieldRecordType();
    MET_InitWriteField( mF, "NRidgeScales", MET_INT, m_RidgeScales.size() );
    m_Fields.push_back( mF );
    int nRidgeScales = m_RidgeScales.size();
    mF = new MET_FieldRecordType();
    MET_InitWriteField( mF, "RidgeScales", MET_FLOAT_ARRAY, nRidgeScales,
                        m_RidgeScalesTmp.data_block() );
    m_Fields.push_back( mF );
    }
}

bool MetaNJetLDA
::M_Read( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: M_Read: Loading header."
                        << METAIO_STREAM::endl;
    }

  if( !MetaLDA::M_Read() )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: M_Read: Error parsing file."
                        << METAIO_STREAM::endl;
    return false;
    }

  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: M_Read: Parsing header."
                        << METAIO_STREAM::endl;
    }

  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaNJetLDA: M_Read: num fields = "
                        << m_Fields.size() << METAIO_STREAM::endl;
    for( unsigned int i = 0; i < m_Fields.size(); i++ )
      {
      METAIO_STREAM::cout << " Field " << i << " = " << m_Fields[i]->name
                          << METAIO_STREAM::endl;
      }
    }

  MET_FieldRecordType * mF = MET_GetFieldRecord( "NZeroScales", &m_Fields );
  if( mF && mF->defined )
    {
    const unsigned int nZeroScales = ( unsigned int )mF->value[0];
    if( META_DEBUG )
      {
      METAIO_STREAM::cout << "MetaNJetLDA: M_Read: ZeroScales"
                          << METAIO_STREAM::endl;
      METAIO_STREAM::cout << " size = " << nZeroScales << METAIO_STREAM::endl;
      }

    m_ZeroScales.resize( nZeroScales, 0 );
    mF = MET_GetFieldRecord( "ZeroScales", &m_Fields );
    if( mF && mF->defined )
      {
      for( unsigned int i = 0; i < nZeroScales; i++ )
        {
        m_ZeroScales[i] = (double)mF->value[i];
        }
      }
    }
  else
    {
    m_ZeroScales.clear();
    }

  mF = MET_GetFieldRecord( "NFirstScales", &m_Fields );
  if( mF && mF->defined )
    {
    const unsigned int nFirstScales = ( unsigned int )mF->value[0];
    m_FirstScales.resize( nFirstScales, 0 );
    mF = MET_GetFieldRecord( "FirstScales", &m_Fields );
    if( mF && mF->defined )
      {
      for( unsigned int i = 0; i < nFirstScales; i++ )
        {
        m_FirstScales[i] = (double)mF->value[i];
        }
      }
    }
  else
    {
    m_FirstScales.clear();
    }

  mF = MET_GetFieldRecord( "NSecondScales", &m_Fields );
  if( mF && mF->defined )
    {
    const unsigned int nSecondScales = ( unsigned int )mF->value[0];
    m_SecondScales.resize( nSecondScales, 0 );
    mF = MET_GetFieldRecord( "SecondScales", &m_Fields );
    if( mF && mF->defined )
      {
      for( unsigned int i = 0; i < nSecondScales; i++ )
        {
        m_SecondScales[i] = (double)mF->value[i];
        }
      }
    }
  else
    {
    m_SecondScales.clear();
    }

  mF = MET_GetFieldRecord( "NRidgeScales", &m_Fields );
  if( mF && mF->defined )
    {
    const unsigned int nRidgeScales = ( unsigned int )mF->value[0];
    m_RidgeScales.resize( nRidgeScales, 0 );
    mF = MET_GetFieldRecord( "RidgeScales", &m_Fields );
    if( mF && mF->defined )
      {
      for( unsigned int i = 0; i < nRidgeScales; i++ )
        {
        m_RidgeScales[i] = (double)mF->value[i];
        }
      }
    }
  else
    {
    m_RidgeScales.clear();
    }

  return true;
}

} // End namespace tube

} // End namespace itk
