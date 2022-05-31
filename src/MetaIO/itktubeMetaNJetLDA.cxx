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
    std::cout << "MetaNJetLDA()" << std::endl;
    }

  this->Clear();
}

MetaNJetLDA
::MetaNJetLDA( const char * _headerName )
{
  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA()" << std::endl;
    }

  this->Clear();

  MetaNJetLDA::Read( _headerName );
}

MetaNJetLDA
::MetaNJetLDA( const MetaNJetLDA & metaNJetLDA ) : MetaLDA()
{
  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA()" << std::endl;
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
    std::cout << "MetaNJetLDA()" << std::endl;
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

  std::cout << "ZeroScales = " << m_ZeroScales.size()
                      << std::endl;
  std::cout << "FirstScales = " << m_FirstScales.size()
                      << std::endl;
  std::cout << "SecondScales = " << m_SecondScales.size()
                      << std::endl;
  std::cout << "RidgeScales = " << m_RidgeScales.size()
                      << std::endl;

  std::cout << "ZeroScalesTmp = " << m_ZeroScalesTmp.size()
                      << std::endl;
  std::cout << "FirstScalesTmp = " << m_FirstScalesTmp.size()
                      << std::endl;
  std::cout << "SecondScalesTmp = " << m_SecondScalesTmp.size()
                      << std::endl;
  std::cout << "RidgeScalesTmp = " << m_RidgeScalesTmp.size()
                      << std::endl;
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
    std::cout << "MetaNJetLDA: Clear" << std::endl;
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
    std::cout << "MetaNJetLDA: Initialize" << std::endl;
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
    std::cout << "MetaNJetLDA: SetZeroScales" << std::endl;
    std::cout << " Size = " << zeroScales.size()
                        << std::endl;
    for( unsigned int i = 0; i < zeroScales.size(); i++ )
      {
      std::cout << " Scale " << i << " = " << zeroScales[i]
                          << std::endl;
      }
    }

  m_ZeroScales = zeroScales;
}

const MetaNJetLDA::NJetScalesType & MetaNJetLDA
::GetZeroScales( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA: GetZeroScales" << std::endl;
    }

  return m_ZeroScales;
}

void MetaNJetLDA
::SetFirstScales( const NJetScalesType & firstScales )
{
  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA: SetFirstScales" << std::endl;
    }

  m_FirstScales = firstScales;
}

const MetaNJetLDA::NJetScalesType & MetaNJetLDA
::GetFirstScales( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA: GetFirstScales" << std::endl;
    }

  return m_FirstScales;
}

void MetaNJetLDA
::SetSecondScales( const NJetScalesType & secondScales )
{
  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA: SetSecondScales"
                        << std::endl;
    }

  m_SecondScales = secondScales;
}

const MetaNJetLDA::NJetScalesType & MetaNJetLDA
::GetSecondScales( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA: GetSecondScales"
                        << std::endl;
    }

  return m_SecondScales;
}

void MetaNJetLDA
::SetRidgeScales( const NJetScalesType & ridgeScales )
{
  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA: SetRidgeScales" << std::endl;
    }

  m_RidgeScales = ridgeScales;
}

const MetaNJetLDA::NJetScalesType & MetaNJetLDA
::GetRidgeScales( void ) const
{
  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA: GetRidgeScales" << std::endl;
    }

  return m_RidgeScales;
}

bool MetaNJetLDA
::CanRead( const char * headerName ) const
{
  // First check the extension.
  std::string fname = headerName;
  if( fname == "" )
    {
    return false;
    }

  bool extensionFound = false;

  std::string::size_type stringPos = fname.rfind( ".mnda" );
  if( ( stringPos != std::string::npos )
      && ( stringPos == fname.length() - 5 ) )
    {
    extensionFound = true;
    }

  if( !extensionFound )
    {
    return false;
    }

  // Now check the file content.
  std::ifstream inputStream;

  inputStream.open( headerName,
                    std::ios::in | std::ios::binary );

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

  std::ifstream * const tmpStream = new std::ifstream();

  tmpStream->open( m_FileName,
                   std::ios::in | std::ios::binary );

  if( !tmpStream->rdbuf()->is_open() )
    {
    std::cout << "MetaNJetLDA: Read: Cannot open file _"
                        << m_FileName << "_" << std::endl;
    delete tmpStream;
    return false;
    }

  const bool result = ReadStream( tmpStream );

  tmpStream->close();

  delete tmpStream;

  return result;
}

bool MetaNJetLDA
::CanReadStream( std::ifstream * stream ) const
{
  if( !std::strncmp( MET_ReadForm( *stream ).c_str(), "NJetLDA", 7 ) )
    {
    return true;
    }

  return false;
}

bool MetaNJetLDA
::ReadStream( std::ifstream * stream )
{
  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA: ReadStream" << std::endl;
    }

  this->M_Destroy();
  this->Clear();
  this->M_SetupReadFields();

  if( m_ReadStream )
    {
    std::cout << "MetaNJetLDA: ReadStream: Are two files open?"
                        << std::endl;
    delete m_ReadStream;
    }

  m_ReadStream = stream;

  if( !this->M_Read() )
    {
    std::cout << "MetaNJetLDA: Read: Cannot parse file."
                        << std::endl;
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

  std::ofstream * const tmpWriteStream = new
    std::ofstream();

  tmpWriteStream->open( m_FileName,
    std::ios::binary | std::ios::out );

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
::WriteStream( std::ofstream * stream )
{
  if( m_WriteStream != NULL )
    {
    std::cout << "MetaNJetLDA: WriteStream: Are two files open?"
                        << std::endl;
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
    std::cout << "MetaNJetLDA: M_Destroy" << std::endl;
    }

  MetaLDA::M_Destroy();
}

void MetaNJetLDA
::M_SetupReadFields( void )
{
  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA: M_SetupReadFields"
                        << std::endl;
    }

  MetaLDA::M_SetupReadFields();

  MET_FieldRecordType * mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "NZeroScales", MET_INT, false );
  m_Fields.push_back( mF );
  int nScalesRecNum = MET_GetFieldRecordNumber( "NZeroScales", &m_Fields );

  if( META_DEBUG )
    {
    std::cout << " nZeroScalesRecNum = " << nScalesRecNum
                        << std::endl;
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
    std::cout << "MetaNJetLDA: M_Read: Loading header."
                        << std::endl;
    }

  if( !MetaLDA::M_Read() )
    {
    std::cout << "MetaNJetLDA: M_Read: Error parsing file."
                        << std::endl;
    return false;
    }

  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA: M_Read: Parsing header."
                        << std::endl;
    }

  if( META_DEBUG )
    {
    std::cout << "MetaNJetLDA: M_Read: num fields = "
                        << m_Fields.size() << std::endl;
    for( unsigned int i = 0; i < m_Fields.size(); i++ )
      {
      std::cout << " Field " << i << " = " << m_Fields[i]->name
                          << std::endl;
      }
    }

  MET_FieldRecordType * mF = MET_GetFieldRecord( "NZeroScales", &m_Fields );
  if( mF && mF->defined )
    {
    const unsigned int nZeroScales = ( unsigned int )mF->value[0];
    if( META_DEBUG )
      {
      std::cout << "MetaNJetLDA: M_Read: ZeroScales"
                          << std::endl;
      std::cout << " size = " << nZeroScales << std::endl;
      }

    m_ZeroScales.resize( nZeroScales, 0 );
    mF = MET_GetFieldRecord( "ZeroScales", &m_Fields );
    if( mF && mF->defined )
      {
      for( unsigned int i = 0; i < nZeroScales; i++ )
        {
        m_ZeroScales[i] = ( double )mF->value[i];
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
        m_FirstScales[i] = ( double )mF->value[i];
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
        m_SecondScales[i] = ( double )mF->value[i];
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
        m_RidgeScales[i] = ( double )mF->value[i];
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
