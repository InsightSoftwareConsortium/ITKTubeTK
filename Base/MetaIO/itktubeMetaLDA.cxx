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

#include <cstring>

#include "itktubeMetaLDA.h"

#include "metaUtilsTemp.h"

namespace itk
{

namespace tube
{

MetaLDA
::MetaLDA( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA()" << METAIO_STREAM::endl;
    }

  this->Clear();
}

MetaLDA
::MetaLDA( const char * headerName )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA()" << METAIO_STREAM::endl;
    }

  this->Clear();
  m_ReadStream = NULL;

  MetaLDA::Read( headerName );
}

MetaLDA
::MetaLDA( const MetaLDA & metaLDA ) : MetaForm()
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA()" << METAIO_STREAM::endl;
    }

  this->Clear();
  this->CopyInfo( metaLDA );
}

MetaLDA
::MetaLDA( unsigned int numberOfPCAToUseAsFeatures,
  unsigned int numberOfLDAToUseAsFeatures,
  const LDAValuesType & ldaValues,
  const LDAMatrixType & ldaMatrix,
  const ValueListType & inputWhitenMeans,
  const ValueListType & inputWhitenStdDevs,
  const ValueListType & outputWhitenMeans,
  const ValueListType & outputWhitenStdDevs )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA()" << METAIO_STREAM::endl;
    }

  this->Clear();
  this->InitializeEssential( numberOfPCAToUseAsFeatures,
    numberOfLDAToUseAsFeatures, ldaValues, ldaMatrix, inputWhitenMeans,
    inputWhitenStdDevs, outputWhitenMeans, outputWhitenStdDevs );
}

MetaLDA
::~MetaLDA( void )
{
  this->M_Destroy();
}

void MetaLDA
::PrintInfo( void ) const
{
  MetaForm::PrintInfo();

  METAIO_STREAM::cout << "NumberOfPCABasisToUseAsFeatures = "
    << m_NumberOfPCABasisToUseAsFeatures << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "NumberOfLDABasisToUseAsFeatures = "
    << m_NumberOfLDABasisToUseAsFeatures << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "LDAValues = " << m_LDAValues
    << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "LDAMatrix = " << m_LDAMatrix
    << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "InputWhitenMeans = " << METAIO_STREAM::endl;
  for( unsigned int i = 0; i < m_InputWhitenMeans.size(); i++ )
    {
    METAIO_STREAM::cout << m_InputWhitenMeans[i] << " ";
    }
  METAIO_STREAM::cout << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "InputWhitenStdDevs = " << METAIO_STREAM::endl;
  for( unsigned int i = 0; i < m_InputWhitenStdDevs.size(); i++ )
    {
    METAIO_STREAM::cout << m_InputWhitenStdDevs[i] << " ";
    }

  METAIO_STREAM::cout << "OutputWhitenMeans = " << METAIO_STREAM::endl;
  for( unsigned int i = 0; i < m_OutputWhitenMeans.size(); i++ )
    {
    METAIO_STREAM::cout << m_OutputWhitenMeans[i] << " ";
    }
  METAIO_STREAM::cout << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "OutputWhitenStdDevs = " << METAIO_STREAM::endl;
  for( unsigned int i = 0; i < m_OutputWhitenStdDevs.size(); i++ )
    {
    METAIO_STREAM::cout << m_OutputWhitenStdDevs[i] << " ";
    }

  METAIO_STREAM::cout << METAIO_STREAM::endl;
}

void MetaLDA
::CopyInfo( const MetaLDA & lda )
{
  MetaForm::CopyInfo( dynamic_cast< const MetaForm * >( &lda ) );
  this->SetNumberOfPCABasisToUseAsFeatures(
    lda.GetNumberOfPCABasisToUseAsFeatures() );
  this->SetNumberOfLDABasisToUseAsFeatures(
    lda.GetNumberOfLDABasisToUseAsFeatures() );
  this->SetLDAValues( lda.GetLDAValues() );
  this->SetLDAMatrix( lda.GetLDAMatrix() );
  this->SetInputWhitenMeans( lda.GetInputWhitenMeans() );
  this->SetInputWhitenStdDevs( lda.GetInputWhitenStdDevs() );
  this->SetOutputWhitenMeans( lda.GetOutputWhitenMeans() );
  this->SetOutputWhitenStdDevs( lda.GetOutputWhitenStdDevs() );
}

void MetaLDA
::Clear( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: Clear" << METAIO_STREAM::endl;
    }

  MetaForm::Clear();

  strcpy( m_FormTypeName, "LDA" );

  m_NumberOfPCABasisToUseAsFeatures = 0;
  m_NumberOfLDABasisToUseAsFeatures = 0;
  m_LDAValues.set_size( 0 );
  m_LDAMatrix.set_size( 0, 0 );
  m_InputWhitenMeans.clear();
  m_InputWhitenStdDevs.clear();
  m_OutputWhitenMeans.clear();
  m_OutputWhitenStdDevs.clear();
}

bool MetaLDA
::InitializeEssential( unsigned int numberOfPCAToUseAsFeatures,
  unsigned int numberOfLDAToUseAsFeatures,
  const LDAValuesType & ldaValues,
  const LDAMatrixType & ldaMatrix,
  const ValueListType & inputWhitenMeans,
  const ValueListType & inputWhitenStdDevs,
  const ValueListType & outputWhitenMeans,
  const ValueListType & outputWhitenStdDevs )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: Initialize" << METAIO_STREAM::endl;
    }

  MetaForm::InitializeEssential();
  this->SetNumberOfPCABasisToUseAsFeatures( numberOfPCAToUseAsFeatures );
  this->SetNumberOfLDABasisToUseAsFeatures( numberOfLDAToUseAsFeatures );
  this->SetLDAValues( ldaValues );
  this->SetLDAMatrix( ldaMatrix );
  this->SetInputWhitenMeans( inputWhitenMeans );
  this->SetInputWhitenStdDevs( inputWhitenStdDevs );
  this->SetOutputWhitenMeans( outputWhitenMeans );
  this->SetOutputWhitenStdDevs( outputWhitenStdDevs );
  return true;
}

void MetaLDA
::SetNumberOfPCABasisToUseAsFeatures( unsigned int
  numberOfPCABasisToUseAsFeatures )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: SetNumberOfPCABasisToUseAsFeatures"
      << METAIO_STREAM::endl;
    }

  m_NumberOfPCABasisToUseAsFeatures = numberOfPCABasisToUseAsFeatures;
}

unsigned int MetaLDA
::GetNumberOfPCABasisToUseAsFeatures( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: GetNumberOfPCABasisToUseAsFeatures"
      << METAIO_STREAM::endl;
    }

  return m_NumberOfPCABasisToUseAsFeatures;
}

void MetaLDA
::SetNumberOfLDABasisToUseAsFeatures( unsigned int
  numberOfLDABasisToUseAsFeatures )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: SetNumberOfLDABasisToUseAsFeatures"
      << METAIO_STREAM::endl;
    }

  m_NumberOfLDABasisToUseAsFeatures = numberOfLDABasisToUseAsFeatures;
}

unsigned int MetaLDA
::GetNumberOfLDABasisToUseAsFeatures( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: GetNumberOfLDABasisToUseAsFeatures"
      << METAIO_STREAM::endl;
    }

  return m_NumberOfLDABasisToUseAsFeatures;
}

void MetaLDA
::SetLDAValues( const LDAValuesType & ldaValues )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: SetLDAValues" << METAIO_STREAM::endl;
    }

  m_LDAValues = ldaValues;
}

const MetaLDA::LDAValuesType & MetaLDA
::GetLDAValues( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: GetLDAValues" << METAIO_STREAM::endl;
    }

  return m_LDAValues;
}

void MetaLDA
::SetLDAMatrix( const LDAMatrixType & ldaMatrix )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: SetLDAMatrix" << METAIO_STREAM::endl;
    }

  m_LDAMatrix = ldaMatrix;
}

const MetaLDA::LDAMatrixType & MetaLDA
::GetLDAMatrix( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: GetLDAMatrix" << METAIO_STREAM::endl;
    }

  return m_LDAMatrix;
}

void MetaLDA
::SetInputWhitenMeans( const ValueListType & whitenMeans )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: SetInputWhitenMeans"
      << METAIO_STREAM::endl;
    }

  m_InputWhitenMeans = whitenMeans;
}

const MetaLDA::ValueListType & MetaLDA
::GetInputWhitenMeans( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: GetInputWhitenMeans"
      << METAIO_STREAM::endl;
    }

  return m_InputWhitenMeans;
}

void MetaLDA
::SetInputWhitenStdDevs( const ValueListType & whitenStdDevs )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: SetInputWhitenStdDevs"
      << METAIO_STREAM::endl;
    }

  m_InputWhitenStdDevs = whitenStdDevs;
}

const MetaLDA::ValueListType & MetaLDA
::GetInputWhitenStdDevs( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: GetInputWhitenStdDevs"
      << METAIO_STREAM::endl;
    }

  return m_InputWhitenStdDevs;
}

void MetaLDA
::SetOutputWhitenMeans( const ValueListType & whitenMeans )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: SetOutputWhitenMeans"
      << METAIO_STREAM::endl;
    }

  m_OutputWhitenMeans = whitenMeans;
}

const MetaLDA::ValueListType & MetaLDA
::GetOutputWhitenMeans( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: GetOutputWhitenMeans"
      << METAIO_STREAM::endl;
    }

  return m_OutputWhitenMeans;
}

void MetaLDA
::SetOutputWhitenStdDevs( const ValueListType & whitenStdDevs )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: SetOutputWhitenStdDevs"
      << METAIO_STREAM::endl;
    }

  m_OutputWhitenStdDevs = whitenStdDevs;
}

const MetaLDA::ValueListType & MetaLDA
::GetOutputWhitenStdDevs( void ) const
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: GetOutputWhitenStdDevs"
      << METAIO_STREAM::endl;
    }

  return m_OutputWhitenStdDevs;
}

bool MetaLDA
::CanRead( const char * headerName ) const
{
  // First check the extension.
  METAIO_STL::string fname = headerName;
  if( fname == "" )
    {
    std::cout << "CanRead: Error: Empty file name." << std::endl;
    return false;
    }

  bool extensionFound = false;

  METAIO_STL::string::size_type stringPos = fname.rfind( ".mlda" );
  if( ( stringPos != METAIO_STL::string::npos )
      && ( stringPos == fname.length() - 5 ) )
    {
    extensionFound = true;
    }

  if( !extensionFound )
    {
    std::cout << "CanRead: Error: Extension not supported." << std::endl;
    return false;
    }

  // Now check the file content.
  METAIO_STREAM::ifstream inputStream;

  inputStream.open( headerName,
                    METAIO_STREAM::ios::in | METAIO_STREAM::ios::binary );

  if( !inputStream.rdbuf()->is_open() )
    {
    std::cout << "CanRead: Error: Cannot open file." << std::endl;
    return false;
    }

  const bool result = !std::strncmp( MET_ReadForm(
                                       inputStream ).c_str(), "LDA", 3 );

  if( !result )
    {
    std::cout << "CanRead: Error: Read form failed." << std::endl;
    }

  inputStream.close();

  return result;
}

bool MetaLDA
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
    METAIO_STREAM::cout << "MetaLDA: Read: Cannot open file _" << m_FileName
                        << "_" << METAIO_STREAM::endl;
    delete tmpStream;
    return false;
    }

  const bool result = ReadStream( tmpStream );

  tmpStream->close();
  delete tmpStream;

  return result;
}

bool MetaLDA
::CanReadStream( METAIO_STREAM::ifstream * stream ) const
{
  if( !std::strncmp( MET_ReadForm( *stream ).c_str(), "LDA", 3 ) )
    {
    return true;
    }

  return false;
}

bool MetaLDA
::ReadStream( METAIO_STREAM::ifstream * stream )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: ReadStream" << METAIO_STREAM::endl;
    }

  this->M_Destroy();
  this->Clear();
  this->M_SetupReadFields();

  if( m_ReadStream )
    {
    METAIO_STREAM::cout << "MetaLDA: ReadStream: Are two files open?"
                        << METAIO_STREAM::endl;
    delete m_ReadStream;
    }

  m_ReadStream = stream;

  if( !this->M_Read() )
    {
    METAIO_STREAM::cout << "MetaLDA: Read: Cannot parse file."
                        << METAIO_STREAM::endl;
    m_ReadStream = NULL;
    return false;
    }

  m_ReadStream = NULL;

  this->InitializeEssential( m_NumberOfPCABasisToUseAsFeatures,
    m_NumberOfLDABasisToUseAsFeatures, m_LDAValues, m_LDAMatrix,
    m_InputWhitenMeans, m_InputWhitenStdDevs,
    m_OutputWhitenMeans, m_OutputWhitenStdDevs );

  return true;
}

bool MetaLDA
::Write( const char * headerName )
{
  if( headerName != NULL && std::strlen( headerName ) > 1 )
    {
    this->FileName( headerName );
    }

  MET_SetFileSuffix( m_FileName, "mlda" );

  METAIO_STREAM::ofstream * const tmpWriteStream = new
    METAIO_STREAM::ofstream();

  tmpWriteStream->open( m_FileName, METAIO_STREAM::ios::binary
    | METAIO_STREAM::ios::out );

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

bool MetaLDA
::WriteStream( METAIO_STREAM::ofstream * stream )
{
  if( m_WriteStream != NULL )
    {
    METAIO_STREAM::cout << "MetaLDA: WriteStream: Are two files open?"
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

void MetaLDA
::M_Destroy( void )
{
  MetaForm::M_Destroy();
}

void MetaLDA
::M_SetupReadFields( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: M_SetupReadFields" << METAIO_STREAM::endl;
    }

  MetaForm::M_SetupReadFields();

  MET_FieldRecordType * mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "NDims", MET_INT, true );
  m_Fields.push_back( mF );

  const int nDimsRecNum = MET_GetFieldRecordNumber( "NDims", &m_Fields );

  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "NPCABasis", MET_INT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "NLDABasis", MET_INT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "Values", MET_FLOAT_ARRAY, true, nDimsRecNum );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "Matrix", MET_FLOAT_MATRIX, true, nDimsRecNum );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "InputWhitenMeans", MET_FLOAT_ARRAY, false,
    nDimsRecNum );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "InputWhitenStdDevs", MET_FLOAT_ARRAY, false,
    nDimsRecNum );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "OutputWhitenMeans", MET_FLOAT_ARRAY, false,
    nDimsRecNum );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  MET_InitReadField( mF, "OutputWhitenStdDevs", MET_FLOAT_ARRAY, false,
    nDimsRecNum );
  m_Fields.push_back( mF );
}

void MetaLDA
::M_SetupWriteFields( void )
{
  std::strcpy( m_FormTypeName, "LDA" );
  MetaForm::M_SetupWriteFields();

  MET_FieldRecordType * mF = new MET_FieldRecordType();
  MET_InitWriteField( mF, "NDims", MET_INT, m_LDAValues.size() );
  m_Fields.push_back( mF );

  const int nDims = m_LDAValues.size();

  mF = new MET_FieldRecordType();
  MET_InitWriteField( mF, "NPCABasis", MET_INT,
    m_NumberOfPCABasisToUseAsFeatures );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  MET_InitWriteField( mF, "NLDABasis", MET_INT,
    m_NumberOfLDABasisToUseAsFeatures );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  ::tube::MET_InitWriteField_Temp( mF, "Values", MET_FLOAT_ARRAY, nDims,
                      m_LDAValues.data_block() );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType();
  ::tube::MET_InitWriteField_Temp( mF, "Matrix", MET_FLOAT_MATRIX, nDims,
                      m_LDAMatrix.data_block() );
  m_Fields.push_back( mF );

  int tfCount = m_InputWhitenMeans.size();

  if( tfCount > 0 )
    {
    double tf[4096];
    mF = new MET_FieldRecordType();
    for( int i = 0; i < tfCount; i++ )
      {
      tf[i] = m_InputWhitenMeans[i];
      }
    for( int i = tfCount; i < nDims; i++ )
      {
      tf[i] = 0;
      }

    MET_InitWriteField( mF, "InputWhitenMeans", MET_FLOAT_ARRAY, nDims,
      tf );
    m_Fields.push_back( mF );

    mF = new MET_FieldRecordType();
    for( int i = 0; i < tfCount; i++ )
      {
      tf[i] = m_InputWhitenStdDevs[i];
      }

    MET_InitWriteField( mF, "InputWhitenStdDevs", MET_FLOAT_ARRAY, tfCount,
      tf );
    m_Fields.push_back( mF );
    }

  tfCount = m_OutputWhitenMeans.size();

  if( tfCount > 0 )
    {
    double tf[4096];
    mF = new MET_FieldRecordType();
    for( int i = 0; i < tfCount; i++ )
      {
      tf[i] = m_OutputWhitenMeans[i];
      }
    for( int i = tfCount; i < nDims; i++ )
      {
      tf[i] = 0;
      }

    MET_InitWriteField( mF, "OutputWhitenMeans", MET_FLOAT_ARRAY, nDims,
      tf );
    m_Fields.push_back( mF );

    mF = new MET_FieldRecordType();
    for( int i = 0; i < tfCount; i++ )
      {
      tf[i] = m_OutputWhitenStdDevs[i];
      }

    MET_InitWriteField( mF, "OutputWhitenStdDevs", MET_FLOAT_ARRAY, nDims,
      tf );
    m_Fields.push_back( mF );
    }
}

bool MetaLDA
::M_Read( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: M_Read: Loading header."
                        << METAIO_STREAM::endl;
    }

  if( !MetaForm::M_Read() )
    {
    METAIO_STREAM::cout << "MetaLDA: M_Read: Error parsing file."
                        << METAIO_STREAM::endl;
    return false;
    }

  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaLDA: M_Read: Parsing header."
                        << METAIO_STREAM::endl;
    }

  unsigned int nDims = 0;
  MET_FieldRecordType * mF = MET_GetFieldRecord( "NDims", &m_Fields );
  if( mF && mF->defined )
    {
    nDims = ( unsigned int )mF->value[0];
    m_LDAValues.set_size( nDims );
    m_LDAValues.fill( 0 );
    m_LDAMatrix.set_size( nDims, nDims );
    m_LDAMatrix.fill( 0 );
    m_InputWhitenMeans.resize( nDims );
    m_InputWhitenStdDevs.resize( nDims );
    m_OutputWhitenMeans.resize( nDims );
    m_OutputWhitenStdDevs.resize( nDims );
    }
  else
    {
    METAIO_STREAM::cout << "MetaLDA: M_Read: Error: NDims required."
                        << METAIO_STREAM::endl;
    return false;
    }

  mF = MET_GetFieldRecord( "NPCABasis", &m_Fields );
  if( mF && mF->defined )
    {
    m_NumberOfPCABasisToUseAsFeatures = ( unsigned int )mF->value[0];
    }

  mF = MET_GetFieldRecord( "NLDABasis", &m_Fields );
  if( mF && mF->defined )
    {
    m_NumberOfLDABasisToUseAsFeatures = ( unsigned int )mF->value[0];
    }

  mF = MET_GetFieldRecord( "Values", &m_Fields );
  if( mF && mF->defined )
    {
    for( unsigned int i = 0; i < nDims; i++ )
      {
      m_LDAValues[i] = (double)mF->value[i];
      }
    }
  else
    {
    METAIO_STREAM::cout << "MetaLDA: M_Read: Error: Values required."
                        << METAIO_STREAM::endl;
    return false;
    }

  mF = MET_GetFieldRecord( "Matrix", &m_Fields );
  if( mF && mF->defined )
    {
    for( unsigned int i = 0; i < nDims; i++ )
      {
      for( unsigned int j = 0; j < nDims; j++ )
        {
        m_LDAMatrix[i][j] = (double)mF->value[i * nDims + j];
        }
      }
    }
  else
    {
    METAIO_STREAM::cout << "MetaLDA: M_Read: Error: Matrix required."
                        << METAIO_STREAM::endl;
    return false;
    }

  mF = MET_GetFieldRecord( "InputWhitenMeans", &m_Fields );
  if( mF && mF->defined )
    {
    for( unsigned int i = 0; i < nDims; i++ )
      {
      m_InputWhitenMeans[i] = (double)mF->value[i];
      }
    }
  else
    {
    m_InputWhitenMeans.clear();
    }

  mF = MET_GetFieldRecord( "InputWhitenStdDevs", &m_Fields );
  if( mF && mF->defined )
    {
    for( unsigned int i = 0; i < nDims; i++ )
      {
      m_InputWhitenStdDevs[i] = (double)mF->value[i];
      }
    }
  else
    {
    m_InputWhitenStdDevs.clear();
    }

  mF = MET_GetFieldRecord( "OutputWhitenMeans", &m_Fields );
  if( mF && mF->defined )
    {
    for( unsigned int i = 0; i < nDims; i++ )
      {
      m_OutputWhitenMeans[i] = (double)mF->value[i];
      }
    }
  else
    {
    m_OutputWhitenMeans.clear();
    }

  mF = MET_GetFieldRecord( "OutputWhitenStdDevs", &m_Fields );
  if( mF && mF->defined )
    {
    for( unsigned int i = 0; i < nDims; i++ )
      {
      m_OutputWhitenStdDevs[i] = (double)mF->value[i];
      }
    }
  else
    {
    m_OutputWhitenStdDevs.clear();
    }

  return true;
}

} // End namespace tube

} // End namespace itk
