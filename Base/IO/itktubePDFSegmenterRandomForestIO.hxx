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
#ifndef __itktubePDFSegmenterRandomForestIO_hxx
#define __itktubePDFSegmenterRandomForestIO_hxx

#include "itktubePDFSegmenterRandomForestIO.h"
#include "metaUtils.h"

#include "andres/ml/decision-trees.hxx"

namespace itk
{

namespace tube
{

template< class TImage, class TLabelMap >
PDFSegmenterRandomForestIO< TImage, TLabelMap >::
PDFSegmenterRandomForestIO( void )
{
  Clear();
}

template< class TImage, class TLabelMap >
PDFSegmenterRandomForestIO< TImage, TLabelMap >::
PDFSegmenterRandomForestIO( const char * _headerName )
{
  Clear();

  PDFSegmenterRandomForestIO::Read( _headerName );
}

template< class TImage, class TLabelMap >
PDFSegmenterRandomForestIO< TImage, TLabelMap >::
PDFSegmenterRandomForestIO( const typename
  PDFSegmenterType::Pointer & _filter )
{
  Clear();

  InitializeEssential( _filter );
}

template< class TImage, class TLabelMap >
PDFSegmenterRandomForestIO< TImage, TLabelMap >::
~PDFSegmenterRandomForestIO()
{
}

template< class TImage, class TLabelMap >
void
PDFSegmenterRandomForestIO< TImage, TLabelMap >::
PrintInfo() const
{
  if( m_PDFSegmenter.IsNotNull() )
    {
    std::cout << m_PDFSegmenter << std::endl;
    }
  else
    {
    std::cout << "PDFSegmenter = NULL" << std::endl;
    }
}

template< class TImage, class TLabelMap >
void
PDFSegmenterRandomForestIO< TImage, TLabelMap >::
CopyInfo( const PDFSegmenterIOType & _filterIO )
{
  Clear();

  InitializeEssential( _filterIO.GetPDFSegmenter() );
}

template< class TImage, class TLabelMap >
void
PDFSegmenterRandomForestIO< TImage, TLabelMap >::
Clear( void )
{
  m_PDFSegmenter = NULL;
}


template< class TImage, class TLabelMap >
bool PDFSegmenterRandomForestIO< TImage, TLabelMap >::
InitializeEssential( const typename
  PDFSegmenterType::Pointer & _filter )
{
  m_PDFSegmenter = _filter;

  return true;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterRandomForestIO< TImage, TLabelMap >::
SetPDFSegmenter( const typename
  PDFSegmenterType::Pointer & _filter )
{
  m_PDFSegmenter = _filter;
}

template< class TImage, class TLabelMap >
const typename PDFSegmenterRandomForest< TImage, TLabelMap >::Pointer
PDFSegmenterRandomForestIO< TImage, TLabelMap >::
GetPDFSegmenter( void ) const
{
  return m_PDFSegmenter;
}

template< class TImage, class TLabelMap >
bool PDFSegmenterRandomForestIO< TImage, TLabelMap >::
CanRead( const char * _headerName ) const
{
  // First check the extension
  METAIO_STL::string fname = _headerName;
  if( fname == "" )
    {
    return false;
    }

  bool extensionFound = false;

  METAIO_STL::string::size_type stringPos = fname.rfind( ".mpd" );
  if ( ( stringPos != METAIO_STL::string::npos )
      && ( stringPos == fname.length() - 4 ) )
    {
    extensionFound = true;
    }

  if( !extensionFound )
    {
    return false;
    }

  // Now check the file content
  METAIO_STREAM::ifstream inputStream;

  inputStream.open( fname.c_str(), METAIO_STREAM::ios::in |
    METAIO_STREAM::ios::binary );

  if( inputStream.fail() )
    {
    return false;
    }

  char* buf = new char[8001];
  inputStream.read( buf,8000 );
  unsigned long fileSize = inputStream.gcount();
  buf[fileSize] = 0;
  METAIO_STL::string header( buf );
  header.resize( fileSize );
  delete [] buf;
  inputStream.close();

  stringPos = header.find( "NDims" );
  if( stringPos == METAIO_STL::string::npos )
    {
    return false;
    }

  stringPos = header.find( "ObjectPDFFile" );
  if( stringPos == METAIO_STL::string::npos )
    {
    return false;
    }

  return true;
}

template< class TImage, class TLabelMap >
bool PDFSegmenterRandomForestIO< TImage, TLabelMap >::
Read( const char * _headerName )
{
  if( m_PDFSegmenter.IsNull() )
    {
    m_PDFSegmenter = PDFSegmenterType::New();
    }

  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaClassPDF: M_SetupReadFields"
                        << METAIO_STREAM::endl;
    }

  std::vector< MET_FieldRecordType * > metaFields;

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "NFeatures", MET_INT, true );
  metaFields.push_back( mF );
  MET_GetFieldRecordNumber( "NFeatures", &metaFields );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "NObjects", MET_INT, true );
  metaFields.push_back( mF );
  int nObjectsRec = MET_GetFieldRecordNumber( "NObjects", &metaFields );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ObjectId", MET_INT_ARRAY, true, nObjectsRec );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ObjectPDFWeight", MET_FLOAT_ARRAY, true,
    nObjectsRec );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "VoidId", MET_INT, true );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ErodeRadius", MET_INT, true );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "HoleFillIterations", MET_INT, true );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ProbabilityImageSmoothingStandardDeviation",
    MET_FLOAT, true );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ReclassifyObjectLabels", MET_STRING, true );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ReclassifyNotObjectLabels", MET_STRING, true );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ForceClassification", MET_STRING, true );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ObjectPDFFile", MET_STRING, true );
  metaFields.push_back( mF );

  // READ
  METAIO_STREAM::ifstream tmpReadStream;

  tmpReadStream.open( _headerName, METAIO_STREAM::ios::binary |
    METAIO_STREAM::ios::in );

  if( !tmpReadStream.rdbuf()->is_open() )
    {
    std::cout << "PDF::Read Could not open file." << std::endl;
    for( unsigned int i=0; i<metaFields.size(); ++i )
      {
      delete metaFields[i];
      }
    metaFields.clear();
    return false;
    }

  if( !MET_Read( tmpReadStream, &metaFields ) )
    {
    METAIO_STREAM::cerr
      << "PDFSegmenterRandomForestIO: Read: MET_Read Failed"
      << METAIO_STREAM::endl;
    for( unsigned int i=0; i<metaFields.size(); ++i )
      {
      delete metaFields[i];
      }
    metaFields.clear();
    return false;
    }

  mF = MET_GetFieldRecord( "NFeatures", &metaFields );
  if( static_cast< unsigned int >( mF->value[0] ) !=
    m_PDFSegmenter->GetNumberOfFeatures() )
    {
    std::cout << "NDims don't match: " << static_cast< int >( mF->value[0] )
      << " != " << m_PDFSegmenter->GetNumberOfFeatures() << std::endl;
    for( unsigned int i=0; i<metaFields.size(); ++i )
      {
      delete metaFields[i];
      }
    metaFields.clear();
    throw( "Expected features and features in PDF file do not match" );
    }

  std::vector< unsigned int > tmpVectorUInt;
  std::vector< int > tmpVectorInt;
  std::vector< double > tmpVectorDouble;

  mF = MET_GetFieldRecord( "NObjects", &metaFields );
  unsigned int nObjects = static_cast< unsigned int >( mF->value[0] );

  mF = MET_GetFieldRecord( "ObjectId", &metaFields );
  tmpVectorInt.resize( nObjects );
  for( unsigned int i = 0; i < nObjects; ++i )
    {
    tmpVectorInt[i] = static_cast< int >( mF->value[i] );
    }
  m_PDFSegmenter->SetObjectId( tmpVectorInt );

  mF = MET_GetFieldRecord( "ObjectPDFWeight", &metaFields );
  tmpVectorDouble.resize( nObjects );
  for( unsigned int i = 0; i < nObjects; ++i )
    {
    tmpVectorDouble[i] = static_cast< double >( mF->value[i] );
    }
  m_PDFSegmenter->SetObjectPDFWeight( tmpVectorDouble );

  mF = MET_GetFieldRecord( "VoidId", &metaFields );
  m_PDFSegmenter->SetVoidId( static_cast< int >( mF->value[0] ) );

  mF = MET_GetFieldRecord( "ErodeRadius", &metaFields );
  m_PDFSegmenter->SetErodeRadius( static_cast< int >( mF->value[0] ) );

  mF = MET_GetFieldRecord( "HoleFillIterations", &metaFields );
  m_PDFSegmenter->SetHoleFillIterations( static_cast< int >( 
    mF->value[0] ) );

  mF = MET_GetFieldRecord( "ProbabilityImageSmoothingStandardDeviation",
    &metaFields );
  m_PDFSegmenter->SetProbabilityImageSmoothingStandardDeviation( 
    static_cast< double >( mF->value[0] ) );

  mF = MET_GetFieldRecord( "ReclassifyObjectLabels", &metaFields );
  if( ( ( char * )( mF->value ) )[0] == 'T' || ( ( char * )( mF->value ) )[0] == 't' )
    {
    m_PDFSegmenter->SetReclassifyObjectLabels( true );
    }
  else
    {
    m_PDFSegmenter->SetReclassifyObjectLabels( false );
    }

  mF = MET_GetFieldRecord( "ReclassifyNotObjectLabels", &metaFields );
  if( ( ( char * )( mF->value ) )[0] == 'T' || ( ( char * )( mF->value ) )[0] == 't' )
    {
    m_PDFSegmenter->SetReclassifyNotObjectLabels( true );
    }
  else
    {
    m_PDFSegmenter->SetReclassifyNotObjectLabels( false );
    }

  mF = MET_GetFieldRecord( "ForceClassification", &metaFields );
  if( ( ( char * )( mF->value ) )[0] == 'T' || ( ( char * )( mF->value ) )[0] == 't' )
    {
    m_PDFSegmenter->SetForceClassification( true );
    }
  else
    {
    m_PDFSegmenter->SetForceClassification( false );
    }

  mF = MET_GetFieldRecord( "ObjectPDFFile", &metaFields );

  METAIO_STREAM::ifstream sFile;
  sFile.open( ( char * )( mF->value ) );
  m_PDFSegmenter->GetModel().deserialize( sFile );
  sFile.close();


  for( unsigned int i=0; i<metaFields.size(); ++i )
    {
    delete metaFields[i];
    }
  metaFields.clear();
  return true;
}

template< class TImage, class TLabelMap >
bool PDFSegmenterRandomForestIO< TImage, TLabelMap >::
Write( const char * _headerName )
{
  if( m_PDFSegmenter.IsNull() )
    {
    return false;
    }

  std::vector< MET_FieldRecordType * > metaFields;

  unsigned int numFeatures = m_PDFSegmenter->GetNumberOfFeatures();

  MET_FieldRecordType * mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "NFeatures", MET_INT, numFeatures );
  metaFields.push_back( mF );

  int tmpI[4096];
  float tmpF[4096];

  unsigned int nObjects = m_PDFSegmenter->GetNumberOfObjectIds();

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "NObjects", MET_INT, nObjects );
  metaFields.push_back( mF );

  for( unsigned int i = 0; i < nObjects; ++i )
    {
    tmpI[i] = m_PDFSegmenter->GetObjectId()[i];
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField< int >( mF, "ObjectId", MET_INT_ARRAY, nObjects, tmpI );
  metaFields.push_back( mF );

  for( unsigned int i = 0; i < nObjects; ++i )
    {
    tmpF[i] = m_PDFSegmenter->GetObjectPDFWeight()[i];
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField< float >( mF, "ObjectPDFWeight", MET_FLOAT_ARRAY,
    nObjects, tmpF );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "VoidId", MET_INT, m_PDFSegmenter->GetVoidId() );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ErodeRadius", MET_INT,
    m_PDFSegmenter->GetErodeRadius() );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "HoleFillIterations", MET_INT,
    m_PDFSegmenter->GetHoleFillIterations() );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ProbabilityImageSmoothingStandardDeviation",
    MET_FLOAT, m_PDFSegmenter->
    GetProbabilityImageSmoothingStandardDeviation() );
  metaFields.push_back( mF );

  char tmpC[4096];
  if( m_PDFSegmenter->GetReclassifyObjectLabels() )
    {
    strcpy( tmpC, "True" );
    }
  else
    {
    strcpy( tmpC, "False" );
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField< const char >( mF, "ReclassifyObjectLabels",
    MET_STRING, strlen( tmpC ), tmpC );
  metaFields.push_back( mF );

  if( m_PDFSegmenter->GetReclassifyNotObjectLabels() )
    {
    strcpy( tmpC, "True" );
    }
  else
    {
    strcpy( tmpC, "False" );
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField< const char >( mF, "ReclassifyNotObjectLabels",
    MET_STRING, strlen( tmpC ), tmpC );
  metaFields.push_back( mF );

  if( m_PDFSegmenter->GetForceClassification() )
    {
    strcpy( tmpC, "True" );
    }
  else
    {
    strcpy( tmpC, "False" );
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField< const char >( mF, "ForceClassification", MET_STRING,
    strlen( tmpC ), tmpC );
  metaFields.push_back( mF );

  char filePath[255];
  MET_GetFilePath( _headerName, filePath );
  int skip = strlen( filePath );
  char shortFileName[255];
  sprintf( shortFileName, "%s", &( _headerName[skip] ) );
  MET_SetFileSuffix( shortFileName, "mrf" );
  std::string fullFileName = filePath;
  fullFileName = fullFileName + shortFileName;

  std::string modelFileName = shortFileName;
  modelFileName = modelFileName + ".rff";
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ObjectPDFFile", MET_STRING,
    modelFileName.size(), modelFileName.c_str() );
  metaFields.push_back( mF );

  METAIO_STREAM::ofstream writeStream;

  writeStream.open( fullFileName.c_str(), METAIO_STREAM::ios::binary |
    METAIO_STREAM::ios::out );

  writeStream.precision( 10 );

  if( !MET_Write( writeStream, & metaFields ) )
    {
    METAIO_STREAM::cerr << "MetaObject: Write: MET_Write Failed"
                        << METAIO_STREAM::endl;
    for( unsigned int i=0; i<metaFields.size(); ++i )
      {
      delete metaFields[i];
      }
    metaFields.clear();
    return false;
    }

  METAIO_STREAM::ofstream sFile;
  sFile.open( modelFileName.c_str() );
  m_PDFSegmenter->GetModel().serialize( sFile );
  sFile.close();

  for( unsigned int i=0; i<metaFields.size(); ++i )
    {
    delete metaFields[i];
    }
  metaFields.clear();
  return true;
}

} // End namespace tube

} // End namespace itk

#endif // __itktubePDFSegmenterRandomForestIO_hxx
