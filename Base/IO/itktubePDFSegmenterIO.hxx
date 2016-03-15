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
#ifndef __itktubePDFSegmenterIO_hxx
#define __itktubePDFSegmenterIO_hxx

#include "itktubePDFSegmenterIO.h"
#include "itktubeMetaClassPDF.h"
#include "metaUtils.h"

namespace itk
{

namespace tube
{

template< class TImage, unsigned int N, class TLabelMap >
PDFSegmenterIO< TImage, N, TLabelMap >::
PDFSegmenterIO( void )
{
  Clear();
}

template< class TImage, unsigned int N, class TLabelMap >
PDFSegmenterIO< TImage, N, TLabelMap >::
PDFSegmenterIO( const char * _headerName )
{
  Clear();

  PDFSegmenterIO::Read( _headerName );
}

template< class TImage, unsigned int N, class TLabelMap >
PDFSegmenterIO< TImage, N, TLabelMap >::
PDFSegmenterIO( const typename
  PDFSegmenterType::Pointer & _filter )
{
  Clear();

  InitializeEssential( _filter );
}

template< class TImage, unsigned int N, class TLabelMap >
PDFSegmenterIO< TImage, N, TLabelMap >::
~PDFSegmenterIO()
{
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenterIO< TImage, N, TLabelMap >::
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

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenterIO< TImage, N, TLabelMap >::
CopyInfo( const PDFSegmenterIOType & _filterIO )
{
  Clear();

  InitializeEssential( _filterIO.GetPDFSegmenter() );
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenterIO< TImage, N, TLabelMap >::
Clear( void )
{
  m_PDFSegmenter = NULL;
}


template< class TImage, unsigned int N, class TLabelMap >
bool PDFSegmenterIO< TImage, N, TLabelMap >::
InitializeEssential( const typename
  PDFSegmenterType::Pointer & _filter )
{
  m_PDFSegmenter = _filter;

  return true;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenterIO< TImage, N, TLabelMap >::
SetPDFSegmenter( const typename
  PDFSegmenterType::Pointer & _filter )
{
  m_PDFSegmenter = _filter;
}

template< class TImage, unsigned int N, class TLabelMap >
const typename PDFSegmenter< TImage, N, TLabelMap >::Pointer
PDFSegmenterIO< TImage, N, TLabelMap >::
GetPDFSegmenter( void ) const
{
  return m_PDFSegmenter;
}

template< class TImage, unsigned int N, class TLabelMap >
bool PDFSegmenterIO< TImage, N, TLabelMap >::
CanRead( const char * _headerName ) const
{
  // First check the extension
  METAIO_STL::string fname = _headerName;
  if( fname == "" )
    {
    return false;
    }

  bool extensionFound = false;

  METAIO_STL::string::size_type stringPos = fname.rfind(".mpd");
  if ((stringPos != METAIO_STL::string::npos)
      && (stringPos == fname.length() - 4))
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
  inputStream.read(buf,8000);
  unsigned long fileSize = inputStream.gcount();
  buf[fileSize] = 0;
  METAIO_STL::string header(buf);
  header.resize(fileSize);
  delete [] buf;
  inputStream.close();

  stringPos = header.find("NDims");
  if( stringPos == METAIO_STL::string::npos )
    {
    return false;
    }

  stringPos = header.find("ObjectPDFFile");
  if( stringPos == METAIO_STL::string::npos )
    {
    return false;
    }

  return true;
}

template< class TImage, unsigned int N, class TLabelMap >
bool PDFSegmenterIO< TImage, N, TLabelMap >::
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

  MetaClassPDF classPDFReader;

  std::vector< MET_FieldRecordType * > metaFields;

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "NFeatures", MET_INT, true );
  metaFields.push_back( mF );
  int nFeaturesRec = MET_GetFieldRecordNumber( "NFeatures", &metaFields );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "BinsPerFeature", MET_INT_ARRAY, true,
    nFeaturesRec );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "BinMin", MET_FLOAT_ARRAY, true, nFeaturesRec );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "BinSize", MET_FLOAT_ARRAY, true, nFeaturesRec );
  metaFields.push_back( mF );

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
  MET_InitReadField( mF, "HistogramSmoothingStandardDeviation",
    MET_FLOAT, true );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "OutlierRejectPortion", MET_FLOAT, true );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "Draft", MET_STRING, true );
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

  tmpReadStream.open(_headerName, METAIO_STREAM::ios::binary |
    METAIO_STREAM::ios::in);

  if(!tmpReadStream.rdbuf()->is_open())
    {
    std::cout << "PDF::Read Could not open file." << std::endl;
    return false;
    }

  if(!MET_Read( tmpReadStream, &metaFields ) )
    {
    METAIO_STREAM::cerr << "PDFSegmenterIO: Read: MET_Read Failed"
      << METAIO_STREAM::endl;
    return false;
    }

  mF = MET_GetFieldRecord( "NFeatures", &metaFields );
  if( static_cast< unsigned int >( mF->value[0] ) !=
    m_PDFSegmenter->GetNumberOfFeatures() )
    {
    std::cout << "NDims don't match: " << static_cast< int >( mF->value[0] )
      << " != " << m_PDFSegmenter->GetNumberOfFeatures() << std::endl;
    throw( "Expected features and features in PDF file do not match" );
    }

  std::vector< unsigned int > tmpVectorUInt;
  std::vector< int > tmpVectorInt;
  std::vector< double > tmpVectorDouble;

  mF = MET_GetFieldRecord( "BinsPerFeature", &metaFields );
  tmpVectorUInt.resize( N );
  for( unsigned int i = 0; i < N; ++i )
    {
    tmpVectorUInt[i] = static_cast< int >( mF->value[i] );
    }
  m_PDFSegmenter->SetNumberOfBinsPerFeature( tmpVectorUInt );

  mF = MET_GetFieldRecord( "BinMin", &metaFields );
  tmpVectorDouble.resize( N );
  for( unsigned int i = 0; i < N; ++i )
    {
    tmpVectorDouble[i] = static_cast< double >( mF->value[i] );
    }
  m_PDFSegmenter->SetBinMin( tmpVectorDouble );

  mF = MET_GetFieldRecord( "BinSize", &metaFields );
  tmpVectorDouble.resize( N );
  for( unsigned int i = 0; i < N; ++i )
    {
    tmpVectorDouble[i] = static_cast< double >( mF->value[i] );
    }
  m_PDFSegmenter->SetBinSize( tmpVectorDouble );

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

  mF = MET_GetFieldRecord( "HistogramSmoothingStandardDeviation",
    &metaFields );
  m_PDFSegmenter->SetHistogramSmoothingStandardDeviation(
    static_cast< double >( mF->value[0] ) );

  mF = MET_GetFieldRecord( "OutlierRejectPortion", &metaFields );
  m_PDFSegmenter->SetOutlierRejectPortion( static_cast< double >(
    mF->value[0] ) );

  mF = MET_GetFieldRecord( "Draft", &metaFields );
  if( ((char *)( mF->value))[0] == 'T' || ((char *)( mF->value))[0] == 't' )
    {
    m_PDFSegmenter->SetDraft( true );
    }
  else
    {
    m_PDFSegmenter->SetDraft( false );
    }

  mF = MET_GetFieldRecord( "ReclassifyObjectLabels", &metaFields );
  if( ((char *)( mF->value))[0] == 'T' || ((char *)( mF->value))[0] == 't' )
    {
    m_PDFSegmenter->SetReclassifyObjectLabels( true );
    }
  else
    {
    m_PDFSegmenter->SetReclassifyObjectLabels( false );
    }

  mF = MET_GetFieldRecord( "ReclassifyNotObjectLabels", &metaFields );
  if( ((char *)( mF->value))[0] == 'T' || ((char *)( mF->value))[0] == 't' )
    {
    m_PDFSegmenter->SetReclassifyNotObjectLabels( true );
    }
  else
    {
    m_PDFSegmenter->SetReclassifyNotObjectLabels( false );
    }

  mF = MET_GetFieldRecord( "ForceClassification", &metaFields );
  if( ((char *)( mF->value))[0] == 'T' || ((char *)( mF->value))[0] == 't' )
    {
    m_PDFSegmenter->SetForceClassification( true );
    }
  else
    {
    m_PDFSegmenter->SetForceClassification( false );
    }

  mF = MET_GetFieldRecord( "ObjectPDFFile", &metaFields );
  std::string str = (char *)(mF->value);
  std::vector< std::string > fileName;
  MET_StringToVector( str, fileName );
  if( fileName.size() != nObjects )
    {
    std::cerr << "Number of PDFFiles != number of objects" << std::endl;
    m_PDFSegmenter = NULL;
    return false;
    }

  for( unsigned int i = 0; i < nObjects; ++i )
    {
    typedef Image< float, N >  pdfImageType;

    char filePath[255];
    MET_GetFilePath( _headerName, filePath );
    std::string fullFileName = filePath + fileName[i];

    MetaClassPDF pdfClassReader( fullFileName.c_str() );

    typename pdfImageType::Pointer img = pdfImageType::New();
    typename pdfImageType::RegionType region;
    typename pdfImageType::SizeType size;
    typename pdfImageType::PointType origin;
    typename pdfImageType::SpacingType spacing;
    for( unsigned int j = 0; j < N; ++j )
      {
      size[j] = pdfClassReader.GetNumberOfBinsPerFeature()[j];
      if( size[j] != m_PDFSegmenter->GetNumberOfBinsPerFeature()[j] )
        {
        std::cout << "ERROR: N mismatch" << std::endl;
        return false;
        }
      spacing[j] = pdfClassReader.GetBinSize()[j];
      if( vnl_math_abs( spacing[j] - m_PDFSegmenter->GetBinSize()[j] ) >
        0.005 * spacing[j] )
        {
        std::cout << "ERROR: Spacing mismatch" << std::endl;
        return false;
        }
      origin[j] = pdfClassReader.GetBinMin()[j];
      if( vnl_math_abs( origin[j] - m_PDFSegmenter->GetBinMin()[j] ) >
        0.005 * spacing[j] )
        {
        std::cout << "ERROR: Min mismatch" << std::endl;
        std::cout << "   Origin[" << j << "] = " << origin[j] << std::endl;
        std::cout << "   PDFMin[" << j << "] = "
          << m_PDFSegmenter->GetBinMin()[j] << std::endl;
        std::cout << "      spacing[" << j << "] = "
          << spacing[j] << std::endl;
        return false;
        }
      }
    region.SetSize( size );
    img->SetRegions( region );
    img->SetOrigin( origin );
    img->SetSpacing( spacing );

    int nPixels = size[0];
    for( unsigned int j = 1; j < N; ++j )
      {
      nPixels *= size[j];
      }
    img->GetPixelContainer()->SetImportPointer( pdfClassReader.ExportPDF(),
      nPixels, true );

    m_PDFSegmenter->SetClassPDFImage( i, img );
    }

  return true;
}

template< class TImage, unsigned int N, class TLabelMap >
bool PDFSegmenterIO< TImage, N, TLabelMap >::
Write( const char * _headerName )
{
  if( m_PDFSegmenter.IsNull() )
    {
    return false;
    }

  std::vector< MET_FieldRecordType * > metaFields;

  MET_FieldRecordType * mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "NFeatures", MET_INT, N );
  metaFields.push_back( mF );

  int tmpI[4096];
  for( unsigned int i = 0; i < N; ++i )
    {
    tmpI[i] = m_PDFSegmenter->GetNumberOfBinsPerFeature()[i];
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "BinsPerFeature", MET_INT_ARRAY, N, tmpI );
  metaFields.push_back( mF );

  float tmpF[4096];
  for( unsigned int i = 0; i < N; ++i )
    {
    tmpF[i] = m_PDFSegmenter->GetBinMin()[i];
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "BinMin", MET_FLOAT_ARRAY, N, tmpF );
  metaFields.push_back( mF );

  for( unsigned int i = 0; i < N; ++i )
    {
    tmpF[i] = m_PDFSegmenter->GetBinSize()[i];
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "BinSize", MET_FLOAT_ARRAY, N, tmpF );
  metaFields.push_back( mF );

  unsigned int nObjects = m_PDFSegmenter->GetNumberOfObjectIds();

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "NObjects", MET_INT, nObjects );
  metaFields.push_back( mF );

  for( unsigned int i = 0; i < nObjects; ++i )
    {
    tmpI[i] = m_PDFSegmenter->GetObjectId()[i];
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ObjectId", MET_INT_ARRAY, nObjects, tmpI );
  metaFields.push_back( mF );

  for( unsigned int i = 0; i < nObjects; ++i )
    {
    tmpF[i] = m_PDFSegmenter->GetObjectPDFWeight()[i];
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ObjectPDFWeight", MET_FLOAT_ARRAY,
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

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "HistogramSmoothingStandardDeviation",
    MET_FLOAT, m_PDFSegmenter->GetHistogramSmoothingStandardDeviation() );
  metaFields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "OutlierRejectPortion", MET_FLOAT,
    m_PDFSegmenter->GetOutlierRejectPortion() );
  metaFields.push_back( mF );

  char tmpC[4096];
  if( m_PDFSegmenter->GetDraft() )
    {
    strcpy( tmpC, "True" );
    }
  else
    {
    strcpy( tmpC, "False" );
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "Draft", MET_STRING, strlen(tmpC), tmpC );
  metaFields.push_back( mF );

  if( m_PDFSegmenter->GetReclassifyObjectLabels() )
    {
    strcpy( tmpC, "True" );
    }
  else
    {
    strcpy( tmpC, "False" );
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ReclassifyObjectLabels", MET_STRING,
    strlen(tmpC), tmpC );
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
  MET_InitWriteField( mF, "ReclassifyNotObjectLabels", MET_STRING,
    strlen(tmpC), tmpC );
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
  MET_InitWriteField( mF, "ForceClassification", MET_STRING,
    strlen(tmpC), tmpC );
  metaFields.push_back( mF );

  char filePath[255];
  MET_GetFilePath( _headerName, filePath );
  int skip = strlen( filePath );
  char shortFileName[255];
  sprintf( shortFileName, "%s", &(_headerName[skip]) );
  MET_SetFileSuffix( shortFileName, "mpd" );
  std::string fullFileName = filePath;
  fullFileName = fullFileName + shortFileName;

  std::string tmpString;
  for( unsigned int i = 0; i < nObjects; ++i )
    {
    char objectFileName[4096];
    sprintf( objectFileName, "%s.%02d.mha", shortFileName, i );
    tmpString = tmpString + objectFileName;
    if( i < nObjects-1 )
      {
      tmpString = tmpString + ",";
      }
    }

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ObjectPDFFile", MET_STRING,
    tmpString.size(), tmpString.c_str() );
  metaFields.push_back( mF );

  METAIO_STREAM::ofstream writeStream;

  writeStream.open( fullFileName.c_str(), METAIO_STREAM::ios::binary |
    METAIO_STREAM::ios::out );

  writeStream.precision( 10 );

  if(!MET_Write(writeStream, & metaFields))
    {
    METAIO_STREAM::cerr << "MetaObject: Write: MET_Write Failed"
                        << METAIO_STREAM::endl;
    return false;
    }

  for( unsigned int i = 0; i < nObjects; ++i )
    {
    MetaClassPDF pdfClassWriter( N,
      m_PDFSegmenter->GetNumberOfBinsPerFeature(),
      m_PDFSegmenter->GetBinMin(),
      m_PDFSegmenter->GetBinSize(),
      m_PDFSegmenter->GetClassPDFImage( i )->GetPixelContainer()
        ->GetBufferPointer() );

    pdfClassWriter.SetObjectId( m_PDFSegmenter->GetObjectId() );
    pdfClassWriter.SetObjectPDFWeight(
      m_PDFSegmenter->GetObjectPDFWeight() );
    pdfClassWriter.SetVoidId( m_PDFSegmenter->GetVoidId() );
    pdfClassWriter.SetErodeRadius( m_PDFSegmenter->GetErodeRadius() );
    pdfClassWriter.SetHoleFillIterations(
      m_PDFSegmenter->GetHoleFillIterations() );
    pdfClassWriter.SetHistogramSmoothingStandardDeviation(
      m_PDFSegmenter->GetHistogramSmoothingStandardDeviation() );
    pdfClassWriter.SetProbabilityImageSmoothingStandardDeviation(
      m_PDFSegmenter->GetProbabilityImageSmoothingStandardDeviation() );
    pdfClassWriter.SetOutlierRejectPortion(
      m_PDFSegmenter->GetOutlierRejectPortion() );
    pdfClassWriter.SetReclassifyObjectLabels(
      m_PDFSegmenter->GetReclassifyObjectLabels() );
    pdfClassWriter.SetReclassifyNotObjectLabels(
      m_PDFSegmenter->GetReclassifyNotObjectLabels() );
    pdfClassWriter.SetForceClassification(
      m_PDFSegmenter->GetForceClassification() );
    pdfClassWriter.SetDraft( m_PDFSegmenter->GetDraft() );

    char objectFileName[4096];
    sprintf( objectFileName, "%s.%02d.mha", fullFileName.c_str(), i );

    pdfClassWriter.Write( objectFileName );
    }

  return true;
}

} // End namespace tube

} // End namespace itk

#endif // __itktubePDFSegmenterIO_hxx
