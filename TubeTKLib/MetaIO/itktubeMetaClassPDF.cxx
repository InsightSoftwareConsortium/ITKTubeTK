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

#include <cstdlib>
#include <cstring>

#include "itktubeMetaClassPDF.h"

namespace itk
{

namespace tube
{

MetaClassPDF::
MetaClassPDF( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaClassPDF()" << METAIO_STREAM::endl;
    }

  Clear();
}

MetaClassPDF::
MetaClassPDF( const char * _headerName )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaClassPDF()" << METAIO_STREAM::endl;
    }

  Clear();

  MetaClassPDF::Read( _headerName );
}

MetaClassPDF::
MetaClassPDF( const MetaClassPDF & _metaPDF )
: MetaImage()
{
  if( META_DEBUG )
   {
   METAIO_STREAM::cout << "MetaClassPDF()" << METAIO_STREAM::endl;
   }

  Clear();

  CopyInfo( & _metaPDF );
}

MetaClassPDF::
MetaClassPDF(
  unsigned int _nFeatures,
  const VectorUIntType & _nBinsPerFeature,
  const VectorDoubleType & _binMin,
  const VectorDoubleType & _binSize,
  float * _elementData )
{
  if( META_DEBUG )
   {
   METAIO_STREAM::cout << "MetaClassPDF()" << METAIO_STREAM::endl;
   }

  Clear();

  InitializeEssential( _nFeatures, _nBinsPerFeature, _binMin, _binSize,
    _elementData );
}

MetaClassPDF::
MetaClassPDF( unsigned int _x,
  unsigned int _y,
  double _binMinX,
  double _binMinY,
  double _binSizeX,
  double _binSizeY,
  float * _elementData )
{
  VectorUIntType nBinsPerFeature;
  nBinsPerFeature.resize( 2 );
  nBinsPerFeature[0] = _x;
  nBinsPerFeature[1] = _y;

  VectorDoubleType binMin;
  binMin.resize( 2 );
  binMin[0] = _binMinX;
  binMin[1] = _binMinY;

  VectorDoubleType binSize;
  binSize.resize( 2 );
  binSize[0] = _binSizeX;
  binSize[1] = _binSizeY;

  Clear();

  InitializeEssential( 2, nBinsPerFeature, binMin, binSize, _elementData );
}


MetaClassPDF::
MetaClassPDF( unsigned int _x,
  unsigned int _y,
  unsigned int _z,
  double _binMinX,
  double _binMinY,
  double _binMinZ,
  double _binSizeX,
  double _binSizeY,
  double _binSizeZ,
  float * _elementData )
{
  VectorUIntType nBinsPerFeature;
  nBinsPerFeature.resize( 3 );
  nBinsPerFeature[0] = _x;
  nBinsPerFeature[1] = _y;
  nBinsPerFeature[2] = _z;

  VectorDoubleType binMin;
  binMin.resize( 3 );
  binMin[0] = _binMinX;
  binMin[1] = _binMinY;
  binMin[2] = _binMinZ;

  VectorDoubleType binSize;
  binSize.resize( 3 );
  binSize[0] = _binSizeX;
  binSize[1] = _binSizeY;
  binSize[2] = _binSizeZ;

  Clear();

  InitializeEssential( 3, nBinsPerFeature, binMin, binSize, _elementData );
}


MetaClassPDF::
~MetaClassPDF()
{
}

void MetaClassPDF::
PrintInfo( void ) const
{
  MetaImage::PrintInfo();

  std::cout << "NumberOfBinsPerFeature : ";
  if( ! m_NumberOfBinsPerFeature.empty() )
    {
    std::cout << m_NumberOfBinsPerFeature[0];
    for( int i = 1; i < MetaImage::NDims(); i++ )
      {
      std::cout << ", " << m_NumberOfBinsPerFeature[i];
      }
    }
  std::cout << std::endl;

  std::cout << "BinMin : ";
  if( ! m_BinMin.empty() )
    {
    std::cout << m_BinMin[0];
    for( int i = 1; i < MetaImage::NDims(); i++ )
      {
      std::cout << ", " << m_BinMin[i];
      }
    }
  std::cout << std::endl;

  std::cout << "BinSize : ";
  if( ! m_BinSize.empty() )
    {
    std::cout << m_BinSize[0];
    for( int i = 1; i < MetaImage::NDims(); i++ )
      {
      std::cout << ", " << m_BinSize[i];
      }
    }
  std::cout << std::endl;

  std::cout << "ObjectId : ";
  if( ! m_ObjectId.empty() )
    {
    std::cout << m_ObjectId[0];
    for( int i = 1; i < MetaImage::NDims(); i++ )
      {
      std::cout << ", " << m_ObjectId[i];
      }
    }
  std::cout << std::endl;

  std::cout << "ObjectPDFWeight : ";
  if( ! m_ObjectPDFWeight.empty() )
    {
    std::cout << m_ObjectPDFWeight[0];
    for( int i = 1; i < MetaImage::NDims(); i++ )
      {
      std::cout << ", " << m_ObjectPDFWeight[i];
      }
    }
  std::cout << std::endl;

  std::cout << "VoidId : " << m_VoidId << std::endl;
  std::cout << "ErodeRadius : " << m_ErodeRadius << std::endl;
  std::cout << "HoleFillIterations : " << m_HoleFillIterations << std::endl;
  std::cout << "ProbabilityImageSmoothingStandardDeviation : "
    << m_ProbabilityImageSmoothingStandardDeviation << std::endl;
  std::cout << "HistogramSmoothingStandardDeviation : "
    << m_HistogramSmoothingStandardDeviation << std::endl;
  std::cout << "OutlierRejectPortion : " << m_OutlierRejectPortion
    << std::endl;
  if( m_Draft )
    {
    std::cout << "Draft : True"  << std::endl;
    }
  else
    {
    std::cout << "Draft : False"  << std::endl;
    }
  if( m_ReclassifyObjectLabels )
    {
    std::cout << "ReclassifyObjectLabels : True"  << std::endl;
    }
  else
    {
    std::cout << "ReclassifyObjectLabels : False"  << std::endl;
    }
  if( m_ReclassifyNotObjectLabels )
    {
    std::cout << "ReclassifyNotObjectLabels : True"  << std::endl;
    }
  else
    {
    std::cout << "ReclassifyNotObjectLabels : False"  << std::endl;
    }
  if( m_ForceClassification )
    {
    std::cout << "ForceClassification : True"  << std::endl;
    }
  else
    {
    std::cout << "ForceClassification : False"  << std::endl;
    }
}

void MetaClassPDF::
CopyInfo( const MetaObject * _obj )
{
  Clear();

  MetaImage::CopyInfo( _obj );

  if( _obj == NULL )
    {
    return;
    }

  const MetaClassPDF * tmpPDF;
  try
    {
    tmpPDF = static_cast< const MetaClassPDF * >( _obj );
    }
  catch( ... )
    {
    return;
    }

  if( tmpPDF )
    {
    InitializeEssential( tmpPDF->GetNumberOfFeatures(),
      tmpPDF->GetNumberOfBinsPerFeature(), tmpPDF->GetBinMin(),
      tmpPDF->GetBinSize(), NULL );
 
    this->SetObjectId( tmpPDF->GetObjectId() );
    this->SetObjectPDFWeight( tmpPDF->GetObjectPDFWeight() );
    this->SetVoidId( tmpPDF->GetVoidId() );
    this->SetErodeRadius( tmpPDF->GetErodeRadius() );
    this->SetHoleFillIterations( tmpPDF->GetHoleFillIterations() );
    this->SetProbabilityImageSmoothingStandardDeviation(
      tmpPDF->GetProbabilityImageSmoothingStandardDeviation() );
    this->SetHistogramSmoothingStandardDeviation(
      tmpPDF->GetHistogramSmoothingStandardDeviation() );
    this->SetOutlierRejectPortion( tmpPDF->GetOutlierRejectPortion() );
    this->SetDraft( tmpPDF->GetDraft() );
    this->SetReclassifyObjectLabels( tmpPDF->GetReclassifyObjectLabels() );
    this->SetReclassifyNotObjectLabels(
      tmpPDF->GetReclassifyNotObjectLabels() );
    this->SetForceClassification( tmpPDF->GetForceClassification() );
    }
 
}

void MetaClassPDF::
Clear( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaClassPDF: Clear" << METAIO_STREAM::endl;
    }

  MetaImage::Clear();

  m_BinMin.clear();
  m_BinSize.clear();

  m_ObjectId.resize( 2 );
  m_ObjectId[0] = 255;
  m_ObjectId[1] = 127;
  m_ObjectPDFWeight.resize( 2 );
  m_ObjectPDFWeight[0] = 1;
  m_ObjectPDFWeight[1] = 1;
  m_VoidId = 0;

  m_ErodeRadius = 1;
  m_HoleFillIterations = 5;
  m_ProbabilityImageSmoothingStandardDeviation = 0.5;
  m_HistogramSmoothingStandardDeviation = 2;
  m_OutlierRejectPortion = 0.01;
  m_Draft = false;
  m_ReclassifyObjectLabels = false;
  m_ReclassifyNotObjectLabels = false;
  m_ForceClassification = false;
}


bool MetaClassPDF::
InitializeEssential( unsigned int _nFeatures,
  const VectorUIntType & _nBinsPerFeature,
  const VectorDoubleType & _binMin,
  const VectorDoubleType & _binSize,
  float * _elementData )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaClassPDF: Initialize" << METAIO_STREAM::endl;
    }

  m_NumberOfBinsPerFeature = _nBinsPerFeature;
  m_BinMin = _binMin;
  m_BinSize = _binSize;

  int nBins[10];
  float binSize[10];
  double minD[10];
  for( unsigned int i = 0; i < _nFeatures; ++i )
    {
    nBins[i] = _nBinsPerFeature[i];
    binSize[i] = _binSize[i];
    minD[i] = _binMin[i];
    }

  MetaImage::InitializeEssential( _nFeatures, nBins, binSize,
    MET_FLOAT, 1, ( void * )_elementData, true );
  MetaImage::Origin( minD );

  MetaImage::CompressedData( true );

  return true;
}

unsigned int MetaClassPDF::
GetNumberOfFeatures( void ) const
{
  return MetaImage::NDims();
}

void MetaClassPDF::
SetNumberOfBinsPerFeature( const VectorUIntType & _nBinsPerFeature )
{
  m_NumberOfBinsPerFeature = _nBinsPerFeature;

  int nBins[10];
  for( unsigned int i = 0; i < _nBinsPerFeature.size(); ++i )
    {
    nBins[i] = _nBinsPerFeature[i];
    }
  MetaImage::InitializeEssential( MetaImage::NDims(), nBins,
    MetaImage::ElementSpacing(), MET_FLOAT, 1,
    ( float * )( MetaImage::ElementData() ) );
}

const MetaClassPDF::VectorUIntType & MetaClassPDF::
GetNumberOfBinsPerFeature( void ) const
{
  for( int i = 0; i < MetaImage::NDims(); i++ )
    {
    m_NumberOfBinsPerFeature[i] = MetaImage::DimSize()[i];
    }
  return m_NumberOfBinsPerFeature;
}

void MetaClassPDF::
SetBinMin( const VectorDoubleType & _binMin )
{
  m_BinMin = _binMin;

  double binMinTemp[10];
  for( int i = 0; i < MetaImage::NDims(); i++ )
    {
    binMinTemp[i] = _binMin[i];
    }
  MetaImage::Origin( binMinTemp );
}

const MetaClassPDF::VectorDoubleType & MetaClassPDF::
GetBinMin( void ) const
{
  for( int i = 0; i < MetaImage::NDims(); i++ )
    {
    m_BinMin[i] = MetaImage::Origin()[i];
    }
  return m_BinMin;
}

void MetaClassPDF::
SetBinSize( const VectorDoubleType & _binSize )
{
  m_BinSize = _binSize;

  float binSizeTemp[10];
  for( int i = 0; i < MetaImage::NDims(); i++ )
    {
    binSizeTemp[i] = _binSize[i];
    }
  MetaImage::ElementSpacing( binSizeTemp );
}

const MetaClassPDF::VectorDoubleType & MetaClassPDF::
GetBinSize( void ) const
{
  for( int i = 0; i < MetaImage::NDims(); i++ )
    {
    m_BinSize[i] = MetaImage::ElementSpacing()[i];
    }
  return m_BinSize;
}

void MetaClassPDF::
SetPDF( float * _pdf )
{
  MetaImage::ElementData( ( float * )( _pdf ) );
}

float * MetaClassPDF::
GetPDF( void )
{
  return ( float * )( MetaImage::ElementData() );
}

float * MetaClassPDF::
ExportPDF( void )
{
  MetaImage::AutoFreeElementData( false );
  return ( float * )( MetaImage::ElementData() );
}

void MetaClassPDF::
SetObjectId( const VectorIntType & _objectId )
{
  m_ObjectId = _objectId;
}

const MetaClassPDF::VectorIntType & MetaClassPDF::
GetObjectId( void ) const
{
  return m_ObjectId;
}

void MetaClassPDF::
SetObjectPDFWeight( const VectorDoubleType & _objectPDFWeight )
{
  m_ObjectPDFWeight = _objectPDFWeight;
}

const MetaClassPDF::VectorDoubleType & MetaClassPDF::
GetObjectPDFWeight( void ) const
{
  return m_ObjectPDFWeight;
}

void MetaClassPDF::
SetVoidId( int _VoidId )
{
  m_VoidId = _VoidId;
}

int MetaClassPDF::
GetVoidId( void ) const
{
  return m_VoidId;
}

void MetaClassPDF::
SetErodeRadius( unsigned int _ErodeRadius )
{
  m_ErodeRadius = _ErodeRadius;
}

unsigned int MetaClassPDF::
GetErodeRadius( void ) const
{
  return m_ErodeRadius;
}

void MetaClassPDF::
SetHoleFillIterations( unsigned int _HoleFillIterations )
{
  m_HoleFillIterations = _HoleFillIterations;
}

unsigned int MetaClassPDF::
GetHoleFillIterations( void ) const
{
  return m_HoleFillIterations;
}

void MetaClassPDF::
SetProbabilityImageSmoothingStandardDeviation( double
  _ProbabilityImageSmoothingStandardDeviation )
{
  m_ProbabilityImageSmoothingStandardDeviation =
    _ProbabilityImageSmoothingStandardDeviation;
}

double MetaClassPDF::
GetProbabilityImageSmoothingStandardDeviation( void ) const
{
  return m_ProbabilityImageSmoothingStandardDeviation;
}

void MetaClassPDF::
SetHistogramSmoothingStandardDeviation( double
  _HistogramSmoothingStandardDeviation )
{
  m_HistogramSmoothingStandardDeviation =
    _HistogramSmoothingStandardDeviation;
}

double MetaClassPDF::
GetHistogramSmoothingStandardDeviation( void ) const
{
  return m_HistogramSmoothingStandardDeviation;
}

void MetaClassPDF::
SetOutlierRejectPortion( double _OutlierRejectPortion )
{
  m_OutlierRejectPortion = _OutlierRejectPortion;
}

double MetaClassPDF::
GetOutlierRejectPortion( void ) const
{
  return m_OutlierRejectPortion;
}

void MetaClassPDF::
SetDraft( bool _Draft )
{
  m_Draft = _Draft;
}

bool MetaClassPDF::
GetDraft( void ) const
{
  return m_Draft;
}

void MetaClassPDF::
SetReclassifyObjectLabels( bool _ReclassifyObjectLabels )
{
  m_ReclassifyObjectLabels = _ReclassifyObjectLabels;
}

bool MetaClassPDF::
GetReclassifyObjectLabels( void ) const
{
  return m_ReclassifyObjectLabels;
}

void MetaClassPDF::
SetReclassifyNotObjectLabels( bool _ReclassifyNotObjectLabels )
{
  m_ReclassifyNotObjectLabels = _ReclassifyNotObjectLabels;
}

bool MetaClassPDF::
GetReclassifyNotObjectLabels( void ) const
{
  return m_ReclassifyNotObjectLabels;
}

void MetaClassPDF::
SetForceClassification( bool _ForceClassification )
{
  m_ForceClassification = _ForceClassification;
}

bool MetaClassPDF::
GetForceClassification( void ) const
{
  return m_ForceClassification;
}

bool MetaClassPDF::
CanRead( const char * _headerName ) const
{
  return MetaImage::CanRead( _headerName );
}

bool MetaClassPDF::
Read( const char * _headerName )
{
  return MetaImage::Read( _headerName, true, NULL );
}

bool MetaClassPDF::
CanReadStream( METAIO_STREAM::ifstream * _stream ) const
{
  return MetaImage::CanReadStream( _stream );
}

bool MetaClassPDF::
ReadStream( METAIO_STREAM::ifstream * _stream )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaClassPDF: ReadStream"
      << METAIO_STREAM::endl;
    }

  M_Destroy();

  Clear();

  M_SetupReadFields();

  if( m_ReadStream )
    {
    METAIO_STREAM::cout << "MetaClassPDF: ReadStream: two files open?"
                        << METAIO_STREAM::endl;
    delete m_ReadStream;
    }

  m_ReadStream = _stream;

  if( !M_Read() )
    {
    METAIO_STREAM::cout << "MetaClassPDF: Read: Cannot parse file"
                        << METAIO_STREAM::endl;
    m_ReadStream = NULL;
    return false;
    }

  m_ReadStream = NULL;

  InitializeEssential( MetaImage::NDims(), m_NumberOfBinsPerFeature,
    m_BinMin, m_BinSize, ( float * )( m_ElementData ) );

  return true;
}

bool MetaClassPDF::
Write( const char * _headerName )
{
  return MetaImage::Write( _headerName );
}

bool MetaClassPDF::
WriteStream( METAIO_STREAM::ofstream * _stream )
{
  return MetaImage::WriteStream( _stream );
}

void MetaClassPDF::
M_SetupReadFields( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaClassPDF: M_SetupReadFields"
                        << METAIO_STREAM::endl;
    }

  MetaImage::Clear();

  MetaImage::M_SetupReadFields();

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "NObjects", MET_INT, true );
  m_Fields.push_back( mF );

  unsigned int nObjects = MET_GetFieldRecordNumber( "NObjects", &m_Fields );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ObjectId", MET_INT_ARRAY, true, nObjects );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ObjectPDFWeight", MET_FLOAT_ARRAY, true,
    nObjects );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "VoidId", MET_INT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ErodeRadius", MET_INT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "HoleFillIterations", MET_INT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ProbabilityImageSmoothingStandardDeviation",
    MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "HistogramSmoothingStandardDeviation",
    MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "OutlierRejectPortion", MET_FLOAT, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "Draft", MET_STRING, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ReclassifyObjectLabels", MET_STRING, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ReclassifyNotObjectLabels", MET_STRING, true );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ForceClassification", MET_STRING, true );
  m_Fields.push_back( mF );
}

void MetaClassPDF::
M_SetupWriteFields( void )
{
  MetaImage::M_SetupWriteFields();

  double binMinTemp[10];
  for( int i = 0; i < MetaImage::NDims(); i++ )
    {
    binMinTemp[i] = m_BinMin[i];
    }
  MetaImage::Origin( binMinTemp );

  float binSizeTemp[10];
  for( int i = 0; i < MetaImage::NDims(); i++ )
    {
    binSizeTemp[i] = m_BinSize[i];
    }
  MetaImage::ElementSpacing( binSizeTemp );

  MET_FieldRecordType * mF_LastField; // ElementDataFileName
  mF_LastField = m_Fields.back();
  m_Fields.erase( m_Fields.end()-1 );

  unsigned int nObjects = static_cast< unsigned int >( m_ObjectId.size() );

  MET_FieldRecordType * mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "NObjects", MET_INT, nObjects );
  m_Fields.push_back( mF );

  int tmpI[4096];
  for( unsigned int i = 0; i < nObjects; ++i )
    {
    tmpI[i] = m_ObjectId[i];
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ObjectId", MET_INT_ARRAY, nObjects, tmpI );
  m_Fields.push_back( mF );

  float tmpF[4096];
  for( unsigned int i = 0; i < nObjects; ++i )
    {
    tmpF[i] = m_ObjectPDFWeight[i];
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ObjectPDFWeight", MET_FLOAT_ARRAY,
    nObjects, tmpF );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "VoidId", MET_INT, m_VoidId );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ErodeRadius", MET_INT, m_ErodeRadius );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "HoleFillIterations", MET_INT,
    m_HoleFillIterations );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ProbabilityImageSmoothingStandardDeviation",
    MET_FLOAT, m_ProbabilityImageSmoothingStandardDeviation );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "HistogramSmoothingStandardDeviation",
    MET_FLOAT, m_HistogramSmoothingStandardDeviation );
  m_Fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "OutlierRejectPortion", MET_FLOAT,
    m_OutlierRejectPortion );
  m_Fields.push_back( mF );

  char tmpC[4096];
  if( m_Draft )
    {
    strcpy( tmpC, "True" );
    }
  else
    {
    strcpy( tmpC, "False" );
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "Draft", MET_STRING, std::strlen( tmpC ), tmpC );
  m_Fields.push_back( mF );

  if( m_ReclassifyObjectLabels )
    {
    strcpy( tmpC, "True" );
    }
  else
    {
    strcpy( tmpC, "False" );
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ReclassifyObjectLabels", MET_STRING,
    std::strlen( tmpC ), tmpC );
  m_Fields.push_back( mF );

  if( m_ReclassifyNotObjectLabels )
    {
    strcpy( tmpC, "True" );
    }
  else
    {
    strcpy( tmpC, "False" );
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ReclassifyNotObjectLabels", MET_STRING,
    std::strlen( tmpC ), tmpC );
  m_Fields.push_back( mF );

  if( m_ForceClassification )
    {
    strcpy( tmpC, "True" );
    }
  else
    {
    strcpy( tmpC, "False" );
    }
  mF = new MET_FieldRecordType;
  MET_InitWriteField( mF, "ForceClassification", MET_STRING,
    std::strlen( tmpC ), tmpC );
  m_Fields.push_back( mF );

  m_Fields.push_back( mF_LastField );
}

bool MetaClassPDF::
M_Read( void )
{
  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaClassPDF: M_Read: Loading Header"
                        << METAIO_STREAM::endl;
    }
  if( !MetaImage::M_Read() )
    {
    METAIO_STREAM::cout << "MetaClassPDF: M_Read: Error parsing file"
                        << METAIO_STREAM::endl;
    return false;
    }

  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaClassPDF: M_Read: Parsing Header"
                        << METAIO_STREAM::endl;
    }

  if( META_DEBUG )
    {
    METAIO_STREAM::cout << "MetaClassPDF: M_Read: num fields = "
      << m_Fields.size() << METAIO_STREAM::endl;
    for( unsigned int i = 0; i < m_Fields.size(); i++ )
      {
      METAIO_STREAM::cout << "  Field " << i << " = "
        << m_Fields[i]->name << METAIO_STREAM::endl;
      }
    }

  unsigned int nFeatures = static_cast< unsigned int >(
    MetaImage::NDims() );

  m_BinMin.resize( nFeatures );
  m_BinSize.resize( nFeatures );
  m_NumberOfBinsPerFeature.resize( nFeatures );
  for( unsigned int i = 0; i < nFeatures; ++i )
    {
    m_NumberOfBinsPerFeature[i] = MetaImage::DimSize()[i];
    m_BinMin[i] = MetaImage::Origin()[i];
    m_BinSize[i] = MetaImage::ElementSpacing()[i];
    }

  MET_FieldRecordType * mF = MET_GetFieldRecord( "NObjects",
    &m_Fields );
  unsigned int nObjects = ( unsigned int )mF->value[0];

  m_ObjectId.resize( nObjects );
  mF = MET_GetFieldRecord( "ObjectId", &m_Fields );
  for( unsigned int i = 0; i < nObjects; i++ )
    {
    m_ObjectId[i] = static_cast< int >( mF->value[i] );
    }

  m_ObjectPDFWeight.resize( nObjects );
  mF = MET_GetFieldRecord( "ObjectPDFWeight", &m_Fields );
  if( mF && mF->defined )
    {
    for( unsigned int i = 0; i < nObjects; i++ )
      {
      m_ObjectPDFWeight[i] = ( double )mF->value[i];
      }
    }

  mF = MET_GetFieldRecord( "VoidId", &m_Fields );
  if( mF && mF->defined )
    {
    m_VoidId = static_cast< int >( mF->value[0] );
    }

  mF = MET_GetFieldRecord( "ErodeRadius", &m_Fields );
  if( mF && mF->defined )
    {
    m_ErodeRadius = static_cast< unsigned int >( mF->value[0] );
    }

  mF = MET_GetFieldRecord( "HoleFillIterations", &m_Fields );
  if( mF && mF->defined )
    {
    m_HoleFillIterations = static_cast< int >( mF->value[0] );
    }

  mF = MET_GetFieldRecord( "ProbabilityImageSmoothingStandardDeviation",
    &m_Fields );
  if( mF && mF->defined )
    {
    m_ProbabilityImageSmoothingStandardDeviation = ( double )mF->value[0];
    }

  mF = MET_GetFieldRecord( "HistogramSmoothingStandardDeviation",
    &m_Fields );
  if( mF && mF->defined )
    {
    m_HistogramSmoothingStandardDeviation = ( double )mF->value[0];
    }

  mF = MET_GetFieldRecord( "OutlierRejectPortion", &m_Fields );
  if( mF && mF->defined )
    {
    m_OutlierRejectPortion = ( double )mF->value[0];
    }

  mF = MET_GetFieldRecord( "Draft", &m_Fields );
  if( mF && mF->defined )
    {
    m_Draft = false;
    if( ( ( char * )( mF->value ) )[0] == 'T' ||
      ( ( char * )( mF->value ) )[0] == 't' )
      {
      m_Draft = true;
      }
    }

  mF = MET_GetFieldRecord( "ReclassifyObjectLabels", &m_Fields );
  if( mF && mF->defined )
    {
    m_ReclassifyObjectLabels = false;
    if( ( ( char * )( mF->value ) )[0] == 'T' ||
      ( ( char * )( mF->value ) )[0] == 't' )
      {
      m_ReclassifyObjectLabels = true;
      }
    }

  mF = MET_GetFieldRecord( "ReclassifyNotObjectLabels", &m_Fields );
  if( mF && mF->defined )
    {
    m_ReclassifyNotObjectLabels = false;
    if( ( ( char * )( mF->value ) )[0] == 'T' ||
      ( ( char * )( mF->value ) )[0] == 't' )
      {
      m_ReclassifyNotObjectLabels = true;
      }
    }

  mF = MET_GetFieldRecord( "ForceClassification", &m_Fields );
  if( mF && mF->defined )
    {
    m_ForceClassification = false;
    if( ( ( char * )( mF->value ) )[0] == 'T' ||
      ( ( char * )( mF->value ) )[0] == 't' )
      {
      m_ForceClassification = true;
      }
    }

  return true;
}


} // End namespace tube

} // End namespace itk
