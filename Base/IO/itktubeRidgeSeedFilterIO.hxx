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
#ifndef __itktubeRidgeSeedFilterIO_hxx
#define __itktubeRidgeSeedFilterIO_hxx

#include "itktubeRidgeSeedFilterIO.h"
#include "itktubePDFSegmenterIO.h"

namespace itk
{

namespace tube
{

template< class TImage, class TLabelMap >
RidgeSeedFilterIO< TImage, TLabelMap >::
RidgeSeedFilterIO( void )
{
  Clear();
}

template< class TImage, class TLabelMap >
RidgeSeedFilterIO< TImage, TLabelMap >::
RidgeSeedFilterIO( const char * _headerName )
{
  Clear();

  RidgeSeedFilterIO::Read( _headerName );
}

template< class TImage, class TLabelMap >
RidgeSeedFilterIO< TImage, TLabelMap >::
RidgeSeedFilterIO( const typename
  RidgeSeedFilterType::Pointer & _filter )
{
  Clear();

  InitializeEssential( _filter );
}

template< class TImage, class TLabelMap >
RidgeSeedFilterIO< TImage, TLabelMap >::
~RidgeSeedFilterIO()
{
}

template< class TImage, class TLabelMap >
void RidgeSeedFilterIO< TImage,
  TLabelMap >::
PrintInfo() const
{
  if( m_RidgeSeedFilter.IsNotNull() )
    {
    std::cout << m_RidgeSeedFilter << std::endl;
    }
  else
    {
    std::cout << "RidgeSeedFilter = NULL" << std::endl;
    }
}

template< class TImage, class TLabelMap >
void RidgeSeedFilterIO< TImage,
  TLabelMap >::
CopyInfo( const RidgeSeedFilterIOType & _filterIO )
{
  Clear();

  InitializeEssential( _filterIO.GetRidgeSeedFilter() );
}

template< class TImage, class TLabelMap >
void RidgeSeedFilterIO< TImage,
  TLabelMap >::
Clear( void )
{
  m_RidgeSeedFilter = NULL;
}


template< class TImage, class TLabelMap >
bool RidgeSeedFilterIO< TImage, TLabelMap >::
InitializeEssential( const typename
  RidgeSeedFilterType::Pointer & _filter )
{
  m_RidgeSeedFilter = _filter;

  return true;
}

template< class TImage, class TLabelMap >
void RidgeSeedFilterIO< TImage,
  TLabelMap >::
SetRidgeSeedFilter( const typename
  RidgeSeedFilterType::Pointer & _filter )
{
  m_RidgeSeedFilter = _filter;
}

template< class TImage, class TLabelMap >
const typename RidgeSeedFilter< TImage, TLabelMap >::Pointer
RidgeSeedFilterIO< TImage, TLabelMap >::
GetRidgeSeedFilter( void ) const
{
  return m_RidgeSeedFilter;
}

template< class TImage, class TLabelMap >
bool RidgeSeedFilterIO< TImage, TLabelMap >::
CanRead( const char * _headerName ) const
{
  MetaRidgeSeed seedReader;
  MetaClassPDF  pdfReader;

  std::string pdfName = _headerName;
  pdfName = pdfName + ".mpd";
  char pdfPath[255];
  MET_GetFilePath( _headerName, pdfPath );
  pdfName = pdfPath + pdfName;

  bool result = seedReader.CanRead( _headerName );
  bool result2 = pdfReader.CanRead( pdfName.c_str() );

  if( result && result2 )
    {
    return true;
    }
  else
    {
    return false;
    }
}

template< class TImage, class TLabelMap >
bool RidgeSeedFilterIO< TImage, TLabelMap >::
Read( const char * _headerName )
{
  if( m_RidgeSeedFilter.IsNull() )
    {
    m_RidgeSeedFilter = RidgeSeedFilterType::New();
    }

  MetaRidgeSeed seedReader;

  if( !seedReader.Read( _headerName ) )
    {
    m_RidgeSeedFilter = NULL;
    return false;
    }

  m_RidgeSeedFilter->SetScales( seedReader.GetRidgeSeedScales() );
  m_RidgeSeedFilter->SetUseIntensityOnly( seedReader.GetUseIntensityOnly() );
  m_RidgeSeedFilter->SetRidgeId( seedReader.GetRidgeId() );
  m_RidgeSeedFilter->SetBackgroundId( seedReader.GetBackgroundId() );
  m_RidgeSeedFilter->SetUnknownId( seedReader.GetUnknownId() );
  m_RidgeSeedFilter->SetSeedTolerance( seedReader.GetSeedTolerance() );
  m_RidgeSeedFilter->SetSkeletonize( seedReader.GetSkeletonize() );

  m_RidgeSeedFilter->SetNumberOfPCABasisToUseAsFeatures(
    seedReader.GetNumberOfPCABasisToUseAsFeatures() );
  m_RidgeSeedFilter->SetNumberOfLDABasisToUseAsFeatures(
    seedReader.GetNumberOfLDABasisToUseAsFeatures() );

  m_RidgeSeedFilter->SetBasisValues( seedReader.GetLDAValues() );
  m_RidgeSeedFilter->SetBasisMatrix( seedReader.GetLDAMatrix() );
  m_RidgeSeedFilter->SetWhitenMeans( seedReader.GetWhitenMeans() );
  m_RidgeSeedFilter->SetWhitenStdDevs( seedReader.GetWhitenStdDevs() );

  PDFSegmenterIO< Image< float, TImage::ImageDimension >, 4, TLabelMap >
    pdfReader( m_RidgeSeedFilter->GetPDFSegmenter() );
  std::string pdfFileName = seedReader.GetPDFFileName();
  char pdfPath[255];
  MET_GetFilePath( _headerName, pdfPath );
  pdfFileName = pdfPath + pdfFileName;
  if( !pdfReader.Read( pdfFileName.c_str() ) )
    {
    m_RidgeSeedFilter = NULL;
    return false;
    }

  return true;
}

template< class TImage, class TLabelMap >
bool RidgeSeedFilterIO< TImage, TLabelMap >::
Write( const char * _headerName )
{
  if( m_RidgeSeedFilter.IsNull() )
    {
    return false;
    }

  MetaRidgeSeed seedWriter;

  seedWriter.SetRidgeSeedScales( m_RidgeSeedFilter->GetScales() );
  seedWriter.SetUseIntensityOnly( m_RidgeSeedFilter->GetUseIntensityOnly() );
  seedWriter.SetRidgeId( m_RidgeSeedFilter->GetRidgeId() );
  seedWriter.SetBackgroundId( m_RidgeSeedFilter->GetBackgroundId() );
  seedWriter.SetUnknownId( m_RidgeSeedFilter->GetUnknownId() );
  seedWriter.SetSeedTolerance( m_RidgeSeedFilter->GetSeedTolerance() );
  seedWriter.SetSkeletonize( m_RidgeSeedFilter->GetSkeletonize() );

  seedWriter.SetNumberOfPCABasisToUseAsFeatures(
    m_RidgeSeedFilter->GetNumberOfPCABasisToUseAsFeatures() );
  seedWriter.SetNumberOfLDABasisToUseAsFeatures(
    m_RidgeSeedFilter->GetNumberOfLDABasisToUseAsFeatures() );

  seedWriter.SetLDAValues( m_RidgeSeedFilter->GetBasisValues() );
  seedWriter.SetLDAMatrix( m_RidgeSeedFilter->GetBasisMatrix() );
  seedWriter.SetWhitenMeans( m_RidgeSeedFilter->GetWhitenMeans() );
  seedWriter.SetWhitenStdDevs( m_RidgeSeedFilter->GetWhitenStdDevs() );

  char fileName[255];
  MET_GetFilePath( _headerName, fileName );
  int skip = strlen( fileName );
  std::string pdfName = &(_headerName[skip]);

  pdfName = pdfName + ".mpd";

  seedWriter.SetPDFFileName(  pdfName.c_str() );

  char pdfPath[255];
  MET_GetFilePath( _headerName, pdfPath );
  std::string pdfWriteName = pdfPath + pdfName;

  PDFSegmenterIO< Image< float, TImage::ImageDimension >, 4, TLabelMap >
    pdfWriter( m_RidgeSeedFilter->GetPDFSegmenter() );

  bool result = seedWriter.Write( _headerName );
  if( !pdfWriter.Write( pdfWriteName.c_str() ) )
    {
    result = false;
    }

  return result;
}

} // End namespace tube

} // End namespace itk

#endif // __itktubeRidgeSeedFilterIO_hxx
