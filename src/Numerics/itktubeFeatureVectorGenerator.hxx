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

#ifndef __itktubeFeatureVectorGenerator_hxx
#define __itktubeFeatureVectorGenerator_hxx



#include "tubeMatrixMath.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <limits>
#include <iostream>

namespace itk
{

namespace tube
{

template< class TImage >
FeatureVectorGenerator< TImage >
::FeatureVectorGenerator( void )
{
  m_InputImageList.clear();

  m_UpdateWhitenStatisticsOnUpdate = false;
  m_WhitenMean.clear();
  m_WhitenStdDev.clear();
}

template< class TImage >
FeatureVectorGenerator< TImage >
::~FeatureVectorGenerator( void )
{
  m_InputImageList.clear();
}

template< class TImage >
void
FeatureVectorGenerator< TImage >
::SetInput( const ImageType * img )
{
  m_WhitenMean.clear();
  m_WhitenMean.push_back( 0 );
  m_WhitenStdDev.clear();
  m_WhitenStdDev.push_back( 1 );
  m_InputImageList.clear();
  m_InputImageList.push_back( img );
}

template< class TImage >
void
FeatureVectorGenerator< TImage >
::SetInput( unsigned int id, const ImageType * img )
{
  m_WhitenMean[id] = 0;
  m_WhitenStdDev[id] = 1;
  m_InputImageList[id] = img;
}

template< class TImage >
void
FeatureVectorGenerator< TImage >
::AddInput( const ImageType * img )
{
  m_InputImageList.push_back( img );
  m_WhitenMean.push_back( 0 );
  m_WhitenStdDev.push_back( 1 );
}

template< class TImage >
typename FeatureVectorGenerator< TImage >::ImageType::ConstPointer
FeatureVectorGenerator< TImage >
::GetInput( unsigned int imageNum )
{
  return m_InputImageList[ imageNum ];
}

template< class TImage >
unsigned int
FeatureVectorGenerator< TImage >
::GetNumberOfInputImages( void ) const
{
  return m_InputImageList.size();
}

template< class TImage >
void
FeatureVectorGenerator< TImage >
::SetUpdateWhitenStatisticsOnUpdate( bool updateOnUpdate )
{
  m_UpdateWhitenStatisticsOnUpdate = updateOnUpdate;
}

template< class TImage >
bool
FeatureVectorGenerator< TImage >
::GetUpdateWhitenStatisticsOnUpdate( void )
{
  return m_UpdateWhitenStatisticsOnUpdate;
}


template< class TImage >
void
FeatureVectorGenerator< TImage >
::SetWhitenMeans( const ValueListType & means )
{
  m_WhitenMean = means;
}

template< class TImage >
const typename FeatureVectorGenerator< TImage >::ValueListType &
FeatureVectorGenerator< TImage >
::GetWhitenMeans( void ) const
{
  return m_WhitenMean;
}

template< class TImage >
void
FeatureVectorGenerator< TImage >
::SetWhitenStdDevs( const ValueListType & stdDevs )
{
  m_WhitenStdDev = stdDevs;
}

template< class TImage >
const typename FeatureVectorGenerator< TImage >::ValueListType &
FeatureVectorGenerator< TImage >
::GetWhitenStdDevs( void ) const
{
  return m_WhitenStdDev;
}

template< class TImage >
void
FeatureVectorGenerator< TImage >
::SetWhitenMean( unsigned int num, double mean )
{
  if( num < m_WhitenMean.size() )
    {
    m_WhitenMean[num] = mean;
    }
}

template< class TImage >
double
FeatureVectorGenerator< TImage >
::GetWhitenMean( unsigned int num ) const
{
  if( num < m_WhitenMean.size() )
    {
    return m_WhitenMean[num];
    }
  else
    {
    return 0;
    }
}

template< class TImage >
void
FeatureVectorGenerator< TImage >
::SetWhitenStdDev( unsigned int num, double stdDev )
{
  if( num < m_WhitenStdDev.size() )
    {
    m_WhitenStdDev[num] = stdDev;
    }
}

template< class TImage >
double
FeatureVectorGenerator< TImage >
::GetWhitenStdDev( unsigned int num ) const
{
  if( num < m_WhitenStdDev.size() )
    {
    return m_WhitenStdDev[num];
    }

  return 1;
}

template< class TImage >
unsigned int
FeatureVectorGenerator< TImage >
::GetNumberOfFeatures( void ) const
{
  return m_InputImageList.size();
}

template< class TImage >
typename FeatureVectorGenerator< TImage >::FeatureVectorType
FeatureVectorGenerator< TImage >
::GetFeatureVector( const IndexType & indx ) const
{
  const unsigned int numFeatures = this->GetNumberOfFeatures();

  FeatureVectorType featureVector;
  featureVector.set_size( numFeatures );

  for( unsigned int i = 0; i < numFeatures; i++ )
    {
    featureVector[i] = this->GetFeatureVectorValue( indx, i );
    }

  return featureVector;
}

template< class TImage >
typename FeatureVectorGenerator< TImage >::FeatureValueType
FeatureVectorGenerator< TImage >
::GetFeatureVectorValue( const IndexType & indx, unsigned int fNum ) const
{
  if( m_WhitenStdDev.size() > 0 &&
    m_WhitenStdDev[fNum] > 0 )
    {
    return static_cast< FeatureValueType >(
        ( m_InputImageList[fNum]->GetPixel( indx )
          - m_WhitenMean[fNum] )
        / m_WhitenStdDev[fNum] );
    }
  else
    {
    return static_cast< FeatureValueType >(
      m_InputImageList[fNum]->GetPixel( indx ) );
    }
}

template< class TImage >
typename FeatureVectorGenerator< TImage >::FeatureImageType::Pointer
FeatureVectorGenerator< TImage >
::GetFeatureImage( unsigned int featureNum ) const
{
  const unsigned int numFeatures = this->GetNumberOfFeatures();
  if( featureNum < numFeatures )
    {
    typedef itk::ImageRegionIteratorWithIndex< FeatureImageType >
      ImageIteratorType;

    typename FeatureImageType::Pointer fi;

    typename FeatureImageType::RegionType region;
    region = m_InputImageList[ 0 ]->GetLargestPossibleRegion();

    fi = FeatureImageType::New();
    fi->SetRegions( region );
    fi->CopyInformation( m_InputImageList[ 0 ] );
    fi->Allocate();

    ImageIteratorType itFeatureIm( fi, fi->GetLargestPossibleRegion() );

    FeatureValueType v;
    while( !itFeatureIm.IsAtEnd() )
      {
      IndexType indx = itFeatureIm.GetIndex();

      v = this->GetFeatureVectorValue( indx, featureNum );

      itFeatureIm.Set( v );

      ++itFeatureIm;
      }

    return fi;
    }
  else
    {
    throw itk::ExceptionObject( "Feature does not exist." );
    }
}


template< class TImage >
void
FeatureVectorGenerator< TImage >::
Update( void )
{
  if( m_UpdateWhitenStatisticsOnUpdate )
    {
    this->UpdateWhitenStatistics();
    }
}

template< class TImage >
void
FeatureVectorGenerator< TImage >::
UpdateWhitenStatistics( void )
{
  const unsigned int numFeatures = this->GetNumberOfFeatures();

  m_WhitenMean.resize( numFeatures );
  m_WhitenStdDev.resize( numFeatures );
  ValueListType delta;
  delta.resize( numFeatures );
  ValueListType imMean;
  imMean.resize( numFeatures );
  ValueListType imStdDev;
  imStdDev.resize( numFeatures );
  for( unsigned int i = 0; i < numFeatures; i++ )
    {
    m_WhitenMean[i] = 0;
    m_WhitenStdDev[i] = 1;
    delta[i] = 0;
    imMean[i] = 0;
    imStdDev[i] = 0;
    }
  unsigned int imCount = 0;

  typedef itk::ImageRegionConstIteratorWithIndex< TImage >
    ImageConstIteratorType;
  ImageConstIteratorType itIm( m_InputImageList[0],
    m_InputImageList[0]->GetLargestPossibleRegion() );

  IndexType indx;
  double imVal;
  FeatureVectorType fv;
  while( !itIm.IsAtEnd() )
    {
    indx = itIm.GetIndex();
    fv = this->GetFeatureVector( indx );
    ++imCount;
    for( unsigned int i = 0; i < numFeatures; i++ )
      {
      imVal = fv[i];
      delta[i] = imVal - imMean[i];
      imMean[i] += delta[i] / imCount;
      imStdDev[i] += delta[i] * ( imVal - imMean[i] );
      }
    ++itIm;
    }
  if( imCount > 1 )
    {
    for( unsigned int i = 0; i < numFeatures; i++ )
      {
      imStdDev[i] = std::sqrt( imStdDev[i] / ( imCount - 1 ) );
      }
    }
  else
    {
    for( unsigned int i = 0; i < numFeatures; i++ )
      {
      imStdDev[i] = 1;
      }
    }

  for( unsigned int i = 0; i < numFeatures; i++ )
    {
    m_WhitenMean[i] = imMean[i];
    m_WhitenStdDev[i] = imStdDev[i];
    }
}

template< class TImage >
void
FeatureVectorGenerator< TImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "InputImageList.size = " << m_InputImageList.size()
    << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeFeatureVectorGenerator_hxx )
