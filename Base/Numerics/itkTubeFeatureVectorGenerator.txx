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

#ifndef __itkTubeFeatureVectorGenerator_txx
#define __itkTubeFeatureVectorGenerator_txx

#include <limits>
#include <iostream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkTimeProbesCollectorBase.h>

#include "tubeMatrixMath.h"

#include "itkTubeFeatureVectorGenerator.h"

namespace itk
{

namespace tube
{

template< class ImageT >
FeatureVectorGenerator< ImageT >
::FeatureVectorGenerator( void )
{
  m_InputImageList.clear();

  m_WhitenFeatureImageMean.clear();
  m_WhitenFeatureImageStdDev.clear();
}

template< class ImageT >
FeatureVectorGenerator< ImageT >
::~FeatureVectorGenerator( void )
{
}

template < class ImageT >
void
FeatureVectorGenerator< ImageT >
::SetInputImage( typename ImageType::Pointer img )
{
  m_InputImageList.clear();
  m_WhitenFeatureImageMean.clear();
  m_WhitenFeatureImageStdDev.clear();
  m_InputImageList.push_back( img );
}

template < class ImageT >
void
FeatureVectorGenerator< ImageT >
::AddInputImage( typename ImageType::Pointer img )
{
  m_InputImageList.push_back( img );
}

template < class ImageT >
typename ImageT::Pointer
FeatureVectorGenerator< ImageT >
::GetInputImage( unsigned int num ) const
{
  if( num < m_InputImageList.size() )
    {
    return m_InputImageList[num];
    }

  return NULL;
}

template < class ImageT >
unsigned int
FeatureVectorGenerator< ImageT >
::GetNumberOfInputImages( void ) const
{
  return m_InputImageList.size();
}

template < class ImageT >
void
FeatureVectorGenerator< ImageT >
::UpdateWhitenFeatureImageStats( void )
{
  const unsigned int numFeatures = this->GetNumberOfFeatures();

  m_WhitenFeatureImageMean.resize( numFeatures );
  m_WhitenFeatureImageStdDev.resize( numFeatures );
  ValueListType imVal;
  imVal.resize( numFeatures );
  ValueListType delta;
  delta.resize( numFeatures );
  ValueListType imMean;
  imMean.resize( numFeatures );
  ValueListType imStdDev;
  imStdDev.resize( numFeatures );
  for( unsigned int i = 0; i < numFeatures; i++ )
    {
    m_WhitenFeatureImageMean[i] = 0;
    m_WhitenFeatureImageStdDev[i] = 1;
    imVal[i] = 0;
    delta[i] = 0;
    imMean[i] = 0;
    imStdDev[i] = 0;
    }
  unsigned int imCount = 0;

  typedef itk::ImageRegionIteratorWithIndex< ImageT >
    ImageIteratorType;
  ImageIteratorType itIm( m_InputImageList[0],
    m_InputImageList[0]->GetLargestPossibleRegion() );

  IndexType indx;
  while( !itIm.IsAtEnd() )
    {
    indx = itIm.GetIndex();
    imVal = this->GetFeatureVector( indx );
    ++imCount;
    for( unsigned int i = 0; i < numFeatures; i++ )
      {
      delta[i] = imVal[i] - imMean[i];
      imMean[i] += delta[i] / imCount;
      imStdDev[i] += delta[i] * ( imVal[i] - imMean[i] );
      }
    ++itIm;
    }
  if( imCount > 1 )
    {
    for( unsigned int i = 0; i < numFeatures; i++ )
      {
      imStdDev[i] = vcl_sqrt( imStdDev[i] / ( imCount - 1 ) );
      }
    }
  else
    {
    for( unsigned int i = 0; i < numFeatures; i++ )
      {
      imStdDev[i] = 0;
      }
    }

  for( unsigned int i = 0; i < numFeatures; i++ )
    {
    m_WhitenFeatureImageMean[i] = imMean[i];
    m_WhitenFeatureImageStdDev[i] = imStdDev[i];
    }
}


template < class ImageT >
void
FeatureVectorGenerator< ImageT >
::SetWhitenMeans( const ValueListType & means )
{
  m_WhitenFeatureImageMean = means;
}

template < class ImageT >
const typename FeatureVectorGenerator< ImageT >::ValueListType &
FeatureVectorGenerator< ImageT >
::GetWhitenMeans( void ) const
{
  return m_WhitenFeatureImageMean;
}

template < class ImageT >
void
FeatureVectorGenerator< ImageT >
::SetWhitenStdDevs( const ValueListType & stdDevs )
{
  m_WhitenFeatureImageStdDev = stdDevs;
}

template < class ImageT >
const typename FeatureVectorGenerator< ImageT >::ValueListType &
FeatureVectorGenerator< ImageT >
::GetWhitenStdDevs( void ) const
{
  return m_WhitenFeatureImageStdDev;
}

template < class ImageT >
void
FeatureVectorGenerator< ImageT >
::SetWhitenFeatureImageMean( unsigned int num, double mean )
{
  if( num < m_WhitenFeatureImageMean.size() )
    {
    m_WhitenFeatureImageMean[num] = mean;
    }
}

template < class ImageT >
double
FeatureVectorGenerator< ImageT >
::GetWhitenFeatureImageMean( unsigned int num ) const
{
  if( num < m_WhitenFeatureImageMean.size() )
    {
    return m_WhitenFeatureImageMean[num];
    }
  else
    {
    return 0;
    }
}

template < class ImageT >
void
FeatureVectorGenerator< ImageT >
::SetWhitenFeatureImageStdDev( unsigned int num, double stdDev )
{
  if( num < m_WhitenFeatureImageStdDev.size() )
    {
    m_WhitenFeatureImageStdDev[num] = stdDev;
    }
}

template < class ImageT >
double
FeatureVectorGenerator< ImageT >
::GetWhitenFeatureImageStdDev( unsigned int num ) const
{
  if( num < m_WhitenFeatureImageStdDev.size() )
    {
    return m_WhitenFeatureImageStdDev[num];
    }

  return 1;
}

template < class ImageT >
unsigned int
FeatureVectorGenerator< ImageT >
::GetNumberOfFeatures( void ) const
{
  return m_InputImageList.size();
}

template < class ImageT >
typename FeatureVectorGenerator< ImageT >::FeatureVectorType
FeatureVectorGenerator< ImageT >
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

template < class ImageT >
typename FeatureVectorGenerator< ImageT >::FeatureValueType
FeatureVectorGenerator< ImageT >
::GetFeatureVectorValue( const IndexType & indx, unsigned int fNum ) const
{
  if( m_WhitenFeatureImageStdDev.size() > 0 &&
    m_WhitenFeatureImageStdDev[fNum] > 0 )
    {
    return static_cast< FeatureValueType >(
        ( m_InputImageList[fNum]->GetPixel( indx )
          - m_WhitenFeatureImageMean[fNum] )
        / m_WhitenFeatureImageStdDev[fNum] );
    }
  else
    {
    return static_cast< FeatureValueType >(
      m_InputImageList[fNum]->GetPixel( indx ) );
    }
}

template < class ImageT >
typename FeatureVectorGenerator< ImageT >::FeatureImageType::Pointer
FeatureVectorGenerator< ImageT >
::GetFeatureImage( unsigned int featureNum ) const
{
  const unsigned int numFeatures = this->GetNumberOfFeatures();
  if( featureNum < numFeatures )
    {
    itk::TimeProbesCollectorBase timeCollector;

    timeCollector.Start( "GenerateFeatureImage" );

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

    timeCollector.Stop( "GenerateFeatureImage" );
    timeCollector.Report();

    return fi;
    }
  else
    {
    return NULL;
    }
}


template <class ImageT >
void
FeatureVectorGenerator< ImageT >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "InputImageList.size = " << m_InputImageList.size()
    << std::endl;
}

} // tube namespace

} // itk namespace

#endif //FeatureVectorGenerator_txx
