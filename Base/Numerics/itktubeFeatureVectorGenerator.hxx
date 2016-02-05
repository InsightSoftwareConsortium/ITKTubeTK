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

#ifndef __itktubeFeatureVectorGenerator_hxx
#define __itktubeFeatureVectorGenerator_hxx


#include "itktubeFeatureVectorGenerator.h"

#include "tubeMatrixMath.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkTimeProbesCollectorBase.h>

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

  m_WhitenFeatureImageMean.clear();
  m_WhitenFeatureImageStdDev.clear();
}

template< class TImage >
FeatureVectorGenerator< TImage >
::~FeatureVectorGenerator( void )
{
}

template< class TImage >
void
FeatureVectorGenerator< TImage >
::SetInput( typename ImageType::Pointer img )
{
  m_InputImageList.clear();
  m_WhitenFeatureImageMean.clear();
  m_WhitenFeatureImageStdDev.clear();
  m_InputImageList.push_back( img );
}

template< class TImage >
void
FeatureVectorGenerator< TImage >
::AddInput( typename ImageType::Pointer img )
{
  m_InputImageList.push_back( img );
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
::UpdateWhitenFeatureImageStats( void )
{
  const unsigned int numFeatures = m_InputImageList.size();

  m_WhitenFeatureImageMean.resize( numFeatures );
  m_WhitenFeatureImageStdDev.resize( numFeatures );
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
    delta[i] = 0;
    imMean[i] = 0;
    imStdDev[i] = 0;
    }
  unsigned int imCount = 0;

  typedef itk::ImageRegionIteratorWithIndex< TImage >
    ImageIteratorType;
  ImageIteratorType itIm( m_InputImageList[0],
    m_InputImageList[0]->GetLargestPossibleRegion() );

  IndexType indx;
  double imVal;
  while( !itIm.IsAtEnd() )
    {
    indx = itIm.GetIndex();
    ++imCount;
    for( unsigned int i = 0; i < numFeatures; i++ )
      {
      imVal = m_InputImageList[i]->GetPixel( indx );
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


template< class TImage >
void
FeatureVectorGenerator< TImage >
::SetWhitenMeans( const ValueListType & means )
{
  m_WhitenFeatureImageMean = means;
}

template< class TImage >
const typename FeatureVectorGenerator< TImage >::ValueListType &
FeatureVectorGenerator< TImage >
::GetWhitenMeans( void ) const
{
  return m_WhitenFeatureImageMean;
}

template< class TImage >
void
FeatureVectorGenerator< TImage >
::SetWhitenStdDevs( const ValueListType & stdDevs )
{
  m_WhitenFeatureImageStdDev = stdDevs;
}

template< class TImage >
const typename FeatureVectorGenerator< TImage >::ValueListType &
FeatureVectorGenerator< TImage >
::GetWhitenStdDevs( void ) const
{
  return m_WhitenFeatureImageStdDev;
}

template< class TImage >
void
FeatureVectorGenerator< TImage >
::SetWhitenFeatureImageMean( unsigned int num, double mean )
{
  if( num < m_WhitenFeatureImageMean.size() )
    {
    m_WhitenFeatureImageMean[num] = mean;
    }
}

template< class TImage >
double
FeatureVectorGenerator< TImage >
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

template< class TImage >
void
FeatureVectorGenerator< TImage >
::SetWhitenFeatureImageStdDev( unsigned int num, double stdDev )
{
  if( num < m_WhitenFeatureImageStdDev.size() )
    {
    m_WhitenFeatureImageStdDev[num] = stdDev;
    }
}

template< class TImage >
double
FeatureVectorGenerator< TImage >
::GetWhitenFeatureImageStdDev( unsigned int num ) const
{
  if( num < m_WhitenFeatureImageStdDev.size() )
    {
    return m_WhitenFeatureImageStdDev[num];
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

template< class TImage >
typename FeatureVectorGenerator< TImage >::FeatureImageType::Pointer
FeatureVectorGenerator< TImage >
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
    throw;
    return NULL;
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

#endif // End !defined(__itktubeFeatureVectorGenerator_hxx)
