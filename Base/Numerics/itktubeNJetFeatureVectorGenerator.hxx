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

#ifndef __itktubeNJetFeatureVectorGenerator_hxx
#define __itktubeNJetFeatureVectorGenerator_hxx

#include "itktubeNJetFeatureVectorGenerator.h"
#include "itktubeNJetImageFunction.h"
#include "tubeMatrixMath.h"

#include <itkImage.h>
#include <itkTimeProbesCollectorBase.h>

#include <limits>

namespace itk
{

namespace tube
{

template< class TImage >
NJetFeatureVectorGenerator< TImage >
::NJetFeatureVectorGenerator( void )
{
  m_ZeroScales.clear();
  m_FirstScales.clear();
  m_SecondScales.clear();
  m_RidgeScales.clear();
}

template< class TImage >
NJetFeatureVectorGenerator< TImage >
::~NJetFeatureVectorGenerator()
{
  m_ZeroScales.clear();
  m_FirstScales.clear();
  m_SecondScales.clear();
  m_RidgeScales.clear();
}

template< class TImage >
unsigned int
NJetFeatureVectorGenerator< TImage >
::GetNumberOfFeatures( void ) const
{
  const unsigned int featuresPerImage = m_ZeroScales.size()
    + ( m_FirstScales.size()*( ImageDimension+1 ) )
    + ( m_SecondScales.size()*( ImageDimension+1 ) )
    + m_RidgeScales.size() * 4;

  const unsigned int numFeatures = this->GetNumberOfInputImages()
    * featuresPerImage;

  return numFeatures;
}

template< class TImage >
typename NJetFeatureVectorGenerator< TImage >::FeatureVectorType
NJetFeatureVectorGenerator< TImage >
::GetFeatureVector( const IndexType & indx ) const
{
  const unsigned int numFeatures = this->GetNumberOfFeatures();

  const unsigned int numInputImages = this->GetNumberOfInputImages();

  typedef NJetImageFunction< ImageType > NJetFunctionType;
  typename NJetFunctionType::Pointer njet = NJetFunctionType::New();
  typename NJetFunctionType::VectorType v;
  typename NJetFunctionType::MatrixType m;

  double val = 0.0;
  FeatureVectorType featureVector;
  featureVector.set_size( numFeatures );
  unsigned int featureCount = 0;
  for( unsigned int inputImageNum = 0; inputImageNum < numInputImages;
    inputImageNum++ )
    {
    njet->SetInputImage( this->m_InputImageList[inputImageNum] );

    for( unsigned int s = 0; s < m_ZeroScales.size(); s++ )
      {
      featureVector[featureCount++] = njet->EvaluateAtIndex( indx,
        this->m_ZeroScales[s] );
      }

    for( unsigned int s = 0; s < m_FirstScales.size(); s++ )
      {
      val = 0.0;
      njet->DerivativeAtIndex( indx, m_FirstScales[s], v );
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        featureVector[featureCount++] = v[d];
        val += v[d] * v[d];
        }
      featureVector[featureCount++] = std::sqrt( val );
      }

    for( unsigned int s = 0; s < m_SecondScales.size(); s++ )
      {
      val = 0.0;
      njet->HessianAtIndex( indx, m_SecondScales[s], m );
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        featureVector[featureCount++] = m[d][d];
        val += m[d][d]*m[d][d];
        }
      featureVector[featureCount++] = std::sqrt( val );
      }

    for( unsigned int s = 0; s < m_RidgeScales.size(); s++ )
      {
      featureVector[featureCount++] = njet->RidgenessAtIndex( indx,
        m_RidgeScales[s] );
      featureVector[featureCount++] = njet->GetMostRecentRidgeRoundness();
      featureVector[featureCount++] = njet->GetMostRecentRidgeCurvature();
      featureVector[featureCount++] = njet->GetMostRecentRidgeLevelness();
      }
    }

  if( numFeatures != featureCount )
    {
    std::cerr << "BUG: featureCount != Expected number of features"
      << std::endl;
    }

  for( unsigned int i=0; i<numFeatures; ++i )
    {
    if( this->GetWhitenStdDev( i ) > 0 )
      {
      featureVector[i] = ( featureVector[i] - this->GetWhitenMean( i ) )
        / this->GetWhitenStdDev( i );
      }
    }

  return featureVector;
}

template< class TImage >
typename NJetFeatureVectorGenerator< TImage >::FeatureValueType
NJetFeatureVectorGenerator< TImage >
::GetFeatureVectorValue( const IndexType & indx, unsigned int fNum ) const
{
  const unsigned int numInputImages = this->GetNumberOfInputImages();

  typedef NJetImageFunction< ImageType > NJetFunctionType;
  typename NJetFunctionType::Pointer njet = NJetFunctionType::New();
  typename NJetFunctionType::VectorType v;
  typename NJetFunctionType::MatrixType m;

  double fNumMean = this->GetWhitenMean( fNum );
  double fNumStdDev = this->GetWhitenStdDev( fNum );
  if( fNumStdDev <= 0 )
    {
    fNumMean = 0;
    fNumStdDev = 1;
    }
  double val = 0.0;
  unsigned int featureCount = 0;
  for( unsigned int inputImageNum = 0; inputImageNum < numInputImages;
    inputImageNum++ )
    {
    njet->SetInputImage( this->m_InputImageList[inputImageNum] );

    if( fNum < m_ZeroScales.size() )
      {
      for( unsigned int s = 0; s < m_ZeroScales.size(); s++ )
        {
        if( featureCount == fNum )
          {
          return ( njet->EvaluateAtIndex( indx, m_ZeroScales[s] )
            - fNumMean ) / fNumStdDev;
          }
        featureCount++;
        }
      }
    else if( fNum < ( m_ZeroScales.size()
      + m_FirstScales.size() * ( ImageDimension + 1 ) ) )
      {
      featureCount = m_ZeroScales.size();
      for( unsigned int s = 0; s < m_FirstScales.size(); s++ )
        {
        val = 0.0;
        njet->DerivativeAtIndex( indx, m_FirstScales[s], v );
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          if( featureCount == fNum )
            {
            return ( v[d] - fNumMean ) / fNumStdDev;
            }
          featureCount++;
          val += v[d]*v[d];
          }
        if( featureCount == fNum )
          {
          return ( std::sqrt( val ) - fNumMean ) / fNumStdDev;
          }
        featureCount++;
        }
      }
    else if( fNum < m_ZeroScales.size()
      + m_FirstScales.size() * ( ImageDimension + 1 )
      + m_SecondScales.size() * ( ImageDimension + 1 ) )
      {
      featureCount = m_ZeroScales.size()
        + m_FirstScales.size() * ( ImageDimension + 1 );
      for( unsigned int s = 0; s < m_SecondScales.size(); s++ )
        {
        val = 0.0;
        njet->HessianAtIndex( indx, m_SecondScales[s], m );
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          if( featureCount == fNum )
            {
            return ( m[d][d] - fNumMean ) / fNumStdDev;
            }
          featureCount++;
          val += m[d][d]*m[d][d];
          }
        if( featureCount == fNum )
          {
          return ( std::sqrt( val ) - fNumMean ) / fNumStdDev;
          }
        featureCount++;
        }
      }
    else
      {
      featureCount = m_ZeroScales.size()
        + m_FirstScales.size() * ( ImageDimension + 1 )
        + m_SecondScales.size() * ( ImageDimension + 1 );
      for( unsigned int s = 0; s < m_RidgeScales.size(); s++ )
        {
        if( featureCount == fNum )
          {
          njet->RidgenessAtIndex( indx, m_RidgeScales[s] );
          return ( njet->GetMostRecentRidgeness()
            - fNumMean ) / fNumStdDev;
          }
        featureCount++;
        if( featureCount == fNum )
          {
          njet->RidgenessAtIndex( indx, m_RidgeScales[s] );
          return ( njet->GetMostRecentRidgeRoundness()
            - fNumMean ) / fNumStdDev;
          }
        featureCount++;
        if( featureCount == fNum )
          {
          njet->RidgenessAtIndex( indx, m_RidgeScales[s] );
          return ( njet->GetMostRecentRidgeCurvature()
            - fNumMean ) / fNumStdDev;
          }
        featureCount++;
        if( featureCount == fNum )
          {
          njet->RidgenessAtIndex( indx, m_RidgeScales[s] );
          return ( njet->GetMostRecentRidgeLevelness()
            - fNumMean ) / fNumStdDev;
          }
        featureCount++;
        }
      }
    }
  itkExceptionMacro( << "Requested non-existent FeatureVectorValue." );
}

template< class TImage >
void
NJetFeatureVectorGenerator< TImage >
::SetZeroScales( const NJetScalesType & scales )
{
  m_ZeroScales = scales;
}

template< class TImage >
void
NJetFeatureVectorGenerator< TImage >
::SetFirstScales( const NJetScalesType & scales )
{
  m_FirstScales = scales;
}

template< class TImage >
void
NJetFeatureVectorGenerator< TImage >
::SetSecondScales( const NJetScalesType & scales )
{
  m_SecondScales = scales;
}

template< class TImage >
void
NJetFeatureVectorGenerator< TImage >
::SetRidgeScales( const NJetScalesType & scales )
{
  m_RidgeScales = scales;
}

template< class TImage >
const std::vector< double > &
NJetFeatureVectorGenerator< TImage >
::GetZeroScales( void ) const
{
  return m_ZeroScales;
}

template< class TImage >
const std::vector< double > &
NJetFeatureVectorGenerator< TImage >
::GetFirstScales( void ) const
{
  return m_FirstScales;
}

template< class TImage >
const std::vector< double > &
NJetFeatureVectorGenerator< TImage >
::GetSecondScales( void ) const
{
  return m_SecondScales;
}

template< class TImage >
const std::vector< double > &
NJetFeatureVectorGenerator< TImage >
::GetRidgeScales( void ) const
{
  return m_RidgeScales;
}

template< class TImage >
void
NJetFeatureVectorGenerator< TImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "ZeroScales.size() = " << m_ZeroScales.size()
    << std::endl;
  os << indent << "FirstScales.size() = " << m_FirstScales.size()
    << std::endl;
  os << indent << "SecondScales.size() = " << m_SecondScales.size()
    << std::endl;
  os << indent << "RidgeScales.size() = " << m_RidgeScales.size()
    << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeNJetFeatureVectorGenerator_hxx )
