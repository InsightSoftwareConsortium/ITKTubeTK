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
#ifndef __itkTubeNJetFeatureVectorGenerator_txx
#define __itkTubeNJetFeatureVectorGenerator_txx

#include <limits>

#include "itkTubeNJetFeatureVectorGenerator.h"

#include "itkTimeProbesCollectorBase.h"

#include "itkImage.h"

#include "tubeMatrixMath.h"
#include "itkTubeNJetImageFunction.h"

namespace itk
{

namespace tube
{

template< class ImageT >
NJetFeatureVectorGenerator< ImageT >
::NJetFeatureVectorGenerator()
{
  m_ZeroScales.resize( 0 );
  m_FirstScales.resize( 0 );
  m_SecondScales.resize( 0 );
  m_RidgeScales.resize( 0 );

  m_ForceOrientationInsensitivity = true;
}

template< class ImageT >
NJetFeatureVectorGenerator< ImageT >
::~NJetFeatureVectorGenerator()
{
}

template < class ImageT >
unsigned int
NJetFeatureVectorGenerator< ImageT >
::GetNumberOfFeatures( void ) const
{
  unsigned int featuresPerImage = m_ZeroScales.size()
    + (m_FirstScales.size()*(ImageDimension+1))
    + (m_SecondScales.size()*(ImageDimension+1))
    + m_RidgeScales.size() * 4;

  unsigned int numFeatures = this->GetNumberOfInputImages()
    * featuresPerImage;

  return numFeatures;
}

template < class ImageT >
typename NJetFeatureVectorGenerator< ImageT >::FeatureVectorType
NJetFeatureVectorGenerator< ImageT >
::GetFeatureVector( const IndexType & indx ) const
{
  unsigned int numFeatures = this->GetNumberOfFeatures();

  unsigned int numInputImages = this->GetNumberOfInputImages();

  typedef NJetImageFunction< ImageType > NJetFunctionType;
  typename NJetFunctionType::Pointer njet = NJetFunctionType::New();
  typename NJetFunctionType::VectorType v;
  typename NJetFunctionType::MatrixType m;

  double val;
  FeatureVectorType fv;
  fv.set_size( numFeatures );
  unsigned int vCount = 0;
  for( unsigned int inputImageNum=0; inputImageNum<numInputImages;
    inputImageNum++ )
    {
    njet->SetInputImage( this->m_InputImageList[ inputImageNum ] );

    for( unsigned int s=0; s<m_ZeroScales.size(); s++ )
      {
      fv[ vCount++ ] = njet->EvaluateAtIndex( indx,
        this->m_ZeroScales[s] );
      }

    for( unsigned int s=0; s<m_FirstScales.size(); s++ )
      {
      val = 0;
      njet->DerivativeAtIndex( indx, m_FirstScales[s], v );
      for( unsigned int d=0; d<ImageDimension; d++ )
        {
        if( m_ForceOrientationInsensitivity && v[d] < 0 )
          {
          v[d] *= -1;
          }
        fv[ vCount++ ] = v[d];
        val += v[d]*v[d];
        }
      fv[ vCount++ ] = vcl_sqrt( val );
      }

    for( unsigned int s=0; s<m_SecondScales.size(); s++ )
      {
      njet->HessianAtIndex( indx, m_SecondScales[s], m );
      for( unsigned int d=0; d<ImageDimension; d++ )
        {
        if( m_ForceOrientationInsensitivity && m[d][d] < 0 )
          {
          m[d][d] *= -1;
          }
        fv[ vCount++ ] = m[d][d];
        val += m[d][d]*m[d][d];
        }
      fv[ vCount++ ] = vcl_sqrt( val );
      }

    for( unsigned int s=0; s<m_RidgeScales.size(); s++ )
      {
      fv[ vCount++ ] = njet->RidgenessAtIndex( indx, m_RidgeScales[s] );
      fv[ vCount++ ] = njet->GetMostRecentRidgeRoundness();
      fv[ vCount++ ] = njet->GetMostRecentRidgeCurvature();
      fv[ vCount++ ] = njet->GetMostRecentRidgeLevelness();
      }
    }
  return fv;
}

template < class ImageT >
typename NJetFeatureVectorGenerator< ImageT >::FeatureValueType
NJetFeatureVectorGenerator< ImageT >
::GetFeatureVectorValue( const IndexType & indx, unsigned int fNum ) const
{
  unsigned int numFeatures = this->GetNumberOfFeatures();

  unsigned int numInputImages = this->GetNumberOfInputImages();

  typedef NJetImageFunction< ImageType > NJetFunctionType;
  typename NJetFunctionType::Pointer njet = NJetFunctionType::New();
  typename NJetFunctionType::VectorType v;
  typename NJetFunctionType::MatrixType m;

  double val;
  unsigned int vCount = 0;
  for( unsigned int inputImageNum=0; inputImageNum<numInputImages;
    inputImageNum++ )
    {
    njet->SetInputImage( this->m_InputImageList[ inputImageNum ] );

    if( fNum < m_ZeroScales.size() )
      {
      for( unsigned int s=0; s<m_ZeroScales.size(); s++ )
        {
        if( vCount == fNum )
          {
          return njet->EvaluateAtIndex( indx, this->m_ZeroScales[s] );
          }
        vCount++;
        }
      }
    else if( fNum < m_ZeroScales.size()
      + (m_FirstScales.size() * ImageDimension) + 1 )
      {
      vCount = m_ZeroScales.size();
      for( unsigned int s=0; s<m_FirstScales.size(); s++ )
        {
        val = 0;
        njet->DerivativeAtIndex( indx, 4, v );
        for( unsigned int d=0; d<ImageDimension; d++ )
          {
          if( m_ForceOrientationInsensitivity && v[d] < 0 )
            {
            v[d] *= -1;
            }
          if( vCount == fNum )
            {
            return v[d];
            }
          vCount++;
          val += v[d]*v[d];
          }
        if( vCount == fNum )
          {
          return vcl_sqrt( val );
          }
        vCount++;
        }
      }
    else if( fNum < m_ZeroScales.size()
      + (m_FirstScales.size() * ImageDimension) + 1
      + (m_SecondScales.size() * ImageDimension) + 1 )
      {
      vCount = m_ZeroScales.size()
        + (m_FirstScales.size() * ImageDimension) + 1;
      for( unsigned int s=0; s<m_SecondScales.size(); s++ )
        {
        njet->HessianAtIndex( indx, m_SecondScales[s], m );
        for( unsigned int d=0; d<ImageDimension; d++ )
          {
          if( m_ForceOrientationInsensitivity && m[d][d] < 0 )
            {
            m[d][d] *= -1;
            }
          if( vCount == fNum )
            {
            return m[d][d];
            }
          vCount++;
          val += m[d][d]*m[d][d];
          }
        if( vCount == fNum )
          {
          return vcl_sqrt( val );
          }
        vCount++;
        }
      }
    else
      {
      vCount = m_ZeroScales.size()
        + (m_FirstScales.size() * ImageDimension) + 1
        + (m_SecondScales.size() * ImageDimension) + 1;
      for( unsigned int s=0; s<m_RidgeScales.size(); s++ )
        {
        if( vCount == fNum )
          {
          return njet->RidgenessAtIndex( indx, m_RidgeScales[s] );
          }
        vCount++;
        if( vCount == fNum )
          {
          return njet->GetMostRecentRidgeRoundness();
          }
        vCount++;
        if( vCount == fNum )
          {
          return njet->GetMostRecentRidgeCurvature();
          }
        vCount++;
        if( vCount == fNum )
          {
          return njet->GetMostRecentRidgeLevelness();
          }
        vCount++;
        }
      }
    }
}

template < class ImageT >
void
NJetFeatureVectorGenerator< ImageT >
::SetZeroScales( const NJetScalesType & scales )
{
  m_ZeroScales = scales;
}

template < class ImageT >
void
NJetFeatureVectorGenerator< ImageT >
::SetFirstScales( const NJetScalesType & scales )
{
  m_FirstScales = scales;
}

template < class ImageT >
void
NJetFeatureVectorGenerator< ImageT >
::SetSecondScales( const NJetScalesType & scales )
{
  m_SecondScales = scales;
}

template < class ImageT >
void
NJetFeatureVectorGenerator< ImageT >
::SetRidgeScales( const NJetScalesType & scales )
{
  m_RidgeScales = scales;
}

template < class ImageT >
std::vector< double > &
NJetFeatureVectorGenerator< ImageT >
::GetZeroScales( void ) const
{
  return m_ZeroScales;
}

template < class ImageT >
std::vector< double > &
NJetFeatureVectorGenerator< ImageT >
::GetFirstScales( void ) const
{
  return m_FirstScales;
}

template < class ImageT >
std::vector< double > &
NJetFeatureVectorGenerator< ImageT >
::GetSecondScales( void ) const
{
  return m_SecondScales;
}

template < class ImageT >
std::vector< double > &
NJetFeatureVectorGenerator< ImageT >
::GetRidgeScales( void ) const
{
  return m_RidgeScales;
}

template < class ImageT >
void
NJetFeatureVectorGenerator< ImageT >
::SetForceOrientationInsensitivity( bool _forceOrientationInsensitivity )
{
  m_ForceOrientationInsensitivity = _forceOrientationInsensitivity;
}

template < class ImageT >
bool
NJetFeatureVectorGenerator< ImageT >
::GetForceOrientationInsensitivity( void ) const
{
  return m_ForceOrientationInsensitivity;
}

template <class ImageT >
void
NJetFeatureVectorGenerator< ImageT >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_ForceOrientationInsensitivity )
    {
    os << indent << "ForceOrientationInsensitivity = true" << std::endl;
    }
  else
    {
    os << indent << "ForceOrientationInsensitivity = false" << std::endl;
    }

  os << indent << "ZeroScales.size() = " << m_ZeroScales.size()
    << std::endl;
  os << indent << "FirstScales.size() = " << m_FirstScales.size()
    << std::endl;
  os << indent << "SecondScales.size() = " << m_SecondScales.size()
    << std::endl;
  os << indent << "RidgeScales.size() = " << m_RidgeScales.size()
    << std::endl;
}

}

}

#endif //NJetFeatureVectorGenerator_txx
