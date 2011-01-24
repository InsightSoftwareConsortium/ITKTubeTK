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
#ifndef __itkTubeNJetLDAGenerator_txx
#define __itkTubeNJetLDAGenerator_txx

#include <limits>

#include "itkTubeNJetLDAGenerator.h"

#include "itkTimeProbesCollectorBase.h"

#include "itkOrientedImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"


#include "tubeMatrixMath.h"
#include "itkNJetImageFunction.h"

namespace itk
{

namespace tube
{

template< class ImageT, class LabelmapT >
NJetLDAGenerator< ImageT, LabelmapT >
::NJetLDAGenerator()
{
  m_ZeroScales.resize( 0 );
  m_FirstScales.resize( 0 );
  m_SecondScales.resize( 0 );
  m_RidgeScales.resize( 0 );
}

template< class ImageT, class LabelmapT >
NJetLDAGenerator< ImageT, LabelmapT >
::~NJetLDAGenerator()
{
}

template < class ImageT, class LabelmapT >
unsigned int
NJetLDAGenerator< ImageT, LabelmapT >
::GetNumberOfFeatures( void )
{
  unsigned int featuresPerImage = 1 + m_ZeroScales.size()
    + (m_FirstScales.size()*(ImageDimension+1))
    + (m_SecondScales.size()*(ImageDimension+1))
    + m_RidgeScales.size();

  unsigned int numFeatures = this->GetNumberOfFeatureImages()
    * featuresPerImage;

  return numFeatures;
}

template < class ImageT, class LabelmapT >
void
NJetLDAGenerator< ImageT, LabelmapT >
::SetZeroScales( const NJetScalesType & scales )
  {
  m_ZeroScales = scales;
  }

template < class ImageT, class LabelmapT >
void
NJetLDAGenerator< ImageT, LabelmapT >
::SetFirstScales( const NJetScalesType & scales )
  {
  m_FirstScales = scales;
  }

template < class ImageT, class LabelmapT >
void
NJetLDAGenerator< ImageT, LabelmapT >
::SetSecondScales( const NJetScalesType & scales )
  {
  m_SecondScales = scales;
  }

template < class ImageT, class LabelmapT >
void
NJetLDAGenerator< ImageT, LabelmapT >
::SetRidgeScales( const NJetScalesType & scales )
  {
  m_RidgeScales = scales;
  }

template < class ImageT, class LabelmapT >
std::vector< double > &
NJetLDAGenerator< ImageT, LabelmapT >
::GetZeroScales( void )
  {
  return m_ZeroScales;
  }

template < class ImageT, class LabelmapT >
std::vector< double > &
NJetLDAGenerator< ImageT, LabelmapT >
::GetFirstScales( void )
  {
  return m_FirstScales;
  }

template < class ImageT, class LabelmapT >
std::vector< double > &
NJetLDAGenerator< ImageT, LabelmapT >
::GetSecondScales( void )
  {
  return m_SecondScales;
  }

template < class ImageT, class LabelmapT >
std::vector< double > &
NJetLDAGenerator< ImageT, LabelmapT >
::GetRidgeScales( void )
  {
  return m_RidgeScales;
  }

template < class ImageT, class LabelmapT >
vnl_vector< double >
NJetLDAGenerator< ImageT, LabelmapT >
::GetFeatureVector( const ContinuousIndexType & indx )
{
  unsigned int numFeatureImages = this->GetNumberOfFeatureImages();
  unsigned int numFeatures = this->GetNumberOfFeatures();

  LDAValuesType v( numFeatures );

  typedef NJetImageFunction< ImageType > NJetFunctionType;
  typename NJetFunctionType::Pointer njet = NJetFunctionType::New();

  unsigned int vCount = 0;
  typename ImageType::IndexType indxI;
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    indxI[i] = indx[i];
    }
  for( unsigned int i=0; i<numFeatureImages; i++ )
    {
    njet->SetInputImage( this->GetFeatureImage(i) );
    v[vCount++] = static_cast< FeatureType >( this->GetFeatureImage(i)
      ->GetPixel( indxI ) );
    for( unsigned int s=0; s<m_ZeroScales.size(); s++ )
      {
      v[vCount++] = njet->EvaluateAtContinuousIndex( indx,
        m_ZeroScales[s] );
      }
    for( unsigned int s=0; s<m_FirstScales.size(); s++ )
      {
      double t = 1;
      typename NJetFunctionType::VectorType deriv;
      njet->DerivativeAtContinuousIndex( indx, m_FirstScales[s], deriv );
      for( unsigned int d=0; d<ImageDimension; d++ )
        {
        v[vCount++] = deriv[d];
        t *= deriv[d];
        }
      v[vCount++] = t;
      }
    for( unsigned int s=0; s<m_SecondScales.size(); s++ )
      {
      double t = 1;
      typename NJetFunctionType::MatrixType h;
      njet->HessianAtContinuousIndex( indx, m_SecondScales[s], h );
      for( unsigned int d=0; d<ImageDimension; d++ )
        {
        v[vCount++] = h(d, d);
        t *= h(d, d);
        }
      v[vCount++] = t;
      }
    for( unsigned int s=0; s<m_RidgeScales.size(); s++ )
      {
      v[vCount++] = njet->RidgenessAtContinuousIndex( indx,
        m_RidgeScales[s] );
      }
    }

  return v;
}

template <class ImageT, class LabelmapT >
void
NJetLDAGenerator< ImageT, LabelmapT >
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

}

}

#endif //NJetLDAGenerator_txx
