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

#ifndef __itkTubeRidgeFeatureVectorGenerator_txx
#define __itkTubeRidgeFeatureVectorGenerator_txx

#include "itkTubeNJetImageFunction.h"
#include "itkTubeRidgeFeatureVectorGenerator.h"
#include "tubeMatrixMath.h"

#include <itkImage.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkProgressReporter.h>
#include <itkTimeProbesCollectorBase.h>

#include <limits>

namespace itk
{

namespace tube
{

template< class ImageT >
RidgeFeatureVectorGenerator< ImageT >
::RidgeFeatureVectorGenerator( void )
{
  m_Scales.resize( 0 );
}

template< class ImageT >
RidgeFeatureVectorGenerator< ImageT >
::~RidgeFeatureVectorGenerator( void )
{
}

template < class ImageT >
unsigned int
RidgeFeatureVectorGenerator< ImageT >
::GetNumberOfFeatures( void ) const
{
  const unsigned int numFeatures = m_Scales.size() * 5 + 8;

  return numFeatures;
}

template < class ImageT >
typename RidgeFeatureVectorGenerator< ImageT >::FeatureVectorType
RidgeFeatureVectorGenerator< ImageT >
::GetFeatureVector( const IndexType & indx ) const
{
  const unsigned int numFeatures = this->GetNumberOfFeatures();

  FeatureVectorType featureVector;
  featureVector.set_size( numFeatures );

  typedef NJetImageFunction< ImageType > NJetFunctionType;
  typename NJetFunctionType::Pointer njet = NJetFunctionType::New();
  typename NJetFunctionType::VectorType v;
  typename NJetFunctionType::MatrixType m;
  njet->SetInputImage( this->m_InputImageList[0] );

  unsigned int featureCount = 0;
  double extremeScale = 0;
  double extremeIntensity = 0;
  double extremeRidgeness = 0;
  double extremeRoundness = 0;
  double extremeLevelness = 0;
  double extremeCurvature = 0;
  typename NJetFunctionType::VectorType extremeTangent;
  extremeTangent.Fill( 0 );
  for( unsigned int s=0; s<m_Scales.size(); s++ )
    {
    double ridgeness = njet->RidgenessAtIndex( indx, m_Scales[s] );
    featureVector[featureCount++] = njet->GetMostRecentIntensity();
    featureVector[featureCount++] = ridgeness;
    featureVector[featureCount++] = njet->GetMostRecentRidgeRoundness();
    featureVector[featureCount++] = njet->GetMostRecentRidgeLevelness();
    featureVector[featureCount++] = njet->GetMostRecentRidgeCurvature();
    if( s == 0 || ridgeness > extremeRidgeness )
      {
      extremeScale = m_Scales[s];
      extremeIntensity = njet->GetMostRecentIntensity();
      extremeRidgeness = ridgeness;
      extremeRoundness = njet->GetMostRecentRidgeRoundness();
      extremeLevelness = njet->GetMostRecentRidgeLevelness();
      extremeCurvature = njet->GetMostRecentRidgeCurvature();
      extremeTangent = njet->GetMostRecentRidgeTangent();
      }
    }
  featureVector[featureCount++] = extremeScale;
  featureVector[featureCount++] = extremeIntensity;
  featureVector[featureCount++] = extremeRidgeness;
  featureVector[featureCount++] = extremeRoundness;
  featureVector[featureCount++] = extremeLevelness;
  featureVector[featureCount++] = extremeCurvature;
  typename ImageType::IndexType indx2 = indx;
  typename NJetFunctionType::VectorType t;
  typename NJetFunctionType::VectorType t2;
  indx2[0] = indx[0] + vnl_math_rnd( extremeScale *
    extremeTangent[0] );
  indx2[1] = indx[1] + vnl_math_rnd( extremeScale *
    extremeTangent[1] );
  if( ImageDimension > 2 )
    {
    indx2[2] = indx[2] + vnl_math_rnd( extremeScale *
      extremeTangent[2] );
    }
  double intensity = extremeIntensity;
  if( this->m_InputImageList[0]
    ->GetLargestPossibleRegion().IsInside( indx2 ) )
    {
    intensity = njet->EvaluateAtIndex( indx2, extremeScale );
    }
  featureVector[featureCount++] = intensity;

  indx2[0] = indx[0] - vnl_math_rnd( extremeScale *
    extremeTangent[0] );
  indx2[1] = indx[1] - vnl_math_rnd( extremeScale *
    extremeTangent[1] );
  if( ImageDimension > 2 )
    {
    indx2[2] = indx[2] - vnl_math_rnd( extremeScale *
      extremeTangent[2] );
    }
  intensity = extremeIntensity;
  if( this->m_InputImageList[0]
    ->GetLargestPossibleRegion().IsInside( indx2 ) )
    {
    intensity = njet->EvaluateAtIndex( indx2, extremeScale );
    }
  featureVector[featureCount++] = intensity;

  return featureVector;
}

template < class ImageT >
typename RidgeFeatureVectorGenerator< ImageT >::FeatureValueType
RidgeFeatureVectorGenerator< ImageT >
::GetFeatureVectorValue( const IndexType & indx, unsigned int fNum ) const
{
  const unsigned int numFeatures = this->GetNumberOfFeatures();

  if( fNum < numFeatures - 8 )
    {
    typedef NJetImageFunction< ImageType > NJetFunctionType;
    typename NJetFunctionType::Pointer njet = NJetFunctionType::New();
    typename NJetFunctionType::VectorType v;
    typename NJetFunctionType::MatrixType m;
    njet->SetInputImage( this->m_InputImageList[0] );

    unsigned int featureCount = 0;
    for( unsigned int s=0; s<m_Scales.size(); s++ )
      {
      double ridgeness = njet->RidgenessAtIndex( indx, m_Scales[s] );
      if( featureCount++ == fNum )
        {
        return njet->GetMostRecentIntensity();
        }
      if( featureCount++ == fNum )
        {
        return ridgeness;
        }
      if( featureCount++ == fNum )
        {
        return njet->GetMostRecentRidgeRoundness();
        }
      if( featureCount++ == fNum )
        {
        return njet->GetMostRecentRidgeLevelness();
        }
      if( featureCount++ == fNum )
        {
        return njet->GetMostRecentRidgeCurvature();
        }
      }
    }

  return this->GetFeatureVector( indx )[fNum];
}

template < class ImageT >
void
RidgeFeatureVectorGenerator< ImageT >
::SetScales( const RidgeScalesType & scales )
{
  m_Scales = scales;
}

template < class ImageT >
const std::vector< double > &
RidgeFeatureVectorGenerator< ImageT >
::GetScales( void ) const
{
  return m_Scales;
}

template <class ImageT >
void
RidgeFeatureVectorGenerator< ImageT >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Scales.size() = " << m_Scales.size() << std::endl;
}

}  // tube namespace

}  // itk namespace

#endif // RidgeFeatureVectorGenerator_txx
