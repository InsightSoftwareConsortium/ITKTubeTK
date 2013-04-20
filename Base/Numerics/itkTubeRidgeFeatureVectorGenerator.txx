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

#include <limits>

#include "itkTubeRidgeFeatureVectorGenerator.h"

#include "itkProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"


#include "tubeMatrixMath.h"
#include "itkTubeNJetImageFunction.h"

namespace itk
{

namespace tube
{

template< class ImageT, class LabelmapT >
RidgeFeatureVectorGenerator< ImageT, LabelmapT >
::RidgeFeatureVectorGenerator()
{
  m_Scales.resize( 0 );
}

template< class ImageT, class LabelmapT >
RidgeFeatureVectorGenerator< ImageT, LabelmapT >
::~RidgeFeatureVectorGenerator()
{
}

template < class ImageT, class LabelmapT >
unsigned int
RidgeFeatureVectorGenerator< ImageT, LabelmapT >
::GetNumberOfFeatures( void )
{
  unsigned int numFeatures = m_Scales.size() * 5 + 9;

  return numFeatures;
}

template < class ImageT, class LabelmapT >
FeatureVectorType
RidgeFeatureVectorGenerator< ImageT, LabelmapT >
::GetFeatureVector( IndexType indx )
{
  unsigned int numFeatures = this->GetNumberOfFeatures();

  this->m_FeatureVector.resize( numFeatures );

  typedef NJetImageFunction< RidgeImageType > NJetFunctionType;
  typename NJetFunctionType::Pointer njet = NJetFunctionType::New();
  typename NJetFunctionType::VectorType v;
  typename NJetFunctionType::MatrixType m;
  njet->SetInputImage( this->m_InputImage[0] );

  ProgressReporter progress( this, 0,
    m_RidgeImage->GetLargestPossibleRegion().GetNumberOfPixels()*2, 100 );

        unsigned int fcount = 0;
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
          this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
            njet->GetMostRecentIntensity() );
          this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
            ridgeness );
          this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
            njet->GetMostRecentRidgeRoundness() );
          this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
            njet->GetMostRecentRidgeLevelness() );
          this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
            njet->GetMostRecentRidgeCurvature() );
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
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
          extremeScale );
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
          extremeIntensity );
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
          extremeRidgeness );
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
          extremeRoundness );
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
          extremeLevelness );
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
          extremeCurvature );
        typename RidgeImageType::IndexType indx2 = indx;
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
        if( this->m_FeatureImageList[ 0 ]
          ->GetLargestPossibleRegion().IsInside( indx2 ) )
          {
          intensity = njet->EvaluateAtIndex( indx2, extremeScale );
          }
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, intensity );

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
        if( this->m_FeatureImageList[ 0 ]
          ->GetLargestPossibleRegion().IsInside( indx2 ) )
          {
          intensity = njet->EvaluateAtIndex( indx2, extremeScale );
          }
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, intensity );
}

template < class ImageT, class LabelmapT >
void
RidgeFeatureVectorGenerator< ImageT, LabelmapT >
::SetScales( const RidgeScalesType & scales )
{
  m_Scales = scales;
}

template < class ImageT, class LabelmapT >
std::vector< double > &
RidgeFeatureVectorGenerator< ImageT, LabelmapT >
::GetScales( void )
{
  return m_Scales;
}

template <class ImageT, class LabelmapT >
void
RidgeFeatureVectorGenerator< ImageT, LabelmapT >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Scales.size() = " << m_Scales.size() << std::endl;
}

}

}

#endif //RidgeFeatureVectorGenerator_txx
