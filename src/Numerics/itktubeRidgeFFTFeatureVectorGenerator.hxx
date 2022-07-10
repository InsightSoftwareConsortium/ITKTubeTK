/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef __itktubeRidgeFFTFeatureVectorGenerator_hxx
#define __itktubeRidgeFFTFeatureVectorGenerator_hxx

#include "itktubeRidgeFFTFilter.h"
#include "tubeMatrixMath.h"

#include <itkImage.h>
#include <itkSubtractImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkStatisticsImageFilter.h>

#include <limits>

namespace itk
{

namespace tube
{

template< class TImage >
RidgeFFTFeatureVectorGenerator< TImage >
::RidgeFFTFeatureVectorGenerator( void )
{
  m_UseIntensityOnly = false;
  m_Scales.resize( 0 );
  m_FeatureImageList.resize( 0 );
}

template< class TImage >
RidgeFFTFeatureVectorGenerator< TImage >
::~RidgeFFTFeatureVectorGenerator( void )
{
}

template< class TImage >
unsigned int
RidgeFFTFeatureVectorGenerator< TImage >
::GetNumberOfImageFeaturesPerScale( void ) const
{
  if( m_UseIntensityOnly )
    {
    return 2;
    }
  else
    {
    return 5;
    }
}

template< class TImage >
unsigned int
RidgeFFTFeatureVectorGenerator< TImage >
::GetNumberOfImageFeatures( void ) const
{
  unsigned int featuresPerImage =
    m_Scales.size() * this->GetNumberOfImageFeaturesPerScale() +
    this->GetNumberOfImageFeaturesPerScale() + 1;

  unsigned int numFeatures = this->GetNumberOfInputImages() * featuresPerImage;

  return numFeatures;
}

template< class TImage >
unsigned int
RidgeFFTFeatureVectorGenerator< TImage >
::GetNumberOfMathFeatures( void ) const
{
  unsigned int imageFeatures = this->GetNumberOfImageFeatures();

  unsigned int mathFeatures = 0;
  if( this->m_UseFeatureAddition )
    {
    mathFeatures += (imageFeatures*(imageFeatures-1))/2;
    }
  if( this->m_UseFeatureSubtraction )
    {
    mathFeatures += (imageFeatures*(imageFeatures-1))/2;
    }
  if( this->m_UseFeatureMultiplication )
    {
    mathFeatures += (imageFeatures*(imageFeatures-1))/2;
    }
  if( this->m_UseFeatureRatio )
    {
    mathFeatures += (imageFeatures*(imageFeatures-1))/2;
    }

  unsigned int numFeatures = this->GetNumberOfInputImages() * mathFeatures;

  return numFeatures;
}

template< class TImage >
unsigned int
RidgeFFTFeatureVectorGenerator< TImage >
::GetNumberOfFeatures( void ) const
{
  unsigned int numFeatures = this->GetNumberOfImageFeatures() +
    this->GetNumberOfMathFeatures();

  return numFeatures;
}


template< class TImage >
void
RidgeFFTFeatureVectorGenerator< TImage >
::UpdateWhitenStatistics( void )
{
  unsigned int numFeatureImages = m_FeatureImageList.size();

  this->m_WhitenMean.resize(numFeatureImages);
  this->m_WhitenStdDev.resize(numFeatureImages);
  typedef itk::StatisticsImageFilter<FeatureImageType> StatsFilterType;
  for( unsigned int f=0; f<numFeatureImages; ++f )
    {
    typename StatsFilterType::Pointer stats = StatsFilterType::New();
    stats->SetInput( m_FeatureImageList[f] );
    stats->Update();
    this->m_WhitenMean[f] = stats->GetMean();
    this->m_WhitenStdDev[f] = stats->GetSigma();
    }
}

template< class TImage >
void
RidgeFFTFeatureVectorGenerator< TImage >
::Update( void )
{

  unsigned int numFeatures = this->GetNumberOfFeatures();
  unsigned int numImageFeatures = this->GetNumberOfImageFeatures();

  typename FeatureImageType::RegionType region =
    this->m_InputImageList[0]->GetLargestPossibleRegion();
  m_FeatureImageList.resize( numImageFeatures );
  for( unsigned int feat=0; feat < numImageFeatures; ++feat )
    {
    m_FeatureImageList[feat] = FeatureImageType::New();
    m_FeatureImageList[feat]->CopyInformation( this->m_InputImageList[0] );
    m_FeatureImageList[feat]->SetRegions( region );
    m_FeatureImageList[feat]->Allocate();
    }

  unsigned int numFeaturesPerScale = this->GetNumberOfImageFeaturesPerScale();
  unsigned int optScaleFeat = 1;
  unsigned int feat = 0;
  for( unsigned int img=0; img<this->m_InputImageList.size(); ++img )
    {
    unsigned int imageFirstFeat = feat;
    if( m_UseIntensityOnly )
      {
      typedef itk::DiscreteGaussianImageFilter<TImage, FeatureImageType>
        BlurFilterType;
      for( unsigned int s=0; s<m_Scales.size(); ++s )
        {
        typename BlurFilterType::Pointer blurFilter = BlurFilterType::New();
        blurFilter->SetInput( this->m_InputImageList[img] );
        blurFilter->SetSigma(m_Scales[s]);
        blurFilter->SetUseImageSpacing(true);
        blurFilter->Update();
        m_FeatureImageList[feat++] = blurFilter->GetOutput();
        if( s == 0 )
          {
          typedef itk::SubtractImageFilter< FeatureImageType, TImage >
            DiffFilterType;
          typename DiffFilterType::Pointer diffFilter = DiffFilterType::New();
          diffFilter->SetInput1(m_FeatureImageList[feat-1]);
          diffFilter->SetInput2(this->m_InputImageList[img]);
          diffFilter->Update();
          m_FeatureImageList[feat++] = diffFilter->GetOutput();
          }
        else
          {
          typedef itk::SubtractImageFilter< FeatureImageType, FeatureImageType >
            DiffFilterType;
          typename DiffFilterType::Pointer diffFilter = DiffFilterType::New();
          diffFilter->SetInput1(m_FeatureImageList[feat-1]);
          diffFilter->SetInput2(m_FeatureImageList[feat-numFeaturesPerScale-1]);
          diffFilter->Update();
          m_FeatureImageList[feat++] = diffFilter->GetOutput();
          }
        }
      }
    else
      {
      typedef RidgeFFTFilter< TImage > RidgeFilterType;
      typename RidgeFilterType::Pointer ridgeF = RidgeFilterType::New();
      ridgeF->SetInput( this->m_InputImageList[img] );
      ridgeF->SetUseIntensityOnly( false );
      for( unsigned int s=0; s<m_Scales.size(); ++s )
        {
        ridgeF->SetScale( m_Scales[s] );
        ridgeF->Update();
        m_FeatureImageList[feat++] = ridgeF->GetIntensity();
        m_FeatureImageList[feat++] = ridgeF->GetRidgeness();
        m_FeatureImageList[feat++] = ridgeF->GetRoundness();
        m_FeatureImageList[feat++] = ridgeF->GetCurvature();
        m_FeatureImageList[feat++] = ridgeF->GetLevelness();
        }
      }

    typedef ImageRegionIterator< FeatureImageType >  IterType;
    unsigned int numIterators = numFeaturesPerScale * m_Scales.size() +
      numFeaturesPerScale + 1;
    std::vector< IterType > iterF( numIterators );
    for( unsigned int f=0; f<numIterators; ++f )
      {
      iterF[f] = IterType( m_FeatureImageList[f+imageFirstFeat], region );
      }
    unsigned int optScaleDest = numFeaturesPerScale * m_Scales.size();
    unsigned int optFeatDest = optScaleDest + 1;
    while( !iterF[0].IsAtEnd() )
      {
      // Init
      double optScaleV = iterF[optScaleFeat].Get();
      iterF[ optScaleDest ].Set( m_Scales[ 0 ] );
      for( unsigned int f=0; f<numFeaturesPerScale; ++f )
        {
        iterF[ optFeatDest + f ].Set( iterF[ f ].Get() );
        }
      // Calc max
      for( unsigned int s=1; s<m_Scales.size(); ++s )
        {
        unsigned int scaleFirstFeat = s * numFeaturesPerScale;
        if( iterF[scaleFirstFeat + optScaleFeat].Get() > optScaleV )
          {
          optScaleV = iterF[scaleFirstFeat + optScaleFeat].Get();
          iterF[ optScaleDest ].Set( m_Scales[ s ] );
          for( unsigned int f=0; f<numFeaturesPerScale; ++f )
            {
            iterF[ optFeatDest + f ].Set( iterF[ scaleFirstFeat + f ].Get() );
            }
          }
        }
      for( unsigned int f=0; f<numIterators; ++f )
        {
        ++iterF[ f ];
        }
      }
    }

  if( this->GetUpdateWhitenStatisticsOnUpdate() )
    {
    this->UpdateWhitenStatistics();
    }

}

template< class TImage >
typename RidgeFFTFeatureVectorGenerator< TImage >::FeatureVectorType
RidgeFFTFeatureVectorGenerator< TImage >
::GetFeatureVector( const IndexType & indx ) const
{
  unsigned int imageFeatureCount = this->GetNumberOfImageFeatures();

  unsigned int numFeatures = this->GetNumberOfFeatures();

  FeatureVectorType featureVector( numFeatures );

  for( unsigned int f=0; f<imageFeatureCount; ++f )
    {
    featureVector[f] = m_FeatureImageList[f]->GetPixel( indx );
    }

  if( this->m_WhitenStdDev.size() > 0 )
    {
    for( unsigned int f=0; f<imageFeatureCount; ++f )
      {
      if( this->m_WhitenStdDev.size() > f &&
        this->m_WhitenStdDev[f] > 0 )
        {
        featureVector[f] = (featureVector[f]-this->m_WhitenMean[f])
          / this->m_WhitenStdDev[f];
        }
      }
    }
  unsigned int featureCount = imageFeatureCount;
  if( this->m_UseFeatureAddition )
    {
    for( unsigned int f0 = 0; f0 < imageFeatureCount; f0++ )
      {
      for( unsigned int f1 = f0+1; f1 < imageFeatureCount; f1++ )
        {
        featureVector[featureCount++] = (featureVector[f0]+featureVector[f1])/2.0;
        }
      }
    }
  if( this->m_UseFeatureSubtraction )
    {
    for( unsigned int f0 = 0; f0 < imageFeatureCount; f0++ )
      {
      for( unsigned int f1 = f0+1; f1 < imageFeatureCount; f1++ )
        {
        featureVector[featureCount++] = (featureVector[f0]-featureVector[f1])/2.0;
        }
      }
    }
  if( this->m_UseFeatureMultiplication )
    {
    for( unsigned int f0 = 0; f0 < imageFeatureCount; f0++ )
      {
      for( unsigned int f1 = f0+1; f1 < imageFeatureCount; f1++ )
        {
        featureVector[featureCount++] = sqrt(fabs(featureVector[f0]*featureVector[f1]));
        }
      }
    }
  if( this->m_UseFeatureRatio )
    {
    for( unsigned int f0 = 0; f0 < imageFeatureCount; f0++ )
      {
      for( unsigned int f1 = f0+1; f1 < imageFeatureCount; f1++ )
        {
        if( fabs(featureVector[f0])+fabs(featureVector[f1]) != 0 )
          {
          featureVector[featureCount++] = (featureVector[f0]-featureVector[f1])
            / (fabs(featureVector[f0])+fabs(featureVector[f1]));
          }
        else
          {
          featureVector[featureCount++] = 0;
          }
        }
      }
    }

  if( this->m_WhitenStdDev.size() > imageFeatureCount )
    {
    for( unsigned int f=imageFeatureCount; f<featureCount; ++f )
      {
      if( this->m_WhitenStdDev.size() > f &&
        this->m_WhitenStdDev[f] > 0 )
        {
        featureVector[f] = (featureVector[f] - this->m_WhitenMean[f])
          / this->m_WhitenStdDev[f];
        }
      }
    }

  return featureVector;
}

template< class TImage >
typename RidgeFFTFeatureVectorGenerator< TImage >::FeatureValueType
RidgeFFTFeatureVectorGenerator< TImage >
::GetFeatureVectorValue( const IndexType & indx, unsigned int fNum ) const
{
  return this->GetFeatureVector(indx)[fNum];
}

template< class TImage >
typename RidgeFFTFeatureVectorGenerator< TImage >::FeatureImageType::Pointer
RidgeFFTFeatureVectorGenerator< TImage >
::GetFeatureImage( unsigned int i ) const
{
  if( i < this->m_FeatureImageList.size() )
    {
    return this->m_FeatureImageList[ i ];
    }
  else
    {
    return Superclass::GetFeatureImage(i);
    }
}

template< class TImage >
void
RidgeFFTFeatureVectorGenerator< TImage >
::SetScales( const RidgeScalesType & scales )
{
  this->m_Scales = scales;
}

template< class TImage >
const std::vector< double > &
RidgeFFTFeatureVectorGenerator< TImage >
::GetScales( void ) const
{
  return this->m_Scales;
}

template< class TImage >
void
RidgeFFTFeatureVectorGenerator< TImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Scales.size() = " << m_Scales.size() << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeRidgeFFTFeatureVectorGenerator_hxx )
