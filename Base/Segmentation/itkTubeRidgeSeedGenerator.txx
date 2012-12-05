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
#ifndef __itkTubeRidgeSeedGenerator_txx
#define __itkTubeRidgeSeedGenerator_txx

#include <limits>

#include "itkTubeRidgeSeedGenerator.h"

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
RidgeSeedGenerator< ImageT, LabelmapT >
::RidgeSeedGenerator()
{
  m_RidgeImage = NULL;

  m_IntensityMin = 1;
  m_IntensityMax = 0;

  m_Scales.resize( 0 );
}

template< class ImageT, class LabelmapT >
RidgeSeedGenerator< ImageT, LabelmapT >
::~RidgeSeedGenerator()
{
}

template < class ImageT, class LabelmapT >
void
RidgeSeedGenerator< ImageT, LabelmapT >
::SetRidgeImage( typename RidgeImageType::Pointer img )
{
  m_RidgeImage = img;
}

template < class ImageT, class LabelmapT >
typename ImageT::Pointer
RidgeSeedGenerator< ImageT, LabelmapT >
::GetRidgeImage( void )
{
  return m_RidgeImage;
}


template < class ImageT, class LabelmapT >
unsigned int
RidgeSeedGenerator< ImageT, LabelmapT >
::GetNumberOfFeatures( void )
{
  unsigned int numFeatures = m_Scales.size() * (ImageDimension+6);

  return numFeatures;
}

template < class ImageT, class LabelmapT >
void
RidgeSeedGenerator< ImageT, LabelmapT >
::SetIntensityRange( float intensityMin, float intensityMax )
{
  m_IntensityMin = intensityMin;
  m_IntensityMax = intensityMax;
}

template < class ImageT, class LabelmapT >
float
RidgeSeedGenerator< ImageT, LabelmapT >
::GetIntensityMin( void )
{
  return m_IntensityMin;
}

template < class ImageT, class LabelmapT >
float
RidgeSeedGenerator< ImageT, LabelmapT >
::GetIntensityMax( void )
{
  return m_IntensityMax;
}

template < class ImageT, class LabelmapT >
void
RidgeSeedGenerator< ImageT, LabelmapT >
::SetIntensityRangeByPercentile( float percentile,
  bool findBrightPoints )
{
  // HERE
}

template < class ImageT, class LabelmapT >
void
RidgeSeedGenerator< ImageT, LabelmapT >
::GenerateFeatureImages( void )
{
  unsigned int numFeatures = this->GetNumberOfFeatures();

  this->m_FeatureImageList.resize( numFeatures );

  for( unsigned int i=0; i<numFeatures; i++ )
    {
    this->m_FeatureImageList[ i ] = LDAImageType::New();
    this->m_FeatureImageList[ i ]->SetRegions(
      m_RidgeImage->GetLargestPossibleRegion() );
    this->m_FeatureImageList[ i ]->Allocate();
    this->m_FeatureImageList[ i ]->CopyInformation( m_RidgeImage );
    this->m_FeatureImageList[ i ]->FillBuffer( 0 );
    }

  typedef NJetImageFunction< RidgeImageType > NJetFunctionType;
  typename NJetFunctionType::Pointer njet = NJetFunctionType::New();
  typename NJetFunctionType::VectorType v;
  typename NJetFunctionType::MatrixType m;
  njet->SetInputImage( m_RidgeImage );
  if( this->m_Labelmap.IsNotNull() )
    {
    unsigned int numClasses = this->GetNumberOfObjectIds();

    itk::ImageRegionIteratorWithIndex< RidgeImageType > iter(
      m_RidgeImage, m_RidgeImage->GetLargestPossibleRegion() );
    typedef itk::ImageRegionConstIteratorWithIndex< MaskImageType >
      ConstMaskImageIteratorType;
    ConstMaskImageIteratorType itInMask( this->m_Labelmap,
      this->m_Labelmap->GetLargestPossibleRegion() );
    while( !iter.IsAtEnd() )
      {
      ObjectIdType val = static_cast<ObjectIdType>( itInMask.Get() );
      bool found = false;
      for( unsigned int c=0; c<numClasses; c++ )
        {
        if( val == this->m_ObjectIdList[c] )
          {
          found = true;
          break;
          }
        }
      if( found )
        {
        float tf = iter.Get();
        if( m_IntensityMin > m_IntensityMax ||
            ( tf >= m_IntensityMin && tf <= m_IntensityMax ) )
          {
          typename RidgeImageType::IndexType indx = iter.GetIndex();
          unsigned int fcount = 0;
          for( unsigned int s=0; s<m_Scales.size(); s++ )
            {
            double ridgeness = njet->RidgenessAtIndex( indx, m_Scales[s] );
            this->m_FeatureImageList[ fcount++ ]->
              SetPixel( indx, njet->GetMostRecentIntensity() );
            v = njet->GetMostRecentDerivative();
            double dMag = 0;
            for( unsigned int d=0; d<ImageDimension; d++ )
              {
              dMag += v[d]*v[d];
              }
            this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, dMag );
            m = njet->GetMostRecentHessian();
            vnl_symmetric_eigensystem< double > eigSys( m.GetVnlMatrix() );
            for( unsigned int d=0; d<ImageDimension; d++ )
              {
              this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
                eigSys.get_eigenvalue( d ) );
              }
            this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
              ridgeness );
            this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
              njet->GetMostRecentRidgeRoundness() );
            this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
              njet->GetMostRecentRidgeLevelness() );
            this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
              njet->GetMostRecentRidgeCurvature() );
            }
          }
        }
      ++iter;
      ++itInMask;
      }
    }
  else
    {
    itk::ImageRegionIteratorWithIndex< LDAImageType > iter(
      m_RidgeImage, m_RidgeImage->GetLargestPossibleRegion() );
    while( !iter.IsAtEnd() )
      {
      float tf = iter.Get();
      if( m_IntensityMin > m_IntensityMax ||
          ( tf >= m_IntensityMin && tf <= m_IntensityMax ) )
        {
        typename RidgeImageType::IndexType indx = iter.GetIndex();
        unsigned int fcount = 0;
        for( unsigned int s=0; s<m_Scales.size(); s++ )
          {
          double ridgeness = njet->RidgenessAtIndex( indx, m_Scales[s] );
          this->m_FeatureImageList[ fcount++ ]->
            SetPixel( indx, njet->GetMostRecentIntensity() );
          v = njet->GetMostRecentDerivative();
          double dMag = 0;
          for( unsigned int d=0; d<ImageDimension; d++ )
            {
            dMag += v[d]*v[d];
            }
          this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, dMag );
          m = njet->GetMostRecentHessian();
          LDAMatrixType eVects;
          LDAVectorType eVals;
          ::tube::ComputeEigen<double>( m.GetVnlMatrix().as_ref(),
             eVects, eVals, false, true );
          for( unsigned int d=0; d<ImageDimension; d++ )
            {
            this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
              eVals[d] );
            }
          this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
            ridgeness );
          this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
            njet->GetMostRecentRidgeRoundness() );
          this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
            njet->GetMostRecentRidgeLevelness() );
          this->m_FeatureImageList[ fcount++ ]->SetPixel( indx,
            njet->GetMostRecentRidgeCurvature() );
          }
        }
      ++iter;
      }
    }

}

template < class ImageT, class LabelmapT >
void
RidgeSeedGenerator< ImageT, LabelmapT >
::SetScales( const RidgeScalesType & scales )
{
  m_Scales = scales;
}

template < class ImageT, class LabelmapT >
std::vector< double > &
RidgeSeedGenerator< ImageT, LabelmapT >
::GetScales( void )
{
  return m_Scales;
}

template < class ImageT, class LabelmapT >
void
RidgeSeedGenerator< ImageT, LabelmapT >
::GenerateLDA()
{
  Superclass::GenerateLDA();
}

template <class ImageT, class LabelmapT >
void
RidgeSeedGenerator< ImageT, LabelmapT >
::Update( void )
{
  this->GenerateFeatureImages();

  unsigned int numFeatures = this->GetNumberOfFeatures();
  for( unsigned int i=0; i<numFeatures; i++ )
    {
    this->UpdateWhitenFeatureImageStats( i );
    this->WhitenFeatureImage( i );
    }

  Superclass::Update();
}

template <class ImageT, class LabelmapT >
void
RidgeSeedGenerator< ImageT, LabelmapT >
::UpdateLDAImages( void )
{
  this->GenerateFeatureImages();

  unsigned int numFeatures = this->GetNumberOfFeatures();
  for( unsigned int i=0; i<numFeatures; i++ )
    {
    this->WhitenFeatureImage( i );
    }

  Superclass::UpdateLDAImages();
}

template <class ImageT, class LabelmapT >
void
RidgeSeedGenerator< ImageT, LabelmapT >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "IntensityMin = " << m_IntensityMin << std::endl;
  os << indent << "IntensityMax = " << m_IntensityMax << std::endl;

  os << indent << "Scales.size() = " << m_Scales.size() << std::endl;
}

}

}

#endif //RidgeSeedGenerator_txx
