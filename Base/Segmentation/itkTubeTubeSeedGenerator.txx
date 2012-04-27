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
#ifndef __itkTubeTubeSeedGenerator_txx
#define __itkTubeTubeSeedGenerator_txx

#include <limits>

#include "itkTubeTubeSeedGenerator.h"

#include "itkTimeProbesCollectorBase.h"

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkRecursiveGaussianImageFilter.h"


#include "tubeMatrixMath.h"
#include "itkTubeNJetImageFunction.h"

namespace itk
{

namespace tube
{

template< class ImageT, class LabelmapT >
TubeSeedGenerator< ImageT, LabelmapT >
::TubeSeedGenerator()
{
  m_IntensityMin = 0;
  m_IntensityMax = 0;

  m_Scales.resize( 0 );

  m_NJetFeatureImageList.clear();
}

template< class ImageT, class LabelmapT >
TubeSeedGenerator< ImageT, LabelmapT >
::~TubeSeedGenerator()
{
}

template < class ImageT, class LabelmapT >
unsigned int
TubeSeedGenerator< ImageT, LabelmapT >
::GetNumberOfFeatures( void )
{
  unsigned int featuresPerImage = m_Scales.size() * (ImageDimension+3);

  unsigned int numFeatures = this->GetNumberOfFeatureImages()
    * featuresPerImage;

  return numFeatures;
}

template < class ImageT, class LabelmapT >
const typename TubeSeedGenerator< ImageT, LabelmapT >::LDAImageType::Pointer &
TubeSeedGenerator< ImageT, LabelmapT >
::GetNJetFeatureImage( unsigned int num )
{
  return m_NJetFeatureImageList[ num ];
}

template < class ImageT, class LabelmapT >
void
TubeSeedGenerator< ImageT, LabelmapT >
::GenerateNJetFeatureImages( void )
{
  unsigned int numFeatures = this->GetNumberOfFeatures();
  unsigned int numFeatureImages = this->GetNumberOfFeatureImages();

  m_NJetFeatureImageList.resize( numFeatures );

  typedef NJetImageFunction< ImageType > NJetFunctionType;
  typename NJetFunctionType::Pointer njet = NJetFunctionType::New();
  typename NJetFunctionType::VectorType v;
  unsigned int vCount = 0;
  for( unsigned int i=0; i<numFeatureImages; i++ )
    {
    for( unsigned int s=0; s<m_Scales.size(); s++ )
      {
      m_NJetFeatureImageList[vCount] = LDAImageType::New();
      m_NJetFeatureImageList[vCount]->CopyInformation(
        this->GetFeatureImage(i) );
      m_NJetFeatureImageList[vCount]->SetRegions(
        this->GetFeatureImage(i)->GetLargestPossibleRegion() );
      m_NJetFeatureImageList[vCount]->Allocate();
      itk::ImageRegionIteratorWithIndex< LDAImageType > iter(
        this->GetFeatureImage(i),
        this->GetFeatureImage(i)->GetLargestPossibleRegion() );
      itk::ImageRegionIteratorWithIndex< LDAImageType > iterNJet(
        m_NJetFeatureImageList[vCount],
        m_NJetFeatureImageList[vCount]->GetLargestPossibleRegion() );
      njet->SetInputImage( this->GetFeatureImage(i) );
      while( !iter.IsAtEnd() )
        {
        float tf = iter.Get();
        if( m_IntensityMin > m_IntensityMax ||
            ( tf >= m_IntensityMin && tf <= m_IntensityMax ) )
          {
          ridgeness = njet->RidgenessAndDerivativeAtIndex( iter.GetIndex(),
            m_Scales[s], v );
          }
        ++iter;
        }
      vCount++;
      }
    }
}

template < class ImageT, class LabelmapT >
void
TubeSeedGenerator< ImageT, LabelmapT >
::SetZeroScales( const NJetScalesType & scales )
{
  m_ZeroScales = scales;
}

template < class ImageT, class LabelmapT >
void
TubeSeedGenerator< ImageT, LabelmapT >
::SetFirstScales( const NJetScalesType & scales )
{
  m_FirstScales = scales;
}

template < class ImageT, class LabelmapT >
void
TubeSeedGenerator< ImageT, LabelmapT >
::SetSecondScales( const NJetScalesType & scales )
{
  m_SecondScales = scales;
}

template < class ImageT, class LabelmapT >
void
TubeSeedGenerator< ImageT, LabelmapT >
::SetRidgeScales( const NJetScalesType & scales )
{
  m_RidgeScales = scales;
}

template < class ImageT, class LabelmapT >
std::vector< double > &
TubeSeedGenerator< ImageT, LabelmapT >
::GetZeroScales( void )
{
  return m_ZeroScales;
}

template < class ImageT, class LabelmapT >
std::vector< double > &
TubeSeedGenerator< ImageT, LabelmapT >
::GetFirstScales( void )
{
  return m_FirstScales;
}

template < class ImageT, class LabelmapT >
std::vector< double > &
TubeSeedGenerator< ImageT, LabelmapT >
::GetSecondScales( void )
{
  return m_SecondScales;
}

template < class ImageT, class LabelmapT >
std::vector< double > &
TubeSeedGenerator< ImageT, LabelmapT >
::GetRidgeScales( void )
{
  return m_RidgeScales;
}

template < class ImageT, class LabelmapT >
void
TubeSeedGenerator< ImageT, LabelmapT >
::SetForceIntensityConsistency( bool _forceIntensity )
{
  m_ForceIntensityConsistency = _forceIntensity;
}

template < class ImageT, class LabelmapT >
bool
TubeSeedGenerator< ImageT, LabelmapT >
::GetForceIntensityConsistency( void )
{
  return m_ForceIntensityConsistency;
}

template < class ImageT, class LabelmapT >
void
TubeSeedGenerator< ImageT, LabelmapT >
::SetForceOrientationInsensitivity( bool _forceOrientationInsensitivity )
{
  m_ForceOrientationInsensitivity = _forceOrientationInsensitivity;
}

template < class ImageT, class LabelmapT >
bool
TubeSeedGenerator< ImageT, LabelmapT >
::GetForceOrientationInsensitivity( void )
{
  return m_ForceOrientationInsensitivity;
}

template < class ImageT, class LabelmapT >
vnl_vector< double >
TubeSeedGenerator< ImageT, LabelmapT >
::GetFeatureVector( const ContinuousIndexType & indx )
{
  unsigned int numFeatures = this->GetNumberOfFeatures();

  m_NJetFeatureVector.set_size( numFeatures );

  typename ImageType::IndexType indxI;
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    indxI[i] = indx[i];
    }
  for( unsigned int i=0; i<numFeatures; i++ )
    {
    m_NJetFeatureVector[i] = m_NJetFeatureImageList[i]->GetPixel( indxI );
    }

  return m_NJetFeatureVector;
}

template < class ImageT, class LabelmapT >
void
TubeSeedGenerator< ImageT, LabelmapT >
::GenerateLDA()
{
  Superclass::GenerateLDA();

  if( m_ForceIntensityConsistency || m_ForceOrientationInsensitivity )
    {
    unsigned int vCount = 0;
    std::vector< int > orientationNum( this->GetNumberOfFeatures(), 0 );
    for( unsigned int i=0; i<this->GetNumberOfFeatureImages(); i++ )
      {
      orientationNum[ vCount ] = -1;
      vCount++;
      int orientationBase = vCount;
      for( unsigned int s=0; s<m_ZeroScales.size(); s++ )
        {
        orientationNum[ vCount ] = -1;
        vCount++;
        }
      for( unsigned int s=0; s<m_FirstScales.size(); s++ )
        {
        orientationBase = vCount;
        for( unsigned int d=0; d<ImageDimension; d++ )
          {
          orientationNum[ vCount ] = orientationBase;
          vCount++;
          }
        orientationNum[ vCount ] = -1;
        vCount++;
        }
      orientationBase = vCount;
      for( unsigned int s=0; s<m_SecondScales.size(); s++ )
        {
        orientationBase = vCount;
        for( unsigned int d=0; d<ImageDimension; d++ )
          {
          orientationNum[ vCount ] = orientationBase;
          vCount++;
          }
        orientationNum[ vCount ] = -1;
        vCount++;
        }
      for( unsigned int s=0; s<m_RidgeScales.size(); s++ )
        {
        orientationNum[ vCount ] = -1;
        vCount++;
        }
      }

    for( unsigned int i=0; i<this->GetNumberOfLDA(); i++ )
      {
      LDAVectorType v;

      v = this->GetLDAVector( i );

      if( m_ForceIntensityConsistency )
        {
        itk::ImageRegionIterator< LDAImageType > iterF(
          this->GetFeatureImage(0),
          this->GetFeatureImage(0)->GetLargestPossibleRegion() );
        itk::ImageRegionIterator< LDAImageType > iterM(
          this->GetLDAImage(i),
          this->GetLDAImage(i)->GetLargestPossibleRegion() );
        double fVal;
        double mVal;
        double sff = 0;
        double smm = 0;
        double sfm = 0;
        double fMean = 0;
        double mMean = 0;
        unsigned int count = 0;
        while( !iterF.IsAtEnd() )
          {
          fVal = iterF.Get();
          mVal = iterM.Get();
          fMean += fVal;
          mMean += mVal;
          ++count;
          ++iterF;
          ++iterM;
          }
        fMean /= count;
        mMean /= count;
        iterF.GoToBegin();
        iterM.GoToBegin();
        while( !iterF.IsAtEnd() )
          {
          fVal = iterF.Get() - fMean;
          mVal = iterM.Get() - mMean;
          sff += fVal * fVal;
          smm += mVal * mVal;
          sfm += fVal * mVal;
          ++iterF;
          ++iterM;
          }
        double denom = 1.0 * vcl_sqrt( sff * smm );
        double measure = 0;
        if( denom != 0 )
          {
          measure = sfm / denom;
          }
        if( measure < 0 )
          {
          for( unsigned int f=0; f<this->GetNumberOfFeatures(); f++ )
            {
            v[f] *= -1;
            }
          }
        }

      if( m_ForceOrientationInsensitivity )
        {
        for( unsigned int f=0; f<this->GetNumberOfFeatures(); f++ )
          {
          if( orientationNum[f] != -1 )
            {
            double fSumS = 0;
            int fStart = f;
            while( orientationNum[f] == fStart )
              {
              fSumS += v[f] * v[f];
              ++f;
              }
            double orientVal = vcl_sqrt( fSumS / (f-fStart) );
            if( fSumS > 0 )
              {
              f = fStart;
              while( orientationNum[f] == fStart )
                {
                v[f] = orientVal;
                ++f;
                }
              }
            --f;
            }
          }
        }
      this->SetLDAVector( i, v );
      }
    }
}

template <class ImageT, class LabelmapT >
void
TubeSeedGenerator< ImageT, LabelmapT >
::Update( void )
{
  this->GenerateNJetFeatureImages();

  Superclass::Update();
}

template <class ImageT, class LabelmapT >
void
TubeSeedGenerator< ImageT, LabelmapT >
::UpdateLDAImages( void )
{
  this->GenerateNJetFeatureImages();

  Superclass::UpdateLDAImages();
}

template <class ImageT, class LabelmapT >
void
TubeSeedGenerator< ImageT, LabelmapT >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_ForceIntensityConsistency )
    {
    os << indent << "ForceIntensityConsistency = true" << std::endl;
    }
  else
    {
    os << indent << "ForceIntensityConsistency = false" << std::endl;
    }

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

  os << indent << "NJetFeatureImageList.size() = "
    << m_NJetFeatureImageList.size() << std::endl;
}

}

}

#endif //TubeSeedGenerator_txx
