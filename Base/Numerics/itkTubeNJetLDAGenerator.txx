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
NJetLDAGenerator< ImageT, LabelmapT >
::NJetLDAGenerator()
{
  m_ZeroScales.resize( 0 );
  m_FirstScales.resize( 0 );
  m_SecondScales.resize( 0 );
  m_RidgeScales.resize( 0 );

  m_ForceIntensityConsistency = true;
  m_ForceOrientationInsensitivity = true;

  m_NJetImageList.clear();
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
  unsigned int featuresPerImage = m_ZeroScales.size()
    + (m_FirstScales.size()*(ImageDimension+1))
    + (m_SecondScales.size()*(ImageDimension+1))
    + m_RidgeScales.size();

  unsigned int numFeatures = this->GetNumberOfNJetImages()
    * featuresPerImage;

  return numFeatures;
}

template < class ImageT, class LabelmapT >
void
NJetLDAGenerator< ImageT, LabelmapT >
::SetNJetImage( typename NJetImageType::Pointer img )
{
  m_NJetImageList.clear();
  m_NJetImageList.push_back( img );
}

template < class ImageT, class LabelmapT >
void
NJetLDAGenerator< ImageT, LabelmapT >
::AddNJetImage( typename NJetImageType::Pointer img )
{
  m_NJetImageList.push_back( img );
}

template < class ImageT, class LabelmapT >
typename ImageT::Pointer
NJetLDAGenerator< ImageT, LabelmapT >
::GetNJetImage( unsigned int num )
{
  if( num < m_NJetImageList.size() )
    {
    return m_NJetImageList[num];
    }
  else
    {
    return NULL;
    }
}

template < class ImageT, class LabelmapT >
unsigned int
NJetLDAGenerator< ImageT, LabelmapT >
::GetNumberOfNJetImages( void )
{
  return m_NJetImageList.size();
}


template < class ImageT, class LabelmapT >
void
NJetLDAGenerator< ImageT, LabelmapT >
::GenerateFeatureImages( void )
{
  unsigned int numFeatures = this->GetNumberOfFeatures();

  unsigned int numNJetImages = this->GetNumberOfNJetImages();

  this->m_FeatureImageList.resize( numFeatures );

  typedef itk::RecursiveGaussianImageFilter< NJetImageType, LDAImageType >
    FirstGaussFilterType;
  typedef itk::RecursiveGaussianImageFilter< LDAImageType, LDAImageType >
    GaussFilterType;
  unsigned int vCount = 0;
  for( unsigned int njetImageNum=0; njetImageNum<numNJetImages;
    njetImageNum++ )
    {
    for( unsigned int s=0; s<m_ZeroScales.size(); s++ )
      {
      typename FirstGaussFilterType::Pointer firstFilter =
        FirstGaussFilterType::New();
      firstFilter->SetInput( m_NJetImageList[ njetImageNum ] );
      firstFilter->SetSigma( m_ZeroScales[s] );
      firstFilter->SetNormalizeAcrossScale( true );
      firstFilter->SetOrder( GaussFilterType::ZeroOrder );
      firstFilter->SetDirection( 0 );
      firstFilter->Update();
      typename LDAImageType::Pointer curImage = firstFilter->GetOutput();
      for( unsigned int d=1; d<ImageDimension; d++ )
        {
        typename GaussFilterType::Pointer filter = GaussFilterType::New();
        filter->SetInput( curImage );
        filter->SetSigma( m_ZeroScales[s] );
        filter->SetNormalizeAcrossScale( true );
        filter->SetOrder( GaussFilterType::ZeroOrder );
        filter->SetDirection( d );
        filter->Update();
        curImage = filter->GetOutput();
        }
      this->m_FeatureImageList[ vCount++ ] = curImage;
      }
    for( unsigned int s=0; s<m_FirstScales.size(); s++ )
      {
      for( unsigned int d=0; d<ImageDimension; d++ )
        {
        typename FirstGaussFilterType::Pointer firstFilter =
          FirstGaussFilterType::New();
        firstFilter->SetInput( m_NJetImageList[ njetImageNum ] );
        firstFilter->SetNormalizeAcrossScale( true );
        if( d == 0 )
          {
          firstFilter->SetSigma( m_FirstScales[s] );
          firstFilter->SetOrder( GaussFilterType::FirstOrder );
          firstFilter->SetDirection( 0 );
          }
        else
          {
          firstFilter->SetSigma( m_FirstScales[s] );
          firstFilter->SetOrder( GaussFilterType::ZeroOrder );
          firstFilter->SetDirection( 0 );
          }
        firstFilter->Update();
        typename LDAImageType::Pointer curImage = firstFilter->GetOutput();
        for( unsigned int d2=1; d2<ImageDimension; d2++ )
          {
          typename GaussFilterType::Pointer filter = GaussFilterType::New();
          filter->SetInput( curImage );
          filter->SetNormalizeAcrossScale( true );
          if( d == d2 )
            {
            filter->SetSigma( m_FirstScales[s] );
            filter->SetOrder( GaussFilterType::FirstOrder );
            filter->SetDirection( d2 );
            }
          else
            {
            filter->SetSigma( m_FirstScales[s] );
            filter->SetOrder( GaussFilterType::ZeroOrder );
            filter->SetDirection( d2 );
            }
          filter->Update();
          curImage = filter->GetOutput();
          }
        if( m_ForceOrientationInsensitivity )
          {
          itk::ImageRegionIteratorWithIndex< LDAImageType > iter( curImage,
            curImage->GetLargestPossibleRegion() );
          while( !iter.IsAtEnd() )
            {
            if( iter.Get() < 0 )
              {
              iter.Set( -iter.Get() );
              }
            ++iter;
            }
          }
        this->m_FeatureImageList[ vCount++ ] = curImage;
        }
      typename FirstGaussFilterType::Pointer firstFilter =
        FirstGaussFilterType::New();
      firstFilter->SetInput( m_NJetImageList[ njetImageNum ] );
      firstFilter->SetNormalizeAcrossScale( true );
      firstFilter->SetSigma( m_FirstScales[s] );
      firstFilter->SetOrder( GaussFilterType::FirstOrder );
      firstFilter->SetDirection( 0 );
      firstFilter->Update();
      typename LDAImageType::Pointer curImage = firstFilter->GetOutput();
      for( unsigned int d=1; d<ImageDimension; d++ )
        {
        typename GaussFilterType::Pointer filter = GaussFilterType::New();
        filter->SetInput( curImage );
        filter->SetNormalizeAcrossScale( true );
        filter->SetSigma( m_FirstScales[s] );
        filter->SetOrder( GaussFilterType::FirstOrder );
        filter->SetDirection( d );
        filter->Update();
        curImage = filter->GetOutput();
        }
      if( m_ForceOrientationInsensitivity )
        {
        itk::ImageRegionIteratorWithIndex< LDAImageType > iter( curImage,
          curImage->GetLargestPossibleRegion() );
        while( !iter.IsAtEnd() )
          {
          if( iter.Get() < 0 )
            {
            iter.Set( -iter.Get() );
            }
          ++iter;
          }
        }
      this->m_FeatureImageList[ vCount++ ] = curImage;
      }
    for( unsigned int s=0; s<m_SecondScales.size(); s++ )
      {
      for( unsigned int d=0; d<ImageDimension; d++ )
        {
        typename FirstGaussFilterType::Pointer firstFilter =
          FirstGaussFilterType::New();
        firstFilter->SetInput( m_NJetImageList[ njetImageNum ] );
        firstFilter->SetNormalizeAcrossScale( true );
        if( d == 0 )
          {
          firstFilter->SetSigma( m_SecondScales[s] );
          firstFilter->SetOrder( GaussFilterType::SecondOrder );
          firstFilter->SetDirection( 0 );
          }
        else
          {
          firstFilter->SetSigma( m_SecondScales[s] );
          firstFilter->SetOrder( GaussFilterType::ZeroOrder );
          firstFilter->SetDirection( 0 );
          }
        firstFilter->Update();
        typename LDAImageType::Pointer curImage = firstFilter->GetOutput();
        for( unsigned int d2=1; d2<ImageDimension; d2++ )
          {
          typename GaussFilterType::Pointer filter = GaussFilterType::New();
          filter->SetInput( curImage );
          filter->SetNormalizeAcrossScale( true );
          if( d == d2 )
            {
            filter->SetSigma( m_SecondScales[s] );
            filter->SetOrder( GaussFilterType::SecondOrder );
            filter->SetDirection( d2 );
            }
          else
            {
            filter->SetSigma( m_SecondScales[s] );
            filter->SetOrder( GaussFilterType::ZeroOrder );
            filter->SetDirection( d2 );
            }
          filter->Update();
          curImage = filter->GetOutput();
          }
        this->m_FeatureImageList[ vCount++ ] = curImage;
        }
      typename FirstGaussFilterType::Pointer firstFilter =
        FirstGaussFilterType::New();
      firstFilter->SetInput( m_NJetImageList[ njetImageNum ] );
      firstFilter->SetNormalizeAcrossScale( true );
      firstFilter->SetSigma( m_SecondScales[s] );
      firstFilter->SetOrder( GaussFilterType::SecondOrder );
      firstFilter->SetDirection( 0 );
      firstFilter->Update();
      typename LDAImageType::Pointer curImage = firstFilter->GetOutput();
      for( unsigned int d=1; d<ImageDimension; d++ )
        {
        typename GaussFilterType::Pointer filter = GaussFilterType::New();
        filter->SetInput( curImage );
        filter->SetNormalizeAcrossScale( true );
        filter->SetSigma( m_SecondScales[s] );
        filter->SetOrder( GaussFilterType::SecondOrder );
        filter->SetDirection( d );
        filter->Update();
        curImage = filter->GetOutput();
        }
      this->m_FeatureImageList[ vCount++ ] = curImage;
      }

    typedef NJetImageFunction< NJetImageType > NJetFunctionType;
    typename NJetFunctionType::Pointer njet = NJetFunctionType::New();

    for( unsigned int s=0; s<m_RidgeScales.size(); s++ )
      {
      this->m_FeatureImageList[ vCount ] = LDAImageType::New();
      this->m_FeatureImageList[ vCount ]->SetRegions(
        m_NJetImageList[ njetImageNum ]->GetLargestPossibleRegion() );
      this->m_FeatureImageList[ vCount ]->Allocate();
      this->m_FeatureImageList[ vCount ]->CopyInformation(
        m_NJetImageList[ njetImageNum ] );
      itk::ImageRegionIteratorWithIndex< LDAImageType > iter(
        this->m_FeatureImageList[ vCount ],
        this->m_FeatureImageList[ vCount ]->GetLargestPossibleRegion() );
      njet->SetInputImage( m_NJetImageList[ njetImageNum ] );
      while( !iter.IsAtEnd() )
        {
        iter.Set( njet->RidgenessAtIndex( iter.GetIndex(),
            m_RidgeScales[s] ) );
        ++iter;
        }
      vCount++;
      }
    }
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
void
NJetLDAGenerator< ImageT, LabelmapT >
::SetForceIntensityConsistency( bool _forceIntensity )
{
  m_ForceIntensityConsistency = _forceIntensity;
}

template < class ImageT, class LabelmapT >
bool
NJetLDAGenerator< ImageT, LabelmapT >
::GetForceIntensityConsistency( void )
{
  return m_ForceIntensityConsistency;
}

template < class ImageT, class LabelmapT >
void
NJetLDAGenerator< ImageT, LabelmapT >
::SetForceOrientationInsensitivity( bool _forceOrientationInsensitivity )
{
  m_ForceOrientationInsensitivity = _forceOrientationInsensitivity;
}

template < class ImageT, class LabelmapT >
bool
NJetLDAGenerator< ImageT, LabelmapT >
::GetForceOrientationInsensitivity( void )
{
  return m_ForceOrientationInsensitivity;
}

template < class ImageT, class LabelmapT >
void
NJetLDAGenerator< ImageT, LabelmapT >
::GenerateLDA()
{
  Superclass::GenerateLDA();

  if( m_ForceIntensityConsistency || m_ForceOrientationInsensitivity )
    {
    unsigned int vCount = 0;
    std::vector< int > orientationNum( this->GetNumberOfFeatures(), 0 );
    for( unsigned int i=0; i<this->GetNumberOfNJetImages(); i++ )
      {
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
          m_NJetImageList[ 0 ],
          m_NJetImageList[ 0 ]->GetLargestPossibleRegion() );
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
NJetLDAGenerator< ImageT, LabelmapT >
::Update( void )
{
  this->GenerateFeatureImages();

  Superclass::Update();
}

template <class ImageT, class LabelmapT >
void
NJetLDAGenerator< ImageT, LabelmapT >
::UpdateLDAImages( void )
{
  this->GenerateFeatureImages();

  Superclass::UpdateLDAImages();
}

template <class ImageT, class LabelmapT >
void
NJetLDAGenerator< ImageT, LabelmapT >
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

  os << indent << "NJetImageList.size() = " << m_NJetImageList.size()
    << std::endl;
}

}

}

#endif //NJetLDAGenerator_txx
