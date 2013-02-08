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
RidgeSeedGenerator< ImageT, LabelmapT >
::RidgeSeedGenerator()
{
  m_RidgeImage = NULL;

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
  unsigned int numFeatures = m_Scales.size() * 5 + 9;

  return numFeatures;
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

  ProgressReporter progress( this, 0,
    m_RidgeImage->GetLargestPossibleRegion().GetNumberOfPixels()*2, 100 );

  if( this->m_Labelmap.IsNotNull() )
    {
    unsigned int numClasses = this->GetNumberOfObjectIds();

    ImageRegionIteratorWithIndex< RidgeImageType > iter(
      m_RidgeImage, m_RidgeImage->GetLargestPossibleRegion() );
    typedef ImageRegionConstIteratorWithIndex< MaskImageType >
      ConstMaskImageIteratorType;
    ConstMaskImageIteratorType itInMask( this->m_Labelmap,
      this->m_Labelmap->GetLargestPossibleRegion() );
    bool found = false;
    ObjectIdType prevObjVal = static_cast<ObjectIdType>( itInMask.Get() )+1;
    while( !iter.IsAtEnd() )
      {
      progress.CompletedPixel();

      ObjectIdType val = static_cast<ObjectIdType>( itInMask.Get() );
      if( val != prevObjVal )
        {
        found = false;
        prevObjVal = val;
        for( unsigned int c=0; c<numClasses; c++ )
          {
          if( val == this->m_ObjectIdList[c] )
            {
            found = true;
            break;
            }
          }
        }
      if( found )
        {
        typename RidgeImageType::IndexType indx = iter.GetIndex();
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
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeScale );
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeIntensity );
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeRidgeness );
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeRoundness );
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeLevelness );
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeCurvature );
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeTangent[0] );
        this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeTangent[1] );
        }
      ++iter;
      ++itInMask;
      }
    int xPos = numFeatures-3;
    int yPos = numFeatures-2;
    int dotPos = numFeatures-1;
    iter.GoToBegin();
    itInMask.GoToBegin();
    while( !iter.IsAtEnd() )
      {
      progress.CompletedPixel();
      ObjectIdType val = static_cast<ObjectIdType>( itInMask.Get() );
      if( val != prevObjVal )
        {
        found = false;
        prevObjVal = val;
        for( unsigned int c=0; c<numClasses; c++ )
          {
          if( val == this->m_ObjectIdList[c] )
            {
            found = true;
            break;
            }
          }
        }
      if( found )
        {
        typename RidgeImageType::IndexType indx = iter.GetIndex();
        typename RidgeImageType::IndexType indx2 = iter.GetIndex();
        typename NJetFunctionType::VectorType t;
        typename NJetFunctionType::VectorType t2;
        t[0] = this->m_FeatureImageList[ xPos ]->GetPixel( indx );
        t[1] = this->m_FeatureImageList[ yPos ]->GetPixel( indx );
        indx2[0] = indx[0] + vnl_math_rnd( 2*t[0] );
        indx2[1] = indx[1] + vnl_math_rnd( 2*t[1] );
        if( ImageDimension > 2 )
          {
          t[2] = vcl_sqrt( vnl_math_abs( 1 - ( t[0]*t[0] + t[1]*t[1] ) ) );
          indx2[2] = indx[2] + vnl_math_rnd( 2*t[2] );
          }
        double dot = 0;
        if( this->m_FeatureImageList[ xPos ]->GetLargestPossibleRegion().IsInside( indx2 ) )
          {
          t2[0] = this->m_FeatureImageList[ xPos ]->GetPixel( indx2 );
          t2[1] = this->m_FeatureImageList[ yPos ]->GetPixel( indx2 );
          dot = t[0]*t2[0] + t[1]*t2[1];
          if( ImageDimension>2 )
            {
            t[2] = vcl_sqrt( vnl_math_abs( 1 - ( t[0]*t[0] + t[1]*t[1] ) ) );
            t2[2] = vcl_sqrt( vnl_math_abs( 1 - ( t2[0]*t2[0] + t2[1]*t2[1] ) ) );
            dot += t[2]*t2[2];
            }
          }
        this->m_FeatureImageList[ dotPos ]->SetPixel( indx, vnl_math_abs( dot ) );
        }
      ++iter;
      ++itInMask;
      }
    }
  else
    {
    ImageRegionIteratorWithIndex< LDAImageType > iter(
      m_RidgeImage, m_RidgeImage->GetLargestPossibleRegion() );
    while( !iter.IsAtEnd() )
      {
      progress.CompletedPixel();

      typename RidgeImageType::IndexType indx = iter.GetIndex();
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
      this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeScale );
      this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeIntensity );
      this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeRidgeness );
      this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeRoundness );
      this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeLevelness );
      this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeCurvature );
      this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeTangent[0] );
      this->m_FeatureImageList[ fcount++ ]->SetPixel( indx, extremeTangent[1] );
      ++iter;
      }
    int xPos = numFeatures-3;
    int yPos = numFeatures-2;
    int dotPos = numFeatures-1;
    iter.GoToBegin();
    while( !iter.IsAtEnd() )
      {
      progress.CompletedPixel();

      typename RidgeImageType::IndexType indx = iter.GetIndex();
      typename RidgeImageType::IndexType indx2 = iter.GetIndex();
      typename NJetFunctionType::VectorType t;
      typename NJetFunctionType::VectorType t2;
      t[0] = this->m_FeatureImageList[ xPos ]->GetPixel( indx );
      t[1] = this->m_FeatureImageList[ yPos ]->GetPixel( indx );
      indx2[0] = indx[0] + vnl_math_rnd( 2*t[0] );
      indx2[1] = indx[1] + vnl_math_rnd( 2*t[1] );
      if( ImageDimension > 2 )
        {
        t[2] = vcl_sqrt( vnl_math_abs( 1 - ( t[0]*t[0] + t[1]*t[1] ) ) );
        indx2[2] = indx[2] + vnl_math_rnd( 2*t[2] );
        }
      double dot = 0;
      if( this->m_FeatureImageList[ xPos ]->GetLargestPossibleRegion().IsInside( indx2 ) )
        {
        t2[0] = this->m_FeatureImageList[ xPos ]->GetPixel( indx2 );
        t2[1] = this->m_FeatureImageList[ yPos ]->GetPixel( indx2 );
        dot = t[0]*t2[0] + t[1]*t2[1];
        if( ImageDimension>2 )
          {
          t[2] = vcl_sqrt( vnl_math_abs( 1 - ( t[0]*t[0] + t[1]*t[1] ) ) );
          t2[2] = vcl_sqrt( vnl_math_abs( 1 - ( t2[0]*t2[0] + t2[1]*t2[1] ) ) );
          dot += t[2]*t2[2];
          }
        }
      this->m_FeatureImageList[ dotPos ]->SetPixel( indx, vnl_math_abs( dot ) );
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

  os << indent << "Scales.size() = " << m_Scales.size() << std::endl;
}

}

}

#endif //RidgeSeedGenerator_txx
