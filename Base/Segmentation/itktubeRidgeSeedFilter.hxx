/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeRidgeSeedFilter_hxx
#define __itktubeRidgeSeedFilter_hxx

#include "itktubeRidgeSeedFilter.h"

#include "tubeMatrixMath.h"

#include <itkImage.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkProgressReporter.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkBinaryThinningImageFilter.h>

#include <limits>

namespace itk
{

namespace tube
{

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::RidgeSeedFilter( void )
{
  m_SeedFeatureGenerator = SeedFeatureGeneratorType::New();
  m_RidgeFeatureGenerator = RidgeFeatureGeneratorType::New();
  m_SeedFeatureGenerator->SetInputFeatureVectorGenerator(
    m_RidgeFeatureGenerator );

  m_SeedFeatureGenerator->SetNumberOfLDABasisToUseAsFeatures( 1 );
  m_SeedFeatureGenerator->SetNumberOfPCABasisToUseAsFeatures(
    TNumberOfFeatures-1 );

  m_PDFSegmenter = PDFSegmenterType::New();
  m_PDFSegmenter->SetReclassifyObjectLabels( true );
  m_PDFSegmenter->SetReclassifyNotObjectLabels( true );
  m_PDFSegmenter->SetForceClassification( true );
  m_PDFSegmenter->SetErodeRadius( 0 );
  m_PDFSegmenter->SetHoleFillIterations( 0 );
  m_PDFSegmenter->SetOutlierRejectPortion( 0.01 );
  m_PDFSegmenter->SetProbabilityImageSmoothingStandardDeviation( 0.3 );
  m_PDFSegmenter->SetHistogramSmoothingStandardDeviation( 4 );

  m_RidgeId = 255;
  m_BackgroundId = 127;
  m_UnknownId = 0;

  m_SeedTolerance = 1.0;

  m_LabelMap = NULL;

  m_Skeletonize = true;

  m_UseIntensityOnly = false;

  m_TrainClassifier = true;
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::~RidgeSeedFilter( void )
{
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetInput( typename ImageType::Pointer img )
{
  m_SeedFeatureGenerator->SetInput( img );
  m_RidgeFeatureGenerator->SetInput( img );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::AddInput( typename ImageType::Pointer img )
{
  m_SeedFeatureGenerator->AddInput( img );
  m_RidgeFeatureGenerator->AddInput( img );
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetLabelMap( typename LabelMapType::Pointer img )
{
  m_SeedFeatureGenerator->SetLabelMap( img );
  m_PDFSegmenter->SetLabelMap( img );
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::SeedFeatureGeneratorType::
Pointer
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetSeedFeatureGenerator( void )
{
  return m_SeedFeatureGenerator;
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::RidgeFeatureGeneratorType::
Pointer
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetRidgeFeatureGenerator( void )
{
  return m_RidgeFeatureGenerator;
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::PDFSegmenterType::Pointer
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetPDFSegmenter( void )
{
  return m_PDFSegmenter;
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetScales( const RidgeScalesType & scales )
{
  m_RidgeFeatureGenerator->SetScales( scales );
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::RidgeScalesType
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetScales( void ) const
{
  return m_RidgeFeatureGenerator->GetScales();
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetNumberOfPCABasisToUseAsFeatures( unsigned int num )
{
  m_SeedFeatureGenerator->SetNumberOfPCABasisToUseAsFeatures( num );
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
unsigned int
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetNumberOfPCABasisToUseAsFeatures( void ) const
{
  return m_SeedFeatureGenerator->GetNumberOfPCABasisToUseAsFeatures();
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetNumberOfLDABasisToUseAsFeatures( unsigned int num )
{
  m_SeedFeatureGenerator->SetNumberOfLDABasisToUseAsFeatures( num );
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
unsigned int
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetNumberOfLDABasisToUseAsFeatures( void ) const
{
  return m_SeedFeatureGenerator->GetNumberOfLDABasisToUseAsFeatures();
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetInputWhitenMeans( const WhitenMeansType & means )
{
  m_RidgeFeatureGenerator->SetWhitenMeans( means );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
const typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::WhitenMeansType &
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetInputWhitenMeans( void ) const
{
  return m_RidgeFeatureGenerator->GetWhitenMeans();
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetInputWhitenStdDevs( const WhitenStdDevsType & stdDevs )
{
  m_RidgeFeatureGenerator->SetWhitenStdDevs( stdDevs );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
const typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::WhitenStdDevsType &
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetInputWhitenStdDevs( void ) const
{
  return m_RidgeFeatureGenerator->GetWhitenStdDevs();
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetOutputWhitenMeans( const WhitenMeansType & means )
{
  m_SeedFeatureGenerator->SetWhitenMeans( means );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
const typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::WhitenMeansType &
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetOutputWhitenMeans( void ) const
{
  return m_SeedFeatureGenerator->GetWhitenMeans();
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetOutputWhitenStdDevs( const WhitenStdDevsType & stdDevs )
{
  m_SeedFeatureGenerator->SetWhitenStdDevs( stdDevs );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
const typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::WhitenStdDevsType &
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetOutputWhitenStdDevs( void ) const
{
  return m_SeedFeatureGenerator->GetWhitenStdDevs();
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
unsigned int
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetNumberOfBasis( void ) const
{
  return m_SeedFeatureGenerator->GetNumberOfFeatures();
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
double
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetBasisValue( unsigned int num ) const
{
  return m_SeedFeatureGenerator->GetBasisValue( num );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >::VectorType
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetBasisVector( unsigned int num ) const
{
  return m_SeedFeatureGenerator->GetBasisVector( num );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >::MatrixType
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetBasisMatrix( void ) const
{
  return m_SeedFeatureGenerator->GetBasisMatrix();
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >::VectorType
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetBasisValues( void ) const
{
  return m_SeedFeatureGenerator->GetBasisValues();
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::BasisImageType::Pointer
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetBasisImage( unsigned int num ) const
{
  return m_SeedFeatureGenerator->GetFeatureImage( num );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetBasisValue( unsigned int num, double value )
{
  return m_SeedFeatureGenerator->SetBasisValue( num, value );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetBasisVector( unsigned int num, const VectorType & vector )
{
  return m_SeedFeatureGenerator->SetBasisVector( num, vector );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetBasisMatrix( const MatrixType & matrix )
{
  return m_SeedFeatureGenerator->SetBasisMatrix( matrix );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::SetBasisValues( const VectorType & values )
{
  return m_SeedFeatureGenerator->SetBasisValues( values );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::ProbabilityImageType::Pointer
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetClassProbabilityForInput( unsigned int objectNum ) const
{
  return m_PDFSegmenter->GetClassProbabilityForInput( objectNum );
}

template < class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::ProbabilityImageType::Pointer
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetClassProbabilityDifferenceForInput( unsigned int objectNum ) const
{
  typename ProbabilityImageType::Pointer classImage = m_PDFSegmenter->
    GetClassProbabilityForInput( objectNum );

  typename ProbabilityImageType::RegionType region = classImage->
    GetLargestPossibleRegion();

  typename ProbabilityImageType::Pointer resultImage =
    ProbabilityImageType::New();
  resultImage->SetRegions( region );
  resultImage->CopyInformation( classImage );
  resultImage->Allocate();
  resultImage->FillBuffer( 0 );

  itk::ImageRegionIterator< ProbabilityImageType > resultIter(
    resultImage, region );
  for( unsigned int c = 0; c < 2; c++ )
    {
    if( c != objectNum )
      {
      itk::ImageRegionConstIterator< ProbabilityImageType >
        classIter( m_PDFSegmenter->GetClassProbabilityForInput( c ),
          region );
      resultIter.GoToBegin();
      while( ! resultIter.IsAtEnd() )
        {
        if( classIter.Get() > resultIter.Get() )
          {
          resultIter.Set( classIter.Get() );
          }
        ++resultIter;
        ++classIter;
        }
      }
    }
  itk::ImageRegionConstIterator< ProbabilityImageType > classIter(
    m_PDFSegmenter->GetClassProbabilityForInput( objectNum ), region );
  resultIter.GoToBegin();
  while( ! resultIter.IsAtEnd() )
    {
    double denum = classIter.Get() + resultIter.Get();
    if( denum == denum && denum != 0 )
      {
      resultIter.Set( classIter.Get()
        / ( classIter.Get() + resultIter.Get() ) );
      }
    ++resultIter;
    ++classIter;
    }

  return resultImage;
}


template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::Update( void )
{
  m_RidgeFeatureGenerator->SetUseIntensityOnly( m_UseIntensityOnly );
  m_RidgeFeatureGenerator->Update();

  m_SeedFeatureGenerator->SetObjectId( m_RidgeId );
  m_SeedFeatureGenerator->AddObjectId( m_BackgroundId );
  m_PDFSegmenter->SetObjectId( m_RidgeId );
  m_PDFSegmenter->AddObjectId( m_BackgroundId );
  m_PDFSegmenter->SetVoidId( m_UnknownId );

  m_PDFSegmenter->SetObjectPDFWeight( 0, m_SeedTolerance );

  if( m_TrainClassifier )
    {
    m_RidgeFeatureGenerator->SetUpdateWhitenStatisticsOnUpdate( true );
    m_RidgeFeatureGenerator->Update();

    m_SeedFeatureGenerator->SetUpdateWhitenStatisticsOnUpdate( true );
    m_SeedFeatureGenerator->Update();

    for( unsigned int i=0; i<m_SeedFeatureGenerator->GetNumberOfFeatures();
      ++i )
      {
      m_PDFSegmenter->SetInput( i,
        m_SeedFeatureGenerator->GetFeatureImage( i ) );
      }

    m_PDFSegmenter->Update();
    }
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::ClassifyImages( void )
{
  typename LabelMapType::Pointer tmpLabelMap =
    m_SeedFeatureGenerator->GetLabelMap();
  m_SeedFeatureGenerator->SetLabelMap( NULL );

  for( unsigned int i=0; i<m_SeedFeatureGenerator->GetNumberOfFeatures();
    ++i )
    {
    m_PDFSegmenter->SetInput( i,
      m_SeedFeatureGenerator->GetFeatureImage( i ) );
    }

  m_PDFSegmenter->ClassifyImages();

  m_SeedFeatureGenerator->SetLabelMap( tmpLabelMap );

  m_LabelMap = m_PDFSegmenter->GetLabelMap();

  itk::ImageRegionIterator< LabelMapType > resultIter(
    m_LabelMap, m_LabelMap->GetLargestPossibleRegion() );
  while( !resultIter.IsAtEnd() )
    {
    if( resultIter.Get() == m_RidgeId )
      {
      resultIter.Set( 1 );
      }
    else
      {
      resultIter.Set( 0 );
      }
    ++resultIter;
    }

  if( m_Skeletonize )
    {
    typedef itk::BinaryThinningImageFilter< LabelMapType, LabelMapType >
      FilterType;

    typename FilterType::Pointer filter;

    filter = FilterType::New();
    filter->SetInput( m_LabelMap );
    filter->Update();
    m_LabelMap = filter->GetOutput();
    }
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::LabelMapType::Pointer
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetOutput( void )
{
  return m_LabelMap;
}

template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
typename RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures
  >::OutputImageType::Pointer
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::GetOutputSeedScales( void )
{
  int num = m_SeedFeatureGenerator->GetInputFeatureVectorGenerator()->
    GetNumberOfFeatures();
  if( !m_UseIntensityOnly )
    {
    num = num - 6;
    }
  else
    {
    num = num - 3;
    }
  return m_SeedFeatureGenerator->GetInputFeatureVectorGenerator()->
    GetFeatureImage( num );
}


template< class TImage, class TLabelMap, unsigned int TNumberOfFeatures >
void
RidgeSeedFilter< TImage, TLabelMap, TNumberOfFeatures >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "PDFSegmenter = " << m_PDFSegmenter << std::endl;
  os << indent << "RidgeFeatureGenerator = " << m_RidgeFeatureGenerator
    << std::endl;
  os << indent << "SeedFeatureGenerator = " << m_SeedFeatureGenerator
    << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeRidgeSeedFilter_hxx)
