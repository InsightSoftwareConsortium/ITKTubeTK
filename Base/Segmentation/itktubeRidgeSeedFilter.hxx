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

template< class TImage, class TLabelMap >
RidgeSeedFilter< TImage, TLabelMap >
::RidgeSeedFilter( void )
{
  m_SeedFeatureGenerator = SeedFeatureGeneratorType::New();
  m_RidgeFeatureGenerator = RidgeFeatureGeneratorType::New();
  m_SeedFeatureGenerator->SetInputFeatureVectorGenerator(
    m_RidgeFeatureGenerator );

  m_PDFSegmenter = PDFSegmenterType::New();
  m_PDFSegmenter->SetReclassifyObjectLabels( true );
  m_PDFSegmenter->SetReclassifyNotObjectLabels( true );
  m_PDFSegmenter->SetForceClassification( true );
  m_PDFSegmenter->SetErodeRadius( 0 );
  m_PDFSegmenter->SetHoleFillIterations( 0 );
  m_PDFSegmenter->SetOutlierRejectPortion( 0.01 );
  m_PDFSegmenter->SetProbabilityImageSmoothingStandardDeviation( 0.1 );
  m_PDFSegmenter->SetHistogramSmoothingStandardDeviation( 0.5 );

  m_RidgeId = 255;
  m_BackgroundId = 127;
  m_UnknownId = 0;

  m_SeedTolerance = 1.0;

  m_LabelMap = NULL;

  m_Skeletonize = true;
}

template< class TImage, class TLabelMap >
RidgeSeedFilter< TImage, TLabelMap >
::~RidgeSeedFilter( void )
{
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetInput( typename ImageType::Pointer img )
{
  m_SeedFeatureGenerator->SetInput( img );
  m_RidgeFeatureGenerator->SetInput( img );
}

template < class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::AddInput( typename ImageType::Pointer img )
{
  m_SeedFeatureGenerator->AddInput( img );
  m_RidgeFeatureGenerator->AddInput( img );
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetLabelMap( typename LabelMapType::Pointer img )
{
  m_SeedFeatureGenerator->SetLabelMap( img );
  m_PDFSegmenter->SetLabelMap( img );
}

template< class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::SeedFeatureGeneratorType::
Pointer
RidgeSeedFilter< TImage, TLabelMap >
::GetSeedFeatureGenerator( void )
{
  return m_SeedFeatureGenerator;
}

template< class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::RidgeFeatureGeneratorType::
Pointer
RidgeSeedFilter< TImage, TLabelMap >
::GetRidgeFeatureGenerator( void )
{
  return m_RidgeFeatureGenerator;
}

template< class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::PDFSegmenterType::Pointer
RidgeSeedFilter< TImage, TLabelMap >
::GetPDFSegmenter( void )
{
  return m_PDFSegmenter;
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetIntensityRange( float intensityMin, float intensityMax )
{
  m_RidgeFeatureGenerator->SetIntensityRange( intensityMin, intensityMax );
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetIntensityMin( float intensityMin )
{
  m_RidgeFeatureGenerator->SetIntensityMin( intensityMin );
}

template< class TImage, class TLabelMap >
float
RidgeSeedFilter< TImage, TLabelMap >
::GetIntensityMin( void ) const
{
  return m_RidgeFeatureGenerator->GetIntensityMin();
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetIntensityMax( float intensityMax )
{
  m_RidgeFeatureGenerator->SetIntensityMax( intensityMax );
}

template< class TImage, class TLabelMap >
float
RidgeSeedFilter< TImage, TLabelMap >
::GetIntensityMax( void ) const
{
  return m_RidgeFeatureGenerator->GetIntensityMax();
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetScales( const RidgeScalesType & scales )
{
  m_RidgeFeatureGenerator->SetScales( scales );
}

template< class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::RidgeScalesType
RidgeSeedFilter< TImage, TLabelMap >
::GetScales( void ) const
{
  return m_RidgeFeatureGenerator->GetScales();
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetWhitenMeans( const WhitenMeansType & means )
{
  m_RidgeFeatureGenerator->SetWhitenMeans( means );
}

template < class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetWhitenStdDevs( const WhitenStdDevsType & stdDevs )
{
  m_RidgeFeatureGenerator->SetWhitenStdDevs( stdDevs );
}

template < class TImage, class TLabelMap >
const typename RidgeSeedFilter< TImage, TLabelMap >::WhitenMeansType &
RidgeSeedFilter< TImage, TLabelMap >
::GetWhitenMeans( void ) const
{
  return m_RidgeFeatureGenerator->GetWhitenMeans();
}

template < class TImage, class TLabelMap >
const typename RidgeSeedFilter< TImage, TLabelMap >::WhitenStdDevsType &
RidgeSeedFilter< TImage, TLabelMap >
::GetWhitenStdDevs( void ) const
{
  return m_RidgeFeatureGenerator->GetWhitenStdDevs();
}

template < class TImage, class TLabelMap >
unsigned int
RidgeSeedFilter< TImage, TLabelMap >
::GetNumberOfBasis( void ) const
{
  return m_SeedFeatureGenerator->GetNumberOfBasis();
}

template < class TImage, class TLabelMap >
double
RidgeSeedFilter< TImage, TLabelMap >
::GetBasisValue( unsigned int num ) const
{
  return m_SeedFeatureGenerator->GetBasisValue( num );
}

template < class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::VectorType
RidgeSeedFilter< TImage, TLabelMap >
::GetBasisVector( unsigned int num ) const
{
  return m_SeedFeatureGenerator->GetBasisVector( num );
}

template < class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::MatrixType
RidgeSeedFilter< TImage, TLabelMap >
::GetBasisMatrix( void ) const
{
  return m_SeedFeatureGenerator->GetBasisMatrix();
}

template < class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::VectorType
RidgeSeedFilter< TImage, TLabelMap >
::GetBasisValues( void ) const
{
  return m_SeedFeatureGenerator->GetBasisValues();
}

template < class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::BasisImageType::Pointer
RidgeSeedFilter< TImage, TLabelMap >
::GetBasisImage( unsigned int num ) const
{
  return m_SeedFeatureGenerator->GetFeatureImage( num );
}

template < class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetBasisValue( unsigned int num, double value )
{
  return m_SeedFeatureGenerator->SetBasisValue( num, value );
}

template < class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetBasisVector( unsigned int num, const VectorType & vector )
{
  return m_SeedFeatureGenerator->SetBasisVector( num, vector );
}

template < class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetBasisMatrix( const MatrixType & matrix )
{
  return m_SeedFeatureGenerator->SetBasisMatrix( matrix );
}

template < class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetBasisValues( const VectorType & values )
{
  return m_SeedFeatureGenerator->SetBasisValues( values );
}

template < class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::ProbabilityImageType::Pointer
RidgeSeedFilter< TImage, TLabelMap >
::GetClassProbabilityForInput( unsigned int objectNum ) const
{
  return m_PDFSegmenter->GetClassProbabilityForInput( objectNum );
}

template < class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::ProbabilityImageType::Pointer
RidgeSeedFilter< TImage, TLabelMap >
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
  for( unsigned int c = 0; c < this->GetNumberOfObjectIds(); c++ )
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
    resultIter.Set( classIter.Get()
      * ( classIter.Get() - resultIter.Get() ) / ( classIter.Get() ) );
    ++resultIter;
    ++classIter;
    }

  return resultImage;
}


template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::Update( void )
{
  m_RidgeFeatureGenerator->GenerateData();

  m_SeedFeatureGenerator->SetObjectId( m_RidgeId );
  m_SeedFeatureGenerator->AddObjectId( m_BackgroundId );
  m_PDFSegmenter->SetObjectId( m_RidgeId );
  m_PDFSegmenter->AddObjectId( m_BackgroundId );
  m_PDFSegmenter->SetVoidId( m_UnknownId );

  m_PDFSegmenter->SetObjectPDFWeight( 0, m_SeedTolerance * 0.25 );

  m_SeedFeatureGenerator->SetNumberOfBasisToUseAsFeatures( 3 );

  m_SeedFeatureGenerator->GenerateBasis();
  m_PDFSegmenter->SetInput( 0,
    m_SeedFeatureGenerator->GetFeatureImage( 0 ) );
  m_PDFSegmenter->SetInput( 1,
    m_SeedFeatureGenerator->GetFeatureImage( 1 ) );
  m_PDFSegmenter->SetInput( 2,
    m_SeedFeatureGenerator->GetFeatureImage( 2 ) );
  m_PDFSegmenter->Update();
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::ClassifyImages( void )
{
  m_RidgeFeatureGenerator->GenerateData();

  m_SeedFeatureGenerator->SetObjectId( m_RidgeId );
  m_SeedFeatureGenerator->AddObjectId( m_BackgroundId );
  m_PDFSegmenter->SetObjectId( m_RidgeId );
  m_PDFSegmenter->AddObjectId( m_BackgroundId );
  m_PDFSegmenter->SetVoidId( m_UnknownId );

  m_SeedFeatureGenerator->SetNumberOfBasisToUseAsFeatures( 3 );

  typename LabelMapType::Pointer tmpLabelMap =
    m_SeedFeatureGenerator->GetLabelMap();
  m_SeedFeatureGenerator->SetLabelMap( NULL );
  m_PDFSegmenter->SetInput( 0,
    m_SeedFeatureGenerator->GetFeatureImage( 0 ) );
  m_PDFSegmenter->SetInput( 1,
    m_SeedFeatureGenerator->GetFeatureImage( 1 ) );
  m_PDFSegmenter->SetInput( 2,
    m_SeedFeatureGenerator->GetFeatureImage( 2 ) );
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
    typedef itk::BinaryBallStructuringElement< LabelMapPixelType,
      ImageDimension> SEType;

    typename FilterType::Pointer filter;

    filter = FilterType::New();
    filter->SetInput( m_LabelMap );
    filter->Update();
    m_LabelMap = filter->GetOutput();
    }
}

template< class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::LabelMapType::Pointer
RidgeSeedFilter< TImage, TLabelMap >
::GetOutput( void )
{
  return m_LabelMap;
}

template< class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::OutputImageType::Pointer
RidgeSeedFilter< TImage, TLabelMap >
::GetOutputSeedScales( void )
{
  int num = m_SeedFeatureGenerator->GetInputFeatureVectorGenerator()->
    GetNumberOfFeatures();
  num = num - 6;
  return m_SeedFeatureGenerator->GetInputFeatureVectorGenerator()->
    GetFeatureImage( num );
}


template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
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
