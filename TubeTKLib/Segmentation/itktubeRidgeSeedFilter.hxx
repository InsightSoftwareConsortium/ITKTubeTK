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

#ifndef __itktubeRidgeSeedFilter_hxx
#define __itktubeRidgeSeedFilter_hxx

#include "itktubeRidgeSeedFilter.h"

#include "tubeMatrixMath.h"
#include "itktubePDFSegmenterParzen.h"

#include <itkImage.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkProgressReporter.h>
#include <itkBinaryThinningImageFilter.h>

#include <limits>

namespace itk
{

namespace tube
{

template< class TImage, class TLabelMap >
RidgeSeedFilter< TImage, TLabelMap >
::RidgeSeedFilter()
{
  m_RidgeFeatureGenerator = RidgeFeatureGeneratorType::New();

  m_SeedFeatureGenerator = SeedFeatureGeneratorType::New();
  m_SeedFeatureGenerator->SetInputFeatureVectorGenerator(
    m_RidgeFeatureGenerator );

  m_SeedFeatureGenerator->SetNumberOfLDABasisToUseAsFeatures( 1 );
  m_SeedFeatureGenerator->SetNumberOfPCABasisToUseAsFeatures( 3 );

  m_PDFSegmenter = nullptr;

  m_RidgeId = 255;
  m_BackgroundId = 127;
  m_UnknownId = 0;

  m_SeedTolerance = 1.0;

  m_LabelMap = nullptr;

  m_Skeletonize = true;

  m_UseIntensityOnly = false;

  m_TrainClassifier = true;

  m_RatioImageVector.resize(0);
}

template< class TImage, class TLabelMap >
RidgeSeedFilter< TImage, TLabelMap >
::~RidgeSeedFilter( void )
{
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetInput( const InputImageType * img )
{
  m_RidgeFeatureGenerator->SetInput( img );
  m_SeedFeatureGenerator->SetInput( img );
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetInput( unsigned int id, const InputImageType * img )
{
  m_RidgeFeatureGenerator->SetInput( id, img );
  m_SeedFeatureGenerator->SetInput( id, img );
}

template < class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::AddInput( const InputImageType * img )
{
  m_RidgeFeatureGenerator->AddInput( img );
  m_SeedFeatureGenerator->AddInput( img );
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetLabelMap( LabelMapType * img )
{
  m_SeedFeatureGenerator->SetLabelMap( img );
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
::SetPDFSegmenter( typename RidgeSeedFilter< TImage, TLabelMap >
  ::PDFSegmenterType * pdfSegmenter )
{
  m_PDFSegmenter = pdfSegmenter;
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
::SetInputWhitenMeans( const WhitenMeansType & means )
{
  m_RidgeFeatureGenerator->SetWhitenMeans( means );
}

template < class TImage, class TLabelMap >
const typename RidgeSeedFilter< TImage, TLabelMap >::WhitenMeansType &
RidgeSeedFilter< TImage, TLabelMap >
::GetInputWhitenMeans( void ) const
{
  return m_RidgeFeatureGenerator->GetWhitenMeans();
}

template < class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetInputWhitenStdDevs( const WhitenStdDevsType & stdDevs )
{
  m_RidgeFeatureGenerator->SetWhitenStdDevs( stdDevs );
}

template < class TImage, class TLabelMap >
const typename RidgeSeedFilter< TImage, TLabelMap >::WhitenStdDevsType &
RidgeSeedFilter< TImage, TLabelMap >
::GetInputWhitenStdDevs( void ) const
{
  return m_RidgeFeatureGenerator->GetWhitenStdDevs();
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetOutputWhitenMeans( const WhitenMeansType & means )
{
  m_SeedFeatureGenerator->SetWhitenMeans( means );
}

template < class TImage, class TLabelMap >
const typename RidgeSeedFilter< TImage, TLabelMap >::WhitenMeansType &
RidgeSeedFilter< TImage, TLabelMap >
::GetOutputWhitenMeans( void ) const
{
  return m_SeedFeatureGenerator->GetWhitenMeans();
}

template < class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetOutputWhitenStdDevs( const WhitenStdDevsType & stdDevs )
{
  m_SeedFeatureGenerator->SetWhitenStdDevs( stdDevs );
}

template < class TImage, class TLabelMap >
const typename RidgeSeedFilter< TImage, TLabelMap >::WhitenStdDevsType &
RidgeSeedFilter< TImage, TLabelMap >
::GetOutputWhitenStdDevs( void ) const
{
  return m_SeedFeatureGenerator->GetWhitenStdDevs();
}

template < class TImage, class TLabelMap >
unsigned int
RidgeSeedFilter< TImage, TLabelMap >
::GetNumberOfBasis( void ) const
{
  return m_SeedFeatureGenerator->GetNumberOfFeatures();
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
typename RidgeSeedFilter< TImage, TLabelMap >::FeatureImageType::Pointer
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
::GetClassProbabilityImage( unsigned int objectNum ) const
{
  return m_PDFSegmenter->GetClassProbabilityImage( objectNum );
}

template < class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::ProbabilityImageType::Pointer
RidgeSeedFilter< TImage, TLabelMap >
::GetClassLikelihoodRatioImage( unsigned int objectNum )
{
  unsigned int numClasses = m_PDFSegmenter->GetNumberOfClasses();
  if( m_RatioImageVector.size() != numClasses )
    {
    m_RatioImageVector.resize(numClasses);
    }

  typename ProbabilityImageType::Pointer classImage = m_PDFSegmenter->
    GetClassProbabilityImage( objectNum );

  typename ProbabilityImageType::RegionType region = classImage->
    GetLargestPossibleRegion();

  m_RatioImageVector[objectNum] = ProbabilityImageType::New();
  m_RatioImageVector[objectNum]->SetRegions( region );
  m_RatioImageVector[objectNum]->CopyInformation( classImage );
  m_RatioImageVector[objectNum]->Allocate();
  m_RatioImageVector[objectNum]->FillBuffer( 0 );

  itk::ImageRegionIterator< ProbabilityImageType > resultIter(
    m_RatioImageVector[objectNum], region );
  double backgroundMax = 0;
  for( unsigned int c = 0; c < numClasses; c++ )
    {
    if( c != objectNum )
      {
      itk::ImageRegionConstIterator< ProbabilityImageType >
        classIter( m_PDFSegmenter->GetClassProbabilityImage( c ),
          region );
      resultIter.GoToBegin();
      while( ! resultIter.IsAtEnd() )
        {
        double tf = classIter.Get();
        if( tf > resultIter.Get() )
          {
          resultIter.Set( tf );
          if( tf > backgroundMax )
            {
            backgroundMax = tf;
            }
          }
        ++resultIter;
        ++classIter;
        }
      }
    }
  double tubeMax = 0;
  itk::ImageRegionConstIterator< ProbabilityImageType > classIter(
    m_PDFSegmenter->GetClassProbabilityImage( objectNum ), region );
  resultIter.GoToBegin();
  while( ! resultIter.IsAtEnd() )
    {
    double tf = classIter.Get();
    if( tf > tubeMax )
      {
      tubeMax = tf;
      }
    double denum = tf + resultIter.Get();
    if( denum == denum && denum > 0 )
      {
      resultIter.Set( tf / denum );
      }
    else
      {
      resultIter.Set( 0 );
      }
    ++resultIter;
    ++classIter;
    }

  return m_RatioImageVector[objectNum];
}


template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::Update( void )
{
  //itk::TimeProbesCollectorBase timeCollector;

  //timeCollector.Start( "RidgeSeedFilter Update" );

  if( m_PDFSegmenter.IsNull() )
    {
    m_PDFSegmenter = PDFSegmenterParzenType::New();
    typename PDFSegmenterParzenType::Pointer tmpPDF =
      static_cast< PDFSegmenterParzenType  * >( m_PDFSegmenter.GetPointer() );
    tmpPDF->SetHistogramSmoothingStandardDeviation( 2 );
    tmpPDF->SetOutlierRejectPortion( 0.001 );
    }

  m_PDFSegmenter->SetFeatureVectorGenerator(
    m_SeedFeatureGenerator.GetPointer() );
  m_PDFSegmenter->SetReclassifyObjectLabels( true );
  m_PDFSegmenter->SetReclassifyNotObjectLabels( true );
  m_PDFSegmenter->SetForceClassification( true );
  m_PDFSegmenter->SetErodeDilateRadius( 0 );
  m_PDFSegmenter->SetHoleFillIterations( 5 );

  m_PDFSegmenter->SetInputLabelMap( m_SeedFeatureGenerator->GetLabelMap() );

  m_RidgeFeatureGenerator->SetUseIntensityOnly( m_UseIntensityOnly );

  //timeCollector.Start( "RidgeSeedFilter FeatureGenerator" );
  m_RidgeFeatureGenerator->Update();
  //timeCollector.Stop( "RidgeSeedFilter FeatureGenerator" );

  m_SeedFeatureGenerator->SetObjectId( m_RidgeId );
  m_SeedFeatureGenerator->AddObjectId( m_BackgroundId );
  m_PDFSegmenter->SetObjectId( m_RidgeId );
  m_PDFSegmenter->AddObjectId( m_BackgroundId );
  m_PDFSegmenter->SetVoidId( m_UnknownId );

  m_PDFSegmenter->SetObjectPDFWeight( 0, m_SeedTolerance );

  if( m_TrainClassifier )
    {
    //timeCollector.Start( "RidgeSeedFilter RidgeFeatureGenerator Update" );
    m_RidgeFeatureGenerator->SetUpdateWhitenStatisticsOnUpdate( true );
    m_RidgeFeatureGenerator->Update();
    //timeCollector.Stop( "RidgeSeedFilter RidgeFeatureGenerator Update" );

    //timeCollector.Start( "RidgeSeedFilter SeedFeatureGenerator Update" );
    m_SeedFeatureGenerator->SetUpdateWhitenStatisticsOnUpdate( true );
    m_SeedFeatureGenerator->Update();
    //timeCollector.Stop( "RidgeSeedFilter SeedFeatureGenerator Update" );

    //timeCollector.Start( "RidgeSeedFilter PDFSegmenter Update" );
    m_PDFSegmenter->Update();
    //timeCollector.Start( "RidgeSeedFilter PDFSegmenter Update" );
    }

  //timeCollector.Stop( "RidgeSeedFilter Update" );
  //timeCollector.Report();
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::ClassifyImages( void )
{
  typename LabelMapType::Pointer tmpLabelMap =
    m_SeedFeatureGenerator->GetLabelMap();
  m_SeedFeatureGenerator->SetLabelMap( nullptr );

  m_PDFSegmenter->ClassifyImages();

  m_SeedFeatureGenerator->SetLabelMap( tmpLabelMap );

  m_LabelMap = m_PDFSegmenter->GetOutputLabelMap();

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
    m_LabelMap = filter->GetOutput();
    filter->Update();
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

#endif // End !defined( __itktubeRidgeSeedFilter_hxx )
