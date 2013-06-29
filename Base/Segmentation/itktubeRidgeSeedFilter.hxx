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

template< class TImage, class TLabelMap >
RidgeSeedFilter< TImage, TLabelMap >
::RidgeSeedFilter( void )
{
  m_SeedFeatureGenerator = SeedFeatureGeneratorType::New();
  m_RidgeFeatureGenerator = RidgeFeatureGeneratorType::New();
  m_SeedFeatureGenerator->SetInputFeatureVectorGenerator(
    m_RidgeFeatureGenerator );

  m_PDFSegmenter = PDFSegmenterType::New();
  m_PDFSegmenter->SetVoidId( 10 );
  m_PDFSegmenter->SetReclassifyObjectLabels( true );
  m_PDFSegmenter->SetReclassifyNotObjectLabels( true );
  m_PDFSegmenter->SetForceClassification( true );
  m_PDFSegmenter->SetErodeRadius( 0 );
  m_PDFSegmenter->SetHoleFillIterations( 0 );
  m_PDFSegmenter->SetOutlierRejectPortion( 0.01 );
  m_PDFSegmenter->SetProbabilityImageSmoothingStandardDeviation( 0.1 );
  m_PDFSegmenter->SetHistogramSmoothingStandardDeviation( 0.5 );
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
  m_SeedFeatureGenerator->SetInputImage( img );
  m_RidgeFeatureGenerator->SetInputImage( img );
}

template < class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::AddInput( typename ImageType::Pointer img )
{
  m_SeedFeatureGenerator->AddInputImage( img );
  m_RidgeFeatureGenerator->AddInputImage( img );
}

template < class TImage, class TLabelMap >
typename TImage::Pointer
RidgeSeedFilter< TImage, TLabelMap >
::GetInput( unsigned int num )
{
  return m_SeedFeatureGenerator->GetInputImage( num );
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
typename RidgeSeedFilter< TImage, TLabelMap >::LabelMapType::Pointer
RidgeSeedFilter< TImage, TLabelMap >
::GetLabelMap( void )
{
  return m_SeedFeatureGenerator->GetLabelMap();
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
::GetIntensityMin( void )
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
::GetIntensityMax( void )
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
::GetScales( void )
{
  return m_RidgeFeatureGenerator->GetScales();
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::SetObjectId( ObjectIdType id )
{
  m_SeedFeatureGenerator->SetObjectId( id );
  m_PDFSegmenter->SetObjectId( id );
  m_PDFSegmenter->SetObjectPDFWeight( 0, 0.5 );
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::AddObjectId( ObjectIdType id )
{
  m_SeedFeatureGenerator->AddObjectId( id );
  m_PDFSegmenter->AddObjectId( id );
}

template< class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::ObjectIdType
RidgeSeedFilter< TImage, TLabelMap >
::GetObjectId( unsigned int num ) const
{
  return m_PDFSegmenter->GetObjectId( num );
}

template< class TImage, class TLabelMap >
unsigned int
RidgeSeedFilter< TImage, TLabelMap >
::GetNumberOfObjectIds( void ) const
{
  return m_PDFSegmenter->GetNumberOfObjectIds();
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::Update( void )
{
  m_SeedFeatureGenerator->GenerateBasis();
  m_PDFSegmenter->SetInputVolume( 0,
    m_SeedFeatureGenerator->GetFeatureImage( 0 ) );
  m_PDFSegmenter->SetInputVolume( 1,
    m_SeedFeatureGenerator->GetFeatureImage( 1 ) );
  m_PDFSegmenter->SetInputVolume( 2,
    m_SeedFeatureGenerator->GetFeatureImage( 2 ) );
  m_PDFSegmenter->Update();
}

template< class TImage, class TLabelMap >
void
RidgeSeedFilter< TImage, TLabelMap >
::ClassifyImages( void )
{
  typename LabelMapType::Pointer tmpLabelMap =
    m_SeedFeatureGenerator->GetLabelMap();
  m_SeedFeatureGenerator->SetLabelMap( NULL );
  m_PDFSegmenter->SetInputVolume( 0,
    m_SeedFeatureGenerator->GetFeatureImage( 0 ) );
  m_PDFSegmenter->SetInputVolume( 1,
    m_SeedFeatureGenerator->GetFeatureImage( 1 ) );
  m_PDFSegmenter->SetInputVolume( 2,
    m_SeedFeatureGenerator->GetFeatureImage( 2 ) );
  m_PDFSegmenter->ClassifyImages();
  m_SeedFeatureGenerator->SetLabelMap( tmpLabelMap );
}

template< class TImage, class TLabelMap >
typename RidgeSeedFilter< TImage, TLabelMap >::LabelMapType::Pointer
RidgeSeedFilter< TImage, TLabelMap >
::GetOutput( void )
{
  return m_PDFSegmenter->GetLabelMap();
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
