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

#ifndef __itktubeStructureTensorRecursiveGaussianImageFilter_hxx
#define __itktubeStructureTensorRecursiveGaussianImageFilter_hxx

#include "itktubeStructureTensorRecursiveGaussianImageFilter.h"

#include <itkImageRegionIteratorWithIndex.h>

namespace itk
{

namespace tube
{

/**
 * Constructor
 */
template< class TInputImage, class TOutputImage >
StructureTensorRecursiveGaussianImageFilter<TInputImage, TOutputImage>
::StructureTensorRecursiveGaussianImageFilter( void )
{
  m_NormalizeAcrossScale = false;

  unsigned int imageDimensionMinus1 = static_cast<int>( ImageDimension )-1;
  if( ImageDimension > 1 )
    {
    m_SmoothingFilters.resize( imageDimensionMinus1 );
    }

  for( unsigned int i = 0; i < imageDimensionMinus1; i++ )
    {
    m_SmoothingFilters[ i ] = GaussianFilterType::New();
    m_SmoothingFilters[ i ]->SetOrder( GaussianOrderEnum::ZeroOrder );
    m_SmoothingFilters[ i ]->SetNormalizeAcrossScale( m_NormalizeAcrossScale );
    m_SmoothingFilters[ i ]->ReleaseDataFlagOn();
    }

  // Outer Gaussian smoothing filter
  m_TensorComponentSmoothingFilter = GaussianFilterType::New();
  m_TensorComponentSmoothingFilter->SetOrder(
    GaussianOrderEnum::ZeroOrder );
  m_TensorComponentSmoothingFilter->SetNormalizeAcrossScale(
    m_NormalizeAcrossScale );
  //m_TensorComponentSmoothingFilter->ReleaseDataFlagOn();

  m_DerivativeFilter = DerivativeFilterType::New();
  m_DerivativeFilter->SetOrder( GaussianOrderEnum::FirstOrder );
  m_DerivativeFilter->SetNormalizeAcrossScale( m_NormalizeAcrossScale );
  m_DerivativeFilter->SetInput( this->GetInput() );

  if( ImageDimension > 1 )
    {
    m_SmoothingFilters[0]->SetInput( m_DerivativeFilter->GetOutput() );
    }

  for( unsigned int i = 1; i < imageDimensionMinus1; i++ )
    {
    m_SmoothingFilters[ i ]->SetInput( m_SmoothingFilters[i-1]->GetOutput() );
    }

  m_ImageAdaptor = OutputImageAdaptorType::New();

  this->SetSigma( 1.0 );
  this->SetSigmaOuter( 1.0 );

}

/**
 * Set value of Sigma
 */
template< class TInputImage, class TOutputImage >
void
StructureTensorRecursiveGaussianImageFilter<TInputImage, TOutputImage>
::SetSigma( RealType sigma )
{

  m_Sigma = sigma;
  for( unsigned int i = 0; i < ImageDimension - 1; i++ )
    {
    m_SmoothingFilters[ i ]->SetSigma( sigma );
    }
  m_DerivativeFilter->SetSigma( sigma );

  this->Modified();

}


/**
 * Set value of SigmaOuter
 */
template< class TInputImage, class TOutputImage >
void
StructureTensorRecursiveGaussianImageFilter<TInputImage, TOutputImage>
::SetSigmaOuter( RealType sigma )
{
  m_SigmaOuter = sigma;
  m_TensorComponentSmoothingFilter->SetSigma( sigma );
  this->Modified();
}

/**
 * Set Normalize Across Scale Space
 */
template< class TInputImage, class TOutputImage >
void
StructureTensorRecursiveGaussianImageFilter<TInputImage, TOutputImage>
::SetNormalizeAcrossScale( bool normalize )
{

  m_NormalizeAcrossScale = normalize;

  for( unsigned int i = 0; i<ImageDimension-1; i++ )
    {
    m_SmoothingFilters[ i ]->SetNormalizeAcrossScale( normalize );
    }
  m_DerivativeFilter->SetNormalizeAcrossScale( normalize );

  this->Modified();

}

template< class TInputImage, class TOutputImage >
void
StructureTensorRecursiveGaussianImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method. this should
  // copy the output requested region to the input requested region
  Superclass::GenerateInputRequestedRegion();

  // This filter needs all of the input
  typename
      StructureTensorRecursiveGaussianImageFilter< TInputImage, TOutputImage >
      ::InputImagePointer image
      = const_cast< InputImageType * >( this->GetInput() );
  image->SetRequestedRegion( this->GetInput()->GetLargestPossibleRegion() );
}

template< class TInputImage, class TOutputImage >
void
StructureTensorRecursiveGaussianImageFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion( DataObject *output )
{
  TOutputImage *out = dynamic_cast< TOutputImage* >( output );

  if( out )
    {
    out->SetRequestedRegion( out->GetLargestPossibleRegion() );
    }
}

/**
 * Compute filter for Gaussian kernel
 */
template< class TInputImage, class TOutputImage >
void
StructureTensorRecursiveGaussianImageFilter<TInputImage, TOutputImage >
::GenerateData( void )
{
  // Create a process accumulator for tracking the progress of this
  // mini-pipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter( this );

  // Compute the contribution of each filter to the total progress.
  const double weight = 1.0 / ( ImageDimension * ImageDimension );
  for( unsigned int i = 0; i<ImageDimension-1; i++ )
    {
    progress->RegisterInternalFilter( m_SmoothingFilters[i], weight );
    }
  progress->RegisterInternalFilter( m_DerivativeFilter, weight );

  const typename TInputImage::ConstPointer   inputImage( this->GetInput() );

  m_ImageAdaptor->SetImage( this->GetOutput() );
  m_ImageAdaptor->SetLargestPossibleRegion(
      inputImage->GetLargestPossibleRegion() );
  m_ImageAdaptor->SetBufferedRegion( inputImage->GetBufferedRegion() );
  m_ImageAdaptor->SetRequestedRegion( inputImage->GetRequestedRegion() );
  m_ImageAdaptor->Allocate();

  m_DerivativeFilter->SetInput( inputImage );

  unsigned int imageDimensionMinus1 = static_cast<int>( ImageDimension )-1;
  for( unsigned int dim=0; dim < ImageDimension; dim++ )
    {
    unsigned int i=0;
    unsigned int j=0;
    while( i< imageDimensionMinus1 )
      {
      if( i == dim )
        {
        j++;
        }
      m_SmoothingFilters[ i ]->SetDirection( j );
      i++;
      j++;
      }
    m_DerivativeFilter->SetDirection( dim );

    GaussianFilterPointer lastFilter;
    if( ImageDimension > 1 )
      {
      int imageDimensionMinus2 = static_cast<int>( ImageDimension )-2;
      lastFilter = m_SmoothingFilters[imageDimensionMinus2];
      lastFilter->Update();
      }
    else
      {
      m_DerivativeFilter->Update();
      }

    //progress->ResetFilterProgressAndKeepAccumulatedProgress();

    // Copy the results to the corresponding component
    // on the output image of vectors
    m_ImageAdaptor->SelectNthElement( dim );

    typename RealImageType::Pointer derivativeImage = lastFilter->GetOutput();

    ImageRegionIteratorWithIndex< RealImageType > it(
      derivativeImage,
      derivativeImage->GetRequestedRegion() );

    ImageRegionIteratorWithIndex< OutputImageAdaptorType > ot(
      m_ImageAdaptor,
      m_ImageAdaptor->GetRequestedRegion() );

    const RealType spacing = inputImage->GetSpacing()[ dim ];

    it.GoToBegin();
    ot.GoToBegin();
    while( !it.IsAtEnd() )
      {
      ot.Set( it.Get() / spacing );
      ++it;
      ++ot;
      }

    }

  //Calculate the outer ( diadic ) product of the gradient.
  ImageRegionIteratorWithIndex< OutputImageType > ottensor(
    this->GetOutput(),
    this->GetOutput()->GetRequestedRegion() );

  ImageRegionIteratorWithIndex< OutputImageType > itgradient(
    this->GetOutput(),
    this->GetOutput()->GetRequestedRegion() );

  const unsigned int numberTensorElements
      = ( ImageDimension*( ImageDimension+1 ) )/2;
  std::vector<InternalRealType> tmp( numberTensorElements );

  ottensor.GoToBegin();
  itgradient.GoToBegin();
  while( !itgradient.IsAtEnd() )
    {
    unsigned int count = 0;
    for( unsigned int j = 0; j < ImageDimension; ++j )
      {
      for( unsigned int k = j; k < ImageDimension; ++k )
        {
        tmp[count++] = itgradient.Get()[j]*itgradient.Get()[k];
        }
      }
    for( unsigned int j = 0; j < numberTensorElements; ++j )
      {
      ottensor.Value()[j] = tmp[j];
      }

    ++itgradient;
    ++ottensor;
    }

  //Finally, smooth the outer product components
  typedef typename itk::Image<InternalRealType, ImageDimension>
    ComponentImageType;

  for( unsigned int i =0; i < numberTensorElements; i++ )
    {
    typename ComponentImageType::Pointer componentImage =
      ComponentImageType::New();
    componentImage->SetLargestPossibleRegion(
      inputImage->GetLargestPossibleRegion() );
    componentImage->SetBufferedRegion( inputImage->GetBufferedRegion() );
    componentImage->SetRequestedRegion( inputImage->GetRequestedRegion() );
    componentImage->Allocate();

    ImageRegionIteratorWithIndex< ComponentImageType > compit(
      componentImage, componentImage->GetRequestedRegion() );

    ottensor.GoToBegin();
    compit.GoToBegin();
    while( !compit.IsAtEnd() )
      {
      compit.Value() = ottensor.Get()[i];
      ++compit;
      ++ottensor;
      }

    m_TensorComponentSmoothingFilter->SetInput( componentImage );
    m_TensorComponentSmoothingFilter->Update();

    ImageRegionIteratorWithIndex< ComponentImageType >
    smoothedCompIt( m_TensorComponentSmoothingFilter->GetOutput(),
      m_TensorComponentSmoothingFilter->GetOutput()->GetRequestedRegion() );

    ottensor.GoToBegin();
    smoothedCompIt.GoToBegin();

    while( !ottensor.IsAtEnd() )
      {
      ottensor.Value()[i] = smoothedCompIt.Get();
      ++smoothedCompIt;
      ++ottensor;
      }

    }
}


template< class TInputImage, class TOutputImage >
void
StructureTensorRecursiveGaussianImageFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << "NormalizeAcrossScale: " << m_NormalizeAcrossScale << std::endl;
  os << "Sigma: " << m_Sigma << std::endl;
  os << "SigmaOuter: " << m_SigmaOuter << std::endl;
}

} // End namespace tube

} // End namespace itk

// End !defined( __itktubeStructureTensorRecursiveGaussianImageFilter_hxx )
#endif
