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
#ifndef __itkGradientBasedAngleOfIncidenceImageFilter_txx
#define __itkGradientBasedAngleOfIncidenceImageFilter_txx

#include "itkGradientBasedAngleOfIncidenceImageFilter.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkGradientImageFilter.h"
#include "itkProgressReporter.h"

namespace itk
{

template< class TInputImage, class TOutputImage, class TOperatorValue >
GradientBasedAngleOfIncidenceImageFilter< TInputImage,
  TOutputImage,
  TOperatorValue >
::GradientBasedAngleOfIncidenceImageFilter()
{
  typedef GradientImageFilter<
      OperatorImageType, OperatorValueType, OperatorValueType >
    DefaultGradientFilterType;
  this->m_GradientFilter = DefaultGradientFilterType::New().GetPointer();

  this->m_CastImageFilter = CastImageFilterType::New();

  this->m_UltrasoundProbeOrigin.Fill( 0.0 );

  this->m_GradientMagnitudeTolerance = 1.0e-7;
}


template< class TInputImage, class TOutputImage, class TOperatorValue >
void
GradientBasedAngleOfIncidenceImageFilter< TInputImage,
  TOutputImage,
  TOperatorValue >
::BeforeThreadedGenerateData( void )
{
  this->m_CastImageFilter->SetInput( this->GetInput() );
  this->m_GradientFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->m_GradientFilter->SetInput( this->m_CastImageFilter->GetOutput() );
  this->m_GradientFilter->Update();
}


template< class TInputImage, class TOutputImage, class TOperatorValue >
void
GradientBasedAngleOfIncidenceImageFilter< TInputImage,
  TOutputImage,
  TOperatorValue >
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
  ThreadIdType threadId )
{
  const InputImageType * input = this->GetInput();
  const GradientOutputImageType * gradient =
    this->m_GradientFilter->GetOutput();
  OutputImageType * output = this->GetOutput();

  const OriginType origin = this->m_UltrasoundProbeOrigin;
  const double gradientMagnitudeTolerance = this->m_GradientMagnitudeTolerance;

  typedef ImageRegionConstIteratorWithIndex< InputImageType > InputIteratorType;
  typedef ImageRegionConstIterator< GradientOutputImageType > GradientIteratorType;
  typedef ImageRegionIterator< OutputImageType >              OutputIteratorType;

  InputIteratorType inputIt( input, outputRegionForThread );
  GradientIteratorType gradientIt( gradient, outputRegionForThread );
  OutputIteratorType outputIt( output, outputRegionForThread );

  ProgressReporter progress( this, threadId,
    outputRegionForThread.GetNumberOfPixels() );

  for( inputIt.GoToBegin(), gradientIt.GoToBegin(), outputIt.GoToBegin();
       !outputIt.IsAtEnd();
       ++inputIt, ++gradientIt, ++outputIt )
    {
    typename InputImageType::IndexType index = inputIt.GetIndex();
    typename InputImageType::PointType point;
    input->TransformIndexToPhysicalPoint( index, point );

    typedef Vector< OperatorValueType, ImageDimension > VectorType;
    VectorType beamDirection = point - origin;
    beamDirection.Normalize();

    GradientOutputPixelType gradientPixel = gradientIt.Get();
    const typename GradientOutputPixelType::RealValueType gradientNorm =
      gradientPixel.GetNorm();
    gradientPixel /= gradientNorm;

    // output  scalar product of the two normalized vectors
    typedef typename OutputImageType::PixelType OutputPixelType;
    const OutputPixelType outputPixel = gradientPixel * beamDirection;
    if( vnl_math_isnan( outputPixel ) || gradientNorm < gradientMagnitudeTolerance )
      {
      outputIt.Set( NumericTraits< OutputPixelType >::Zero );
      }
    else
      {
      outputIt.Set( gradientPixel * beamDirection );
      }
    progress.CompletedPixel();
    }
}


template< class TInputImage, class TOutputImage, class TOperatorValue >
void
GradientBasedAngleOfIncidenceImageFilter< TInputImage,
  TOutputImage,
  TOperatorValue >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "UltrasoundProbeOrigin: "
     << this->m_UltrasoundProbeOrigin
     << std::endl;
}

} // end namespace itk

#endif
