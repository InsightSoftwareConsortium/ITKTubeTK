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

#ifndef __itkAcousticImpulseResponseImageFilter_hxx
#define __itkAcousticImpulseResponseImageFilter_hxx

#include "itkAcousticImpulseResponseImageFilter.h"

#include <itkGradientMagnitudeImageFilter.h>

namespace itk
{

template< class TInputImage, class TOutputImage, class TOperatorValue >
AcousticImpulseResponseImageFilter< TInputImage, TOutputImage, TOperatorValue >
::AcousticImpulseResponseImageFilter( void )
  : m_AngleDependence( 1.0 )
{
  this->SetNumberOfRequiredInputs( 2 );

  typedef GradientMagnitudeImageFilter< OperatorImageType, OperatorImageType >
    DefaultGradientMagnitudeFilterType;
  this->m_GradientMagnitudeFilter = DefaultGradientMagnitudeFilterType::New();

  this->m_CastImageFilter = CastImageFilterType::New();
}


template< class TInputImage, class TOutputImage, class TOperatorValue >
void
AcousticImpulseResponseImageFilter< TInputImage, TOutputImage, TOperatorValue >
::BeforeThreadedGenerateData( void )
{
  this->m_CastImageFilter->SetInput( this->GetInput() );
  this->m_GradientMagnitudeFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->m_GradientMagnitudeFilter->SetInput( this->m_CastImageFilter->GetOutput() );
  this->m_GradientMagnitudeFilter->Update();
}


template< class TInputImage, class TOutputImage, class TOperatorValue >
void
AcousticImpulseResponseImageFilter< TInputImage, TOutputImage, TOperatorValue >
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
  ThreadIdType threadId )
{
  const InputImageType * input = this->GetInput( 0 );
  const OperatorImageType * angleOfIncidence = this->GetInput( 1 );
  const OperatorImageType * gradientMagnitude =
    this->m_GradientMagnitudeFilter->GetOutput();
  OutputImageType * output = this->GetOutput();

  typedef ImageRegionConstIterator< InputImageType > InputIteratorType;
  typedef ImageRegionIterator< OutputImageType >     OutputIteratorType;

  typedef ImageRegionConstIterator< OperatorImageType >
    AngleOfIncidenceIteratorType;
  typedef ImageRegionConstIterator< OperatorImageType >
    GradientMagnitudeIteratorType;

  InputIteratorType inputIt( input, outputRegionForThread );
  AngleOfIncidenceIteratorType angleOfIncidenceIt( angleOfIncidence,
    outputRegionForThread );
  GradientMagnitudeIteratorType gradientMagnitudeIt( gradientMagnitude,
    outputRegionForThread );
  OutputIteratorType outputIt( output, outputRegionForThread );

  ProgressReporter progress( this, threadId,
    outputRegionForThread.GetNumberOfPixels() );

  if( this->m_AngleDependence == 1.0 )  // avoid the pow call
    {
    for( inputIt.GoToBegin(), angleOfIncidenceIt.GoToBegin(),
          gradientMagnitudeIt.GoToBegin(), outputIt.GoToBegin();
         !outputIt.IsAtEnd();
         ++inputIt, ++angleOfIncidenceIt, ++gradientMagnitudeIt, ++outputIt )
      {
      outputIt.Set( angleOfIncidenceIt.Get() *
        gradientMagnitudeIt.Get() / ( 2.0 * inputIt.Get() ) );
      }
    }
  else
    {
    for( inputIt.GoToBegin(), angleOfIncidenceIt.GoToBegin(),
          gradientMagnitudeIt.GoToBegin(), outputIt.GoToBegin();
         !outputIt.IsAtEnd();
         ++inputIt, ++angleOfIncidenceIt, ++gradientMagnitudeIt, ++outputIt )
      {
      outputIt.Set( static_cast< typename OutputImageType::PixelType >( 
        std::pow( static_cast< OperatorValueType >( angleOfIncidenceIt.Get() ),
        static_cast< OperatorValueType >( this->m_AngleDependence *
        gradientMagnitudeIt.Get() / ( 2.0 * inputIt.Get() ) ) ) ) );
      }
    }
}


template< class TInputImage, class TOutputImage, class TOperatorValue >
void
AcousticImpulseResponseImageFilter< TInputImage, TOutputImage, TOperatorValue >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "AngleDependence: "
     << this->m_AngleDependence
     << std::endl;
}

} // End namespace itk

#endif // End !defined( __itkAcousticImpulseResponseImageFilter_hxx )
