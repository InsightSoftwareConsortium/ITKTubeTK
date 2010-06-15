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
#ifndef __itkJointHistogramImageFunction_txx
#define __itkJointHistogramImageFunction_txx

#include "itkMinimumMaximumImageCalculator.h"

#include "itkJointHistogramImageFunction.h"

namespace itk
{

/**
 * Set the input Image
 */
template <class TInputImage, class TCoordRep>
JointHistogramImageFunction<TInputImage,TCoordRep>
::JointHistogramImageFunction()
{
  m_InputMask = 0;
  m_FeatureWidth = 20;
  m_SumHistogram = 0;
  m_SumOfSquaresHistogram = 0;
  m_NumberOfSamples = 0;
  m_NumberOfComputedSamples = 0;
  this->SetHistogramSize( 20 );
}

template <class TInputImage, class TCoordRep>
void 
JointHistogramImageFunction<TInputImage,TCoordRep>
::SetHistogramSize( const unsigned int& size )
{
  m_HistogramSize = size;
  typename HistogramType::IndexType start;
  typename HistogramType::SizeType histSize;
  start[0] = start[1] = 0;
  histSize[0] = histSize[1] = size;
  typename HistogramType::RegionType region;
  region.SetIndex( start );
  region.SetSize( histSize );

  m_SumHistogram = HistogramType::New();
  m_SumHistogram->SetRegions( region );
  m_SumHistogram->Allocate();
  m_SumHistogram->FillBuffer( 0 );

  m_SumOfSquaresHistogram = HistogramType::New();
  m_SumOfSquaresHistogram->SetRegions( region );
  m_SumOfSquaresHistogram->Allocate();
  m_SumOfSquaresHistogram->FillBuffer( 0 );

  m_NumberOfComputedSamples = m_NumberOfSamples = 0;
}

template <class TInputImage, class TCoordRep>
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::SetInputImage( const InputImageType * ptr )
{
  this->Superclass::SetInputImage( ptr );

  typedef itk::MinimumMaximumImageCalculator<InputImageType> MinMaxType;
  
  
}

template <class TInputImage, class TCoordRep>
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::SetInputMask( const typename InputImageType::Pointer mask )
{
    m_InputMask = mask;
}

template <class TInputImage, class TCoordRep>
double 
JointHistogramImageFunction<TInputImage,TCoordRep>
::EvaluateAtIndex( const IndexType & index ) const
{
  if( m_NumberOfComputedSamples < m_NumberOfSamples )
    {

    m_NumberOfComputedSamples = m_NumberOfSamples;
    }
  PixelType value = this->GetInputImage()->GetPixel( index );
  return static_cast<double>( value );
}

template <class TInputImage, class TCoordRep>
double 
JointHistogramImageFunction<TInputImage,TCoordRep>
::PrecomputeAtIndex( const IndexType & index ) const
{
  PixelType value = this->GetInputImage()->GetPixel( index );
  return static_cast<double>( value );
  ++m_NumberOfSamples;
}

template <class TInputImage, class TCoordRep>
void 
JointHistogramImageFunction<TInputImage,TCoordRep>
::ComputeMeanAndStandardDeviation()
{
  
}

template <class TInputImage, class TCoordRep>
void 
JointHistogramImageFunction<TInputImage,TCoordRep>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "m_FeatureWidth = " << m_FeatureWidth << std::endl;
}

template <class TInputImage, class TCoordRep>
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::ComputeHistogramAtIndex( const IndexType& index, typename HistogramType::Pointer hist ) const
{
}

} // end namespace itk

#endif
