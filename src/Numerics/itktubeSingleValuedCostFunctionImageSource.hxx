/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeSingleValuedCostFunctionImageSource_hxx
#define __itktubeSingleValuedCostFunctionImageSource_hxx


#include <itkImageRegionIteratorWithIndex.h>
#include <itkProgressReporter.h>

namespace itk
{

namespace tube
{

template< class TCostFunction, unsigned int VNumberOfParameters >
SingleValuedCostFunctionImageSource< TCostFunction, VNumberOfParameters >
::SingleValuedCostFunctionImageSource():
  m_ParametersLowerBound( NumberOfParameters ),
  m_ParametersUpperBound( NumberOfParameters ),
  m_ParametersStep( NumberOfParameters )
{
  this->m_ParametersLowerBound.Fill( 0.0 );
  this->m_ParametersUpperBound.Fill( 10.0 );
  this->m_ParametersStep.Fill( 1.0 );

  //Use the ITKv4 Threading Model (call ThreadedGenerateData instead of DynamicThreadedGenerateData)
  this->DynamicMultiThreadingOff();
}

template< class TCostFunction, unsigned int VNumberOfParameters >
void
SingleValuedCostFunctionImageSource< TCostFunction, VNumberOfParameters >
::SetParametersLowerBound( const ParametersType & lowerBound )
{
  const SizeValueType lowerBoundSize = lowerBound.GetSize();
  if( lowerBoundSize != NumberOfParameters )
    {
    itkExceptionMacro( << "The ParametersLowerBound size: " << lowerBoundSize
      << " does not equal the expected number of parameters: "
      << NumberOfParameters << std::endl );
    }
  for( SizeValueType ii = 0; ii < lowerBoundSize; ++ii )
    {
    if( this->m_ParametersLowerBound[ii] != lowerBound[ii] )
      {
      this->m_ParametersLowerBound = lowerBound;
      this->Modified();
      break;
      }
    }
}


template< class TCostFunction, unsigned int VNumberOfParameters >
typename SingleValuedCostFunctionImageSource< TCostFunction,
  VNumberOfParameters >::ParametersType
SingleValuedCostFunctionImageSource< TCostFunction, VNumberOfParameters >
::GetParametersLowerBound() const
{
  return this->m_ParametersLowerBound;
}


template< class TCostFunction, unsigned int VNumberOfParameters >
void
SingleValuedCostFunctionImageSource< TCostFunction, VNumberOfParameters >
::SetParametersUpperBound( const ParametersType & upperBound )
{
  const SizeValueType upperBoundSize = upperBound.GetSize();
  if( upperBoundSize != NumberOfParameters )
    {
    itkExceptionMacro( << "The ParametersUpperBound size: " << upperBoundSize
      << " does not equal the expected number of parameters: "
      << NumberOfParameters << std::endl );
    }
  for( SizeValueType ii = 0; ii < upperBoundSize; ++ii )
    {
    if( this->m_ParametersUpperBound[ii] != upperBound[ii] )
      {
      this->m_ParametersUpperBound = upperBound;
      this->Modified();
      break;
      }
    }
}


template< class TCostFunction, unsigned int VNumberOfParameters >
typename SingleValuedCostFunctionImageSource< TCostFunction,
  VNumberOfParameters >::ParametersType
SingleValuedCostFunctionImageSource< TCostFunction, VNumberOfParameters >
::GetParametersUpperBound() const
{
  return this->m_ParametersUpperBound;
}


template< class TCostFunction, unsigned int VNumberOfParameters >
void
SingleValuedCostFunctionImageSource< TCostFunction, VNumberOfParameters >
::SetParametersStep( const ParametersType & step )
{
  const SizeValueType stepSize = step.GetSize();
  if( stepSize != NumberOfParameters )
    {
    itkExceptionMacro( << "The ParametersStep size: " << stepSize
      << " does not equal the expected number of parameters: "
      << NumberOfParameters << std::endl );
    }
  for( SizeValueType ii = 0; ii < stepSize; ++ii )
    {
    if( this->m_ParametersStep[ii] != step[ii] )
      {
      this->m_ParametersStep = step;
      this->Modified();
      break;
      }
    }
}


template< class TCostFunction, unsigned int VNumberOfParameters >
typename SingleValuedCostFunctionImageSource< TCostFunction,
  VNumberOfParameters >::ParametersType
SingleValuedCostFunctionImageSource< TCostFunction, VNumberOfParameters >
::GetParametersStep() const
{
  return this->m_ParametersStep;
}


template< class TCostFunction, unsigned int VNumberOfParameters >
void
SingleValuedCostFunctionImageSource< TCostFunction, VNumberOfParameters >
::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();

  OutputImageType * outputImage = this->GetOutput();

  typename OutputImageType::IndexType index;
  index.Fill( 0 );
  typename OutputImageType::SizeType size;
  typename OutputImageType::PointType origin;
  typename OutputImageType::SpacingType spacing;
  for( unsigned int ii = 0; ii < NumberOfParameters; ++ii )
    {
    origin[ii] = this->m_ParametersLowerBound[ii];
    spacing[ii] = this->m_ParametersStep[ii];
    size[ii] = static_cast< SizeValueType >(
      ( this->m_ParametersUpperBound[ii] -
      this->m_ParametersLowerBound[ii] ) / this->m_ParametersStep[ii] ) + 1;
    }
  typename OutputImageType::RegionType region;
  region.SetIndex( index );
  region.SetSize( size );

  outputImage->SetLargestPossibleRegion( region );
  outputImage->SetSpacing( spacing );
  outputImage->SetOrigin( origin );
  typename OutputImageType::DirectionType direction;
  direction.SetIdentity();
  outputImage->SetDirection( direction );
}


template< class TCostFunction, unsigned int VNumberOfParameters >
void
SingleValuedCostFunctionImageSource< TCostFunction, VNumberOfParameters >
::BeforeThreadedGenerateData()
{
  if( this->m_CostFunction.IsNull() )
    {
    itkExceptionMacro( << "Cost function must be assigned!" );
    }
}


template< class TCostFunction, unsigned int VNumberOfParameters >
void
SingleValuedCostFunctionImageSource< TCostFunction, VNumberOfParameters >
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
  ThreadIdType threadId )
{
  OutputImageType * outputImage = this->GetOutput( 0 );

  ProgressReporter progress( this, threadId,
    outputRegionForThread.GetNumberOfPixels() );

  typedef ImageRegionIteratorWithIndex< OutputImageType > ImageIteratorType;
  ImageIteratorType imageIt( outputImage, outputRegionForThread );
  for( imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt )
    {
    const typename OutputImageType::IndexType index = imageIt.GetIndex();
    typename OutputImageType::PointType point;
    outputImage->TransformIndexToPhysicalPoint( index, point );
    ParametersType parameters( NumberOfParameters );
    parameters.SetDataSameSize( point.GetDataPointer() );
    const MeasureType measure = this->m_CostFunction->GetValue(
      parameters );
    imageIt.Set( static_cast< typename OutputImageType::PixelType >(
      measure ) );
    progress.CompletedPixel();
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeSingleValuedCostFunctionImageSource_hxx )
