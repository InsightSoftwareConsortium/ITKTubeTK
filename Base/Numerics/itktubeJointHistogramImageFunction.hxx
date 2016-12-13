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

#ifndef __itktubeJointHistogramImageFunction_hxx
#define __itktubeJointHistogramImageFunction_hxx

#include "itktubeJointHistogramImageFunction.h"

#include <itkAddImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkMultiplyImageFilter.h>
#include <itkSqrtImageFilter.h>
#include <itkSquareImageFilter.h>
#include <itkSubtractImageFilter.h>

namespace itk
{

namespace tube
{

/**
 * Set the input Image
 */
template< class TInputImage, class TCoordRep >
JointHistogramImageFunction<TInputImage,TCoordRep>
::JointHistogramImageFunction( void )
{
  m_InputMask = NULL;
  m_FeatureWidth = 40;
  m_StdevBase = 0.001;
  m_Histogram = NULL;
  m_SumHistogram = NULL;
  m_SumOfSquaresHistogram = NULL;
  m_MeanHistogram = NULL;
  m_StandardDeviationHistogram = NULL;
  m_NumberOfSamples = 0;
  m_NumberOfComputedSamples = 0;
  m_ImageMin = 0;
  m_ImageMax = 0;
  m_ImageStep = 0;
  m_MaskMin = 0;
  m_MaskMax = 0;
  m_MaskStep = 0;
  m_ForceDiagonalHistogram = false;
  m_HistogramSize = 40;

  this->SetHistogramSize( m_HistogramSize );
}

template< class TInputImage, class TCoordRep >
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::SetHistogramSize( const unsigned int & size )
{
  m_HistogramSize = size;
  typename HistogramType::IndexType start;
  typename HistogramType::SizeType histSize;
  start[0] = start[1] = 0;
  histSize[0] = histSize[1] = size;
  typename HistogramType::RegionType region;
  region.SetIndex( start );
  region.SetSize( histSize );

  m_Histogram = HistogramType::New();
  m_Histogram->SetRegions( region );
  m_Histogram->Allocate();
  m_Histogram->FillBuffer( 0 );

  m_SumHistogram = HistogramType::New();
  m_SumHistogram->SetRegions( region );
  m_SumHistogram->Allocate();
  m_SumHistogram->FillBuffer( 0 );

  m_SumOfSquaresHistogram = HistogramType::New();
  m_SumOfSquaresHistogram->SetRegions( region );
  m_SumOfSquaresHistogram->Allocate();
  m_SumOfSquaresHistogram->FillBuffer( 0 );

  m_MeanHistogram = HistogramType::New();
  m_MeanHistogram->SetRegions( region );
  m_MeanHistogram->Allocate();
  m_MeanHistogram->FillBuffer( 0 );

  m_StandardDeviationHistogram = HistogramType::New();
  m_StandardDeviationHistogram->SetRegions( region );
  m_StandardDeviationHistogram->Allocate();
  m_StandardDeviationHistogram->FillBuffer( 0 );


  m_NumberOfComputedSamples = m_NumberOfSamples = 0;

  m_ImageStep = ( m_ImageMax - m_ImageMin ) / m_HistogramSize;
  m_MaskStep = ( m_MaskMax - m_MaskMin ) / m_HistogramSize;
}

template< class TInputImage, class TCoordRep >
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::SetInputImage( const InputImageType * ptr )
{
  this->Superclass::SetInputImage( ptr );

  typedef itk::MinimumMaximumImageCalculator<InputImageType> MinMaxType;
  typename MinMaxType::Pointer calculator = MinMaxType::New();
  calculator->SetImage( this->GetInputImage() );
  calculator->Compute();
  m_ImageMin = calculator->GetMinimum();
  m_ImageMax = calculator->GetMaximum();

  m_ImageStep = ( m_ImageMax - m_ImageMin ) / m_HistogramSize;
}

template< class TInputImage, class TCoordRep >
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::SetInputMask( const typename InputImageType::Pointer mask )
{
  m_InputMask = mask;

  typedef itk::MinimumMaximumImageCalculator<InputImageType> MinMaxType;
  typename MinMaxType::Pointer calculator = MinMaxType::New();
  calculator->SetImage( m_InputMask );
  calculator->Compute();
  m_MaskMin = calculator->GetMinimum();
  m_MaskMax = calculator->GetMaximum();

  m_MaskStep = ( m_MaskMax - m_MaskMin ) / m_HistogramSize;
}

template< class TInputImage, class TCoordRep >
double
JointHistogramImageFunction<TInputImage,TCoordRep>
::EvaluateAtIndex( const IndexType & index ) const
{
  try
    {
    if( m_NumberOfComputedSamples < m_NumberOfSamples )
      {
      this->ComputeMeanAndStandardDeviation();
      m_NumberOfComputedSamples = m_NumberOfSamples;
      }
    return this->ComputeZScoreAtIndex( index );
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception thrown: " << e << std::endl;
    return 0;
    }
}

template< class TInputImage, class TCoordRep >
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::PrecomputeAtIndex( const IndexType & index )
{
  typename HistogramType::Pointer hist;
  hist = this->ComputeHistogramAtIndex( index, false );

  itk::ImageRegionIterator< HistogramType > iterHist( hist,
    hist->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< HistogramType > iterSum( m_SumHistogram,
    m_SumHistogram->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< HistogramType > iterSumSquare( 
    m_SumOfSquaresHistogram,
    m_SumOfSquaresHistogram->GetLargestPossibleRegion() );
  while( !iterHist.IsAtEnd() )
    {
    double tf = iterHist.Get();
    iterSum.Set( iterSum.Get() + tf );
    iterSumSquare.Set( iterSumSquare.Get() + tf*tf );
    ++iterHist;
    ++iterSum;
    ++iterSumSquare;
    }

  ++m_NumberOfSamples;
}

template< class TInputImage, class TCoordRep >
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::ComputeMeanAndStandardDeviation( void ) const
{
  typedef itk::DiscreteGaussianImageFilter< HistogramType,
          HistogramType > SmootherType;

  SmootherType::Pointer sumSmoother = SmootherType::New();
  sumSmoother->SetInput( m_SumHistogram );
  sumSmoother->SetVariance( 2 );
  sumSmoother->SetUseImageSpacing( false );
  sumSmoother->Update();
  m_SumHistogram = sumSmoother->GetOutput();

  SmootherType::Pointer sumSquaresSmoother = SmootherType::New();
  sumSquaresSmoother->SetInput( m_SumOfSquaresHistogram );
  sumSquaresSmoother->SetVariance( 2 );
  sumSquaresSmoother->SetUseImageSpacing( false );
  sumSquaresSmoother->Update();
  m_SumOfSquaresHistogram = sumSquaresSmoother->GetOutput();

    // Calculate the mean
  typedef itk::ImageRegionIterator<HistogramType>  HistIteratorType;
  HistIteratorType sumItr( m_SumHistogram,
    m_SumHistogram->GetLargestPossibleRegion() );
  HistIteratorType sumOfSquaresItr( m_SumOfSquaresHistogram,
    m_SumOfSquaresHistogram->GetLargestPossibleRegion() );
  HistIteratorType meanItr( m_MeanHistogram,
    m_MeanHistogram->GetLargestPossibleRegion() );
  HistIteratorType stdItr( m_StandardDeviationHistogram,
    m_StandardDeviationHistogram->GetLargestPossibleRegion() );
  if( m_NumberOfSamples != 0 )
    {
    while( !meanItr.IsAtEnd() )
      {
      meanItr.Set( sumItr.Get() / m_NumberOfSamples );
      stdItr.Set( std::sqrt( vnl_math_abs( 
        sumOfSquaresItr.Get() / m_NumberOfSamples -
        meanItr.Get() * meanItr.Get() ) ) );
      ++sumItr;
      ++sumOfSquaresItr;
      ++meanItr;
      ++stdItr;
      }
    }
}

template< class TInputImage, class TCoordRep >
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  if( m_InputMask.IsNotNull() )
    {
    os << indent << "m_InputMask = " << m_InputMask << std::endl;
    }
  else
    {
    os << indent << "m_InputMask = NULL" << std::endl;
    }
  os << indent << "m_Histogram = " << m_Histogram << std::endl;
  os << indent << "m_SumHistogram = " << m_SumHistogram << std::endl;
  os << indent << "m_SumOfSquaresHistogram = "
     << m_SumOfSquaresHistogram << std::endl;
  os << indent << "m_MeanHistogram = " << m_MeanHistogram << std::endl;
  os << indent << "m_StandardDeviationHistogram = "
     << m_StandardDeviationHistogram << std::endl;
  os << indent << "m_FeatureWidth = " << m_FeatureWidth << std::endl;
  os << indent << "m_StdevBase = " << m_StdevBase << std::endl;
  os << indent << "m_HistogramSize = " << m_HistogramSize << std::endl;
  os << indent << "m_NumberOfSamples = " << m_NumberOfSamples << std::endl;
  os << indent << "m_NumberOfComputedSamples = "
     << m_NumberOfComputedSamples << std::endl;
  os << indent << "m_ImageMin = " << m_ImageMin << std::endl;
  os << indent << "m_ImageMax = " << m_ImageMax << std::endl;
  os << indent << "m_ImageStep = " << m_ImageStep << std::endl;
  os << indent << "m_MaskMin = " << m_MaskMin << std::endl;
  os << indent << "m_MaskMax = " << m_MaskMax << std::endl;
  os << indent << "m_MaskStep = " << m_MaskStep << std::endl;
  os << indent << "m_ForceDiagonalHistogram = " << m_MaskStep << std::endl;
}

template< class TInputImage, class TCoordRep >
double
JointHistogramImageFunction<TInputImage,TCoordRep>
::ComputeZScoreAtIndex( const IndexType & index ) const
{
  typename HistogramType::Pointer hist;
  hist = this->ComputeHistogramAtIndex( index, true );

  typedef itk::ImageRegionConstIterator<HistogramType>  HistIteratorType;
  HistIteratorType histItr( hist, hist->GetLargestPossibleRegion() );
  HistIteratorType meanItr( m_MeanHistogram,
                             m_MeanHistogram->GetLargestPossibleRegion() );
  HistIteratorType stdItr( m_StandardDeviationHistogram,
                            m_StandardDeviationHistogram->
                            GetLargestPossibleRegion() );
  double val = 0;
  unsigned int histCount = 0;
  while( !histItr.IsAtEnd() )
    {
    double t = histItr.Get();
    double m = meanItr.Get();
    double s = stdItr.Get();
    s += m_StdevBase;
    if( s != 0 )
      {
      val += vnl_math_abs( t - m ) / s;
      ++histCount;
      }
    ++histItr;
    ++meanItr;
    ++stdItr;
    }

  if( histCount > 0 )
    {
    return val / histCount;
    }
  else
    {
    return 0;
    }
}

template< class TInputImage, class TCoordRep >
itk::Image<float,2>::Pointer &
JointHistogramImageFunction<TInputImage,TCoordRep>
::ComputeHistogramAtIndex( const IndexType & index, bool blur ) const
{
  m_Histogram->FillBuffer( 0 );

  typedef itk::ImageRegionConstIterator<InputImageType> ConstIteratorType;

  typename InputImageType::IndexType minIndex;
  typename InputImageType::IndexType maxIndex;
  minIndex = m_InputMask->GetLargestPossibleRegion().GetIndex();
  maxIndex = minIndex + m_InputMask->GetLargestPossibleRegion().GetSize();
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    maxIndex[i] -= 1;
    }

  typename InputImageType::RegionType region;
  IndexType origin;
  typename InputImageType::SizeType size;
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    origin[i] = index[i] - ( m_FeatureWidth * 0.5 );
    size[i] = m_FeatureWidth;
    if( origin[i] < minIndex[i] )
      {
      size[i] -= ( minIndex[i] - origin[i] );
      origin[i] = minIndex[i];
      }
    if( origin[i] > maxIndex[i] )
      {
      origin[i] = maxIndex[i];
      size[i] = 1;
      }
    if( ( int )( origin[i]+size[i]-1 ) > maxIndex[i] )
      {
      size[i] = ( maxIndex[i] - origin[i] ) + 1;
      }
    }
  region.SetSize( size );
  region.SetIndex( origin );

  ConstIteratorType inputItr( this->GetInputImage(), region );
  ConstIteratorType maskItr( m_InputMask, region );
  while( !inputItr.IsAtEnd() )
    {
    typename HistogramType::IndexType cur;
    cur[0] = ( ( inputItr.Get() - m_ImageMin ) / m_ImageStep );
    cur[1] = ( ( maskItr.Get() - m_MaskMin ) / m_MaskStep );
    if( cur[0] > ( int )( m_HistogramSize ) - 1 )
      {
      cur[0] = ( int )( m_HistogramSize ) - 1;
      }
    if( cur[1] > ( int )( m_HistogramSize ) - 1 )
      {
      cur[1] = ( int )( m_HistogramSize ) - 1;
      }
    if( cur[0] < 0 )
      {
      cur[0] = 0;
      }
    if( cur[1] < 0 )
      {
      cur[1] = 0;
      }

    m_Histogram->SetPixel( cur, m_Histogram->GetPixel( cur ) + 1 );

    ++inputItr;
    ++maskItr;
    }

  if( blur )
    {
    typedef itk::DiscreteGaussianImageFilter< HistogramType,
            HistogramType > SmootherType;
    SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetInput( m_Histogram );
    smoother->SetVariance( 2 );
    smoother->SetUseImageSpacing( false );
    smoother->Update();
    m_Histogram = smoother->GetOutput();
    }

  if( m_ForceDiagonalHistogram )
    {
    typename HistogramType::IndexType cur;
    for( unsigned int i=0; i<m_HistogramSize; i++ )
      {
      cur[0] = i;
      unsigned int maxJ = 0;
      double maxJV = 0;
      for( unsigned int j=0; j<m_HistogramSize; j++ )
        {
        cur[1] = j;
        if( m_Histogram->GetPixel( cur ) > maxJV )
          {
          maxJV = m_Histogram->GetPixel( cur );
          maxJ = j;
          }
        }
      if( maxJV > 0 )
        {
        if( ( int )maxJ > ( int )m_HistogramSize/2 )
          {
          typename HistogramType::IndexType src;
          cur[1] = 0;
          src[0] = cur[0];
          src[1] = ( int )maxJ - ( int )m_HistogramSize/2;
          src[1] = cur[1] + src[1];
          while( cur[1] >= 0 && cur[1] < ( int )( m_HistogramSize ) )
            {
            if( src[1] < 0 || src[1] >= ( int )( m_HistogramSize ) )
              {
              m_Histogram->SetPixel( cur, 0 );
              }
            else
              {
              m_Histogram->SetPixel( cur, m_Histogram->GetPixel( src ) );
              }
            ++cur[1];
            ++src[1];
            }
          }
        else
          {
          typename HistogramType::IndexType src;
          cur[1] = m_HistogramSize-1;
          src[0] = cur[0];
          src[1] = ( int )maxJ - ( int )m_HistogramSize/2;
          src[1] = cur[1] + src[1];
          while( cur[1] >= 0 && cur[1] < ( int )( m_HistogramSize ) )
            {
            if( src[1] < 0 || src[1] >= ( int )( m_HistogramSize ) )
              {
              m_Histogram->SetPixel( cur, 0 );
              }
            else
              {
              m_Histogram->SetPixel( cur, m_Histogram->GetPixel( src ) );
              }
            --cur[1];
            --src[1];
            }
          }
        }
      }
    }

  return m_Histogram;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeJointHistogramImageFunction_hxx )
