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
/*=========================================================================
*
*  Copyright Insight Software Consortium
*
*  Licensed under the Apache License, Version 2.0 ( the "License" );
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*         http://www.apache.org/licenses/LICENSE-2.0.txt
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*=========================================================================*/
/*=========================================================================
*
*  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
*
*  Copyright ( c ) Ken Martin, Will Schroeder, Bill Lorensen
*
*  For complete copyright, license and disclaimer of warranty information
*  please refer to the NOTICE file at the top of the ITK source tree.
*
*=========================================================================*/
#ifndef __itktubeShrinkWithBlendingImageFilter_hxx
#define __itktubeShrinkWithBlendingImageFilter_hxx

#include "itktubeShrinkWithBlendingImageFilter.h"

#include "itkConstNeighborhoodIterator.h"

namespace itk {

namespace tube {

/**
 *
 */
template< class TInputImage, class TOutputImage >
ShrinkWithBlendingImageFilter< TInputImage, TOutputImage >
::ShrinkWithBlendingImageFilter( void )
{
  m_PointImage = NULL;

  m_Overlap.Fill( 0 );

  m_UseLog = false;

  m_BlendWithMax = true;
  m_BlendWithMean = false;
  m_BlendWithGaussianWeighting = false;

  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    m_ShrinkFactors[j] = 1;
    }
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImageFilter< TInputImage, TOutputImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Overlap" << m_Overlap << std::endl;

  if( m_PointImage.IsNotNull() )
    {
    os << indent << "Point Image: " << m_PointImage << std::endl;
    }
  else
    {
    os << indent << "Index Image: NULL" << std::endl;
    }
}

template< typename TInputImage, typename TOutputImage >
void
ShrinkWithBlendingImageFilter< TInputImage, TOutputImage >
::SetShrinkFactor(unsigned int i, unsigned int factor)
{
  if( m_ShrinkFactors[i] == factor )
    {
    return;
    }

  this->Modified();
  m_ShrinkFactors[i] = factor;
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
  ThreadIdType threadId )
{
  // Get the input and output pointers
  InputImageConstPointer inputPtr = this->GetInput();
  OutputImagePointer     outputPtr = this->GetOutput();

  typename TOutputImage::SizeType factorSize;
  for( unsigned int i = 0; i < TInputImage::ImageDimension; i++ )
    {
    factorSize[ i ] = this->GetShrinkFactors()[ i ];
    }

  // Define a few indices that will be used to transform from an input pixel
  // to an output pixel
  OutputIndexType  outputIndex;
  InputIndexType   inputIndex;
  InputIndexType                   inputWindowStartIndex;
  typename TInputImage::SizeType   inputWindowSize;

  typename TOutputImage::PointType tempPoint;

  // Use this index to compute the offset everywhere in this class
  outputIndex = outputPtr->GetLargestPossibleRegion().GetIndex();

  // We wish to perform the following mapping of outputIndex to
  // inputIndex on all points in our region
  outputPtr->TransformIndexToPhysicalPoint( outputIndex, tempPoint );
  inputPtr->TransformPhysicalPointToIndex( tempPoint, inputIndex );

  // Support progress methods/callbacks
  ProgressReporter progress( this, threadId,
    outputRegionForThread.GetNumberOfPixels() );

  // Define/declare an iterator that will walk the output region for this
  // thread.
  typedef ImageRegionIteratorWithIndex< TOutputImage > OutputIteratorType;
  OutputIteratorType outIt( outputPtr, outputRegionForThread );

  typedef ImageRegionIteratorWithIndex< PointImageType >
    PointIteratorType;
  PointIteratorType pointIt( m_PointImage, outputRegionForThread );

  typedef ImageRegionConstIteratorWithIndex< TInputImage >
    InputIteratorType;
  typename TInputImage::RegionType inputRegion;

  while( !outIt.IsAtEnd() )
    {
    // Determine the index and physical location of the output pixel
    outputIndex = outIt.GetIndex();

    outputPtr->TransformIndexToPhysicalPoint( outputIndex, tempPoint );
    inputPtr->TransformPhysicalPointToIndex( tempPoint, inputIndex );

    for( unsigned int i = 0; i < ImageDimension; ++i )
      {
      inputWindowStartIndex[ i ] = inputIndex[ i ] - factorSize[ i ] / 2
        - m_Overlap[ i ];
      inputWindowSize[ i ] = factorSize[ i ] + m_Overlap[ i ] * 2;
      }
    inputRegion.SetIndex( inputWindowStartIndex );
    inputRegion.SetSize( inputWindowSize );
    inputRegion.Crop( this->GetInput()->GetLargestPossibleRegion() );
    InputIteratorType it( this->GetInput(), inputRegion );

    // Walk the neighborhood
    typename TInputImage::PixelType value;
    if( m_BlendWithMax )
      {
      typename TInputImage::PixelType maxValue = it.Get();
      typename TInputImage::IndexType maxValueIndex = it.GetIndex();
      ++it;
      while( !it.IsAtEnd() )
        {
        value = it.Get();
        if( value > maxValue )
          {
          maxValue = value;
          maxValueIndex = it.GetIndex();
          }
        ++it;
        }

      // Copy the input pixel to the output
      outIt.Set( maxValue );
      ++outIt;

      typename TInputImage::PointType point;
      this->GetInput()->TransformIndexToPhysicalPoint( maxValueIndex,
        point );

      typename PointImageType::PixelType pointVector;
      for( unsigned int i = 0; i < ImageDimension; ++i )
        {
        pointVector[i] = point[i];
        }
      pointIt.Set( pointVector );
      ++pointIt;
      }
    else if( m_BlendWithMean )
      {
      double averageValue = 0;
      unsigned long int count = 0;
      if( m_UseLog )
        {
        while( !it.IsAtEnd() )
          {
          value = it.Get();
          averageValue += value * value;
          ++count;
          ++it;
          }
        if( count > 0 )
          {
          averageValue = std::sqrt( averageValue / count );
          }
        }
      else
        {
        while( !it.IsAtEnd() )
          {
          value = it.Get();
          averageValue += value;
          ++count;
          ++it;
          }
        if( count > 0 )
          {
          averageValue = averageValue / count;
          }
        }
      outIt.Set( averageValue );
      ++outIt;
      }
    else if( m_BlendWithGaussianWeighting )
      {
      typename TInputImage::IndexType valueIndex;
      double weight = 0;
      double averageValue = 0;
      double weightSum = 0;
      while( !it.IsAtEnd() )
        {
        value = it.Get();
        valueIndex = it.GetIndex();
        weight = 0;
        for( unsigned int i = 0; i < ImageDimension; ++i )
          {
          double dist = (valueIndex[i] - inputIndex[i] ) / factorSize[i];
          weight += (1.0 / (factorSize[i] * std::sqrt( 2 * vnl_math::pi )))
            * std::exp( -0.5 * dist * dist );
          }
        if( m_UseLog )
          {
          averageValue += weight * value * value;
          }
        else
          {
          averageValue += weight * value;
          }
        weightSum += weight;
        ++it;
        }

      if( weightSum > 0 )
        {
        if( m_UseLog )
          {
          outIt.Set( std::sqrt( averageValue / weightSum ) );
          }
        else
          {
          outIt.Set( averageValue / weightSum );
          }
        }
      else
        {
        outIt.Set( 0 );
        }
      ++outIt;
      }

    progress.CompletedPixel();
    }
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImageFilter< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
{
  // Call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // Get pointers to the input and output
  InputImagePointer  inputPtr = const_cast< TInputImage * >(
    this->GetInput() );
  OutputImagePointer outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // Compute the input requested region (size and start index)
  // Use the image transformations to insure an input requested region
  // that will provide the proper range
  const typename TOutputImage::SizeType & outputRequestedRegionSize =
    outputPtr->GetRequestedRegion().GetSize();
  const typename TOutputImage::IndexType & outputRequestedRegionStartIndex =
    outputPtr->GetRequestedRegion().GetIndex();

  // Convert the factor for convenient multiplication
  typename TOutputImage::SizeType factorSize;
  for( unsigned int i = 0; i < TInputImage::ImageDimension; i++ )
    {
    factorSize[ i ] = this->GetShrinkFactors()[ i ];
    }

  InputIndexType                   inputRequestedRegionStartIndex;
  typename TInputImage::SizeType   inputRequestedRegionSize;
  typename TOutputImage::PointType tempPoint;

  outputPtr->TransformIndexToPhysicalPoint( outputRequestedRegionStartIndex,
    tempPoint);
  inputPtr->TransformPhysicalPointToIndex( tempPoint,
    inputRequestedRegionStartIndex );
  for( unsigned int i = 0; i < TInputImage::ImageDimension; i++ )
    {
    inputRequestedRegionStartIndex[ i ] =
      inputRequestedRegionStartIndex[ i ] - factorSize[i]
      - m_Overlap[ i ];
    }

  // Modified from base claas to request an expanded input region
  for( unsigned int i = 0; i < TInputImage::ImageDimension; i++ )
    {
    inputRequestedRegionSize[ i ] = ( outputRequestedRegionSize[ i ] + 2 ) *
      factorSize[ i ] + 2 * m_Overlap[ i ];
    }

  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion.SetIndex( inputRequestedRegionStartIndex );
  inputRequestedRegion.SetSize( inputRequestedRegionSize );
  inputRequestedRegion.Crop( inputPtr->GetLargestPossibleRegion() );

  inputPtr->SetRequestedRegion( inputRequestedRegion );
}


/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImageFilter< TInputImage, TOutputImage >
::GenerateOutputInformation( void )
{
  Superclass::GenerateOutputInformation();

  InputImageConstPointer inputPtr = this->GetInput();
  OutputImagePointer     outputPtr = this->GetOutput();
  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  const typename TInputImage::SpacingType & inputSpacing =
    inputPtr->GetSpacing();
  const typename TInputImage::SizeType & inputSize =
    inputPtr->GetLargestPossibleRegion().GetSize();
  const typename TInputImage::IndexType & inputStartIndex =
    inputPtr->GetLargestPossibleRegion().GetIndex();

  typename TOutputImage::SpacingType outputSpacing;
  typename TOutputImage::SizeType outputSize;
  typename TOutputImage::IndexType outputStartIndex;

  for( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
    {
    outputSpacing[i] = inputSpacing[i]
      * (double)this->GetShrinkFactors()[i];

    // Round down so that all output pixels fit input input region
    outputSize[i] = static_cast<SizeValueType>( std::floor(
      (double)inputSize[i] / (double)this->GetShrinkFactors()[i] ) );

    if( outputSize[i] < 1 )
      {
      outputSize[i] = 1;
      }

    outputStartIndex[i] = inputStartIndex[i];
    }

  outputPtr->SetSpacing( outputSpacing );
  outputPtr->SetDirection( inputPtr->GetDirection() );

  // Compute origin offset
  // The physical center's of the input and output should be the same
  ContinuousIndex< SpacePrecisionType, TOutputImage::ImageDimension >
    inputCenterIndex;
  ContinuousIndex< SpacePrecisionType, TOutputImage::ImageDimension >
    outputCenterIndex;
  for( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
    {
    inputCenterIndex[i] = inputStartIndex[i] + ( inputSize[i] - 1 ) / 2.0;
    outputCenterIndex[i] = outputStartIndex[i] + ( outputSize[i] - 1 )
      / 2.0;
    }

  typename TOutputImage::PointType inputCenterPoint;
  typename TOutputImage::PointType outputCenterPoint;
  inputPtr->TransformContinuousIndexToPhysicalPoint(inputCenterIndex,
    inputCenterPoint);
  outputPtr->TransformContinuousIndexToPhysicalPoint(outputCenterIndex,
    outputCenterPoint);

  const typename TOutputImage::PointType & inputOrigin =
    inputPtr->GetOrigin();
  typename TOutputImage::PointType outputOrigin;
  outputOrigin = inputOrigin + (inputCenterPoint - outputCenterPoint);
  outputPtr->SetOrigin(outputOrigin);

  // Set region
  typename TOutputImage::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize(outputSize);
  outputLargestPossibleRegion.SetIndex(outputStartIndex);

  outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);

  m_PointImage = PointImageType::New();
  m_PointImage->SetRegions(
    outputPtr->GetLargestPossibleRegion() );
  m_PointImage->CopyInformation( outputPtr );
  m_PointImage->Allocate();

}

} // end namespace tube

} // end namespace itk

#endif
