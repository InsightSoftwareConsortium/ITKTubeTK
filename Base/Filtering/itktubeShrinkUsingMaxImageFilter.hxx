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
#ifndef __itkShrinkUsingMaxImageFilter_hxx
#define __itkShrinkUsingMaxImageFilter_hxx

#include "itktubeShrinkUsingMaxImageFilter.h"
#include "itkConstNeighborhoodIterator.h"

namespace itk {

namespace tube {

/**
 *
 */
template< class TInputImage, class TOutputImage >
ShrinkUsingMaxImageFilter< TInputImage, TOutputImage >
::ShrinkUsingMaxImageFilter( void )
{
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkUsingMaxImageFilter< TInputImage, TOutputImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_IndexImage.IsNotNull() )
    {
    os << indent << "Index Image: " << m_IndexImage << std::endl;
    }
  else
    {
    os << indent << "Index Image: NULL" << std::endl;
    }
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkUsingMaxImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
  ThreadIdType threadId )
{
  // Get the input and output pointers
  InputImageConstPointer inputPtr = this->GetInput();
  OutputImagePointer     outputPtr = this->GetOutput();

  typename TOutputImage::SizeType factorSize;
  for( unsigned int i = 0; i < TInputImage::ImageDimension; i++ )
    {
    factorSize[i] = this->GetShrinkFactors()[i];
    }

  // Define a few indices that will be used to transform from an input pixel
  // to an output pixel
  OutputIndexType  outputIndex;
  InputIndexType   inputIndex;
  OutputOffsetType offsetIndex;

  typename TOutputImage::PointType tempPoint;

  // Use this index to compute the offset everywhere in this class
  outputIndex = outputPtr->GetLargestPossibleRegion().GetIndex();

  // We wish to perform the following mapping of outputIndex to
  // inputIndex on all points in our region
  outputPtr->TransformIndexToPhysicalPoint( outputIndex, tempPoint );
  inputPtr->TransformPhysicalPointToIndex( tempPoint, inputIndex );

  // Given that the size is scaled by a constant factor eq:
  // inputIndex = outputIndex * factorSize
  // is equivalent up to a fixed offset which we now compute
  OffsetValueType zeroOffset = 0;
  for( unsigned int i = 0; i < TInputImage::ImageDimension; i++ )
    {
    offsetIndex[i] = inputIndex[i] - outputIndex[i] *
      this->GetShrinkFactors()[i];
    // It is plausible that due to small amounts of loss of numerical
    // precision that the offset it negaive, this would cause sampling
    // out of out region, this is insurance against that possibility
    offsetIndex[i] = vnl_math_max( zeroOffset, offsetIndex[i] );
    }

  // Support progress methods/callbacks
  ProgressReporter progress( this, threadId,
    outputRegionForThread.GetNumberOfPixels() );

  // Define/declare an iterator that will walk the output region for this
  // thread.
  typedef ImageRegionIteratorWithIndex< TOutputImage > OutputIterator;
  OutputIterator outIt( outputPtr, outputRegionForThread );

  typedef ImageRegionIteratorWithIndex< IndexImageType > IndexIterator;
  IndexIterator indexIt( m_IndexImage, outputRegionForThread );

  while( !outIt.IsAtEnd() )
    {
    // Determine the index and physical location of the output pixel
    outputIndex = outIt.GetIndex();

    // An optimized version of
    // outputPtr->TransformIndexToPhysicalPoint( outputIndex, tempPoint );
    // inputPtr->TransformPhysicalPointToIndex( tempPoint, inputIndex );
    // but without the rounding and precision issues
    inputIndex = outputIndex * factorSize + offsetIndex;

    ConstNeighborhoodIterator< InputImageType > it( factorSize,
      this->GetInput(), this->GetInput()->GetBufferedRegion() );

    // Set the iterator at the desired location
    it.SetLocation( inputIndex );

    // Walk the neighborhood
    typename TInputImage::PixelType maxValue = it.GetCenterPixel();
    typename TInputImage::IndexType maxValueIndex = it.GetIndex();
    typename itk::Vector< unsigned int, ImageDimension > maxValueIndexVector;
    for( unsigned int i = 0; i < ImageDimension; ++i )
      {
      maxValueIndexVector[ i ] = maxValueIndex[ i ];
      }
    typename TInputImage::PixelType value;
    const unsigned int size = it.Size();
    bool isInBounds = true;
    for ( unsigned int i = 0; i < size; ++i )
      {
      value = it.GetPixel( i, isInBounds );
      if( isInBounds && value > maxValue )
        {
        maxValue = value;
        maxValueIndex = it.GetIndex( i );
        for( unsigned int j = 0; j < ImageDimension; ++j )
          {
          maxValueIndexVector[ j ] = (int)maxValueIndex[ j ];
          }
        }
      }

    // Copy the input pixel to the output
    outIt.Set( maxValue );
    ++outIt;

    indexIt.Set( maxValueIndexVector );
    ++indexIt;

    progress.CompletedPixel();
    }
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkUsingMaxImageFilter< TInputImage, TOutputImage >
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
    factorSize[i] = this->GetShrinkFactors()[i];
    }

  OutputIndexType  outputIndex;
  InputIndexType   inputIndex, inputRequestedRegionIndex;
  OutputOffsetType offsetIndex;

  typename TInputImage::SizeType inputRequestedRegionSize;
  typename TOutputImage::PointType tempPoint;

  // Use this index to compute the offset everywhere in this class
  outputIndex = outputPtr->GetLargestPossibleRegion().GetIndex();

  // We wish to perform the following mapping of outputIndex to
  // inputIndex on all points in our region
  outputPtr->TransformIndexToPhysicalPoint(outputIndex, tempPoint);
  inputPtr->TransformPhysicalPointToIndex(tempPoint, inputIndex);

  // Given that the size is scaled by a constant factor eq:
  // inputIndex = outputIndex * factorSize
  // is equivalent up to a fixed offset which we now compute
  OffsetValueType zeroOffset = 0;
  for( unsigned int i = 0; i < TInputImage::ImageDimension; i++ )
    {
    // Modified from base class to request an expanded input region
    offsetIndex[i] = inputIndex[i] - outputIndex[i] *
      this->GetShrinkFactors()[i];
    offsetIndex[i] = vnl_math_max(zeroOffset, offsetIndex[i]);
    }

  for( unsigned int i = 0; i < TInputImage::ImageDimension; i++ )
    {
    inputRequestedRegionIndex[i] = ( outputRequestedRegionStartIndex[i] -
      1 ) * factorSize[i] + offsetIndex[i];
    }

  // Modified from base claas to request an expanded input region
  for( unsigned int i = 0; i < TInputImage::ImageDimension; i++ )
    {
    inputRequestedRegionSize[i] = ( outputRequestedRegionSize[i] + 2 ) *
      factorSize[i];
    }

  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion.SetIndex( inputRequestedRegionIndex );
  inputRequestedRegion.SetSize( inputRequestedRegionSize );
  inputRequestedRegion.Crop( inputPtr->GetLargestPossibleRegion() );

  inputPtr->SetRequestedRegion( inputRequestedRegion );
}


/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkUsingMaxImageFilter< TInputImage, TOutputImage >
::GenerateOutputInformation( void )
{
  Superclass::GenerateOutputInformation();

  OutputImagePointer     outputPtr = this->GetOutput();
  if ( !outputPtr )
    {
    return;
    }

  m_IndexImage = IndexImageType::New();
  m_IndexImage->SetRegions(
    outputPtr->GetLargestPossibleRegion() );
  m_IndexImage->CopyInformation( outputPtr );
  m_IndexImage->Allocate();
}

} // end namespace tube

} // end namespace itk

#endif
