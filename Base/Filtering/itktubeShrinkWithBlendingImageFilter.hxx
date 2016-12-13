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
#ifndef __itktubeShrinkWithBlendingImageFilter_hxx
#define __itktubeShrinkWithBlendingImageFilter_hxx

#include "itktubeShrinkWithBlendingImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk {

namespace tube {

/**
 *
 */
template< class TInputImage, class TOutputImage >
ShrinkWithBlendingImageFilter< TInputImage, TOutputImage >
::ShrinkWithBlendingImageFilter( void )
{
  m_InputMipPointImage = ITK_NULLPTR;
  m_OutputMipPointImage = ITK_NULLPTR;

  m_Overlap.Fill( 0 );

  m_UseLog = false;

  m_BlendWithMax = true;
  m_BlendWithMean = false;
  m_BlendWithGaussianWeighting = false;

  m_DefaultShrinkFactor = 1.0;
  m_DefaultNewSize = 0.0;
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    m_ShrinkFactors[j] = 1;
    m_NewSize[j] = 0;
    m_InternalShrinkFactors = 1;
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

  os << indent << "Overlap:" << m_Overlap << std::endl;
  os << indent << "ShrinkFactors:" << m_ShrinkFactors << std::endl;
  os << indent << "NewSize:" << m_NewSize << std::endl;
  os << indent << "BlendWithMean:" << m_BlendWithMean << std::endl;
  os << indent << "BlendWithMax:"<< m_BlendWithMax << std::endl;
  os << indent << "BlendWithGaussianWeighting:"
     << m_BlendWithGaussianWeighting << std::endl;
  os << indent << "UseLog:"<< m_UseLog << std::endl;

  if( m_InputMipPointImage.IsNotNull() )
    {
    os << indent << "Input MIP Point Image: "
       << m_InputMipPointImage << std::endl;
    }
  else
    {
    os << indent << "Input MIP Point Image: NULL" << std::endl;
    }


  if( m_OutputMipPointImage.IsNotNull() )
    {
    os << indent << "Output MIP Point Image: "
       << m_OutputMipPointImage << std::endl;
    }
  else
    {
    os << indent << "Output MIP Point Image: NULL" << std::endl;
    }

}

template< typename TInputImage, typename TOutputImage >
void
ShrinkWithBlendingImageFilter< TInputImage, TOutputImage >
::SetShrinkFactor( unsigned int i, unsigned int factor )
{
  if( m_ShrinkFactors[i] == factor )
    {
    return;
    }

  this->Modified();
  m_ShrinkFactors[i] = factor;
}

template< typename TInputImage, typename TOutputImage >
unsigned int
ShrinkWithBlendingImageFilter< TInputImage, TOutputImage >
::GetShrinkFactor( unsigned int i )
{
  return m_ShrinkFactors[i];
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
    factorSize[ i ] = m_InternalShrinkFactors[ i ];
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

  typedef ImageRegionConstIteratorWithIndex<
    PointImageType > PointImageConstIteratorType;

  PointImageConstIteratorType * inputMipPointItPtr = NULL;

  if( m_InputMipPointImage.IsNotNull() )
    {
    inputMipPointItPtr = new PointImageConstIteratorType( 
      m_InputMipPointImage, outputRegionForThread );
    }

  typedef ImageRegionIteratorWithIndex<
    PointImageType > PointImageIteratorType;

  PointImageIteratorType outMipPointIt( m_OutputMipPointImage,
                                        outputRegionForThread );

  typedef ImageRegionConstIteratorWithIndex< TInputImage > InputIteratorType;

  typename TInputImage::RegionType inputRegion;

  while( !outIt.IsAtEnd() )
    {
    if( m_InputMipPointImage.IsNotNull() )
      {
      // get index of mip point from the input MIP point image
      typename TInputImage::PointType curMipPoint;

      for( unsigned int i = 0; i < TInputImage::ImageDimension; i++ )
        {
        curMipPoint[i] = inputMipPointItPtr->Get()[i];
        }

      inputPtr->TransformPhysicalPointToIndex( curMipPoint, inputIndex );

      // set output pixel intensity as intensity of input pixel at mip point
      outIt.Set( inputPtr->GetPixel( inputIndex ) );

      // set output mip point to be same as input mi point
      outMipPointIt.Set( inputMipPointItPtr->Get() );

      // move to next output pixel
      ++outIt;
      ++outMipPointIt;
      ++( *inputMipPointItPtr );

      continue;
      }

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
      outMipPointIt.Set( pointVector );
      ++outMipPointIt;
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
          double dist = ( valueIndex[i] - inputIndex[i] ) / factorSize[i];
          weight += ( 1.0 / ( factorSize[i] * std::sqrt( 2 * vnl_math::pi ) ) )
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

  if( inputMipPointItPtr != ITK_NULLPTR )
    {
    delete inputMipPointItPtr;
    }
}

template< class TInputImage, class TOutputImage >
template <typename ArrayType>
bool
ShrinkWithBlendingImageFilter< TInputImage, TOutputImage >
::NotValue( ArrayType array, double val, double tolerance )
{
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if( fabs( array[i]-val ) > tolerance )
      {
      return true;
      }
    }
  return false;
}

template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImageFilter< TInputImage, TOutputImage >
::VerifyInputInformation()
{
  bool useNewSize;
  bool useShrinkFactors;
  useNewSize = this->NotValue( m_NewSize, m_DefaultNewSize );
  useShrinkFactors = this->NotValue( m_ShrinkFactors, m_DefaultShrinkFactor );
  if( useNewSize && useShrinkFactors )
    {
    itkExceptionMacro( << "Only set one of new size or shrink factors." );
    }
  if( !useNewSize && !useShrinkFactors )
    {
    itkExceptionMacro( << "Set either a new size or shrink factors." );
    }
  Superclass::VerifyInputInformation();
}


/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImageFilter< TInputImage, TOutputImage >
::UpdateInternalShrinkFactors()
{
  if( this->NotValue( m_ShrinkFactors, m_DefaultShrinkFactor ) )
    {
    m_InternalShrinkFactors = m_ShrinkFactors;
    return;
    }
  bool warnSize = false;
  InputImagePointer  inputPtr = const_cast< TInputImage * >( 
    this->GetInput() );
  const typename TOutputImage::SizeType & inputSize =
    inputPtr->GetLargestPossibleRegion().GetSize();
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    m_InternalShrinkFactors[ i ] =
      static_cast< unsigned int >( inputSize[ i ] / m_NewSize[ i ] );
    if( static_cast< unsigned int >
      ( inputSize[ i ] / m_InternalShrinkFactors[ i ] )
        != m_NewSize[ i ] )
      {
      warnSize = true;
      }
    }
  if( warnSize )
    {
    itkWarningMacro( << "Warning: Need for integer resampling factor causes "
                        "output size to not match target m_NewSize given." );
    for( unsigned int i = 0; i < ImageDimension; ++i )
      {
      itkWarningMacro( << "   m_NewSize [" << i << "] = " << m_NewSize[ i ] );
      itkWarningMacro( << "   outSize [" << i << "] = " << static_cast< int >( 
      inputSize[ i ] / m_InternalShrinkFactors[ i ] ) );
      }
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

  // Compute the input requested region ( size and start index )
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
    factorSize[ i ] = m_InternalShrinkFactors[ i ];
    }

  InputIndexType                   inputRequestedRegionStartIndex;
  typename TInputImage::SizeType   inputRequestedRegionSize;
  typename TOutputImage::PointType tempPoint;

  outputPtr->TransformIndexToPhysicalPoint( outputRequestedRegionStartIndex,
    tempPoint );
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

  // Check that only one of the two parameter to compute the
  // new size is given ( ShrinkFactors,NewSize )

  // Update shrink factors if new size given
  this->UpdateInternalShrinkFactors();

  for( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
    {
    outputSpacing[i] = inputSpacing[i] * ( double ) m_InternalShrinkFactors[i];

    // Round down so that all output pixels fit input input region
    outputSize[i] = static_cast<SizeValueType>( std::floor( 
      ( double )inputSize[i] / ( double )m_InternalShrinkFactors[i] ) );

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
  inputPtr->TransformContinuousIndexToPhysicalPoint( inputCenterIndex,
    inputCenterPoint );
  outputPtr->TransformContinuousIndexToPhysicalPoint( outputCenterIndex,
    outputCenterPoint );

  const typename TOutputImage::PointType & inputOrigin =
    inputPtr->GetOrigin();
  typename TOutputImage::PointType outputOrigin;
  outputOrigin = inputOrigin + ( inputCenterPoint - outputCenterPoint );
  outputPtr->SetOrigin( outputOrigin );

  // make sure size of output image is same as the Input MIP Point Image
  if( m_InputMipPointImage )
    {
    bool sizeEqual = true;

    typename PointImageType::SizeType inputMipPointImageSize =
      m_InputMipPointImage->GetLargestPossibleRegion().GetSize();

    for( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
      {
      if( inputMipPointImageSize[i] != outputSize[i] )
        {
        sizeEqual = false;
        break;
        }
      }

    if( !sizeEqual )
      {
      itkExceptionMacro( 
        << "Size of output and input MIP point image do not match. "
           "Make sure you are using the same shrink amount parameters "
           "that were used to generate the input MIP point image." )
      }
    }

  // Set region
  typename TOutputImage::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize( outputSize );
  outputLargestPossibleRegion.SetIndex( outputStartIndex );

  outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );

  m_OutputMipPointImage = PointImageType::New();
  m_OutputMipPointImage->SetRegions( 
    outputPtr->GetLargestPossibleRegion() );
  m_OutputMipPointImage->CopyInformation( outputPtr );
  m_OutputMipPointImage->Allocate();

}

} // end namespace tube

} // end namespace itk

#endif
