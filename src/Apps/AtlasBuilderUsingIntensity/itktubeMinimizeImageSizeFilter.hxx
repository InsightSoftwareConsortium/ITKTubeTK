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

#ifndef __itktubeMinimizeImageSizeFilter_hxx
#define __itktubeMinimizeImageSizeFilter_hxx

#include "itktubeMinimizeImageSizeFilter.h"

namespace itk
{

namespace tube
{

template< class TInputImage >
MinimizeImageSizeFilter< TInputImage >
::MinimizeImageSizeFilter( void )
: m_BufferImage( false ),
  m_ThresholdValue( 0 ),
  m_ThresholdAbove( false ),
  m_DefaultPixelValue( 0 ),
  m_ClipEndIndices( false ),
  m_ClipStartIndices( false )
{
  m_NumberOfBufferPixels.Fill( 0 );
}

template< class TInputImage >
void
MinimizeImageSizeFilter<TInputImage>
::GenerateData( void )
{
  InputImageConstPointer input = this->GetInput();
  RegionType region = input->GetLargestPossibleRegion();

  this->AllocateOutputs();

  if( ImageDimension == 3 )
    {
    // Determine the clipped end of the image
    if( this->GetClipEndIndices() )
      {
      Get3DCroppedEndRegion( input, region );
      }
    // Determine the clipped start of the image
    if( this->GetClipStartIndices() )
      {
      Get3DCroppedStartRegion( input, region );
      }
    }
  else if( ImageDimension == 2 )
    {
    itkExceptionMacro( << "2D image cropping is not implemented yet." );
    return;
    }
  else
    {
    itkExceptionMacro( << "Dimension size is not compatible with filter." );
    return;
    }

  PointType newOrigin;
  input->TransformIndexToPhysicalPoint( region.GetIndex(), newOrigin );

  SizeType size;
  size = region.GetSize();

  // Adjust for the buffer ( if provided )
  if( this->GetBufferImage() )
    {
    if( this->GetClipEndIndices() )
      {
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        size[i] += m_NumberOfBufferPixels[i];
        }
      }
    if( this->GetClipStartIndices() )
      {
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        size[j] += m_NumberOfBufferPixels[j];
        // Get the physical point for the new origin
        newOrigin[j] -= ( m_NumberOfBufferPixels[j] *
          input->GetSpacing()[j] );
        }
      }
    }

  // resample image to given region parameters
  typedef ResampleImageFilter< InputImageType, OutputImageType >
    ResampleFilterType;
  typename ResampleFilterType::Pointer resampleFilter =
    ResampleFilterType::New();

  resampleFilter->SetInput( input );
  resampleFilter->SetOutputSpacing( input->GetSpacing() );
  resampleFilter->SetDefaultPixelValue( this->GetDefaultPixelValue() );
  resampleFilter->SetOutputOrigin( newOrigin );
  resampleFilter->SetSize( size );

  // graft our output to the resample filter
  resampleFilter->GraftOutput( this->GetOutput() );
  // update the resample filter
  resampleFilter->Update();
  // and the output of the resample filter back onto this one
  this->GraftOutput( resampleFilter->GetOutput() );
}

template< class TInputImage >
void
MinimizeImageSizeFilter< TInputImage >
::Get3DCroppedEndRegion( InputImageConstPointer input, RegionType& region )
{
  InputPixelType  threshold = GetThresholdValue();
  bool            isThresAbove = this->GetThresholdAbove();

  typedef ImageSliceConstIteratorWithIndex< InputImageType >
    ImageSliceConstIteratorType;
  SizeType  size = region.GetSize();

  // Run sweep for back end
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    ImageSliceConstIteratorType  it_input( input, region );
    it_input.SetFirstDirection( i );
    it_input.SetSecondDirection( ( ( i+1 )%ImageDimension ) );

    it_input.GoToReverseBegin();

    bool breakPoint = false;
    while( !it_input.IsAtReverseEnd() )
      {
      while( !it_input.IsAtReverseEndOfSlice() )
        {
        while( !it_input.IsAtReverseEndOfLine() )
          {
          if( it_input.Get() > threshold && !isThresAbove )
            {
            breakPoint = true;
            break;
            }
          else if( it_input.Get() < threshold && isThresAbove )
            {
            breakPoint = true;
            break;
            }
          --it_input;
          }
        if( breakPoint )
          {
          break;
          }
        it_input.PreviousLine();
        }
      if( breakPoint )
        {
        break;
        }
      it_input.PreviousSlice();
      }

    // This is the index of the first point found in the plane defined by
    // i & i+1 ( starting from end ). Therefore it is the first point for
    // dimension i+2
    IndexType index = it_input.GetIndex();

    unsigned int dim = ( i+2 )%ImageDimension;  //Dimension that is changed
    if( index[dim] < 0 )
      {
      std::cout << "Index below 0 = " << index[dim] << std::endl;
      index[dim] = 0;
      }
    size[dim] = ( index[dim] - region.GetIndex()[dim] ) + 1;
    region.SetSize( size );
    }
}

template< class TInputImage >
void
MinimizeImageSizeFilter< TInputImage >
::Get3DCroppedStartRegion( InputImageConstPointer input, RegionType& region )
{
  InputPixelType  threshold = GetThresholdValue();
  bool            isThresAbove = this->GetThresholdAbove();

  typedef ImageSliceConstIteratorWithIndex< InputImageType>
    ImageSliceConstIteratorType;
  IndexType index = region.GetIndex();
  SizeType  size = region.GetSize();

  // Run sweep for front end of image
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    ImageSliceConstIteratorType  it_input( input, region );
    it_input.SetFirstDirection( i );
    it_input.SetSecondDirection( ( ( i+1 )%ImageDimension ) );

    it_input.GoToBegin();

    bool breakPoint = false;
    while( !it_input.IsAtEnd() )
      {
      while( !it_input.IsAtEndOfSlice() )
        {
        while( !it_input.IsAtEndOfLine() )
          {
          if( it_input.Get() > threshold && !isThresAbove )
            {
            breakPoint = true;
            break;
            }
          else if( it_input.Get() < threshold && isThresAbove )
            {
            breakPoint = true;
            break;
            }
          ++it_input;
          }
        if( breakPoint )
          {
          break;
          }
        it_input.NextLine();
        }
      if( breakPoint )
        {
        break;
        }
      it_input.NextSlice();
      }

    // This is the index of the first point found in the plane defined by
    // i & i+1 ( starting from beginning ). Therefore it is the first point
    // for dimension i+2
    unsigned int dim = ( i+2 ) % ImageDimension;  //Dimension that is changed
    index[dim] = it_input.GetIndex()[dim];

    size[dim] -= ( index[dim] - region.GetIndex()[dim] );
    region.SetSize( size );
    region.SetIndex( index );
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeMinimizeImageSizeFilter_hxx )
