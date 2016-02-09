/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
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
#ifndef __itktubeCropImageFilter_hxx
#define __itktubeCropImageFilter_hxx

#include "itktubeCropImageFilter.h"

namespace itk
{
namespace tube
{
template< typename TInputImage, typename TOutputImage >
CropImageFilter< TInputImage, TOutputImage >
::CropImageFilter( void )
{
  // itkCropImageFilter
  this->SetDirectionCollapseToSubmatrix();
  m_UpperBoundaryCropSize.Fill(0);
  m_LowerBoundaryCropSize.Fill(0);

  // tubecropROI
  m_ROIMin.Fill( 0 );
  m_UseROIMin = false;
  m_ROIMax.Fill( 0 );
  m_UseROIMax = false;
  m_ROISize.Fill( 0 );
  m_UseROISize = false;
  m_ROICenter.Fill( 0 );
  m_UseROICenter = false;
  m_ROIBoundary.Fill( 0 );
  m_UseROIBoundary = false;

}

template< typename TInputImage, typename TOutputImage >
void
CropImageFilter< TInputImage, TOutputImage >
::SetMin( typename ImageType::IndexType roiMin )
{
  m_ROIMin = roiMin;
  m_UseROIMin = true;
}

template< typename TInputImage, typename TOutputImage >
void
CropImageFilter< TInputImage, TOutputImage >
::SetMax( typename ImageType::IndexType roiMax )
{
  m_ROIMax = roiMax;
  m_UseROIMax = true;
}

template< typename TInputImage, typename TOutputImage >
void
CropImageFilter< TInputImage, TOutputImage >
::SetSize( typename ImageType::SizeType roiSize )
{
  m_ROISize = roiSize;
  m_UseROISize = true;
}

template< typename TInputImage, typename TOutputImage >
void
CropImageFilter< TInputImage, TOutputImage >
::SetCenter( typename ImageType::IndexType roiCenter )
{
  m_ROICenter = roiCenter;
  m_UseROICenter = true;
}

template< typename TInputImage, typename TOutputImage >
void
CropImageFilter< TInputImage, TOutputImage >
::SetBoundary( typename ImageType::IndexType roiBoundary )
{
  m_ROIBoundary = roiBoundary;
  m_UseROIBoundary = true;
}

template< typename TInputImage, typename TOutputImage >
void
CropImageFilter< TInputImage, TOutputImage >
::SetMatchVolume( typename ImageType::ConstPointer matchVolume )
{
  const typename ImageType::RegionType matchRegion =
    matchVolume->GetLargestPossibleRegion();
  const typename ImageType::IndexType matchIndex =
    matchRegion.GetIndex();
  const typename ImageType::SizeType matchSize =
    matchRegion.GetSize();
  const typename ImageType::PointType matchOrigin =
    matchVolume->GetOrigin();
  const typename ImageType::SpacingType matchSpacing =
    matchVolume->GetSpacing();

  const TInputImage *image = this->GetInput();
  const typename ImageType::RegionType imgRegion =
    image->GetLargestPossibleRegion();
  const typename ImageType::IndexType imgIndex =
    imgRegion.GetIndex();
  const typename ImageType::PointType imgOrigin =
    image->GetOrigin();
  const typename ImageType::SpacingType imgSpacing =
    image->GetSpacing();

  typename ImageType::IndexType minI;
  typename ImageType::SizeType sizeI;

  if( imgOrigin != matchOrigin || imgSpacing != matchSpacing )
    {
    for( unsigned int i = 0; i < InputImageDimension; i++ )
      {
      minI[i] = vnl_math_rnd( ( ( matchOrigin[i]
        + matchIndex[i] * matchSpacing[i] ) - ( imgOrigin[i]
        + imgIndex[i] * imgSpacing[i] ) )
        / imgSpacing[i] );
      sizeI[i] = vnl_math_rnd( ( matchSize[i] * matchSpacing[i] )
        / imgSpacing[i] );
      }
    }
  else
    {
    for( unsigned int i = 0; i < InputImageDimension; i++ )
      {
      minI[i] = matchIndex[i];
      sizeI[i] = matchSize[i];
      }
    }
  this->SetMin( minI );
  this->SetSize( sizeI );
}

template< typename TInputImage, typename TOutputImage >
void
CropImageFilter< TInputImage, TOutputImage >
::SetMatchMask( typename ImageType::Pointer maskImage )
{
  typename ImageType::IndexType minI;
  typename ImageType::IndexType maxI;

  itk::ImageRegionConstIterator< ImageType > it( maskImage,
    maskImage->GetLargestPossibleRegion() );
  while( !it.IsAtEnd() && it.Get() == 0 )
    {
    ++it;
    }
  minI = it.GetIndex();
  while( !it.IsAtEnd() && it.Get() != 0 )
    {
    ++it;
    }
  maxI = it.GetIndex();
  while( !it.IsAtEnd() )
    {
    while( !it.IsAtEnd() && it.Get() == 0 )
      {
      ++it;
      }
    if( !it.IsAtEnd() )
      {
      for( unsigned int i=0; i<InputImageDimension; ++i )
        {
        if( it.GetIndex()[i] < minI[i] )
          {
          minI[i] = it.GetIndex()[i];
          }
        }
      }
    while( !it.IsAtEnd() && it.Get() != 0 )
      {
      ++it;
      }
    if( !it.IsAtEnd() )
      {
      for( unsigned int i=0; i<InputImageDimension; ++i )
        {
        if( it.GetIndex()[i] > maxI[i] )
          {
          maxI[i] = it.GetIndex()[i];
          }
        }
      }
    }
  typename ImageType::SizeType sizeI;
  for( unsigned int i=0; i<InputImageDimension; ++i )
    {
    sizeI[i] = maxI[i] - minI[i];
    }
  this->SetMin( minI );
  this->SetSize( sizeI );
}

/**
 * Split the filter input image in every direction according
 * to the splitIndex parameter. Because the filter can only have
 * one input, the roiIndex indicates which one of the split part
 * should be used as the inputImage when calling the filter
 * Update() function.
 * For instance with 2D images, if the splitIndex
 * is [2,2], then a roiIndex equal to [1,1] corresponds to the bottom
 * right part of the original input image divided in 4 pieces.
 */
template< typename TInputImage, typename TOutputImage >
void
CropImageFilter< TInputImage, TOutputImage >
::SetSplitInput( typename ImageType::IndexType splitIndex,
  typename ImageType::IndexType roiIndex )
{
  typename ImageType::SizeType inputImageSize = this->GetInput()->
    GetLargestPossibleRegion().GetSize();

  typename ImageType::IndexType roiSize;
  for( unsigned int i = 0; i < InputImageDimension; i++ )
    {
    roiSize[i] = inputImageSize[i] / splitIndex[i];
    }
  typename ImageType::IndexType roiMin;
  roiMin.Fill( 0 );
  typename ImageType::IndexType roiMax;
  roiMax.Fill( 0 );

  for( unsigned int i = 0; i < InputImageDimension; i++ )
    {
    roiMin[i] = roiIndex[i] * roiSize[i];
    roiMax[i] = roiMin[i] + roiSize[i] - 1;
    if( roiIndex[i] == splitIndex[i]-1 )
      {
      roiMax[i] = inputImageSize[i]-1;
      }
    }

  this->SetMin( roiMin );
  this->SetMax( roiMax );

}

template< typename TInputImage, typename TOutputImage >
void
CropImageFilter< TInputImage, TOutputImage >
::GenerateOutputInformation()
{
  const TInputImage *inputPtr = this->GetInput();
  if ( !inputPtr )
    {
    return;
    }

  if( m_UseROISize || m_UseROIMin || m_UseROIMax )
    {
    if( m_UseROISize && m_UseROIMax )
      {
      //Specify either size or max options.  Not both.
      return;
      }

    if( !m_UseROIMin && !m_UseROICenter )
      {
      for( unsigned int i=0; i<InputImageDimension; i++ )
        {
        m_ROIMin[i] = 0;
        }
      }

    if( m_UseROICenter )
      {
      for( unsigned int i=0; i<InputImageDimension; i++ )
        {
        m_ROIMin[i] = m_ROICenter[i] - m_ROISize[i]/2;
        }
      }

    SizeType imageSize;
    imageSize = inputPtr->GetLargestPossibleRegion().GetSize();
    for( unsigned int i=0; i<InputImageDimension; i++ )
      {
      if( m_ROIMin[i] < 0 )
        {
        //Min is less than 0
        return;
        }
      if( m_ROIMin[i] >= (int)(imageSize[i]) )
        {
        //Min is larger than image size
        return;
        }
      }

    SizeType outputSize;
    outputSize = inputPtr->GetLargestPossibleRegion().GetSize();
    if( m_UseROISize )
      {
      for( unsigned int i=0; i<InputImageDimension; i++ )
        {
        outputSize[i] = m_ROISize[i];
        if( outputSize[i] < 1 )
          {
          outputSize[i] = 1;
          }
        }
      }
    else
      {
      for( unsigned int i=0; i<InputImageDimension; i++ )
        {
        if( m_ROIMin[i] > m_ROIMax[i] )
          {
          int tf = m_ROIMin[i];
          m_ROIMin[i] = m_ROIMax[i];
          m_ROIMax[i] = tf;
          }
        outputSize[i] = m_ROIMax[i]-m_ROIMin[i]+1;
        }
      }

    if( m_UseROIBoundary )
      {
      for( unsigned int i=0; i<InputImageDimension; i++ )
        {
        m_ROIMin[i] -= m_ROIBoundary[i];
        outputSize[i] += 2*m_ROIBoundary[i];
        }
      }

    for( unsigned int i=0; i<InputImageDimension; i++ )
      {
      if( m_ROIMin[i] + outputSize[i] > imageSize[i] )
        {
        outputSize[i] = imageSize[i] - m_ROIMin[i];
        }
      if( m_ROIMin[i] < 0 )
        {
        outputSize[i] += m_ROIMin[i];
        m_ROIMin[i] = 0;
        }
      }

    SizeType lowerCropSize;
    SizeType upperCropSize;
    for( unsigned int i=0; i<InputImageDimension; i++ )
      {
      lowerCropSize[i] = m_ROIMin[i];
      upperCropSize[i] = imageSize[i] - (m_ROIMin[i] + outputSize[i]);
      }

    this->SetLowerBoundaryCropSize( lowerCropSize );
    this->SetUpperBoundaryCropSize( upperCropSize );

    // itkCropImage GenerateOutputInformation()
    // Compute the new region size.
    OutputImageRegionType croppedRegion;
    SizeType              sz;
    OutputImageIndexType  idx;

    InputImageSizeType input_sz =
      inputPtr->GetLargestPossibleRegion().GetSize();
    InputImageIndexType input_idx =
      inputPtr->GetLargestPossibleRegion().GetIndex();

    for( unsigned int i = 0; i < InputImageDimension; ++i )
      {
      idx[i] = input_idx[i] + m_LowerBoundaryCropSize[i];
      sz[i]  = input_sz[i]  - ( m_UpperBoundaryCropSize[i] + m_LowerBoundaryCropSize[i] );
      }

      croppedRegion.SetSize(sz);
      croppedRegion.SetIndex(idx);

    // Set extraction region in the superclass.
    this->SetExtractionRegion(croppedRegion);

    Superclass::GenerateOutputInformation();

    }

}

/**
 *
 */
template< typename TInputImage, typename TOutputImage >
void
CropImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  // itktubeCropImageFilter
  Superclass::PrintSelf(os, indent);

  os << indent << "UpperBoundaryCropSize: " << m_UpperBoundaryCropSize
     << std::endl;
  os << indent << "LowerBoundaryCropSize: " << m_LowerBoundaryCropSize
     << std::endl;
  // tubeCropROI
  os << indent << "ROIMin: " << m_ROIMin << std::endl;
  if( m_UseROIMin )
    {
    os << indent << "Use ROIMin: true" << std::endl;
    }
  else
    {
    os << indent << "Use ROIMin: false" << std::endl;
    }

  os << indent << "ROIMax: " << m_ROIMax << std::endl;
  if( m_UseROIMax )
    {
    os << indent << "Use ROIMax: true" << std::endl;
    }
  else
    {
    os << indent << "Use ROIMax: false" << std::endl;
    }

  os << indent << "ROISize: " << m_ROISize << std::endl;
  if( m_UseROISize )
    {
    os << indent << "Use ROISize: true" << std::endl;
    }
  else
    {
    os << indent << "Use ROISize: false" << std::endl;
    }

  os << indent << "ROICenter: " << m_ROICenter << std::endl;
  if( m_UseROICenter )
    {
    os << indent << "Use ROICenter: true" << std::endl;
    }
  else
    {
    os << indent << "Use ROICenter: false" << std::endl;
    }

  os << indent << "ROIBoundary: " << m_ROIBoundary << std::endl;
  if( m_UseROIBoundary )
    {
    os << indent << "Use ROIBoundary: true" << std::endl;
    }
  else
    {
    os << indent << "Use ROIBoundary: false" << std::endl;
    }

}

} // end namespace tube
} // end namespace itk

#endif
