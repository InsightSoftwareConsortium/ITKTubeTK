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

#ifndef __itktubeConvertShrunkenSeedImageToListFilter_hxx
#define __itktubeConvertShrunkenSeedImageToListFilter_hxx

#include "itktubeConvertShrunkenSeedImageToListFilter.h"

namespace itk
{
namespace tube
{

/** Constructor */
template< class TImage, class TPointsImage >
ConvertShrunkenSeedImageToListFilter< TImage, TPointsImage >
::ConvertShrunkenSeedImageToListFilter( void )
{
  this->ProcessObject::SetNthOutput( 0, OutputType::New().GetPointer() );
  m_Threshold = 0;
  this->ProcessObject::SetNumberOfRequiredInputs( 3 );
}

template< class TImage, class TPointsImage >
void
ConvertShrunkenSeedImageToListFilter< TImage, TPointsImage >
::SetInput( const ImageType* image )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0, const_cast< ImageType* >( image ) );
}

template< class TImage, class TPointsImage >
void
ConvertShrunkenSeedImageToListFilter< TImage, TPointsImage >
::SetScaleImage( const ImageType* image )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 1, const_cast< ImageType* >( image ) );
}

template< class TImage, class TPointsImage >
void
ConvertShrunkenSeedImageToListFilter< TImage, TPointsImage >
::SetPointsImage( const PointsImageType* image )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 2,
    const_cast< PointsImageType* >( image ) );
}

template< class TImage, class TPointsImage >
const TImage*
ConvertShrunkenSeedImageToListFilter< TImage, TPointsImage >
::GetInput( void ) const
{
  return static_cast<const ImageType * >
    ( this->ProcessObject::GetInput( 0 ) );
}

template< class TImage, class TPointsImage >
const TImage*
ConvertShrunkenSeedImageToListFilter< TImage, TPointsImage >
::GetScaleImage( void ) const
{
  return static_cast<const ImageType * >
    ( this->ProcessObject::GetInput( 1 ) );
}

template< class TImage, class TPointsImage >
const TPointsImage*
ConvertShrunkenSeedImageToListFilter< TImage, TPointsImage >
::GetPointsImage( void ) const
{
  return static_cast<const PointsImageType * >
    ( this->ProcessObject::GetInput( 2 ) );
}

template< class TImage, class TPointsImage >
SimpleDataObjectDecorator< vnl_matrix <typename TImage::PixelType> >*
ConvertShrunkenSeedImageToListFilter< TImage, TPointsImage >
::GetOutput()
{
  // we assume that the first output is of the templated type
  return itkDynamicCastInDebugMode< OutputType * >( this->GetPrimaryOutput() );
}

template< class TImage, class TPointsImage >
void
ConvertShrunkenSeedImageToListFilter< TImage, TPointsImage >
::GenerateData( void )
{
  const ImageType* inImage = this->GetInput();
  const ImageType* inScale = const_cast< ImageType * >( this->GetScaleImage() );
  const PointsImageType* inPoint =
    const_cast< PointsImageType * >( this->GetPointsImage() );

  itk::ImageRegionConstIterator< ImageType > itImage( inImage,
    inImage->GetLargestPossibleRegion() );
  itk::ImageRegionConstIterator< ImageType > itScale( inScale,
    inScale->GetLargestPossibleRegion() );
  itk::ImageRegionConstIterator< PointsImageType > itPoint( inPoint,
    inPoint->GetLargestPossibleRegion() );

  const SizeValueType ARows =
    inImage->GetLargestPossibleRegion().GetNumberOfPixels();
  if( ARows > std::numeric_limits<unsigned int>::max() )
    {
    itkExceptionMacro( <<
      "Exception caught ! The image is too big for this filter." );
    }
  m_VnlOutput.set_size( ARows, ImageType::ImageDimension + 1 );

  unsigned int row = 0;
  while( !itImage.IsAtEnd() )
    {
    if( itImage.Get() > m_Threshold )
      {
      for( unsigned int i = 0; i < ImageDimension; ++i )
        {
        m_VnlOutput( row, i ) = itPoint.Get()[ i ];
        }
      m_VnlOutput( row, ImageDimension ) = itScale.Get();
      row++;
      }
    ++itImage;
    ++itScale;
    ++itPoint;
    }
  typename OutputType::Pointer outputPtr = this->GetOutput();
  outputPtr->Set( m_VnlOutput );
}

/** PrintSelf */
template< class TImage, class TPointsImage >
void
ConvertShrunkenSeedImageToListFilter< TImage, TPointsImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "VnlOutput = " << m_VnlOutput << std::endl;
  os << indent << "Threshold = " << m_Threshold << std::endl;
}

template< class TImage, class TPointsImage >
void
ConvertShrunkenSeedImageToListFilter< TImage, TPointsImage >
::VerifyPreconditions()
{
  this->ProcessObject::VerifyPreconditions();
  if( this->GetInput()->GetLargestPossibleRegion().GetSize() !=
    this->GetScaleImage()->GetLargestPossibleRegion().GetSize()
    || this->GetInput()->GetLargestPossibleRegion().GetSize() !=
    this->GetPointsImage()->GetLargestPossibleRegion().GetSize() )
    {
    itkExceptionMacro( << "The input images don't have the same size." );
    }
}

} // End namespace tube

} // End namespace itk
#endif // End !defined( __itktubeConvertShrunkenSeedImageToListFilter_hxx )
