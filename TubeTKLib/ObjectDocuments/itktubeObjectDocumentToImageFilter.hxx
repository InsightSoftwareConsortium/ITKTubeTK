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

#ifndef __itktubeObjectDocumentToImageFilter_hxx
#define __itktubeObjectDocumentToImageFilter_hxx

#include "itktubeObjectDocumentToImageFilter.h"

namespace itk
{

namespace tube
{

template< class TObjectDocument, class TImageType >
ObjectDocumentToImageFilter< TObjectDocument, TImageType >
::ObjectDocumentToImageFilter( void )
{
  this->ProcessObject::SetNthOutput( 0, ImageType::New() );
}

template< class TObjectDocument, class TImageType >
ObjectDocumentToImageFilter< TObjectDocument, TImageType >
::~ObjectDocumentToImageFilter( void )
{
}

template< class TObjectDocument, class TImageType >
typename ObjectDocumentToImageFilter< TObjectDocument, TImageType >::ImageType *
ObjectDocumentToImageFilter< TObjectDocument, TImageType >
::GetOutput( void )
{
  return static_cast< ImageType * >( this->ProcessObject::GetOutput( 0 ) );
}

template< class TObjectDocument, class TImageType >
void
ObjectDocumentToImageFilter< TObjectDocument, TImageType >
::GenerateData( void )
{
  ConstDocumentPointer document = this->GetInput();
  ImagePointer output = this->ReadDocument( document );

  if( this->GetApplyTransforms() )
    {
    TransformPointer transform = this->ComposeTransforms( document,
      this->GetStartTransforms(),
      this->GetEndTransforms() );
    output = this->ResampleImage( output, transform );
    }

  this->ProcessObject::SetNthOutput( 0, output );
}

template< class TObjectDocument, class TImageType >
typename ObjectDocumentToImageFilter< TObjectDocument, TImageType >
::ImagePointer
ObjectDocumentToImageFilter< TObjectDocument, TImageType >
::ReadDocument( ConstDocumentPointer document )
{
  typename ImageFileReaderType::Pointer reader = ImageFileReaderType::New();

  reader->SetFileName( document->GetObjectName().c_str() );
  reader->Update();

  return reader->GetOutput();
}

template< class TObjectDocument, class TImageType >
typename ObjectDocumentToImageFilter< TObjectDocument, TImageType >
::ImagePointer
ObjectDocumentToImageFilter< TObjectDocument, TImageType >
::ResampleImage( ImagePointer image, TransformPointer transform )
{
  // Resample image defaulted to linear interpolation.
  typename ResampleImageFilterType::Pointer filter =
    ResampleImageFilterType::New();

  filter->SetInput( image );
  filter->SetOutputSpacing( image->GetSpacing() );

  PointType outputOrigin;

  /* Adjust the output size and origin so that none of the image is lost
   * when transformed. */
  SizeType size;
  this->GetTransformedBoundingBox( image, transform, size, outputOrigin );

  filter->SetSize( size );
  filter->SetOutputOrigin( outputOrigin );

  /* Filter expects transform from fixed image to moving image, so we must
   * use the inverse. */
  typename TransformType::Pointer inverse = TransformType::New();

  transform->GetInverse( inverse );
  filter->SetTransform( inverse );

  // Use B-spline interpolation for resampling.
  if( m_Interpolator )
    {
    filter->SetInterpolator( m_Interpolator );
    }

  filter->Update();

  return filter->GetOutput();
}

/* Note: Set only for 3D images. Builds a bounding box around the resampled
   image. Means that none of the image is cut off in resampling. */
template< class TObjectDocument, class TImageType >
void
ObjectDocumentToImageFilter< TObjectDocument, TImageType >
::GetTransformedBoundingBox( ImagePointer image,
  TransformPointer transform,
  SizeType & outputSize,
  PointType & origin ) const
{
  SizeType size = image->GetLargestPossibleRegion().GetSize();
  PointType maximum;

  // Fill origin ( minimum value with large number ).
  origin.Fill( 9999999 );

  // Transform all of the image corners and determine the minimum bounding
  // box.
  for( unsigned int x = 0; x <= size[0]; x += size[0] )
    {
    for( unsigned int y = 0; y <= size[1]; y += size[1] )
      {
      for( unsigned int z = 0; z <= size[2]; z += size[2] )
        {
        typename ImageType::IndexType index;

        index[0] = x;
        index[1] = y;
        index[2] = z;

        PointType point;

        image->TransformIndexToPhysicalPoint( index, point );
        point = transform->TransformPoint( point );

        for( unsigned int i = 0; i < Dimension; ++i )
          {
          if( point[i] < origin[i] )
            {
            origin[i] = point[i];
            }

          if( point[i] > maximum[i] )
            {
            maximum[i] = point[i];
            }
          }
        }
      }
    }

  // Origin represents the minimum point on the image.
  typename ImageType::SpacingType spacing = image->GetSpacing();

  for( unsigned int i = 0; i < Dimension; ++i )
    {
    outputSize[i] = ( long int )( ( maximum[i] - origin[i] ) / spacing[i] );
    }
}


template< class TObjectDocument, class TImageType >
void
ObjectDocumentToImageFilter< TObjectDocument, TImageType >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  if( m_Interpolator )
    {
    os << indent << "Interpolator: " << m_Interpolator << std::endl;
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeObjectDocumentToImageFilter_hxx )
