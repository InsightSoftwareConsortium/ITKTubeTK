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

#ifndef __itkObjectDocumentToImageFilter_txx
#define __itkObjectDocumentToImageFilter_txx

#include "itkObjectDocumentToImageFilter.h"

namespace itk
{

namespace tube
{

template< class TInputObjectDocument, class TOutputImageType >
ObjectDocumentToImageFilter<TInputObjectDocument,TOutputImageType>
::ObjectDocumentToImageFilter( void )
{
  OutputImagePointer object = OutputImageType::New();
  this->ProcessObject::SetNthOutput(0, object);
}


template< class TInputObjectDocument, class TOutputImageType >
typename ObjectDocumentToImageFilter<TInputObjectDocument,TOutputImageType>::OutputImageType *
ObjectDocumentToImageFilter<TInputObjectDocument,TOutputImageType>
::GetOutput( void )
{
  return static_cast<OutputImageType *>( (this->ProcessObject::GetOutput(0)) );
}


template< class TInputObjectDocument, class TOutputImageType >
void
ObjectDocumentToImageFilter<TInputObjectDocument,TOutputImageType>
::GenerateData( void )
{
  ConstDocumentPointer document = this->GetInput();

  OutputImagePointer output = ReadDocument( document );
  if( this->m_ApplyTransforms )
    {
    TransformPointer trans = this->ComposeTransforms( document, this->m_StartTransforms, this->m_EndTransforms );
    output = ResampleImage( output, trans );
    }
  this->ProcessObject::SetNthOutput( 0, output );
}


template< class TInputObjectDocument, class TOutputImageType >
typename ObjectDocumentToImageFilter<TInputObjectDocument,TOutputImageType>::OutputImagePointer
ObjectDocumentToImageFilter<TInputObjectDocument,TOutputImageType>
::ReadDocument( ConstDocumentPointer doc )
{

  typename ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
  reader->SetFileName( doc->GetObjectName() );

  reader->Update();
  return reader->GetOutput();
}


template< class TInputObjectDocument, class TOutputImageType >
typename ObjectDocumentToImageFilter<TInputObjectDocument,TOutputImageType>::OutputImagePointer
ObjectDocumentToImageFilter<TInputObjectDocument,TOutputImageType>
::ResampleImage( OutputImagePointer image, TransformPointer trans )
{
  //Resample Image Defaulted to Linear Interpolation
  typename ResampleImageFilterType::Pointer filter = ResampleImageFilterType::New();
  filter->SetInput( image );
  filter->SetOutputSpacing( image->GetSpacing() );

  typename OutputImageType::PointType outputOrigin;
  //Adjust the output size and origin so that none of the image is lost when transformed.
  typename OutputImageType::SizeType size;
  GetTransformedBoundingBox( image, trans, size, outputOrigin);
  filter->SetSize( size );
  filter->SetOutputOrigin( outputOrigin );

  // Filter expects transform from Fixed image to moving image, must use inverse
  typename TransformType::Pointer inverse = TransformType::New();
  trans->GetInverse(inverse);
  filter->SetTransform(inverse);

  //Use BSpline Interpolation for resampling
  if( m_Interpolator.IsNotNull() ) {
    filter->SetInterpolator( m_Interpolator );
    }

  filter->Update();
  return filter->GetOutput();
}

//***************************** NOTE: Set only for 3D images ************************************ //
//Builds a bounding box around the resampled image.  Means that none of the image is cut off in resampling
template< class TInputObjectDocument, class TOutputImageType >
void
ObjectDocumentToImageFilter<TInputObjectDocument,TOutputImageType>
::GetTransformedBoundingBox( OutputImagePointer image, TransformPointer transform, SizeType &outputSize, PointType& origin ) const
  {
  typename OutputImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
  typename OutputImageType::PointType  max;
  origin.Fill( 9999999 ); //Fill origin ( min value with large number )
  //Transform all of the image corners and determine the minimum bounding box
  for( unsigned int x = 0; x <= size[0]; x+= size[0] )
    {
    for( unsigned int y = 0; y <= size[1]; y+= size[1] )
      {
      for( unsigned int z = 0; z <= size[2]; z+= size[2] )
        {
        typename OutputImageType::IndexType index;
        index[0] = x;
        index[1] = y;
        index[2] = z;
        typename OutputImageType::PointType point;
        image->TransformIndexToPhysicalPoint( index, point );
        point = transform->TransformPoint( point );

        for( unsigned int i = 0; i < TDimensions; i++ )
          {
          if( point[i] < origin[i] ) { origin[i] = point[i]; }
          if( point[i] > max[i] ) { max[i] = point[i]; }
          }
        }
      }
    }

  // Origin represents the minimum point on the image
  typename OutputImageType::SpacingType  spacing = image->GetSpacing();
  for( unsigned int i = 0; i < TDimensions; i++ )
    {
    outputSize[i] = (long int)(( max[i] - origin[i] ) / spacing[i] );
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itkObjectDocumentToImageFilter_txx)
