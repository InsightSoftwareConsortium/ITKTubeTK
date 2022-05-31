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

#ifndef __itktubeCompleteImageResampleFilter_hxx
#define __itktubeCompleteImageResampleFilter_hxx


namespace itk
{

namespace tube
{

/**
 * Initialize new instance
 */
template< class TInputImage, class TOutputImage,
  class TNonSingularTransform, class TInterpolatorPrecisionType >
CompleteImageResampleFilter< TInputImage, TOutputImage,
  TNonSingularTransform, TInterpolatorPrecisionType >
::CompleteImageResampleFilter( void )
: m_DefaultPixelValue( 0 )
{
  m_OutputSpacing.Fill( 1.0 );
  m_Transform = TransformType::New();
  m_Interpolator = LinearInterpolateImageFunction< InputImageType,
    TInterpolatorPrecisionType>::New();
}

/**
 * Print out a description of self
 *
 * \todo Add details about this class
 */
template< class TInputImage, class TOutputImage,
  class TNonSingularTransform, class TInterpolatorPrecisionType >
void
CompleteImageResampleFilter< TInputImage, TOutputImage,
  TNonSingularTransform, TInterpolatorPrecisionType>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "DefaultPixelValue: "
     << static_cast< typename NumericTraits<PixelType>::PrintType >(
       m_DefaultPixelValue ) << std::endl;
  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
  os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer()
    << std::endl;

  return;
}

/**
 * Set the output image spacing.
 */
template< class TInputImage, class TOutputImage,
  class TNonSingularTransform, class TInterpolatorPrecisionType >
void
CompleteImageResampleFilter< TInputImage, TOutputImage,
  TNonSingularTransform, TInterpolatorPrecisionType >
::SetOutputSpacing( const double spacing[ImageDimension] )
{
  SpacingType s( spacing );
  this->SetOutputSpacing( s );
}

template< class TInputImage, class TOutputImage,
  class TNonSingularTransform, class TInterpolatorPrecisionType >
void
CompleteImageResampleFilter< TInputImage, TOutputImage,
  TNonSingularTransform, TInterpolatorPrecisionType >
::GenerateData( void )
{
  // Get the input image and output pointer
  InputImageConstPointer  inputImage = this->GetInput();
  OutputImagePointer      outputPtr = this->GetOutput();

  // Allocate
  this->AllocateOutputs();

  // Determine the output size and origin
  SizeType outputSize;
  PointType outputOrigin;
  TransformPointerType    transform = this->GetTransform();

  if( ImageDimension == 3 )
    {
    typename TransformType::Pointer inverse = TransformType::New();
    transform->GetInverse( inverse );
    FindOutput3DParameters( inputImage,
                            inverse.GetPointer(),
                            outputSize,
                            outputOrigin );
    }
  else
    {
    itkExceptionMacro( << "Filter only setup for 3D images." );
    outputSize = inputImage->GetLargestPossibleRegion().GetSize();
    outputOrigin = inputImage->GetOrigin();
    }

  typename ResampleImageFilterType::Pointer resampleFilter =
    ResampleImageFilterType::New();
  resampleFilter->SetInput( inputImage );
  resampleFilter->SetOutputSpacing( this->GetOutputSpacing() );
  resampleFilter->SetOutputOrigin( outputOrigin );
  resampleFilter->SetSize( outputSize );
  resampleFilter->SetTransform( transform );
  resampleFilter->SetDefaultPixelValue( m_DefaultPixelValue );
  resampleFilter->SetInterpolator( m_Interpolator );

  resampleFilter->GraftOutput( this->GetOutput() );
  resampleFilter->Update();
  this->GraftOutput( resampleFilter->GetOutput() );
}

/**
 * Defines the bounding box around the entire resampled image. Means that
 * none of the image is cut off in resampling.
 *
 * Expects the transform maps Moving->Fixed
 */
template< class TInputImage, class TOutputImage,
  class TNonSingularTransform, class TInterpolatorPrecisionType >
void
CompleteImageResampleFilter< TInputImage, TOutputImage,
  TNonSingularTransform, TInterpolatorPrecisionType >
::FindOutput3DParameters( InputImageConstPointer image,
  TransformPointerType transform, SizeType& outputSize,
  PointType& origin ) const
{
  SizeType size = image->GetRequestedRegion().GetSize();
  PointType  max;
  max.Fill( NumericTraits<double>::NonpositiveMin() );
  origin.Fill( 9999999 ); //Fill origin ( min value with large number )

  //Transform all of the image corners and determine the minimum bounding box
  for( int x = 0; x <= ( int )size[0]; x += size[0] )
    {
    for( int y = 0; y <= ( int )size[1]; y += size[1] )
      {
      for( int z = 0; z <= ( int )size[2]; z += size[2] )
        {
        typename InputImageType::IndexType index;
        index[0] = x;
        index[1] = y;
        index[2] = z;
        PointType point;
        image->TransformIndexToPhysicalPoint( index, point );
        point = transform->TransformPoint( point );

        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          if( point[i] < origin[i] )
            {
            origin[i] = point[i];
            }
          if( point[i] > max[i] )
            {
            max[i] = point[i];
            }
          }
        }
      }
    }

  // Origin represents the minimum point on the image
  typename InputImageType::SpacingType  spacing = image->GetSpacing();
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    outputSize[i] = ( long unsigned int )( ( max[i] - origin[i] ) /
      spacing[i] );
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeCompleteImageResampleFilter_hxx )
