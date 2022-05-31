/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef __itktubeDiffusiveRegistrationFilterUtils_hxx
#define __itktubeDiffusiveRegistrationFilterUtils_hxx


#include <itkMinimumMaximumImageCalculator.h>
#include <itkResampleImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkResampleImageFilter.h>

namespace itk
{

namespace tube
{
/**
 * Helper function to allocate space for an image given a template image
 */
template< class TUnallocatedImagePointer, class TemplateImagePointer >
void
DiffusiveRegistrationFilterUtils
::AllocateSpaceForImage( TUnallocatedImagePointer& image,
                         const TemplateImagePointer& templateImage )
{
  assert( image );
  assert( templateImage );
  image->SetOrigin( templateImage->GetOrigin() );
  image->SetSpacing( templateImage->GetSpacing() );
  image->SetDirection( templateImage->GetDirection() );
  image->SetLargestPossibleRegion( templateImage->GetLargestPossibleRegion() );
  image->SetRequestedRegion( templateImage->GetRequestedRegion() );
  image->SetBufferedRegion( templateImage->GetBufferedRegion() );
  image->Allocate();
}

/**
 * Helper function to check whether the attributes of an image matches template
 */
template< class TCheckedImage, class TTemplateImage >
bool
DiffusiveRegistrationFilterUtils
::CompareImageAttributes( const TCheckedImage * image,
                          const TTemplateImage * templateImage )
{
  assert( image );
  assert( templateImage );

  return image->GetOrigin() == templateImage->GetOrigin()
      && image->GetSpacing() == templateImage->GetSpacing()
      && image->GetDirection() == templateImage->GetDirection()
      && image->GetLargestPossibleRegion()
          == templateImage->GetLargestPossibleRegion()
      && image->GetLargestPossibleRegion().GetIndex()
          == templateImage->GetLargestPossibleRegion().GetIndex()
      && image->GetLargestPossibleRegion().GetSize()
          == templateImage->GetLargestPossibleRegion().GetSize();
}

/**
 * Resample an image to match a template
 */
template< class TResampleImagePointer, class TemplateImagePointer >
void
DiffusiveRegistrationFilterUtils
::ResampleImageNearestNeighbor(
    const TResampleImagePointer & highResolutionImage,
    const TemplateImagePointer & templateImage,
    TResampleImagePointer & resampledImage )
{
  // We have to implement nearest neighbors by hand, since we are dealing with
  // pixel types that do not have Numeric Traits
  typedef typename TResampleImagePointer::ObjectType ResampleImageType;

  // Create the resized resampled image
  resampledImage = ResampleImageType::New();
  DiffusiveRegistrationFilterUtils::AllocateSpaceForImage( resampledImage,
    templateImage );

  // Do NN interpolation
  typedef itk::ImageRegionIteratorWithIndex< ResampleImageType >
      ResampleImageRegionType;
  ResampleImageRegionType resampledImageIt = ResampleImageRegionType(
      resampledImage, resampledImage->GetLargestPossibleRegion() );

  typename ResampleImageType::PointType physicalPoint;
  physicalPoint.Fill( 0.0 );
  typename ResampleImageType::IndexType highResolutionIndex;
  highResolutionIndex.Fill( 0.0 );
  typename ResampleImageType::PixelType pixelValue;

  for( resampledImageIt.GoToBegin();
       !resampledImageIt.IsAtEnd();
       ++resampledImageIt )
    {
    resampledImage->TransformIndexToPhysicalPoint(
        resampledImageIt.GetIndex(), physicalPoint );
    highResolutionImage->TransformPhysicalPointToIndex(
        physicalPoint, highResolutionIndex );
    pixelValue = highResolutionImage->GetPixel( highResolutionIndex );
    resampledImageIt.Set( pixelValue );
    }
}

/**
 * Resample an image to match a template
 */
template< class TResampleImagePointer, class TemplateImagePointer >
void
DiffusiveRegistrationFilterUtils
::ResampleImageLinear( const TResampleImagePointer & highResolutionImage,
                       const TemplateImagePointer & templateImage,
                       TResampleImagePointer & resampledImage )
{
  // Do linear interpolation
  typedef itk::ResampleImageFilter
      < typename TResampleImagePointer::ObjectType,
        typename TResampleImagePointer::ObjectType > ResampleFilterType;
  typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput( highResolutionImage );
  resampler->SetOutputParametersFromImage( templateImage );
  resampler->Update();
  resampledImage = resampler->GetOutput();
}

/**
 * Resample a vector image to match a template
 */
template< class TVectorResampleImagePointer, class TemplateImagePointer >
void
DiffusiveRegistrationFilterUtils
::VectorResampleImageLinear(
    const TVectorResampleImagePointer & highResolutionImage,
    const TemplateImagePointer & templateImage,
    TVectorResampleImagePointer & resampledImage,
    bool normalize )
{
  // Do linear interpolation
  typedef itk::ResampleImageFilter
      < typename TVectorResampleImagePointer::ObjectType,
        typename TVectorResampleImagePointer::ObjectType > ResampleFilterType;
  typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput( highResolutionImage );
  resampler->SetOutputOrigin( templateImage->GetOrigin() );
  resampler->SetOutputSpacing( templateImage->GetSpacing() );
  resampler->SetOutputDirection( templateImage->GetDirection() );
  resampler->SetOutputStartIndex(
      templateImage->GetLargestPossibleRegion().GetIndex() );
  resampler->SetSize( templateImage->GetLargestPossibleRegion().GetSize() );
  resampler->Update();
  resampledImage = resampler->GetOutput();

  if( normalize )
    {
    DiffusiveRegistrationFilterUtils::NormalizeVectorField( resampledImage );
    }
}

/**
 * Normalizes a vector field to ensure each vector has length 1
 */
template< class TVectorImagePointer >
void
DiffusiveRegistrationFilterUtils
::NormalizeVectorField( TVectorImagePointer & image )
{
  typedef ImageRegionIterator< typename TVectorImagePointer::ObjectType >
      DeformationVectorImageRegionType;
  DeformationVectorImageRegionType vectorIt(
      image, image->GetLargestPossibleRegion() );
  for( vectorIt.GoToBegin(); !vectorIt.IsAtEnd(); ++vectorIt )
    {
    vectorIt.Value().Normalize();
    }
}


/**
 * Returns whether an image has intensity range between 0 and 1
 */
template< class TImage >
bool
DiffusiveRegistrationFilterUtils
::IsIntensityRangeBetween0And1( TImage * image )
{
  typedef itk::MinimumMaximumImageCalculator< TImage > CalculatorType;
  typename CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetImage( image );
  calculator->Compute();

  if( calculator->GetMinimum() < 0.0 || calculator->GetMaximum() > 1.0 )
    {
    return false;
    }
  return true;
}

/**
 * Update x, y, z components of a deformation field
 */
template< class TDeformationField, class TDeformationComponentImageArray >
void
DiffusiveRegistrationFilterUtils
::ExtractXYZComponentsFromDeformationField(
    const TDeformationField * deformationField,
    TDeformationComponentImageArray& deformationComponentImages )
{
  assert( deformationField );

  typedef itk::VectorIndexSelectionCastImageFilter
      < TDeformationField,
      typename TDeformationComponentImageArray::ValueType::ObjectType >
      VectorIndexSelectionFilterType;
  typename VectorIndexSelectionFilterType::Pointer indexSelector;
  for( unsigned int i = 0; i < TDeformationField::ImageDimension; i++ )
    {
    indexSelector = VectorIndexSelectionFilterType::New();
    indexSelector->SetInput( deformationField );
    indexSelector->SetIndex( i );
    deformationComponentImages[i] = indexSelector->GetOutput();
    indexSelector->Update();
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeDiffusiveRegistrationFilterUtils_hxx )
