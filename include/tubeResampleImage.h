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
#ifndef tubeResampleImage_h
#define tubeResampleImage_h

// ITK Includes
#include "itkProcessObject.h"

// TubeTK Includes
#include "tubeWrappingMacros.h"

#include "itktubeReResampleImageFilter.h"

namespace tube
{
/** \class ResampleImage
 *
 *  \ingroup TubeTK
 */

template< class TImage >
class ResampleImage:
  public itk::ProcessObject
{
public:
  /** Standard class type alias. */
  using Self = ResampleImage;
  using Superclass = itk::ProcessObject;
  using Pointer = itk::SmartPointer< Self >;
  using ConstPointer = itk::SmartPointer< const Self >;

  using FilterType = itk::tube::ReResampleImageFilter< typename TImage::PixelType,
          TImage::ImageDimension >;

  using ImageType = typename FilterType::ImageType;
  using ConstImagePointer = typename ImageType::ConstPointer;
  using ImagePointer = typename ImageType::Pointer;
  using TransformType = typename FilterType::TransformType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ResampleImage, ProcessObject );

  /* Set input image */
  tubeWrapSetObjectMacro( Input, ImageType, Filter );
  tubeWrapGetObjectMacro( Input, ImageType, Filter );

  /** Set/Get input Match Image */
  tubeWrapSetObjectMacro( MatchImage, ImageType, Filter );
  tubeWrapGetObjectMacro( MatchImage, ImageType, Filter );

  /** Set/Get whether Output is isotropic or not */
  tubeWrapSetMacro( MakeIsotropic, bool, Filter );
  tubeWrapGetMacro( MakeIsotropic, bool, Filter );

  /** Set/Get whether Output is high resolution isotropic or not */
  tubeWrapSetMacro( MakeHighResIso, bool, Filter );
  tubeWrapGetMacro( MakeHighResIso, bool, Filter );

  /** Set/Get interpolator */
  tubeWrapSetMacro( Interpolator, std::string, Filter );
  tubeWrapGetMacro( Interpolator, std::string, Filter );

  /** Set/Get whether to load transform or not */
  tubeWrapSetMacro( LoadTransform, bool, Filter );
  tubeWrapGetMacro( LoadTransform, bool, Filter );

  tubeWrapGetObjectMacro( Output, ImageType, Filter );

  /** Set Output Transform */
  tubeWrapForceSetObjectMacro( Transform, TransformType, Filter );

  /** Set Output Spacing */
  tubeWrapForceSetMacro( Spacing, std::vector<double>, Filter );

  /** Set Output Origin */
  tubeWrapForceSetMacro( Origin, std::vector<double>, Filter );

  /** Set Output Index */
  tubeWrapForceSetMacro( Index, std::vector<int>, Filter );

  /** Set Output Index */
  tubeWrapForceSetMacro( Size, std::vector<int>, Filter );

  /** Set Output Resample Factor */
  tubeWrapForceSetMacro( ResampleFactor, std::vector<double>, Filter );

  /* Runs tubes to image conversion */
  tubeWrapUpdateMacro( Filter );

protected:
  ResampleImage( void );
  ~ResampleImage() {}

  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

private:
  /** itktubeResampleImageFilter parameters **/
  ResampleImage( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) override {};

  typename FilterType::Pointer  m_Filter;
};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeResampleImage.hxx"
#endif

#endif // End !defined( __tubeResampleImage_h )
