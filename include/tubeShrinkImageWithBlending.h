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
#ifndef __tubeShrinkImageWithBlending_h
#define __tubeShrinkImageWithBlending_h

#include "itkProcessObject.h"

#include "tubeWrappingMacros.h"

#include "itktubeShrinkWithBlendingImageFilter.h"

namespace tube
{
/** \class ShrinkImageWithBlending
 *
 *  \ingroup TubeTK
 */

template< typename TInputImage, typename TOutputImage >
class ShrinkImageWithBlending:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ShrinkImageWithBlending                         Self;
  typedef itk::ProcessObject                              Superclass;
  typedef itk::SmartPointer< Self >                       Pointer;
  typedef itk::SmartPointer< const Self >                 ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ShrinkImageWithBlending, ProcessObject );


  /** Typedef to images */
  typedef TInputImage                                InputImageType;
  typedef TOutputImage                               OutputImageType;
  typedef typename TInputImage::IndexType            InputIndexType;
  typedef typename TInputImage::SizeType             InputSizeType;

  typedef itk::tube::ShrinkWithBlendingImageFilter< InputImageType,
    OutputImageType >                                FilterType;

  typedef typename FilterType::PointImagePixelType   PointImagePixelType;
  typedef typename FilterType::PointImageType        PointImageType;
  typedef typename FilterType::ShrinkFactorsType     ShrinkFactorsType;

  /** Set the shrink factors. Values are clamped to
   * a minimum value of 1. Default is 1 for all dimensions. */
  tubeWrapSetMacro( ShrinkFactors, ShrinkFactorsType, Filter )
  void SetShrinkFactor( unsigned int i, unsigned int factor );
  unsigned int GetShrinkFactor( unsigned int i );

  tubeWrapSetMacro( NewSize, InputSizeType, Filter );
  tubeWrapGetMacro( NewSize, InputSizeType, Filter );

  /** Get/Set the shrink factors. */
  tubeWrapGetMacro( ShrinkFactors, ShrinkFactorsType, Filter );

  tubeWrapSetMacro( Overlap, InputIndexType, Filter );
  tubeWrapGetMacro( Overlap, InputIndexType, Filter );

  tubeWrapSetMacro( BlendWithMean, bool, Filter );
  tubeWrapGetMacro( BlendWithMean, bool, Filter );

  tubeWrapSetMacro( BlendWithMax, bool, Filter );
  tubeWrapGetMacro( BlendWithMax, bool, Filter );

  tubeWrapSetMacro( BlendWithGaussianWeighting, bool, Filter );
  tubeWrapGetMacro( BlendWithGaussianWeighting, bool, Filter );

  tubeWrapSetMacro( UseLog, bool, Filter );
  tubeWrapGetMacro( UseLog, bool, Filter );

  tubeWrapCallOverrideMacro( GenerateOutputInformation, Filter );
  tubeWrapCallOverrideMacro( GenerateInputRequestedRegion, Filter );

  tubeWrapSetConstObjectMacro( Input, InputImageType, Filter );
  tubeWrapGetConstObjectMacro( Input, InputImageType, Filter );

  tubeWrapSetConstObjectMacro( InputMipPointImage, PointImageType, Filter );
  tubeWrapGetConstObjectMacro( InputMipPointImage, PointImageType, Filter );

  tubeWrapUpdateMacro( Filter );

  tubeWrapGetObjectMacro( Output, OutputImageType, Filter );

  tubeWrapGetObjectMacro( OutputMipPointImage, PointImageType, Filter );

protected:
  ShrinkImageWithBlending( void );
  ~ShrinkImageWithBlending() {}

  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

private:
  /** itkShrinkImageWithBlendingFilter parameters **/
  ShrinkImageWithBlending( const Self & );

  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) override
    {};

  typename FilterType::Pointer m_Filter;

};
} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeShrinkImageWithBlending.hxx"
#endif

#endif // End !defined( __tubeShrinkImageWithBlending_h )
