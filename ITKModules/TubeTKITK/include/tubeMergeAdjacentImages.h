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
#ifndef __tubeMergeAdjacentImages_h
#define __tubeMergeAdjacentImages_h

// standard includes
#include <string>

// ITK includes
#include <itkMacro.h>
#include <itkProcessObject.h>

// TubeTK includes
#include "tubeWrappingMacros.h"

#include "itktubeMergeAdjacentImagesFilter.h"

namespace tube
{
/** \class MergeAdjacentImages
 *  \brief Merge the second image into the space of the first
 *
 *  \ingroup TubeTKITK
 */

template< class TImage >
class MergeAdjacentImages:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef MergeAdjacentImages                                  Self;
  typedef itk::ProcessObject                                   Superclass;
  typedef itk::SmartPointer< Self >                            Pointer;
  typedef itk::SmartPointer< const Self >                      ConstPointer;

  typedef TImage                                               ImageType;
  typedef typename TImage::PixelType                           PixelType;

  typedef itk::tube::MergeAdjacentImagesFilter< ImageType >    FilterType;
  typedef typename FilterType::PaddingType                     PaddingType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MergeAdjacentImages, ProcessObject);

  /** Set input image 1 */
  tubeWrapSetConstObjectMacro( Input1, ImageType, Filter );

  /** Get input image 1 */
  tubeWrapGetConstObjectMacro( Input1, ImageType, Filter );

  /** Set input image 2 */
  tubeWrapSetConstObjectMacro( Input2, ImageType, Filter );

  /** Get input image 2 */
  tubeWrapGetConstObjectMacro( Input2, ImageType, Filter );

  /** Set value used for output pixels that dont intersect with input image */
  tubeWrapSetMacro( Background, PixelType, Filter );

  /** Get value used for output pixels that dont intersect with input image */
  tubeWrapGetMacro( Background, PixelType, Filter );

  /** Set if zero-valued input pixels should be ignored */
  tubeWrapSetMacro( MaskZero, bool, Filter );

  /** Get if zero-valued input pixels should be ignored */
  tubeWrapGetMacro( MaskZero, bool, Filter );

  /** Set number of registration iterations */
  tubeWrapSetMacro( MaxIterations, unsigned int, Filter );

  /** Get number of registration iterations */
  tubeWrapGetMacro( MaxIterations, unsigned int, Filter );

  /** Set padding for second image */
  tubeWrapSetConstReferenceMacro( Padding, PaddingType, Filter );

  /** Get padding for second image */
  tubeWrapGetConstReferenceMacro( Padding, PaddingType, Filter );

  /** Set expected initial misalignment offset */
  tubeWrapSetMacro( ExpectedOffset, double, Filter );

  /** Get expected initial misalignment offset */
  tubeWrapGetMacro( ExpectedOffset, double, Filter );

  /** Set expected initial misalignment rotation */
  tubeWrapSetMacro( ExpectedRotation, double, Filter );

  /** Get expected initial misalignment rotation */
  tubeWrapGetMacro( ExpectedRotation, double, Filter );

  /** Set portion of pixels to use to compute similarity in registration */
  tubeWrapSetMacro( SamplingRatio, double, Filter );

  /** Get portion of pixels to use to compute similarity in registration */
  tubeWrapGetMacro( SamplingRatio, double, Filter );

  /** Set if overlapping pixels should be averaged insted of blending */
  tubeWrapSetMacro( BlendUsingAverage, bool, Filter );

  /** Get if overlapping pixels should be averaged instead of blending */
  tubeWrapGetMacro( BlendUsingAverage, bool, Filter );

  /** Set use of experimental method for fast blending */
  tubeWrapSetMacro( UseFastBlending, bool, Filter );

  /** Get use of experimental method for fast blending */
  tubeWrapGetMacro( UseFastBlending, bool, Filter );

  /** Set filename to load the transform from */
  void LoadTransform( const std::string & filename );

  /** Set filename to save the transform to */
  void SaveTransform( const std::string & filename );

  /** Run algorithm */
  tubeWrapUpdateMacro(Filter);

  /** Get output image */
  tubeWrapGetObjectMacro(Output, ImageType, Filter);

protected:
  MergeAdjacentImages( void );
  ~MergeAdjacentImages() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itkMergeAdjacentImagesFilter parameters **/
  MergeAdjacentImages(const Self &);
  void operator=(const Self &);

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  typename FilterType::Pointer m_Filter;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeMergeAdjacentImages.hxx"
#endif

#endif // End !defined( __tubeMergeAdjacentImages_h )
