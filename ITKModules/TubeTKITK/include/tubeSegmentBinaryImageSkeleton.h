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
#ifndef __tubeSegmentBinaryImageSkeleton_h
#define __tubeSegmentBinaryImageSkeleton_h

// ITK includes
#include "itkProcessObject.h"

// TubeTK includes
#include "tubeWrappingMacros.h"

#include "itktubeSegmentBinaryImageSkeleton.h"

namespace tube
{
/** \class SegmentBinaryImageSkeleton
 *  \brief Computes skeleton of a binary image.
 *  The output skeleton can be dilated if a radius greater than zero is
 *  provided
 *
 *  \ingroup TubeTKITK
 */

template< class TPixel, unsigned int VDimension >
class SegmentBinaryImageSkeleton:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef SegmentBinaryImageSkeleton                 Self;
  typedef itk::ProcessObject                         Superclass;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef itk::tube::SegmentBinaryImageSkeleton<
    TPixel, VDimension >                             FilterType;

  typedef typename FilterType::ImageType             ImageType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SegmentBinaryImageSkeleton, ProcessObject);

  /** Set/Get radius for post-dilatation */
  tubeWrapSetMacro( Radius, unsigned int, Filter );
  tubeWrapGetMacro( Radius, unsigned int, Filter );

  /** Set/Get input image */
  tubeWrapSetConstObjectMacro( Input, ImageType, Filter );
  tubeWrapGetConstObjectMacro( Input, ImageType, Filter );

  /** Compute image similarity */
  tubeWrapUpdateMacro(Filter);

  /** Get image similarity */
  tubeWrapGetObjectMacro(Output, ImageType, Filter);

protected:
  SegmentBinaryImageSkeleton( void );
  ~SegmentBinaryImageSkeleton() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itktubeSegmentBinaryImageSkeletonFilter parameters **/
  SegmentBinaryImageSkeleton(const Self &);
  void operator=(const Self &);

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  typename FilterType::Pointer m_Filter;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeSegmentBinaryImageSkeleton.hxx"
#endif

#endif // End !defined( __tubeSegmentBinaryImageSkeleton_h )
