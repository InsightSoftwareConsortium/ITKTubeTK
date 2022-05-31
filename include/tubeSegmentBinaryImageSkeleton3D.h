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
#ifndef __tubeSegmentBinaryImageSkeleton3D_h
#define __tubeSegmentBinaryImageSkeleton3D_h

// ITK includes
#include "itkProcessObject.h"

// TubeTK includes
#include "tubeWrappingMacros.h"

#include "itktubeBinaryThinningImageFilter3D.h"

namespace tube
{
/** \class SegmentBinaryImageSkeleton3D
 *  \brief Computes skeleton of a binary image.
 *  The output skeleton can be dilated if a radius greater than zero is
 *  provided
 *
 *  \ingroup TubeTK
 */

template< class TImageType >
class SegmentBinaryImageSkeleton3D:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef SegmentBinaryImageSkeleton3D               Self;
  typedef itk::ProcessObject                         Superclass;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef TImageType                                 ImageType;

  typedef itk::tube::BinaryThinningImageFilter3D<
    TImageType, TImageType >                         FilterType;

  typedef typename FilterType::EndPointListType      EndPointListType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( SegmentBinaryImageSkeleton3D, ProcessObject );

  /** Set/Get input image */
  tubeWrapSetConstObjectMacro( Input, ImageType, Filter );
  tubeWrapGetConstObjectMacro( Input, ImageType, Filter );

  /** Compute image similarity */
  tubeWrapUpdateMacro( Filter );

  /** Get image similarity */
  tubeWrapGetObjectMacro( Output, ImageType, Filter );
  tubeWrapGetObjectMacro( Thinning, ImageType, Filter );

  tubeWrapGetConstReferenceMacro( EndPoints, EndPointListType, Filter );

protected:
  SegmentBinaryImageSkeleton3D( void );
  ~SegmentBinaryImageSkeleton3D() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itktubeSegmentBinaryImageSkeleton3DFilter parameters **/
  SegmentBinaryImageSkeleton3D( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  typename FilterType::Pointer m_Filter;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeSegmentBinaryImageSkeleton3D.hxx"
#endif

#endif // End !defined( __tubeSegmentBinaryImageSkeleton3D_h )
