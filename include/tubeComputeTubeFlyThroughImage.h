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
#ifndef __tubeComputeTubeFlyThroughImage_h
#define __tubeComputeTubeFlyThroughImage_h

// ITK includes
#include <itkProcessObject.h>
#include <itkGroupSpatialObject.h>

// TubeTK includes
#include "tubeWrappingMacros.h"

#include "itktubeComputeTubeFlyThroughImageFilter.h"

namespace tube
{
/** \class ComputeTubeFlyThroughImage
 *
 *  \ingroup TubeTK
 */

template< class TPixel, unsigned int Dimension >
class ComputeTubeFlyThroughImage:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ComputeTubeFlyThroughImage                 Self;
  typedef itk::ProcessObject                         Superclass;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef itk::tube::ComputeTubeFlyThroughImageFilter< TPixel,
    Dimension >                                     FilterType;

  typedef typename FilterType::TubeGroupType        TubeGroupType;
  typedef typename FilterType::TubeType             TubeType;
  typedef typename FilterType::InputImageType       InputImageType;
  typedef typename FilterType::OutputImageType      OutputImageType;
  typedef typename FilterType::OutputMaskType       OutputMaskType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ComputeTubeFlyThroughImage, ProcessObject );

  /** Set/Get tube id for which the fly through image is to be generated */
  tubeWrapSetMacro( TubeId, unsigned long, Filter );
  tubeWrapGetMacro( TubeId, unsigned long, Filter );

  /* Set/Get input image from which the tubes were extracted/segmented */
  tubeWrapSetConstObjectMacro( InputImage, InputImageType, Filter );
  tubeWrapGetConstObjectMacro( InputImage, InputImageType, Filter );

  /* Set/Get input tubes */
  tubeWrapSetConstObjectMacro( Input, TubeGroupType, Filter );
  tubeWrapGetConstObjectMacro( Input, TubeGroupType, Filter );

  /* Generates tube fly through image and mask */
  tubeWrapUpdateMacro( Filter );

  /* Get the generated tube fly through image */
  tubeWrapGetObjectMacro( Output, OutputImageType, Filter );

  /* Get the generated tube fly through image */
  tubeWrapGetObjectMacro( OutputMask, OutputMaskType, Filter );

protected:
  ComputeTubeFlyThroughImage( void );
  ~ComputeTubeFlyThroughImage() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itkComputeTubeFlyThroughImageFilter parameters **/
  ComputeTubeFlyThroughImage( const Self & );

  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  typename FilterType::Pointer m_Filter;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeComputeTubeFlyThroughImage.hxx"
#endif

#endif // End !defined( __tubeComputeTubeFlyThroughImage_h )
