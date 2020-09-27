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
#ifndef __tubeComputeTrainingMask_h
#define __tubeComputeTrainingMask_h

// ITK Includes
#include "itkProcessObject.h"

// TubeTK Includes
#include "tubeWrappingMacros.h"

#include "itktubeComputeTrainingMaskFilter.h"

namespace tube
{
/** \class ComputeTrainingMask
 *
 *  \ingroup TubeTK
 */

template< class TImage, class TLabelMap=itk::Image< typename TImage::PixelType,
  TImage::ImageDimension> >
class ComputeTrainingMask : public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ComputeTrainingMask                             Self;
  typedef itk::ProcessObject                              Superclass;
  typedef itk::SmartPointer< Self >                       Pointer;
  typedef itk::SmartPointer< const Self >                 ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ComputeTrainingMask, ProcessObject );


  /** Typedef to images */
  typedef TImage                                          ImageType;
  typedef TLabelMap                                       LabelMapType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    ImageType::ImageDimension );

  typedef typename itk::tube::ComputeTrainingMaskFilter< ImageType, LabelMapType >
                                                          FilterType;

  tubeWrapSetMacro( Gap, double, Filter );
  tubeWrapGetMacro( Gap, double, Filter );
  tubeWrapSetMacro( ObjectWidth, double, Filter );
  tubeWrapGetMacro( ObjectWidth, double, Filter );
  tubeWrapSetMacro( NotObjectWidth, double, Filter );
  tubeWrapGetMacro( NotObjectWidth, double, Filter );
  tubeWrapGetConstObjectMacro( ObjectMask, LabelMapType, Filter );
  tubeWrapGetConstObjectMacro( NotObjectMask, LabelMapType, Filter );

  tubeWrapSetObjectMacro( Input, ImageType, Filter );

  tubeWrapUpdateMacro( Filter );

  tubeWrapGetObjectMacro( Output, LabelMapType, Filter );

protected:
  ComputeTrainingMask( void );
  ~ComputeTrainingMask() {}

  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

private:
  /** itkComputeTrainingMask parameters **/
  ComputeTrainingMask( const Self & );

  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) override
    {};

  typename FilterType::Pointer m_Filter;

};

} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeComputeTrainingMask.hxx"
#endif

#endif // End !defined( __tubeComputeTrainingMask_h )
