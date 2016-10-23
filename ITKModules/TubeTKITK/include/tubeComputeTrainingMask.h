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
 *  \ingroup TubeTKITK
 */

template< typename TImage >
class ComputeTrainingMask:
  public itk::ProcessObject
{
  public:
    /** Standard class typedefs. */
    typedef ComputeTrainingMask                             Self;
    typedef itk::ProcessObject                              Superclass;
    typedef itk::SmartPointer< Self >                       Pointer;
    typedef itk::SmartPointer< const Self >                 ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro( Self );

    /** Run-time type information (and related methods). */
    itkTypeMacro( ComputeTrainingMask, ProcessObject );


    /** Typedef to images */
    typedef TImage                                          ImageType;

    itkStaticConstMacro( ImageDimension, unsigned int,
      ImageType::ImageDimension );

    typedef typename itk::tube::ComputeTrainingMaskFilter< ImageType >
                                                            FilterType;
    typedef typename FilterType::ImageTypeShort             ImageTypeShort;

    tubeWrapSetMacro( Gap, double, ComputeTrainingMaskFilter );
    tubeWrapGetMacro( Gap, double, ComputeTrainingMaskFilter );
    tubeWrapSetMacro( NotVesselWidth, double, ComputeTrainingMaskFilter );
    tubeWrapGetMacro( NotVesselWidth, double, ComputeTrainingMaskFilter );
    tubeWrapGetConstObjectMacro( NotVesselMask, ImageTypeShort,
      ComputeTrainingMaskFilter );

    //tubeWrapSetObjectMacro( Input, ImageType,
      //ComputeTrainingMaskFilter );

    //tubeWrapCallMacro( Update, ComputeTrainingMaskFilter );

    //tubeWrapGetObjectMacro( Output, ImageTypeShort,
      //ComputeTrainingMaskFilter);

  protected:
    ComputeTrainingMask( void );
    ~ComputeTrainingMask() {}
    void PrintSelf( std::ostream & os, itk::Indent indent ) const;
  
  private:
    /** itkComputeTrainingMask parameters **/
    ComputeTrainingMask( const Self & );
    void operator=( const Self & );
  
    typename FilterType::Pointer m_ComputeTrainingMaskFilter;

};
} // End namespace tube
#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeComputeTrainingMask.hxx"
#endif
#endif // End !defined( __tubeComputeTrainingMask_h )
