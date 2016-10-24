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
#ifndef __tubeConvertTubesToImage_h
#define __tubeConvertTubesToImage_h

// ITK includes
#include <itkGroupSpatialObject.h>
#include <itkMacro.h>
#include <itkProcessObject.h>


// TubeTK includes
#include "tubeWrappingMacros.h"

#include "itktubeTubeSpatialObjectToImageFilter.h"

namespace tube
{
/** \class ConvertTubesToImage
 *
 *  \ingroup TubeTKITK
 */

template< unsigned int Dimension, class TOutputPixel >
class ConvertTubesToImage:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ConvertTubesToImage                        Self;
  typedef itk::ProcessObject                         Superclass;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef itk::Image< TOutputPixel, Dimension >      OutputImageType;
  typedef itk::tube::TubeSpatialObjectToImageFilter< Dimension,
    OutputImageType >                                FilterType;
  typedef typename FilterType::SpatialObjectType     TubesType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ConvertTubesToImage, ProcessObject);

  /** Set if the tube should be full inside */
  tubeWrapSetMacro(UseRadius, bool, Filter);
  tubeWrapGetMacro(UseRadius, bool, Filter);

  /* Set template image */
  void SetTemplateImage(const OutputImageType * pTemplateImage);
  itkGetConstObjectMacro(TemplateImage, OutputImageType);

  /* Set input tubes */
  tubeWrapSetConstObjectMacro(Input, TubesType, Filter);
  tubeWrapGetConstObjectMacro(Input, TubesType, Filter);

  /* Runs tubes to image conversion */
  tubeWrapUpdateMacro(Filter);

  /* Get the generated binary tubes image */
  tubeWrapGetObjectMacro(Output, OutputImageType, Filter);

protected:
  ConvertTubesToImage( void );
  ~ConvertTubesToImage() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itkConvertTubesToImageFilter parameters **/
  ConvertTubesToImage(const Self &);
  void operator=(const Self &);

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  typename FilterType::Pointer m_Filter;

  typename OutputImageType::ConstPointer m_TemplateImage;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeConvertTubesToImage.hxx"
#endif

#endif // End !defined( __tubeConvertTubesToImage_h )
