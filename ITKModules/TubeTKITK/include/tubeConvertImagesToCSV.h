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
#ifndef __tubeConvertImagesToCSV_h
#define __tubeConvertImagesToCSV_h

// ITK Includes
#include "itkProcessObject.h"

// TubeTK Includes
#include "tubeWrappingMacros.h"

#include "itktubeConvertImagesToCSVFilter.h"

namespace tube
{
template< class TInputImage, class TInputMask >
class ConvertImagesToCSV:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ConvertImagesToCSV                              Self;
  typedef itk::ProcessObject                              Superclass;
  typedef itk::SmartPointer< Self >                       Pointer;
  typedef itk::SmartPointer< const Self >                 ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ConvertImagesToCSV, ProcessObject );

  typedef TInputMask                                      InputMaskType;
  typedef typename InputMaskType::PixelType               MaskPixelType;

  typedef TInputImage                                     InputImageType;
  typedef typename InputImageType::PixelType              InputPixelType;

  typedef itk::tube::ConvertImagesToCSVFilter< InputImageType, InputMaskType >
    ConvertImagesToCSVFilterType;

  tubeWrapSetObjectMacro( InputMask, InputMaskType, ConvertImagesToCSVFilter );
  tubeWrapGetObjectMacro( InputMask, InputMaskType, ConvertImagesToCSVFilter );
  typename ConvertImagesToCSVFilterType::VnlMatrixType GetOutput();
  tubeWrapGetMacro( Stride, unsigned int, ConvertImagesToCSVFilter );
  tubeWrapSetMacro( Stride, unsigned int, ConvertImagesToCSVFilter );
  tubeWrapSetMacro( NumImages, unsigned int, ConvertImagesToCSVFilter );
  tubeWrapGetMacro( NumImages, unsigned int, ConvertImagesToCSVFilter );
  tubeWrapSetMacro( NumberRows, unsigned int, ConvertImagesToCSVFilter );
  tubeWrapGetMacro( NumberRows, unsigned int, ConvertImagesToCSVFilter );
  /** Set the input image and reinitialize the list of images */
  tubeWrapSetObjectMacro( Input, InputImageType, ConvertImagesToCSVFilter );
  tubeWrapGetConstObjectMacro( Input, InputImageType, ConvertImagesToCSVFilter );
  void AddImage( InputImageType* );

  tubeWrapCallMacro( Update, ConvertImagesToCSVFilter );

protected:
  ConvertImagesToCSV( void );
  ~ConvertImagesToCSV() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itkConvertImagesToCSVFilter parameters **/
  ConvertImagesToCSV( const Self & );
  void operator=( const Self & );

  typename ConvertImagesToCSVFilterType::Pointer m_ConvertImagesToCSVFilter;

};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeConvertImagesToCSV.hxx"
#endif

#endif // End !defined( __tubeCropImage_h )
