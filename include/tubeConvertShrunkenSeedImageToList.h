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
#ifndef __tubeConvertShrunkenSeedImageToList_h
#define __tubeConvertShrunkenSeedImageToList_h

// ITK Includes
#include "itkProcessObject.h"

// TubeTK Includes
#include "tubeWrappingMacros.h"

#include "itktubeConvertShrunkenSeedImageToListFilter.h"

namespace tube
{
/** \class ConvertShrunkenSeedImageToList
 *
 *  \ingroup TubeTK
 */
template< class TImage, class TPointsImage>
class ConvertShrunkenSeedImageToList
  : public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ConvertShrunkenSeedImageToList                        Self;
  typedef itk::ProcessObject                                    Superclass;
  typedef itk::SmartPointer< Self >                             Pointer;
  typedef itk::SmartPointer< const Self >                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ConvertShrunkenSeedImageToList, ProcessObject );

  typedef TImage                                           ImageType;
  typedef typename ImageType::PixelType                    PixelType;

  typedef TPointsImage                                     PointsImageType;
  typedef typename PointsImageType::PixelType              PointsPixelType;

  typedef itk::tube::ConvertShrunkenSeedImageToListFilter< ImageType,
    PointsImageType > ConvertShrunkenSeedImageToListFilterType;

  tubeWrapSetConstObjectMacro( Input, ImageType,
    ConvertShrunkenSeedImageToListFilter );
  tubeWrapGetConstObjectMacro( Input, ImageType,
    ConvertShrunkenSeedImageToListFilter );
  tubeWrapSetConstObjectMacro( ScaleImage, ImageType,
    ConvertShrunkenSeedImageToListFilter );
  tubeWrapGetConstObjectMacro( ScaleImage, ImageType,
    ConvertShrunkenSeedImageToListFilter );
  tubeWrapSetConstObjectMacro( PointsImage, PointsImageType,
    ConvertShrunkenSeedImageToListFilter );
  tubeWrapGetConstObjectMacro( PointsImage, PointsImageType,
    ConvertShrunkenSeedImageToListFilter );

  typename ConvertShrunkenSeedImageToListFilterType::VnlMatrixType GetOutput();

  tubeWrapGetMacro( Threshold, double, ConvertShrunkenSeedImageToListFilter );
  tubeWrapSetMacro( Threshold, double, ConvertShrunkenSeedImageToListFilter );

  tubeWrapUpdateMacro( ConvertShrunkenSeedImageToListFilter );

protected:
  ConvertShrunkenSeedImageToList( void );
  ~ConvertShrunkenSeedImageToList() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

private:
  /** itkConvertShrunkenSeedImageToListFilter parameters **/
  ConvertShrunkenSeedImageToList( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) override
    {};

  typename ConvertShrunkenSeedImageToListFilterType::Pointer
    m_ConvertShrunkenSeedImageToListFilter;
};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeConvertShrunkenSeedImageToList.hxx"
#endif

#endif // End !defined( __tubeConvertShrunkenSeedImageToList_h )
