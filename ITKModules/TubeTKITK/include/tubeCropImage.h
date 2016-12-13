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
#ifndef __tubeCropImage_h
#define __tubeCropImage_h

#include "itkProcessObject.h"

#include "tubeWrappingMacros.h"

#include "itktubeCropImageFilter.h"


namespace tube
{
/** \class CropImage
 *
 *  \ingroup TubeTKITK
 */

template< typename TInputImage, typename TOutputImage >
class CropImage:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef CropImage                       Self;
  typedef itk::ProcessObject              Superclass;
  typedef itk::SmartPointer< Self >       Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( CropImage, ProcessObject );


  /** Typedef to images */
  typedef TInputImage                          InputImageType;
  typedef TOutputImage                         OutputImageType;
  typedef typename InputImageType::IndexType   InputIndexType;
  typedef typename InputImageType::SizeType    InputSizeType;

  typedef itk::tube::CropImageFilter< InputImageType,
    OutputImageType >                          FilterType;

  tubeWrapSetMacro( Min, InputIndexType, Filter );
  tubeWrapGetMacro( Min, InputIndexType, Filter );

  tubeWrapSetMacro( Max, InputIndexType, Filter );
  tubeWrapGetMacro( Max, InputIndexType, Filter );

  tubeWrapSetMacro( Size, InputSizeType, Filter );
  tubeWrapGetMacro( Size, InputSizeType, Filter );

  tubeWrapSetMacro( Center, InputIndexType, Filter );
  tubeWrapGetMacro( Center, InputIndexType, Filter );

  tubeWrapSetMacro( Boundary, InputIndexType, Filter );
  tubeWrapGetMacro( Boundary, InputIndexType, Filter );

  tubeWrapForceSetConstObjectMacro( MatchVolume, InputImageType, Filter );

  tubeWrapForceSetConstObjectMacro( MatchMask, InputImageType, Filter );

  void SetSplitInput( InputIndexType splitIndex, InputIndexType roiIndex );

  tubeWrapSetConstObjectMacro( Input, InputImageType, Filter );
  tubeWrapGetConstObjectMacro( Input, InputImageType, Filter );

  tubeWrapCallMacro( Update, Filter );

  tubeWrapGetObjectMacro( Output, OutputImageType, Filter );

protected:
  CropImage( void );
  ~CropImage() {}

  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itkCropImageFilter parameters **/
  CropImage( const Self & );

  void operator = ( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  typename FilterType::Pointer     m_Filter;

};
} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeCropImage.hxx"
#endif

#endif // End !defined( __tubeCropImage_h )
