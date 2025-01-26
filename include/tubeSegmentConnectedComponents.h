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
#ifndef __tubeSegmentConnectedComponents_h
#define __tubeSegmentConnectedComponents_h

// ITK includes
#include "itkProcessObject.h"
#include "itkConnectedComponentImageFilter.h"

// TubeTK includes
#include "tubeWrappingMacros.h"

namespace tube
{
/** \class SegmentConnectedComponents
 *
 *  \ingroup TubeTK
 */

template< class TImage, class TSeedMask >
class SegmentConnectedComponents:
  public itk::ProcessObject
{
public:
  typedef TImage                                ImageType;

  typedef TSeedMask                             SeedMaskType;

  /** Standard class typedefs. */
  typedef SegmentConnectedComponents            Self;
  typedef itk::ProcessObject                    Superclass;
  typedef itk::SmartPointer< Self >             Pointer;
  typedef itk::SmartPointer< const Self >       ConstPointer;


  typedef itk::ConnectedComponentImageFilter< ImageType,
    ImageType, ImageType >                      FilterType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro( SegmentConnectedComponents);

  static constexpr unsigned int ImageDimension = ImageType::ImageDimension ;

  tubeWrapSetObjectMacro( Input, ImageType, Filter );
  tubeWrapGetConstObjectMacro( Input, ImageType, Filter );

  itkGetMacro( NumberOfComponents, unsigned int );

  itkSetMacro( MinimumVolume, double );
  itkGetMacro( MinimumVolume, double );

  itkSetObjectMacro( SeedMask, SeedMaskType );
  itkGetModifiableObjectMacro( SeedMask, SeedMaskType );

  itkSetMacro( KeepOnlyLargestComponent, bool );
  itkGetMacro( KeepOnlyLargestComponent, bool );

  void Update( void ) override;

  tubeWrapGetConstObjectMacro( Output, ImageType, Filter );

protected:
  SegmentConnectedComponents( void );
  ~SegmentConnectedComponents() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

private:
  /** itkSegmentConnectedComponents parameters **/
  SegmentConnectedComponents( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) override
    {};

  typename FilterType::Pointer                    m_Filter;

  typename SeedMaskType::Pointer                  m_SeedMask;
  double                                          m_MinimumVolume;
  unsigned int                                    m_NumberOfComponents;

  bool                                            m_KeepOnlyLargestComponent;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeSegmentConnectedComponents.hxx"
#endif

#endif // End !defined( __tubeSegmentConnectedComponentsUsingParzenPDFs_h )
