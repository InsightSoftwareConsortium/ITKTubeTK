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

#ifndef __itktubeReResampleImageFilter_h
#define __itktubeReResampleImageFilter_h

#include <itkBSplineInterpolateImageFunction.h>
#include <itkCompensatedSummation.h>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkInterpolateImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkProcessObject.h>
#include <itkResampleImageFilter.h>
#include <itkWindowedSincInterpolateImageFunction.h>

namespace itk
{

namespace tube
{

template< class TPixel, unsigned int VDimension >
class ReResampleImageFilter
  : public itk::ProcessObject
{
public:

  /** Standard class typedefs. */
  typedef Image< TPixel, VDimension >                     ImageType;
  typedef ReResampleImageFilter                           Self;
  typedef itk::ProcessObject                              Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  itkTypeMacro( ReResampleImageFilter, ProcessObject );

  typedef typename ImageType::Pointer                     ImagePointer;
  typedef itk::Transform< double, VDimension, VDimension> TransformType;


  /** Set/Get input Match Image */
  itkSetObjectMacro( MatchImage, ImageType );
  itkGetModifiableObjectMacro( MatchImage, ImageType );

  /** Set/Get whether Output is isotropic or not */
  itkSetMacro( MakeIsotropic, bool );
  itkGetMacro( MakeIsotropic, bool );

  /** Set/Get whether Output is high resolution isotropic or not */
  itkSetMacro( MakeHighResIso, bool );
  itkGetMacro( MakeHighResIso, bool );

  /** Set/Get interpolator */
  itkSetMacro( Interpolator, std::string );
  itkGetMacro( Interpolator, std::string );

  /** Set/Get whether to load transform or not */
  itkSetMacro( LoadTransform, bool );
  itkGetMacro( LoadTransform, bool );

  /** Set Output Transform */
  void SetTransform( TransformType* t );

  /** Set Output Spacing */
  void SetSpacing( std::vector<double> s );

  /** Set Output Origin */
  void SetOrigin( std::vector<double> o );

  /** Set Output Index */
  void SetIndex( std::vector<int> i );

  /** Set Output Resample Factor */
  void SetResampleFactor( std::vector<double> rf );

  itkGetModifiableObjectMacro( Output, ImageType );

  itkSetObjectMacro( Input, ImageType );
  itkGetModifiableObjectMacro( Input, ImageType );

  virtual void Update( void ) override;

protected:
  ReResampleImageFilter( void );
  ~ReResampleImageFilter( void ) {}

  void SetInput( const DataObjectIdentifierType & key, DataObject * input) override
    { Superclass::SetInput( key, input ); };

  void PrintSelf( std::ostream& os, Indent indent ) const override;

private:
  ReResampleImageFilter( const Self& );
  void operator=( const Self& );

  typedef typename itk::ResampleImageFilter< ImageType, ImageType >
    ResampleFilterType;
  typename ResampleFilterType::Pointer   m_Filter;

  typename ImageType::Pointer            m_MatchImage;
  std::vector<double>                    m_Spacing;
  std::vector<double>                    m_Origin;
  std::vector<int>                       m_Index;
  std::vector<double>                    m_ResampleFactor;
  bool                                   m_MakeIsotropic;
  bool                                   m_MakeHighResIso;
  std::string                            m_Interpolator;
  bool                                   m_LoadTransform;
  typename TransformType::Pointer        m_Transform;

  typename ImageType::Pointer            m_Input;
  typename ImageType::Pointer            m_Output;
}; // End class ReResampleImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeReResampleImageFilter.hxx"
#endif

#endif // End !defined( _itktubeReResampleImageFilter_h )
