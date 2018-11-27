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
#ifndef __tubeConvertSpatialGraphToImage_h
#define __tubeConvertSpatialGraphToImage_h

// ITK Includes
#include "itkProcessObject.h"

// TubeTK Includes
#include "tubeWrappingMacros.h"

#include "itktubeConvertSpatialGraphToImageFilter.h"

namespace tube
{
/** \class ConvertSpatialGraphToImage
 *
 *  \ingroup TubeTK
 */

template< typename TInputImage, typename TOutputImage >
class ConvertSpatialGraphToImage:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ConvertSpatialGraphToImage      Self;
  typedef ProcessObject                   Superclass;
  typedef itk::SmartPointer< Self >       Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  typedef itk::tube::ConvertSpatialGraphToImageFilter
    < TInputImage, TOutputImage >         FilterType;

  typedef typename FilterType::InputImageType      InputImageType;
  typedef typename FilterType::InputImagePointer   InputImagePointer;
  typedef typename FilterType::OutputImageType     OutputImageType;
  typedef typename FilterType::OutputImagePointer  OutputImagePointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ConvertSpatialGraphToImage, ProcessObject );

  /** Get Adjacency Matrix Image */
  tubeWrapGetMacro( AdjacencyMatrixImage, OutputImagePointer, Filter );
  tubeWrapGetMacro( BranchnessImage, OutputImagePointer, Filter );
  tubeWrapGetMacro( RadiusImage, OutputImagePointer, Filter );
  tubeWrapGetMacro( CentralityImage, OutputImagePointer, Filter );

  /* Set input tubes */
  tubeWrapSetConstObjectMacro( Input, InputImageType, Filter );

  /* Runs tubes to image conversion */
  tubeWrapUpdateMacro( Filter );

  /* Get the generated binary tubes image */
  tubeWrapGetObjectMacro( Output, OutputImageType, Filter );

  void SetAdjacencyMatrix( vnl_matrix< double > );

  void SetBranchnessVector( vnl_vector< double > );

  void SetRadiusVector( vnl_vector< double > );

  void SetCentralityVector( vnl_vector< double > );

protected:
  ConvertSpatialGraphToImage( void );
  ~ConvertSpatialGraphToImage() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itkConvertSpatialGraphToImageFilter parameters **/
  ConvertSpatialGraphToImage( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  typename FilterType::Pointer  m_Filter;
};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeConvertSpatialGraphToImage.hxx"
#endif

#endif // End !defined( __tubeConvertSpatialGraphToImage_h )
