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

#ifndef __itktubeObjectDocumentToImageFilter_h
#define __itktubeObjectDocumentToImageFilter_h

#include "itktubeObjectDocumentToObjectSource.h"

#include <itkImageFileReader.h>
#include <itkResampleImageFilter.h>
#include <itkSpatialObjectReader.h>

namespace itk
{

namespace tube
{

/**
 * Filter that takes an object document as input and produces an image a
 * output.
 *
 * \ingroup  ObjectDocuments
 */
template< class TObjectDocument, class TImageType >
class ObjectDocumentToImageFilter
: public ObjectDocumentToObjectSource< TObjectDocument,
  TImageType::ImageDimension >
{
public:

  enum { Dimension = 3 };

  typedef ObjectDocumentToImageFilter                     Self;
  typedef ObjectDocumentToObjectSource< TObjectDocument,
    TImageType::ImageDimension >                           Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  typedef TObjectDocument                      DocumentType;
  typedef typename DocumentType::ConstPointer  ConstDocumentPointer;

  typedef TImageType                           ImageType;
  typedef typename ImageType::Pointer          ImagePointer;

  typedef typename Superclass::TransformType   TransformType;
  typedef typename TransformType::Pointer      TransformPointer;

  typedef InterpolateImageFunction< ImageType, double >
    InterpolateImageFunctionType;
  typedef typename InterpolateImageFunctionType::Pointer
    InterpolateImageFunctionPointer;

  itkNewMacro( Self );
  itkTypeMacro( ObjectDocumentToImageFilter, ObjectDocumentToObjectSource );

  /** Return the interpolator. */
  itkGetObjectMacro( Interpolator, InterpolateImageFunctionType );

  /** Set the interpolator. */
  itkSetObjectMacro( Interpolator, InterpolateImageFunctionType );

  /** Return the output. */
  using Superclass::GetOutput;
  ImageType * GetOutput( void ) override;

protected:

  typedef ImageFileReader< ImageType >                ImageFileReaderType;
  typedef ResampleImageFilter< ImageType, ImageType > ResampleImageFilterType;
  typedef typename ImageType::SizeType                SizeType;
  typedef typename ImageType::PointType               PointType;

  /** Constructor. */
  ObjectDocumentToImageFilter( void );

  /** Destructor. */
  virtual ~ObjectDocumentToImageFilter( void );

  /** Generate the output data. */
  virtual void GenerateData( void ) override;

  /** Read the specified object document. */
  virtual ImagePointer ReadDocument( ConstDocumentPointer image );

  /** Resample the specified image. */
  virtual ImagePointer ResampleImage( ImagePointer image,
                                      TransformPointer transform );

  /** Print information about the object. */
  virtual void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  // Copy constructor not implemented.
  ObjectDocumentToImageFilter( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

  /** Return the transformed bounding box. */
  void GetTransformedBoundingBox( ImagePointer image,
    TransformPointer transform, SizeType & size, PointType & origin ) const;

  InterpolateImageFunctionPointer  m_Interpolator;

}; // End class ObjectDocumentToImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeObjectDocumentToImageFilter.hxx"
#endif

#endif // End !defined( __itktubeObjectDocumentToImageFilter_h )
