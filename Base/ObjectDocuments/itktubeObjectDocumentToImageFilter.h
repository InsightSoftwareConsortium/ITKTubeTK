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
 * Builds image from Object Document
 * Reads and composes all transforms from Document for single object
 *
 * Note: Does not hold buffers of objects read (single read function)
 */

template< class TInputObjectDocument, class TOutputImageType >
class ITK_EXPORT ObjectDocumentToImageFilter :
  public ObjectDocumentToObjectSource<TInputObjectDocument, TOutputImageType::ImageDimension>
{
public:

  enum { TDimensions = 3 };

  typedef ObjectDocumentToImageFilter                       Self;
  typedef ObjectDocumentToObjectSource<TInputObjectDocument,
    TOutputImageType::ImageDimension>                       Superclass;

  typedef SmartPointer< Self >                              Pointer;
  typedef SmartPointer< const Self >                        ConstPointer;

  typedef TInputObjectDocument                              DocumentType;
  typedef typename DocumentType::Pointer                    DocumentPointer;
  typedef typename DocumentType::ConstPointer               ConstDocumentPointer;

  typedef TOutputImageType                                  OutputImageType;
  typedef typename OutputImageType::Pointer                 OutputImagePointer;

  typedef typename Superclass::TransformType                TransformType;
  typedef typename TransformType::Pointer                   TransformPointer;

  /** Interpolator function Type def. for image resampling (default is Linear) */
  typedef InterpolateImageFunction<OutputImageType, double> InterpolateImageFunctionType;
  typedef typename InterpolateImageFunctionType::Pointer    InterpolatorImageFunctionPointer;

  /** Set a pointer to the interpolator function. (default is linear) */
  itkSetObjectMacro( Interpolator, InterpolateImageFunctionType );

  /** Get a pointer to the interpolator function. (default is linear) */
  itkGetObjectMacro( Interpolator, InterpolateImageFunctionType );

  itkNewMacro( Self );
  itkTypeMacro( Self, Superclass );

  OutputImageType * GetOutput( void );

protected:

  void GenerateData( void );

  typedef ImageFileReader< OutputImageType >                        ImageFileReaderType;
  typedef ResampleImageFilter<OutputImageType,OutputImageType>      ResampleImageFilterType;
  typedef typename OutputImageType::SizeType                        SizeType;
  typedef typename OutputImageType::PointType                       PointType;

  /** read from Document file the defined image */
  OutputImagePointer ReadDocument( ConstDocumentPointer );

  /** Resample the input image using the inputted transform */
  OutputImagePointer ResampleImage( OutputImagePointer, TransformPointer );

  ObjectDocumentToImageFilter( void );
  ~ObjectDocumentToImageFilter( void ) {}

private:

  /** Get the transformed bounding box that contains all corners of the original image */
  void GetTransformedBoundingBox( OutputImagePointer image, TransformPointer transform, SizeType &, PointType &) const;

  InterpolatorImageFunctionPointer m_Interpolator;

}; // End class ObjectDocumentToImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeObjectDocumentToImageFilter.hxx"
#endif

#endif // End !defined(__itktubeObjectDocumentToImageFilter_h)
