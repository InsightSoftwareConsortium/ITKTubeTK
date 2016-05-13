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
#ifndef __tubeCropImage_h
#define __tubeCropImage_h

#include "itktubeCropImageFilter.h"
#include "itkObject.h"


namespace tube
{
/** \class CropImage
 *
 *  \ingroup TubeTKITK
 */

template< typename TInputImage, typename TOutputImage >
class CropImage:
  public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef CropImage                       Self;
  typedef itk::SmartPointer< Self >       Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CropImage, Object);

  typedef TInputImage ImageType;

  void SetMin( typename ImageType::IndexType roiMin );

  void SetMax( typename ImageType::IndexType roiMax );

  void SetSize( typename ImageType::SizeType roiSize );

  void SetCenter( typename ImageType::IndexType roiCenter );

  void SetBoundary( typename ImageType::IndexType roiBoundary );

  void SetMatchVolume( typename ImageType::ConstPointer matchVolume );

  void SetMatchMask( typename ImageType::Pointer maskImage );

  void SetSplitInput( typename ImageType::IndexType splitIndex,
    typename ImageType::IndexType roiIndex );

  void SetInput( const TInputImage *inputImage );
  void Update();
  typename TOutputImage::Pointer GetOutput();

protected:
  CropImage( void );
  ~CropImage() {}
  void PrintSelf(std::ostream & os, itk::Indent indent) const;

private:
  /** itkCropImageFilter parameters **/
  CropImage(const Self &);
  void operator=(const Self &);

  typedef itk::tube::CropImageFilter< ImageType, ImageType > CropFilterType;
  typename CropFilterType::Pointer m_CropFilter;

};
} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeCropImage.hxx"
#endif

#endif // End !defined( __tubeCropImage_h )
