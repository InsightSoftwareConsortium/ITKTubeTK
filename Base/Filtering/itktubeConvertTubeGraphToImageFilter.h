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

#ifndef __itktubeConvertTubeGraphToImageFilter_h
#define __itktubeConvertTubeGraphToImageFilter_h

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkProcessObject.h>
#include <vector>

namespace itk
{

namespace tube
{

template< class TInputImage, class TOutputImage >
class ConvertTubeGraphToImageFilter
  : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:

  /** Standard class typedefs. */
  typedef ConvertTubeGraphToImageFilter                   Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage>  Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  itkTypeMacro( ConvertTubeGraphToImageFilter, ImageToImageFilter );

  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  typedef TInputImage                                   InputImageType;
  typedef typename InputImageType::PixelType            InputPixelType;
  typedef TOutputImage                                  OutputImageType;

  /** In Graph File Name*/
  itkGetMacro( InGraphFileName, std::string );
  itkSetMacro( InGraphFileName, std::string );

  itkGetMacro( AImage, typename OutputImageType::Pointer);
  itkGetMacro( BImage, typename OutputImageType::Pointer);
  itkGetMacro( RImage, typename OutputImageType::Pointer);
  itkGetMacro( CImage, typename OutputImageType::Pointer);

protected:
  ConvertTubeGraphToImageFilter( void );
  ~ConvertTubeGraphToImageFilter( void ) {}

  void PrintSelf(std::ostream& os, Indent indent) const;
  virtual void GenerateData( void );

private:
  ConvertTubeGraphToImageFilter(const Self&);
  void operator=(const Self&);
  std::string           m_InGraphFileName;

  typename OutputImageType::Pointer            m_AImage;
  typename OutputImageType::Pointer            m_BImage;
  typename OutputImageType::Pointer            m_RImage;
  typename OutputImageType::Pointer            m_CImage;
  typename InputImageType::ConstPointer        m_InputImage;

}; // End class CVTImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeConvertTubeGraphToImageFilter.hxx"
#endif

#endif // End !defined(_itktubeConvertTubeGraphToImageFilter_h)
