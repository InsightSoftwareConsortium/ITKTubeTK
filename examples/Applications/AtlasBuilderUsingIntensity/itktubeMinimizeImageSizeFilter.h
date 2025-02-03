/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#ifndef __itktubeMinimizeImageSizeFilter_h
#define __itktubeMinimizeImageSizeFilter_h

#include <itkImageSliceIteratorWithIndex.h>
#include <itkImageToImageFilter.h>
#include <itkResampleImageFilter.h>

namespace itk
{

namespace tube
{

/** \class MinimizeImageSizeFilter
 * \brief Cuts the edges off of an image ( either from the origin, end
 * extend or both for all dimensions ).
 *
 * Can clip the image to the first valid pixel ( i.e., > threshold ) and/or
 * create a buffer around the image of the first valid pixel
 * ( i.e., pixel + buffer ).
 *
 * Done for all dimensions...Good Filter to follow \sa
 * itkCompleteImageResampleFilter.h
 */
template< class TInputImage >
class MinimizeImageSizeFilter
: public ImageToImageFilter< TInputImage, TInputImage >
{
public:

  typedef MinimizeImageSizeFilter                          Self;
  typedef ImageToImageFilter< TInputImage, TInputImage >   Superclass;

  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;

  itkNewMacro( Self );
  itkOverrideGetNameOfClassMacro( MinimizeImageSizeFilter);

  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::PixelType      InputPixelType;

  typedef TInputImage                             OutputImageType;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename OutputImageType::Pointer       OutputImagePointer;

  typedef typename InputImageType::RegionType     RegionType;
  typedef typename InputImageType::SizeType       SizeType;
  typedef typename InputImageType::IndexType      IndexType;
  typedef typename InputImageType::PointType      PointType;

  /** Number of dimensions. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Get/Set the input */
  itkGetConstObjectMacro( Input, InputImageType );

  /** Set the image pixel buffer on each side of the image ( can be negative
   * to crop inside ) */
  itkGetConstReferenceMacro( NumberOfBufferPixels, SizeType );

  /** Set the number of pixels ( per end ) that will be added before and past
   * the first and last valid pixels */
  void SetNumberOfBufferPixels( SizeType& size )
    {
    m_NumberOfBufferPixels = size;
    this->SetBufferImage( true );
    }

  /** Get the threshold marker for considering a pixel to be kept in image
   * ( value non-inclusive -- threshold value to be removed ) */
  itkGetConstMacro( ThresholdValue, InputPixelType );
  /** Set the threshold marker for considering a pixel to be kept in image
    * ( value non-inclusive -- threshold value to be removed ) */
  itkSetMacro( ThresholdValue, InputPixelType );
  /** Get whether the threshold value is upper or lower boundry
    * ( default is false: i.e., all values below m_ThresholdValue will be
    * excluded ) */
  itkGetConstMacro( ThresholdAbove, bool );
  /** Set whether the threshold value is upper or lower boundry
    *( default is false: i.e., all values below m_ThresholdValue will be
    excluded ) */
  itkSetMacro( ThresholdAbove, bool );

  /** Get the default pixel value when resampling: ( default is 0 ) */
  itkGetConstMacro( DefaultPixelValue, InputPixelType );
  /** Set the default pixel value when resampling: ( default is 0 ) */
  itkSetMacro( DefaultPixelValue, InputPixelType );

  /** Turn on and off the clip the end dimension size function.
    * Will clip the end of each dimension in the image to the first
    * incidence of a pixel > threshold */
  itkBooleanMacro( ClipEndIndices );
  /** Turn on and off the clip the start dimension size function.
    * Will clip the start of each dimension in the image to the first
    * incidence of a pixel > threshold */
  itkBooleanMacro( ClipStartIndices );

  itkGetModifiableObjectMacro( Output, OutputImageType );

protected:

  /** Does the real work! */
  virtual void GenerateData( void ) override;

  MinimizeImageSizeFilter( void );
  ~MinimizeImageSizeFilter( void ) {}

  itkGetConstMacro( ClipEndIndices, bool );
  itkSetMacro( ClipEndIndices, bool );

  itkGetConstMacro( ClipStartIndices, bool );
  itkSetMacro( ClipStartIndices, bool );

  itkGetConstMacro( BufferImage, bool );
  itkSetMacro( BufferImage, bool );

  void Get3DCroppedStartRegion( InputImageConstPointer input,
    RegionType& region );
  void Get3DCroppedEndRegion( InputImageConstPointer input,
    RegionType& region );

private:

  InputImageConstPointer                  m_Input;
  OutputImagePointer                      m_Output;

  SizeType                                m_NumberOfBufferPixels;
  bool                                    m_BufferImage;
  InputPixelType                          m_ThresholdValue;
  bool                                    m_ThresholdAbove;

  InputPixelType                          m_DefaultPixelValue;

  bool                                    m_ClipEndIndices;
  bool                                    m_ClipStartIndices;

}; // End class MinimizeImageSizeFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeMinimizeImageSizeFilter.hxx"
#endif

#endif // End !defined( __itktubeMinimizeImageSizeFilter_h )
