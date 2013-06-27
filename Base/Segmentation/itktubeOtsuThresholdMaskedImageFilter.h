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

#ifndef __itktubeOtsuThresholdMaskedImageFilter_h
#define __itktubeOtsuThresholdMaskedImageFilter_h

#include <itkFixedArray.h>
#include <itkImageToImageFilter.h>

namespace itk
{

namespace tube
{

/** \class OtsuThresholdMaskedImageFilter
 * \brief Threshold an image using the Otsu Threshold
 *
 * This filter creates a binary thresholded image that separates an
 * image into foreground and background components. The filter
 * computes the threshold using the OtsuThresholdImageCalculator and
 * applies that theshold to the input image using the
 * BinaryThresholdImageFilter. The NunberOfHistogram bins can be set
 * for the Calculator. The InsideValue and OutsideValue can be set
 * for the BinaryThresholdImageFilter.
 *
 * \sa OtsuThresholdMaskedImageCalculator
 * \sa BinaryThresholdImageFilter
 * \ingroup IntensityImageFilters  Multithreaded
 */

template< class TInputImage, class TOutputImage >
class OtsuThresholdMaskedImageFilter
  : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard Self typedef */
  typedef OtsuThresholdMaskedImageFilter                   Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage >  Superclass;
  typedef SmartPointer< Self >                             Pointer;
  typedef SmartPointer< const Self >                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Runtime information support. */
  itkTypeMacro( OtsuThresholdMaskedImageFilter, ImageToImageFilter );

  /** Image pixel value typedef. */
  typedef typename TInputImage::PixelType   InputPixelType;
  typedef typename TOutputImage::PixelType  OutputPixelType;

  /** Image related typedefs. */
  typedef typename TInputImage::Pointer  InputImagePointer;
  typedef typename TOutputImage::Pointer OutputImagePointer;

  typedef typename TInputImage::SizeType    InputSizeType;
  typedef typename TInputImage::IndexType   InputIndexType;
  typedef typename TInputImage::RegionType  InputImageRegionType;
  typedef typename TOutputImage::SizeType   OutputSizeType;
  typedef typename TOutputImage::IndexType  OutputIndexType;
  typedef typename TOutputImage::RegionType OutputImageRegionType;


  /** Image related typedefs. */
  itkStaticConstMacro( InputImageDimension, unsigned int,
                      TInputImage::ImageDimension );

  itkStaticConstMacro( OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension );

  itkSetObjectMacro( MaskImage, TInputImage );
  itkGetObjectMacro( MaskImage, TInputImage );

  /** Set the "outside" pixel value. The default value
   * NumericTraits<OutputPixelType>::Zero. */
  itkSetMacro( OutsideValue, OutputPixelType );

  /** Get the "outside" pixel value. */
  itkGetConstMacro( OutsideValue, OutputPixelType );

  /** Set the "inside" pixel value. The default value
   * NumericTraits< OutputPixelType >::max(). */
  itkSetMacro( InsideValue, OutputPixelType );

  /** Get the "inside" pixel value. */
  itkGetConstMacro( InsideValue, OutputPixelType );

  /** Set/Get the number of histogram bins. Defaults is 128. */
  itkSetClampMacro( NumberOfHistogramBins, unsigned long, 1,
                    NumericTraits< unsigned long >::max() );
  itkGetConstMacro( NumberOfHistogramBins, unsigned long );

  /** Get the computed threshold. */
  itkGetConstMacro( Threshold, InputPixelType );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( OutputEqualityComparableCheck,
    ( Concept::EqualityComparable< OutputPixelType > ) );
  itkConceptMacro( InputOStreamWritableCheck,
    ( Concept::OStreamWritable< InputPixelType > ) );
  itkConceptMacro( OutputOStreamWritableCheck,
    ( Concept::OStreamWritable< OutputPixelType > ) );
  /** End concept checking */
#endif
protected:
  OtsuThresholdMaskedImageFilter( void );
  ~OtsuThresholdMaskedImageFilter( void ) {};

  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateInputRequestedRegion( void );
  void GenerateData( void );

private:
  OtsuThresholdMaskedImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  InputPixelType      m_Threshold;
  OutputPixelType     m_InsideValue;
  OutputPixelType     m_OutsideValue;
  unsigned long       m_NumberOfHistogramBins;

  typename TInputImage::Pointer m_MaskImage;

}; // End class OtsuThresholdMaskedImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeOtsuThresholdMaskedImageFilter.hxx"
#endif

#endif // End !defined(__itktubeOtsuThresholdMaskedImageFilter_h)
