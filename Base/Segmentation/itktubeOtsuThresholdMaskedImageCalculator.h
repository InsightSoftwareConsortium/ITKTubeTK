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

#ifndef __itktubeOtsuThresholdMaskedImageCalculator_h
#define __itktubeOtsuThresholdMaskedImageCalculator_h

#include <itkNumericTraits.h>
#include <itkObject.h>
#include <itkObjectFactory.h>

namespace itk
{

namespace tube
{

/** \class OtsuThresholdMaskedImageCalculator
 * \brief Computes the Otsu's threshold for an image.
 *
 * This calculator computes the Otsu's threshold which separates an image
 * into foreground and background components. The method relies on a
 * histogram of image intensities. The basic idea is to maximize the
 * between-class variance.
 *
 * This class is templated over the input image type.
 *
 * \warning This method assumes that the input image consists of scalar pixel
 * types.
 *
 * \ingroup Operators
 */
template <class TInputImage>
class OtsuThresholdMaskedImageCalculator : public Object
{
public:
  /** Standard class typedefs. */
  typedef OtsuThresholdMaskedImageCalculator Self;
  typedef Object                             Superclass;
  typedef SmartPointer<Self>                 Pointer;
  typedef SmartPointer<const Self>           ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(OtsuThresholdMaskedImageCalculator, Object);

  /** Type definition for the input image. */
  typedef TInputImage  ImageType;

  /** Pointer type for the image. */
  typedef typename TInputImage::Pointer  ImagePointer;

  /** Const Pointer type for the image. */
  typedef typename TInputImage::ConstPointer ImageConstPointer;

  /** Type definition for the input image pixel type. */
  typedef typename TInputImage::PixelType PixelType;

  /** Type definition for the input image region type. */
  typedef typename TInputImage::RegionType RegionType;

  /** Set the input image. */
  itkSetConstObjectMacro(Image,ImageType);

  /** Set the input image. */
  itkSetConstObjectMacro(MaskImage,ImageType);

  /** Compute the Otsu's threshold for the input image. */
  void Compute( void );

  /** Return the Otsu's threshold value. */
  itkGetConstMacro(Threshold,PixelType);

  /** Set/Get the number of histogram bins. Default is 128. */
  itkSetClampMacro( NumberOfHistogramBins, unsigned long, 1,
                    NumericTraits<unsigned long>::max() );
  itkGetConstMacro( NumberOfHistogramBins, unsigned long );

  /** Set the region over which the values will be computed */
  void SetRegion( const RegionType & region );

protected:
  OtsuThresholdMaskedImageCalculator( void );
  virtual ~OtsuThresholdMaskedImageCalculator( void ) {}
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  //purposely not implemented
  OtsuThresholdMaskedImageCalculator(const Self&);
  void operator=(const Self&); //purposely not implemented

  PixelType            m_Threshold;
  unsigned long        m_NumberOfHistogramBins;
  ImageConstPointer    m_Image;
  ImageConstPointer    m_MaskImage;
  RegionType           m_Region;
  bool                 m_RegionSetByUser;

}; // End class OtsuThresholdMaskedImageCalculator

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeOtsuThresholdMaskedImageCalculator.hxx"
#endif

#endif // End !defined(__itktubeOtsuThresholdMaskedImageCalculator_h)
