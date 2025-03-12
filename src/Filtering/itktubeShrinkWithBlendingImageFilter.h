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

#ifndef itktubeShrinkWithBlendingImageFilter_h
#define itktubeShrinkWithBlendingImageFilter_h

#include "itkShrinkImageFilter.h"

namespace itk
{

namespace tube
{

/** \class ShrinkWithBlendingImageFilter
 * \brief Reduce the size of an image by an integer factor in each
 * dimension.
 *
 * ShrinkWithBlendingImageFilter reduces the size of an image by an integer
 * factor in each dimension. The algorithm implemented is a max over the
 * subsample. The output image size in each dimension is given by:
 *
 * outputSize[j] = max( std::floor( inputSize[j]/shrinkFactor[j] ), 1 );
 *
 * NOTE: The physical centers of the input and output will be the
 * same. Because of this, the Origin of the output may not be the same
 * as the Origin of the input.
 * Since this filter produces an image which is a different
 * resolution, origin and with different pixel spacing than its input
 * image, it needs to override several of the methods defined
 * in ProcessObject in order to properly manage the pipeline execution
 * model. In particular, this filter overrides
 * ProcessObject::GenerateInputRequestedRegion() and
 * ProcessObject::GenerateOutputInformation().
 *
 * This filter is implemented as a multithreaded filter.  It provides a
 * ThreadedGenerateData() method for its implementation.
 *
 * \ingroup GeometricTransform Streamed
 * \ingroup ITKImageGrid
 *
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT ShrinkWithBlendingImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class type alias. */
  using Self = ShrinkWithBlendingImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro(ShrinkWithBlendingImageFilter);

  /** Typedef to images */
  using OutputImageType = TOutputImage;
  using InputImageType = TInputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageConstPointer = typename InputImageType::ConstPointer;

  using OutputIndexType = typename TOutputImage::IndexType;
  using InputIndexType = typename TInputImage::IndexType;
  using OutputOffsetType = typename TOutputImage::OffsetType;

  using OutputImageRegionType = typename TOutputImage::RegionType;
  using InputSizeType = typename TInputImage::SizeType;

  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  using PointImagePixelType = Vector<float, ImageDimension>;
  using PointImageType = Image<PointImagePixelType, OutputImageDimension>;

  using ShrinkFactorsType = FixedArray<unsigned int, ImageDimension>;

  /** Set the shrink factors. Values are clamped to
   * a minimum value of 1. Default is 1 for all dimensions. */
  itkSetMacro(ShrinkFactors, ShrinkFactorsType);
  void
  SetShrinkFactor(unsigned int i, unsigned int factor);
  unsigned int
  GetShrinkFactor(unsigned int i);
  itkGetMacro(NewSize, InputSizeType);
  itkSetMacro(NewSize, InputSizeType);

  /** Get the shrink factors. */
  itkGetConstReferenceMacro(ShrinkFactors, ShrinkFactorsType);

  itkSetMacro(Overlap, InputIndexType);
  itkGetMacro(Overlap, InputIndexType);

  itkSetMacro(BlendWithMean, bool);
  itkGetMacro(BlendWithMean, bool);

  itkSetMacro(BlendWithMax, bool);
  itkGetMacro(BlendWithMax, bool);

  itkSetMacro(BlendWithGaussianWeighting, bool);
  itkGetMacro(BlendWithGaussianWeighting, bool);

  itkSetMacro(UseLog, bool);
  itkGetMacro(UseLog, bool);

  itkSetConstObjectMacro(InputMipPointImage, PointImageType);
  itkGetConstObjectMacro(InputMipPointImage, PointImageType);

  itkGetModifiableObjectMacro(OutputMipPointImage, PointImageType);

  void
  GenerateOutputInformation(void) override;

  void
  GenerateInputRequestedRegion(void) override;

protected:
  ShrinkWithBlendingImageFilter(void);
  ~ShrinkWithBlendingImageFilter(void) {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;

  void
  UpdateInternalShrinkFactors();

  void
  VerifyInputInformation() const override;

  template <class ArrayType>
  bool
  NotValue(ArrayType array, double val, double tolerance = 0.00001) const;

private:
  ShrinkWithBlendingImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  typename PointImageType::Pointer m_OutputMipPointImage;

  typename PointImageType::ConstPointer m_InputMipPointImage;

  typename TInputImage::IndexType m_Overlap;

  bool m_UseLog;

  bool m_BlendWithMean;
  bool m_BlendWithMax;
  bool m_BlendWithGaussianWeighting;

  ShrinkFactorsType m_ShrinkFactors;
  ShrinkFactorsType m_InternalShrinkFactors;
  double            m_DefaultShrinkFactor;
  double            m_DefaultNewSize;
  InputSizeType     m_NewSize;
};

} // end namespace tube

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeShrinkWithBlendingImageFilter.hxx"
#endif

#endif
