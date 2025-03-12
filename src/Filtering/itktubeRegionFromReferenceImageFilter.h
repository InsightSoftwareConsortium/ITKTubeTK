/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itktubeRegionFromReferenceImageFilter_h
#define itktubeRegionFromReferenceImageFilter_h

#include "itkExtractImageFilter.h"

namespace itk
{

namespace tube
{

/** \class RegionFromReferenceImageFilter
 * \brief Decrease the image size by cropping the image by an itk::Size at
 * both the upper and lower bounds of the largest possible region.
 *
 * RegionFromReferenceImageFilter changes the image boundary of an image by
 * removing pixels outside the target region.   The target region is not
 * specified in advance, but calculated in BeforeThreadedGenerateData().
 *
 * This filter uses ExtractImageFilter to perform the cropping.
 *
 * \ingroup GeometricTransforms
 */
template <class TInputImage, class TOutputImage = TInputImage>
class ITK_EXPORT RegionFromReferenceImageFilter : public ExtractImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class type alias. */
  using Self = RegionFromReferenceImageFilter;
  using Superclass = ExtractImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro(RegionFromReferenceImageFilter);

  /** Typedef to describe the output and input image region types. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;
  using InputImageRegionType = typename Superclass::InputImageRegionType;

  /** Typedef to describe the type of pixel. */
  using OutputImagePixelType = typename Superclass::OutputImagePixelType;
  using InputImagePixelType = typename Superclass::InputImagePixelType;

  /** Typedef to describe the output and input image index and size types. */
  using OutputImageIndexType = typename Superclass::OutputImageIndexType;
  using InputImageIndexType = typename Superclass::InputImageIndexType;
  using OutputImageSizeType = typename Superclass::OutputImageSizeType;
  using InputImageSizeType = typename Superclass::InputImageSizeType;
  using SizeType = InputImageSizeType;

  /** ImageDimension constants */
  static constexpr unsigned int InputImageDimension = Superclass::InputImageDimension;
  static constexpr unsigned int OutputImageDimension = Superclass::OutputImageDimension;
  static constexpr unsigned int ImageDimension = Superclass::OutputImageDimension;

  using ReferenceImageType = ImageBase<Self::ImageDimension>;

  /** Copy the output information from another Image. */
  void
  SetReferenceImage(const ReferenceImageType * image);

  const ReferenceImageType *
  GetReferenceImage() const;

  /** Set the input image */
  void
  SetInput1(const TInputImage * input)
  {
    this->SetInput(input);
  }

  /** Set the reference image */
  void
  SetInput2(const ReferenceImageType * input)
  {
    this->SetReferenceImage(input);
  }


#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToOutputCheck, (Concept::Convertible<InputImagePixelType, OutputImagePixelType>));
  itkConceptMacro(SameDimensionCheck, (Concept::SameDimension<InputImageDimension, OutputImageDimension>));
  /** End concept checking */
#endif

protected:
  RegionFromReferenceImageFilter() { this->SetNumberOfRequiredInputs(2); }
  ~RegionFromReferenceImageFilter() {}

  void
  GenerateOutputInformation() override;

private:
  RegionFromReferenceImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented
};


} // end namespace tube

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeRegionFromReferenceImageFilter.hxx"
#endif

#endif
