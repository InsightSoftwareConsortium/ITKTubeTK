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

#ifndef itktubeStructureTensorRecursiveGaussianImageFilter_h
#define itktubeStructureTensorRecursiveGaussianImageFilter_h

#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkNthElementImageAdaptor.h>
#include <itkPixelTraits.h>
#include <itkProgressAccumulator.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkSymmetricSecondRankTensor.h>

namespace itk
{

namespace tube
{

/** \class StructureTensorRecursiveGaussianImageFilter
 * \brief Computes the structure tensor of a multidimensional image
 * \warning Operates in image ( pixel ) space, not physical space
 * \ingroup GradientFilters
 * \ingroup Singlethreaded
 */
template <class TInputImage,
          class TOutputImage =
            Image<SymmetricSecondRankTensor<typename NumericTraits<typename TInputImage::PixelType>::RealType,
                                            TInputImage::ImageDimension>,
                  TInputImage::ImageDimension>>
class StructureTensorRecursiveGaussianImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class type alias. */
  using Self = StructureTensorRecursiveGaussianImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Pixel Type of the input image */
  using InputImageType = TInputImage;
  using PixelType = typename TInputImage::PixelType;
  using RealType = typename NumericTraits<PixelType>::RealType;

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  /** Define the image type for internal computations
      RealType is usually 'double' in NumericTraits.
      Here we prefer float in order to save memory. */
  using InternalRealType = float;
  using RealImageType = Image<InternalRealType, itkGetStaticConstMacro(ImageDimension)>;

  /**  Output Image Nth Element Adaptor
   *  This adaptor allows to use conventional scalar
   *  smoothing filters to compute each one of the
   *  components of the gradient image pixels. */
  using OutputImageAdaptorType = NthElementImageAdaptor<TOutputImage, InternalRealType>;
  typedef typename OutputImageAdaptorType::Pointer OutputImageAdaptorPointer;

  /**  Smoothing filter type */
  using GaussianFilterType = RecursiveGaussianImageFilter<RealImageType, RealImageType>;
  typedef typename GaussianFilterType::Pointer GaussianFilterPointer;

  /**  Derivative filter type, it will be the first in the pipeline  */
  using DerivativeFilterType = RecursiveGaussianImageFilter<InputImageType, RealImageType>;
  typedef typename DerivativeFilterType::Pointer DerivativeFilterPointer;

  /** Type of the output image */
  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputPixelType = typename OutputImageType::PixelType;
  typedef typename PixelTraits<OutputPixelType>::ValueType OutputComponentType;

  /** Set Sigma value. Sigma is measured in the units of image spacing.  */
  void
  SetSigma(RealType sigma);
  void
  SetSigmaOuter(RealType rho);

  /** Define which normalization factor will be used for the Gaussian */
  void
  SetNormalizeAcrossScale(bool normalizeInScaleSpace);
  itkGetMacro(NormalizeAcrossScale, bool);

  // Sigma value for the Gaussian derivative filters
  itkGetMacro(Sigma, RealType);

  // Sigma value for the outer Gaussian smoothing filter
  itkGetMacro(SigmaOuter, RealType);


protected:
  StructureTensorRecursiveGaussianImageFilter(void);
  virtual ~StructureTensorRecursiveGaussianImageFilter(void) {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** StructureTensorRecursiveGaussianImageFilter needs all of the input
   * to produce an output. Therefore,
   * StructureTensorRecursiveGaussianImageFilter needs to provide
   * an implementation for GenerateInputRequestedRegion in order to inform
   * the pipeline execution model.
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void
  GenerateInputRequestedRegion(void) override;

  /** Generate Data */
  void
  GenerateData(void) override;

  // Override since the filter produces the entire data set
  void
  EnlargeOutputRequestedRegion(DataObject * output) override;

private:
  // purposely not implemented
  StructureTensorRecursiveGaussianImageFilter(const Self &);
  // purposely not implemented
  void
  operator=(const Self &);

  std::vector<GaussianFilterPointer> m_SmoothingFilters;
  DerivativeFilterPointer            m_DerivativeFilter;
  GaussianFilterPointer              m_TensorComponentSmoothingFilter;
  OutputImageAdaptorPointer          m_ImageAdaptor;

  /** Normalize the image across scale space */
  bool m_NormalizeAcrossScale;


  RealType m_Sigma;
  RealType m_SigmaOuter;

}; // End class StructureTensorRecursiveGaussianImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeStructureTensorRecursiveGaussianImageFilter.hxx"
#endif

// End !defined( __itktubeStructureTensorRecursiveGaussianImageFilter_h )
#endif
